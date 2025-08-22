#!/usr/bin/env python3
"""
Build tensors ONLY for nodes listed in merged.json (gz ok), AF encoded as bin8 (0..7).
If a node has NO AF -> channel 6 is 0 everywhere and ONLY the variant column is set to 3.

Memory-friendly:
- Scan .idx once and keep offsets only for wanted node IDs.
- Stream .dat per node (bulk read for that node only).
"""

import argparse
import gzip
import json
import os
import re
import struct
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Optional, Tuple

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Constants & logging
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size

BASE_TO_INDEX = {'A': 20, 'C': 30, 'G': 50, 'T': 70, 'N': 10, '*': 90, '_PADDING_': 0}
PADDING_BASE_INDEX = 0

CIGAR_OP_TO_INDEX = {'M': 10, 'N': 20, 'S': 30, 'I': 40, 'D': 50, 'H': 60, 'P': 70, '=': 80, 'X': 90, '_PADDING_': 0}
CIGAR_PADDING_INDEX = 0

INDEX_TO_BASE_FOR_VIEW = {20: 'A', 30: 'C', 50: 'G', 70: 'T', 10: 'N', 90: '*', 0: '0'}

TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200
DEFAULT_QUALITY_PADDING = 0
DEFAULT_MAPPING_QUALITY_PADDING = -1
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = -1

SENTINEL_MISSING_AF_VARIANT_COL = 3  # 无 AF 时仅把候选列置为 3

def ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def log(msg: str) -> None:
    print(f"[{ts()}] {msg}", file=sys.stderr, flush=True)

def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

# ─────────────────────────────────────────────────────────────────────────────
# Worker globals (read-only after fork)
worker_dat_file = None
worker_base_output_dir = None
worker_need_view = False

GLOBAL_NODE_SEQS: Dict[int, str] = {}
GLOBAL_NODE_AF_BINS: Dict[int, Optional[np.ndarray]] = {}  # int8 0..7 per base, or None

# ─────────────────────────────────────────────────────────────────────────────
# Helpers

def calculate_window_start(variant_pos: int, window_size: int) -> int:
    return variant_pos - (window_size // 2)

def reverse_complement(sequence: str) -> str:
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def decode_cigar_to_int_ops(cigar_string: str):
    if not cigar_string or cigar_string == '*':
        return []
    try:
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        sys.stderr.write(f"Warning: Could not parse CIGAR '{cigar_string}': {e}\n")
        return []

# ── AF bin8 ─────────────────────────────────────────────────────────────────
def af_float_to_bin8(x: float) -> int:
    try:
        xf = float(x)
    except Exception:
        return 0
    if xf <= 0.0: return 0
    if xf < 1e-6: return 0
    if xf < 1e-5: return 1
    if xf < 1e-4: return 2
    if xf < 1e-3: return 3
    if xf < 1e-2: return 4
    if xf < 0.1:  return 5
    if xf < 0.5:  return 6
    return 7

def normalize_af_to_bin8_int8_array(af_field, expected_len: int) -> Optional[np.ndarray]:
    """
    Return int8 array of length expected_len with values 0..7, or None if AF truly missing.
    - list of floats -> bin8
    - string of digits '0'..'7' -> direct digits
    - None/empty/others -> None
    """
    if af_field is None:
        return None

    if isinstance(af_field, list):
        vals = [af_float_to_bin8(x) for x in af_field[:expected_len]]
        if len(vals) < expected_len:
            vals.extend([0] * (expected_len - len(vals)))
        return np.asarray(vals, dtype=np.int8)

    if isinstance(af_field, str):
        if not af_field:
            return None
        digits = [(ord(ch) - 48) for ch in af_field[:expected_len]]
        vals = [d if 0 <= d <= 7 else 0 for d in digits]
        if len(vals) < expected_len:
            vals.extend([0] * (expected_len - len(vals)))
        return np.asarray(vals, dtype=np.int8)

    return None

def load_nodes_from_merged(merged_path: str,
                           id_key="node_id", seq_key="sequence", af_key="genomead_af"
                           ) -> Tuple[Dict[int, str], Dict[int, Optional[np.ndarray]]]:
    """
    Read merged.json[.gz] (list or {'nodes':[...]}).
    Returns:
      seqs {node_id:int -> sequence:str}
      af   {node_id:int -> int8 bin8 array, or None if AF missing}
    """
    log(f"Loading merged JSON: {merged_path}")
    with open_maybe_gzip(merged_path, "rt") as f:
        data = json.load(f)

    nodes = data.get("nodes") if isinstance(data, dict) else data
    if not isinstance(nodes, list):
        raise ValueError("Merged JSON must be a list or a dict with key 'nodes'.")

    seqs: Dict[int, str] = {}
    af_bins: Dict[int, Optional[np.ndarray]] = {}
    skipped = 0

    for n in nodes:
        if not isinstance(n, dict) or id_key not in n or seq_key not in n:
            skipped += 1
            continue
        try:
            nid = int(str(n[id_key]))
        except Exception:
            skipped += 1
            continue
        seq = str(n[seq_key]).upper()
        seqs[nid] = seq
        af_bins[nid] = normalize_af_to_bin8_int8_array(n.get(af_key), len(seq))

    log(f"Merged parsed: nodes={len(seqs):,} (skipped {skipped:,})")
    return seqs, af_bins

def select_idx_offsets(idx_path: str, wanted_ids: set) -> Dict[int, Tuple[int, int]]:
    """
    Scan .idx sequentially and keep (offset, n_records) ONLY for node IDs in wanted_ids.
    Keeps memory proportional to |wanted_ids|, not total node count.
    """
    res: Dict[int, Tuple[int, int]] = {}
    wanted = set(wanted_ids)
    found = 0

    with open(idx_path, 'rb') as f:
        header = f.read(4)
        if len(header) < 4:
            raise RuntimeError("idx too small")
        n = struct.unpack('<I', header)[0]
        for _ in range(n):
            rec = f.read(22)
            if len(rec) < 22:
                break
            node_id, offset, _, n_records, _ = struct.unpack('<I Q I I H', rec)
            if node_id in wanted:
                res[node_id] = (offset, n_records)
                found += 1
                if found == len(wanted):
                    break
    log(f"IDX scan done. matched={found:,}, missing={len(wanted)-found:,} (wanted={len(wanted):,})")
    return res

# ─────────────────────────────────────────────────────────────────────────────
# Variant helpers

def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence, read_quality_values, read_cigar_ops_decoded,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None, expected_ref_allele_for_indel=None):
    current_node_pos = read_offset_on_node
    current_read_pos = 0
    for length, op in read_cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            if current_node_pos <= target_node_pos < current_node_pos + length:
                read_idx = current_read_pos + (target_node_pos - current_node_pos)
                if read_idx < len(read_sequence):
                    allele = read_sequence[read_idx].upper()
                    quality = read_quality_values[read_idx] if read_idx < len(read_quality_values) else 0
                    if expected_var_type in ('I', 'D'):
                        return "REF_STATE_FOR_INDEL", quality
                    return allele, quality
                return None, None
            current_node_pos += length
            current_read_pos += length
        elif op == 'I':
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                if current_read_pos + length <= len(read_sequence):
                    q = read_quality_values[current_read_pos: current_read_pos + length]
                    meanq = sum(q) / len(q) if q else 0.0
                    return read_sequence[current_read_pos: current_read_pos + length].upper(), meanq
                return None, None
            current_read_pos += length
        elif op == 'D':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I':
                    return "OTHER_FOR_INDEL", None
                if expected_var_type == 'D':
                    if 0 <= current_node_pos < len(node_sequence) and current_node_pos + length <= len(node_sequence):
                        refdel = node_sequence[current_node_pos: current_node_pos + length]
                        return ("*" if refdel else "OTHER_FOR_INDEL"), None
                    return "OTHER_FOR_INDEL", None
                return "*", None
            current_node_pos += length
        elif op == 'S':
            current_read_pos += length
        elif op == 'N':
            current_node_pos += length
        if current_node_pos > target_node_pos + 1 and not (expected_var_type == 'I' and (current_node_pos - 1) <= target_node_pos):
            break
    return None, None

def detect_variants_from_cigar(offset_on_node, cigar_ops_decoded, read_sequence, node_sequence):
    variants = []
    node_pos, read_pos = offset_on_node, 0
    node_seq_len, read_seq_len = len(node_sequence), len(read_sequence)
    for length, op in cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            for i in range(length):
                n_aln, r_aln = node_pos + i, read_pos + i
                if n_aln < node_seq_len and r_aln < read_seq_len:
                    nb = node_sequence[n_aln].upper()
                    rb = read_sequence[r_aln].upper()
                    if nb != rb and op != '=':
                        variants.append((n_aln, 'X', rb, nb))
                else:
                    break
            node_pos += length
            read_pos += length
        elif op == 'I':
            ins = read_sequence[read_pos: read_pos + length].upper()
            anchor = node_pos - 1 if node_pos > 0 else 0
            ref_anchor = node_sequence[anchor].upper() if 0 <= anchor < node_seq_len else "*"
            variants.append((anchor, 'I', ins, ref_anchor))
            read_pos += length
        elif op == 'D':
            ref_del = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= node_seq_len else ""
            if ref_del:
                variants.append((node_pos, 'D', "*", ref_del))
            node_pos += length
        elif op == 'S':
            read_pos += length
        elif op == 'N':
            node_pos += length
    return variants

def get_read_representation_in_window_for_view(segment_cigar_ops, segment_offset_on_node, segment_read_sequence,
                                               window_start_node, window_size, node_len):
    window_chars = [' '] * window_size
    node_pos, read_pos = segment_offset_on_node, 0
    read_seq_len = len(segment_read_sequence)
    for L, op in segment_cigar_ops:
        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = node_pos + i, read_pos + i
                idx = n_aln - window_start_node
                if 0 <= idx < window_size and r_aln < read_seq_len:
                    window_chars[idx] = segment_read_sequence[r_aln].upper()
            node_pos += L
            read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                idx = (node_pos + i) - window_start_node
                if 0 <= idx < window_size:
                    window_chars[idx] = '*'
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L
        if node_pos >= window_start_node + window_size and read_pos > 0:
            break
        if read_pos >= read_seq_len:
            break
    return window_chars

def get_read_tensor_rows_in_window(segment_cigar_ops, segment_offset_on_node,
                                   segment_read_sequence, segment_quality_values,
                                   mapping_quality,
                                   window_start_node, tensor_win_size, node_len):
    bases = [PADDING_BASE_INDEX] * tensor_win_size
    quals = [DEFAULT_QUALITY_PADDING] * tensor_win_size
    mapqs = [DEFAULT_MAPPING_QUALITY_PADDING] * tensor_win_size
    cigar_ops_indices = [CIGAR_PADDING_INDEX] * tensor_win_size
    node_pos, read_pos = segment_offset_on_node, 0
    read_seq_len = len(segment_read_sequence)
    qual_len = len(segment_quality_values)

    for L, op in segment_cigar_ops:
        op_idx = CIGAR_OP_TO_INDEX.get(op, CIGAR_PADDING_INDEX)
        if node_pos >= window_start_node + tensor_win_size and read_pos > 0:
            break
        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = node_pos + i, read_pos + i
                idx = n_aln - window_start_node
                if 0 <= idx < tensor_win_size:
                    cigar_ops_indices[idx] = op_idx
                    mapqs[idx] = mapping_quality
                    if r_aln < read_seq_len:
                        base_char = segment_read_sequence[r_aln].upper()
                        bases[idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])
                        if r_aln < qual_len:
                            quals[idx] = segment_quality_values[r_aln]
            node_pos += L
            read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                idx = (node_pos + i) - window_start_node
                if 0 <= idx < tensor_win_size:
                    cigar_ops_indices[idx] = op_idx
                    mapqs[idx] = mapping_quality
                    bases[idx] = BASE_TO_INDEX['*']
                    quals[idx] = DEFAULT_QUALITY_PADDING
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L
        if read_pos >= read_seq_len:
            break
    return bases, quals, mapqs, cigar_ops_indices

# ─────────────────────────────────────────────────────────────────────────────
# Workers

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker, need_view_flag):
    global worker_dat_file, worker_base_output_dir, worker_need_view
    worker_dat_file = open(dat_file_path_for_worker, 'rb')
    worker_base_output_dir = base_output_dir_for_worker
    worker_need_view = bool(need_view_flag)

def process_single_node_for_pileup(task_args):
    (node_id, dat_file_offset, n_records,
     min_af_threshold, min_variants_threshold, min_allele_bq_threshold,
     variant_type_to_process) = task_args

    global worker_dat_file, worker_base_output_dir, worker_need_view
    global GLOBAL_NODE_SEQS, GLOBAL_NODE_AF_BINS

    tensor_files_generated_for_node = 0
    if worker_dat_file is None or worker_base_output_dir is None:
        return node_id, None, tensor_files_generated_for_node

    node_sequence = GLOBAL_NODE_SEQS.get(node_id, "")
    if not node_sequence:
        return node_id, {}, tensor_files_generated_for_node
    node_len = len(node_sequence)

    af_bins = GLOBAL_NODE_AF_BINS.get(node_id, None)  # int8 per-base 0..7, or None

    # Read all records for this node
    aligned_read_segments = []
    try:
        worker_dat_file.seek(dat_file_offset + 10)
        bulk = worker_dat_file.read(n_records * RECORD_SIZE)
        if len(bulk) < n_records * RECORD_SIZE:
            n_records = len(bulk) // RECORD_SIZE
            bulk = bulk[: n_records * RECORD_SIZE]

        for (off_from_file, raw_seq, raw_qual, raw_cigar, mapq_val, strand_byte) in RECORD_STRUCT.iter_unpack(bulk):
            if mapq_val < 10:
                continue
            try:
                seq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                qual_values = list(raw_qual.rstrip(b'\0'))
                cigar_str = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                continue
            if not seq or len(seq) != len(qual_values):
                continue
            cigar_ops = decode_cigar_to_int_ops(cigar_str)
            if not cigar_ops and cigar_str != '*':
                continue

            cur_seq = seq
            cur_qual = qual_values
            cur_ops = cigar_ops
            cur_off = off_from_file
            if strand_char == '-':
                cur_seq = reverse_complement(seq)
                cur_qual = qual_values[::-1]
                cur_ops = [op for op in reversed(cigar_ops)] if cigar_ops else []
                aln_span = len(cur_seq)
                cur_off = node_len - aln_span - off_from_file
                if cur_off < 0:
                    continue

            aligned_read_segments.append({
                "offset_on_node": cur_off,
                "read_sequence": cur_seq,
                "processed_quality_values": cur_qual,
                "cigar_ops": cur_ops,
                "original_cigar_str": cigar_str,
                "strand": strand_char,
                "mapping_quality": mapq_val
            })
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()} Node {node_id}]: {e}\n")
        return node_id, None, tensor_files_generated_for_node

    if not aligned_read_segments:
        return node_id, {}, tensor_files_generated_for_node

    # Collect candidate variants
    candidate_variants = defaultdict(int)
    for seg in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    view_oriented_variant_data = {} if worker_need_view else None
    variant_headers_for_summary = []

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate_variants.items():
        if variant_type_to_process == 'snp' and v_type != 'X':
            continue
        if variant_type_to_process == 'indel' and v_type not in ('I', 'D'):
            continue

        alt_allele_count = ref_allele_count = other_allele_count = locus_coverage = 0
        alt_allele_bq = []

        expected_ref = v_ref_from_cigar
        expected_alt = v_alt_from_cigar
        ref_for_indel_ctx = None

        if v_type == 'D':
            expected_alt = "*"
            if 0 <= v_pos < node_len:
                expected_ref = node_sequence[v_pos]
            ref_for_indel_ctx = v_ref_from_cigar
        elif v_type == 'I':
            expected_ref = node_sequence[v_pos] if 0 <= v_pos < node_len else "*"
            ref_for_indel_ctx = expected_ref

        for seg in aligned_read_segments:
            allele, bq = get_allele_from_read_at_node_pos(
                seg["offset_on_node"], seg["read_sequence"], seg["processed_quality_values"], seg["cigar_ops"],
                v_pos, node_sequence, v_type, ref_for_indel_ctx)
            if allele is not None:
                locus_coverage += 1
                if allele == expected_alt:
                    alt_allele_count += 1
                    if bq is not None:
                        alt_allele_bq.append(bq)
                elif allele == expected_ref or (v_type in ('I', 'D') and allele == "REF_STATE_FOR_INDEL"):
                    ref_allele_count += 1
                else:
                    other_allele_count += 1

        if alt_allele_count < min_variants_threshold:
            continue

        if v_type == 'X':
            tmp_af = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
            if tmp_af < min_af_threshold:
                continue
            mean_bq_tmp = sum(alt_allele_bq) / len(alt_allele_bq) if alt_allele_bq else 0.0
            if mean_bq_tmp < min_allele_bq_threshold:
                continue

        current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
        mean_alt_bq = sum(alt_allele_bq) / len(alt_allele_bq) if alt_allele_bq else 0.0

        variant_key_string = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"
        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)

        # ch1 ref row
        ref_base_indices_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i, pos in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
            if 0 <= pos < node_len:
                ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[pos].upper(), BASE_TO_INDEX['N'])

        # ch6 row (bin8 or missing-AF sentinel)
        ch6_row = np.zeros(TENSOR_WINDOW_SIZE, dtype=np.int8)
        if af_bins is not None:
            for i, pos in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
                if 0 <= pos < node_len and pos < len(af_bins):
                    ch6_row[i] = int(af_bins[pos])
        else:
            idx = v_pos - window_start_pos
            if 0 <= idx < TENSOR_WINDOW_SIZE:
                ch6_row[idx] = SENTINEL_MISSING_AF_VARIANT_COL

        # Pre-alloc channels
        H, W = 1 + TENSOR_MAX_READ_ROWS, TENSOR_WINDOW_SIZE
        ch1 = np.full((H, W), PADDING_BASE_INDEX, dtype=np.int8)
        ch2 = np.full((H, W), DEFAULT_QUALITY_PADDING, dtype=np.int8)
        ch3 = np.full((H, W), MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
        ch4 = np.full((H, W), DEFAULT_MAPPING_QUALITY_PADDING, dtype=np.int8)
        ch5 = np.full((H, W), CIGAR_PADDING_INDEX, dtype=np.int8)
        ch6 = np.zeros((H, W), dtype=np.int8)

        ch1[0, :] = np.asarray(ref_base_indices_row, dtype=np.int8)
        ch2[0, :] = DEFAULT_QUALITY_PADDING
        ch3[0, :] = MISMATCH_CHANNEL_REF_ROW_VALUE
        ch4[0, :] = DEFAULT_MAPPING_QUALITY_PADDING
        ch5[0, :] = CIGAR_PADDING_INDEX
        ch6[0, :] = ch6_row

        # Fill read rows
        reads_added = 0
        variant_window_index = v_pos - window_start_pos
        for seg in aligned_read_segments:
            if reads_added >= TENSOR_MAX_READ_ROWS:
                break
            mapq = max(0, min(int(seg["mapping_quality"]), 127))
            base_row, qual_row, mapq_row, cigar_row = get_read_tensor_rows_in_window(
                seg["cigar_ops"], seg["offset_on_node"],
                seg["read_sequence"], seg["processed_quality_values"],
                mapq, window_start_pos, TENSOR_WINDOW_SIZE, node_len)

            if any(b != PADDING_BASE_INDEX for b in base_row):
                r = 1 + reads_added
                ch1[r, :] = np.asarray(base_row, dtype=np.int8)
                ch2[r, :] = np.asarray(qual_row, dtype=np.int8)

                ref_row = ch1[0, :].astype(np.int16)
                read_row = ch1[r, :].astype(np.int16)
                flags = np.full(W, MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
                mask_valid = (read_row != PADDING_BASE_INDEX) & (ref_row != PADDING_BASE_INDEX)
                flags[mask_valid] = (read_row[mask_valid] != ref_row[mask_valid]).astype(np.int8)
                if 0 <= variant_window_index < W and mask_valid[variant_window_index] and flags[variant_window_index] == 1:
                    flags[variant_window_index] = 5
                ch3[r, :] = flags
                ch4[r, :] = np.asarray(mapq_row, dtype=np.int8)
                ch5[r, :] = np.asarray(cigar_row, dtype=np.int8)
                ch6[r, :] = ch6[0, :]  # 复用 AF 行

                reads_added += 1

        tensor = np.stack([ch1, ch2, ch3, ch4, ch5, ch6], axis=0)

        node_dir = os.path.join(worker_base_output_dir, str(node_id))
        os.makedirs(node_dir, exist_ok=True)
        npy_name = f"{variant_key_string}.npy"
        np.save(os.path.join(node_dir, npy_name), tensor)

        variant_headers_for_summary.append({
            "variant_key": variant_key_string, "tensor_file": npy_name,
            "alt_allele_count": alt_allele_count,
            "ref_allele_count_at_locus": ref_allele_count,
            "other_allele_count_at_locus": other_allele_count,
            "coverage_at_locus": locus_coverage,
            "alt_allele_frequency": round(current_alt_freq, 4),
            "mean_alt_allele_base_quality": round(mean_alt_bq, 2)
        })
        tensor_files_generated_for_node += 1

    if variant_headers_for_summary:
        with open(os.path.join(worker_base_output_dir, str(node_id), "variant_summary.json"), 'w') as f:
            json.dump({"node_id": node_id, "node_length": node_len,
                       "variants_passing_af_filter": variant_headers_for_summary}, f, indent=2)

    return node_id, (view_oriented_variant_data or {}), tensor_files_generated_for_node

# ─────────────────────────────────────────────────────────────────────────────
# Optional viewer
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if max_variants_to_display == 0 or not node_data_for_display_view:
        return
    print(f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))
    for i, variant_key in enumerate(sorted_variant_keys):
        if i >= max_variants_to_display:
            print(f"\n  ... ({len(sorted_variant_keys) - i} more variants not shown)")
            break
        v_pos, v_type = int(variant_key.split('_')[0]), variant_key.split('_')[1]
        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)
        print(f"\n--- Variant: {variant_key} ---")
        ref_chars = []
        for j in range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE):
            ref_chars.append(full_node_sequence[j] if 0 <= j < len(full_node_sequence) else '0')
        print(f"  Node Ref: {''.join(ref_chars)}")
        marker_pos = v_pos - window_start_pos
        marker_line = [' '] * TENSOR_WINDOW_SIZE
        if 0 <= marker_pos < TENSOR_WINDOW_SIZE:
            marker_line[marker_pos] = '^'
        print(f"  Marker  : {''.join(marker_line)}")

# ─────────────────────────────────────────────────────────────────────────────
# Main

def main():
    ap = argparse.ArgumentParser(
        description="Build tensors ONLY for nodes listed in merged.json; AF encoded as bin8 (0..7); memory-friendly idx scan.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("dat", help=".dat alignment file")
    ap.add_argument("idx", help=".idx index file")
    ap.add_argument("output", help="Base output directory")
    ap.add_argument("merged_json", help="Merged node JSON (list or {'nodes':[...]}); .gz ok")
    ap.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Worker processes")
    ap.add_argument("--chunksize", type=int, default=512, help="executor.map chunksize")
    ap.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N',
                    help="Print pileups for top N variants per node (-1 for all)")
    ap.add_argument("--max_view_reads", type=int, default=20, help="Reads per pileup in view mode")
    ap.add_argument("--min_af", type=float, default=0.1, help="Min AF for SNP variant to pass")
    ap.add_argument("--min_variants", type=int, default=3, help="Min alt count to pass")
    ap.add_argument("--min_allele_bq", type=float, default=10.0, help="Min mean base quality for alt")
    ap.add_argument("--variant_type", type=str, default='all', choices=['snp', 'indel', 'all'],
                    help="Which variants to emit")
    args = ap.parse_args()

    for p in (args.dat, args.idx, args.merged_json):
        if not os.path.isfile(p):
            sys.exit(f"Error: file not found: {p}")
    os.makedirs(args.output, exist_ok=True)

    need_view = args.view is not None

    # Load only nodes in merged.json
    node_seqs, node_af_bins = load_nodes_from_merged(args.merged_json)

    # Scan idx only for these nodes
    wanted_ids = set(node_seqs.keys())
    idx_map = select_idx_offsets(args.idx, wanted_ids)
    if not idx_map:
        sys.exit("No nodes from merged.json were found in the idx file.")

    # Prepare tasks
    tasks = []
    missing = 0
    for nid in wanted_ids:
        if nid in idx_map:
            offset, nrec = idx_map[nid]
            tasks.append((nid, offset, nrec, args.min_af, args.min_variants, args.min_allele_bq, args.variant_type))
        else:
            missing += 1
    if missing:
        log(f"Warning: {missing} node IDs from merged.json not found in idx (skipped).")
    if not tasks:
        sys.exit("No valid tasks to run after merging JSON and idx.")

    tasks.sort(key=lambda t: t[1])  # by offset for streaming access

    # Publish read-only globals
    global GLOBAL_NODE_SEQS, GLOBAL_NODE_AF_BINS
    GLOBAL_NODE_SEQS = node_seqs
    GLOBAL_NODE_AF_BINS = node_af_bins

    total_tasks = len(tasks)
    log(f"Submitting {total_tasks:,} nodes to {args.num_workers} workers...")

    total_tensors = 0
    processed = 0
    batch_nodes = 0
    batch_tensors = 0
    t0 = time.time()

    with ProcessPoolExecutor(max_workers=args.num_workers,
                             initializer=init_worker,
                             initargs=(args.dat, args.output, need_view)) as ex:
        for node_id, view_data, tensor_count in ex.map(process_single_node_for_pileup, tasks,
                                                       chunksize=max(1, args.chunksize)):
            processed += 1
            batch_nodes += 1
            total_tensors += tensor_count
            batch_tensors += tensor_count

            if need_view and view_data:
                display_pileup_data(view_data, str(node_id), GLOBAL_NODE_SEQS.get(node_id, ""),
                                    args.max_view_reads,
                                    args.view if args.view != -1 else float('inf'))

            if batch_nodes >= 1000 or processed == total_tasks:
                dt = time.time() - t0
                rate = batch_nodes / dt if dt > 0 else 0.0
                log(f"Processed {processed:,}/{total_tasks:,}  "
                    f"batch_nodes={batch_nodes:,}  tensors={batch_tensors:,}  "
                    f"in {dt:.2f}s ({rate:.1f} nodes/s)")
                batch_nodes = 0
                batch_tensors = 0
                t0 = time.time()

    log(f"Done. Nodes processed: {processed:,}; tensors: {total_tensors:,}")

if __name__ == "__main__":
    main()
