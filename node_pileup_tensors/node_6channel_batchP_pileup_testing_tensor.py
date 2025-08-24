#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import math
import gzip
import numpy as np
from collections import defaultdict
import re
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

# ─────────────────────────────────────────────────────────────────────────────
# Constants & Structs
RECORD_STRUCT = struct.Struct("<h150s150s20shc")  # off_from_file, seq, qual, cigar, mapq(h), strand(c)
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

# Globals for worker process state (ProcessPool only)
worker_dat_file = None
worker_base_output_dir = None

# CIGAR regex/cache
CIGAR_RE = re.compile(r'(\d+)([MIDNSHPX=])')
_CIGAR_CACHE = {}

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def _print_progress(s: str):
    sys.stdout.write("\r\033[K" + s)
    sys.stdout.flush()

def calculate_window_start(variant_pos, window_size):
    return variant_pos - (window_size // 2)

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def iter_idx_entries(idx_path):
    """Stream .idx entries without loading all into RAM."""
    with open(idx_path, 'rb') as f:
        if os.fstat(f.fileno()).st_size < 4:
            raise RuntimeError(f"Index too small: {idx_path}")
        n = struct.unpack('<I', f.read(4))[0]
        print(f"  Index file reports {n} total node entries. Streaming...")
        for i in range(n):
            rec = f.read(22)
            if len(rec) < 22:
                sys.stderr.write(f"Error: Index file ended prematurely at record {i + 1}.\n")
                break
            node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', rec)
            yield (node_id_from_idx, offset, n_records)
            if (i + 1) % 5_000_000 == 0:
                print(f"    Streamed {i + 1}/{n} entries...")

def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*':
        return []
    cached = _CIGAR_CACHE.get(cigar_string)
    if cached is not None:
        return cached
    try:
        ops = [(int(L), op) for L, op in CIGAR_RE.findall(cigar_string)]
        if len(_CIGAR_CACHE) < 200_000:
            _CIGAR_CACHE[cigar_string] = ops
        return ops
    except Exception as e:
        sys.stderr.write(f"Warning: Could not parse CIGAR string '{cigar_string}': {e}\n")
        return []

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
                    quals = read_quality_values[current_read_pos: current_read_pos + length]
                    mean_bq = float(np.mean(quals)) if len(quals) > 0 else 0.0
                    return read_sequence[current_read_pos: current_read_pos + length].upper(), mean_bq
                return None, None
            current_read_pos += length
        elif op == 'D':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I':
                    return "OTHER_FOR_INDEL", None
                if expected_var_type == 'D':
                    if 0 <= current_node_pos < len(node_sequence) and current_node_pos + length <= len(node_sequence):
                        del_seq = node_sequence[current_node_pos: current_node_pos + length]
                        return ("*", None) if del_seq == expected_ref_allele_for_indel else ("OTHER_FOR_INDEL", None)
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
                    node_base = node_sequence[n_aln].upper()
                    read_base = read_sequence[r_aln].upper()
                    if node_base != read_base and op != '=':
                        variants.append((n_aln, 'X', read_base, node_base))
                else:
                    break
            node_pos += length
            read_pos += length
        elif op == 'I':
            ins_seq = read_sequence[read_pos: read_pos + length].upper()
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0
            ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if 0 <= ref_anchor_pos < node_seq_len else "*"
            variants.append((ref_anchor_pos, 'I', ins_seq, ref_base_at_anchor))
            read_pos += length
        elif op == 'D':
            del_seq = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= node_seq_len else ""
            if del_seq:
                variants.append((node_pos, 'D', "*", del_seq))
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
                win_idx = n_aln - window_start_node
                if 0 <= win_idx < window_size and r_aln < read_seq_len:
                    window_chars[win_idx] = segment_read_sequence[r_aln].upper()
            node_pos += L
            read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                win_idx = (node_pos + i) - window_start_node
                if 0 <= win_idx < window_size:
                    window_chars[win_idx] = '*'
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
                win_idx = n_aln - window_start_node
                if 0 <= win_idx < tensor_win_size:
                    cigar_ops_indices[win_idx] = op_idx
                    mapqs[win_idx] = mapping_quality
                    if r_aln < read_seq_len:
                        base_char = segment_read_sequence[r_aln].upper()
                        bases[win_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])
                        if r_aln < qual_len:
                            qv = segment_quality_values[r_aln]
                            quals[win_idx] = int(qv) if isinstance(qv, (np.integer,)) else qv
            node_pos += L
            read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                win_idx = (node_pos + i) - window_start_node
                if 0 <= win_idx < tensor_win_size:
                    cigar_ops_indices[win_idx] = op_idx
                    mapqs[win_idx] = mapping_quality
                    bases[win_idx] = BASE_TO_INDEX['*']
                    quals[win_idx] = DEFAULT_QUALITY_PADDING
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L
        if read_pos >= read_seq_len:
            break
    return bases, quals, mapqs, cigar_ops_indices

# ─────────────────────────────────────────────────────────────────────────────
# AF compaction (float list → uint8 1..127; string digits → 0..7)
def _af_float_to_u8(x: float) -> int:
    try:
        xf = float(x)
    except Exception:
        return 0
    if xf <= 0.0:
        return 0
    # Inverted PHRED-like scaling into 1..127
    val = int(127.0 - (10.0 * math.log10(xf)))
    if val < 1: val = 1
    if val > 127: val = 127
    return val

def af_field_to_u8_array(af_field, expected_len):
    """
    af_field may be:
      - list of floats -> convert to 1..127
      - string of digits '0'..'7' -> copy digits (0..7)
      - None / empty -> zeros
    """
    out = np.zeros(expected_len, dtype=np.uint8)
    if af_field is None:
        return out
    if isinstance(af_field, list):
        m = min(expected_len, len(af_field))
        for i in range(m):
            out[i] = _af_float_to_u8(af_field[i])
        return out
    if isinstance(af_field, str):
        m = min(expected_len, len(af_field))
        for i in range(m):
            ch = af_field[i]
            if '0' <= ch <= '9':
                v = ord(ch) - 48
                if v < 0: v = 0
                if v > 7: v = 7
                out[i] = v
        return out
    return out

# ─────────────────────────────────────────────────────────────────────────────
# Worker
def init_worker(dat_file_path_for_worker, base_output_dir_for_worker):
    global worker_dat_file, worker_base_output_dir
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()}] init: {e}\n"); sys.exit(1)

def process_single_node_for_pileup(task_args):
    """
    task_args:
      (node_id, dat_file_offset, n_records, node_sequence:str, af_u8:np.ndarray,
       min_af, min_variants, min_allele_bq, variant_type, need_view,
       use_process:bool, dat_path_if_thread:str, base_out_if_thread:str)
    """
    (node_id, dat_file_offset, n_records, node_sequence, af_u8,
     min_af_threshold, min_variants_threshold, min_allele_bq_threshold, variant_type_to_process,
     need_view, use_process, dat_path_thr, base_out_thr) = task_args

    global worker_dat_file, worker_base_output_dir
    if use_process:
        fdat = worker_dat_file
        base_out = worker_base_output_dir
    else:
        fdat = open(dat_path_thr, 'rb')
        base_out = base_out_thr

    try:
        tensor_files_generated_for_node = 0
        if fdat is None or base_out is None:
            return node_id, None, tensor_files_generated_for_node
        if not node_sequence:
            return node_id, {}, tensor_files_generated_for_node

        node_specific_output_dir = os.path.join(base_out, str(node_id))
        os.makedirs(node_specific_output_dir, exist_ok=True)

        node_len = len(node_sequence)
        aligned_read_segments = []
        try:
            fdat.seek(dat_file_offset + 10)
            buf = fdat.read(RECORD_SIZE * n_records)
            for off_from_file, raw_seq, raw_qual, raw_cigar, mapq_val, strand_byte in RECORD_STRUCT.iter_unpack(buf):
                if mapq_val < 10:
                    continue
                try:
                    seq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                    qual_values = np.frombuffer(raw_qual.rstrip(b'\0'), dtype=np.uint8)
                    cigar_str_original = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                    strand_char = strand_byte.decode('ascii')
                except Exception:
                    continue

                if not seq or len(seq) != len(qual_values):
                    continue

                original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
                if not original_decoded_cigar_ops and cigar_str_original != '*':
                    continue

                current_read_sequence = seq
                current_quality_values = qual_values
                current_decoded_cigar_ops = original_decoded_cigar_ops
                current_offset_on_node = off_from_file

                if strand_char == '-':
                    current_read_sequence = reverse_complement(seq)
                    current_quality_values = current_quality_values[::-1]
                    current_decoded_cigar_ops = list(reversed(original_decoded_cigar_ops)) if original_decoded_cigar_ops else []
                    alignment_span_on_node = len(current_read_sequence)
                    current_offset_on_node = node_len - alignment_span_on_node - off_from_file
                    if current_offset_on_node < 0:
                        continue

                aligned_read_segments.append({
                    "offset_on_node": current_offset_on_node,
                    "read_sequence": current_read_sequence,
                    "processed_quality_values": current_quality_values,
                    "cigar_ops": current_decoded_cigar_ops,
                    "original_cigar_str": cigar_str_original,
                    "strand": strand_char,
                    "mapping_quality": mapq_val
                })
        except Exception as e:
            sys.stderr.write(f"Error [Worker {os.getpid()} for Node {node_id}]: {e}\n")
            return node_id, None, tensor_files_generated_for_node

        if not aligned_read_segments:
            return node_id, {}, tensor_files_generated_for_node

        # Candidate variants
        candidate_variants = defaultdict(int)
        for seg in aligned_read_segments:
            for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                    seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
                candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

        variant_headers_for_summary = []
        view_oriented_variant_data = {} if need_view else None

        for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate_variants.items():
            if variant_type_to_process == 'snp' and v_type != 'X':
                continue
            if variant_type_to_process == 'indel' and v_type not in ('I', 'D'):
                continue

            alt_allele_count = ref_allele_count = other_allele_count = locus_coverage = 0
            alt_allele_base_qualities = []

            expected_ref_for_af = v_ref_from_cigar
            expected_alt_for_af = v_alt_from_cigar
            ref_allele_for_indel_context = None

            if v_type == 'D':
                expected_alt_for_af = "*"
                if 0 <= v_pos < node_len:
                    expected_ref_for_af = node_sequence[v_pos]
                ref_allele_for_indel_context = v_ref_from_cigar
            elif v_type == 'I':
                expected_ref_for_af = node_sequence[v_pos] if 0 <= v_pos < node_len else "*"
                ref_allele_for_indel_context = expected_ref_for_af

            for seg in aligned_read_segments:
                allele_observed, bq = get_allele_from_read_at_node_pos(
                    seg["offset_on_node"], seg["read_sequence"], seg["processed_quality_values"], seg["cigar_ops"],
                    v_pos, node_sequence, v_type, ref_allele_for_indel_context)

                if allele_observed is not None:
                    locus_coverage += 1
                    if allele_observed == expected_alt_for_af:
                        alt_allele_count += 1
                        if bq is not None:
                            alt_allele_base_qualities.append(bq)
                    elif (allele_observed == expected_ref_for_af) or (v_type in ('I','D') and allele_observed == "REF_STATE_FOR_INDEL"):
                        ref_allele_count += 1
                    else:
                        other_allele_count += 1

            if alt_allele_count < min_variants_threshold:
                continue
            if v_type == 'X':
                current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
                if current_alt_freq < min_af_threshold:
                    continue
                mean_alt_bq = float(np.mean(alt_allele_base_qualities)) if alt_allele_base_qualities else 0.0
                if mean_alt_bq < min_allele_bq_threshold:
                    continue

            current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
            mean_alt_bq = float(np.mean(alt_allele_base_qualities)) if alt_allele_base_qualities else 0.0

            variant_key_string = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"
            window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
            window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)

            # Reference top row
            ref_base_indices_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
            for i, node_pos_in_window in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
                if 0 <= node_pos_in_window < node_len:
                    ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[node_pos_in_window].upper(), BASE_TO_INDEX['N'])

            # Channel 6 AF row from compact uint8 (0..127)
            genomead_af_row = [0] * TENSOR_WINDOW_SIZE
            if af_u8 is not None and af_u8.size > 0:
                for i, node_pos_in_window in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
                    if 0 <= node_pos_in_window < af_u8.size:
                        genomead_af_row[i] = int(af_u8[node_pos_in_window])

            # Allocate channels
            H, W = 1 + TENSOR_MAX_READ_ROWS, TENSOR_WINDOW_SIZE
            ch1 = np.full((H, W), PADDING_BASE_INDEX, dtype=np.int8)
            ch2 = np.full((H, W), DEFAULT_QUALITY_PADDING, dtype=np.int8)
            ch3 = np.full((H, W), MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
            ch4 = np.full((H, W), DEFAULT_MAPPING_QUALITY_PADDING, dtype=np.int8)
            ch5 = np.full((H, W), CIGAR_PADDING_INDEX, dtype=np.int8)
            ch6 = np.zeros((H, W), dtype=np.int8)

            ch1[0, :] = np.asarray(ref_base_indices_row, dtype=np.int8)
            ch6[0, :] = np.asarray(genomead_af_row, dtype=np.int8)

            reads_added = 0
            variant_window_index = v_pos - window_start_pos
            for seg_data in aligned_read_segments:
                if reads_added >= TENSOR_MAX_READ_ROWS:
                    break
                base_idx_row, quality_row, mapq_row, cigar_row = get_read_tensor_rows_in_window(
                    seg_data["cigar_ops"], seg_data["offset_on_node"], seg_data["read_sequence"],
                    seg_data["processed_quality_values"], max(0, min(int(seg_data["mapping_quality"]), 127)),
                    window_start_pos, TENSOR_WINDOW_SIZE, node_len)

                if any(b != PADDING_BASE_INDEX for b in base_idx_row):
                    r = 1 + reads_added
                    ch1[r, :] = np.asarray(base_idx_row, dtype=np.int8)
                    ch2[r, :] = np.asarray(quality_row, dtype=np.int8)

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
                    ch6[r, :] = ch6[0, :]
                    reads_added += 1

            arr = np.empty((6, H, W), dtype=np.int8)
            arr[0] = ch1; arr[1] = ch2; arr[2] = ch3; arr[3] = ch4; arr[4] = ch5; arr[5] = ch6
            np.save(os.path.join(node_specific_output_dir, f"{variant_key_string}.npy"), arr, allow_pickle=False)

            variant_headers_for_summary.append({
                "variant_key": variant_key_string, "tensor_file": f"{variant_key_string}.npy",
                "alt_allele_count": int(alt_allele_count),
                "ref_allele_count_at_locus": int(ref_allele_count),
                "other_allele_count_at_locus": int(other_allele_count),
                "coverage_at_locus": int(locus_coverage),
                "alt_allele_frequency": round(float(current_alt_freq), 4),
                "mean_alt_allele_base_quality": round(float(mean_alt_bq), 2)
            })
            tensor_files_generated_for_node += 1

        if variant_headers_for_summary:
            with open(os.path.join(node_specific_output_dir, "variant_summary.json"), 'w') as f:
                json.dump({"node_id": node_id, "node_length": node_len,
                           "variants_passing_af_filter": variant_headers_for_summary}, f, indent=2)

        return node_id, view_oriented_variant_data, tensor_files_generated_for_node
    finally:
        if not use_process and fdat is not None:
            try: fdat.close()
            except Exception: pass

# ─────────────────────────────────────────────────────────────────────────────
# View
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf'),
                        show_empty_info=False):
    if max_variants_to_display == 0:
        return
    if not node_data_for_display_view:
        if show_empty_info:
            print(f"Info: No pileup data for node {node_id_str_for_display} (no variants met all filters).")
        return
    print(f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))
    for i, variant_key in enumerate(sorted_variant_keys):
        if i >= max_variants_to_display:
            print(f"\n  ... ({len(sorted_variant_keys) - i} more variants not shown due to --view limit)")
            break
        variant_data = node_data_for_display_view[variant_key]
        v_pos, v_type = int(variant_key.split('_')[0]), variant_key.split('_')[1]
        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)
        print(f"\n--- Variant: {variant_key} ---")
        ref_chars = [full_node_sequence[j] if 0 <= j < len(full_node_sequence) else '0'
                     for j in range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)]
        print(f"  Node Ref: {''.join(ref_chars)}")
        marker_pos_in_window = v_pos - window_start_pos
        marker_line = [' '] * TENSOR_WINDOW_SIZE
        if 0 <= marker_pos_in_window < TENSOR_WINDOW_SIZE:
            marker_line[marker_pos_in_window] = '^'
        print(f"  Marker  : {''.join(marker_line)}")
        for j, read_entry in enumerate(variant_data.get("pileup_reads_data", [])):
            if j >= max_reads_to_display_per_variant:
                print(f"  ... ({len(variant_data.get('pileup_reads_data', [])) - j} more reads not shown)")
                break
            bases_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in read_entry["bases"]])
            print(f"  Read {j + 1:3d}: {bases_str} (CIGAR:{read_entry['cigar']})")
        print(f"  Alt Count: {variant_data.get('alt_allele_count', 'N/A')}, "
              f"Ref Count: {variant_data.get('ref_allele_count_at_locus', 'N/A')}, "
              f"Coverage: {variant_data.get('coverage_at_locus', 'N/A')}")
        print(f"  Alt Freq: {variant_data.get('alt_allele_frequency', 0.0):.4f}, "
              f"Mean Alt BQ: {variant_data.get('mean_alt_allele_base_quality', 0.0):.2f}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# JSON loader that accepts list or {"nodes":[...]} (and .gz)
def load_nodes_json(path):
    with open_maybe_gzip(path, "rt") as f:
        data = json.load(f)
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        if "nodes" in data and isinstance(data["nodes"], list):
            return data["nodes"]
    raise ValueError("JSON must be a list or a dict with key 'nodes' (list).")

# ─────────────────────────────────────────────────────────────────────────────
# Main
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered tensors from alignment data (batched submission for low RAM).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("dat", help=".dat alignment file")
    parser.add_argument("idx", help=".idx index file")
    parser.add_argument("output", help="Base output directory")
    parser.add_argument("candidate_variants_json",
                        help="Primary JSON with nodes (subset) including sequences and optional genomead_af (.gz ok).")
    parser.add_argument("--all_nodes_json", default=None,
                        help="Optional JSON with sequences for ALL nodes (node_id, sequence) to process nodes not in primary JSON (their ch6 will be zeros). .gz ok.")
    parser.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Number of workers")
    parser.add_argument("--executor", choices=["process", "thread"], default="process",
                        help="Use process (default) or thread pool. Thread saves RAM (no pickling).")
    parser.add_argument("--wave_size", type=int, default=100000,
                        help="Maximum number of in-flight tasks (futures) at a time.")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N',
                        help="Print pileups for top N variants per node (-1 for all). Omit for speed.")
    parser.add_argument("--max_view_reads", type=int, default=20, help="Max reads to show per pileup in view mode")
    parser.add_argument("--min_af", type=float, default=0.1, help="Minimum allele frequency to process a SNP")
    parser.add_argument("--min_variants", type=int, default=3, help="Alt allele count must be >= this value")
    parser.add_argument("--min_allele_bq", type=float, default=10.0, help="Minimum mean base quality of alt alleles")
    parser.add_argument("--variant_type", type=str, default='all', choices=['snp', 'indel', 'all'],
                        help="Variant types to output tensors for: 'snp', 'indel', or 'all'")
    parser.add_argument("--show_empty_info", action="store_true",
                        help="If set, print 'No pileup data' messages for nodes with no passing variants.")
    args = parser.parse_args()

    if not (os.path.isfile(args.dat) and os.path.isfile(args.idx) and os.path.isfile(args.candidate_variants_json)):
        sys.exit("Error: One or more input files (dat, idx, or json) were not found.")
    if args.all_nodes_json and not os.path.isfile(args.all_nodes_json):
        sys.exit(f"Error: --all_nodes_json file not found: {args.all_nodes_json}")
    os.makedirs(args.output, exist_ok=True)

    # Load primary JSON (list or dict-with-nodes)
    node_sequences_primary, node_af_u8_primary = {}, {}
    print(f"Loading nodes from {args.candidate_variants_json}...")
    try:
        nodes = load_nodes_json(args.candidate_variants_json)
        for node_obj in nodes:
            # accept 'node_id' or 'id'
            node_id_val = node_obj.get('node_id', node_obj.get('id'))
            sequence = node_obj.get('sequence')
            af_field = node_obj.get('genomead_af', None)
            if node_id_val is None or sequence is None:
                continue
            try:
                nid = int(str(node_id_val))
            except Exception:
                continue
            seq_up = str(sequence).upper()
            node_sequences_primary[nid] = seq_up
            node_af_u8_primary[nid] = af_field_to_u8_array(af_field, len(seq_up))
    except Exception as e:
        sys.exit(f"Error reading or parsing JSON file: {e}")
    print(f"Primary JSON provides sequences for {len(node_sequences_primary)} node IDs.")

    # Optional fallback JSON (only sequences)
    node_sequences_fallback = {}
    if args.all_nodes_json:
        print(f"Loading fallback sequences from {args.all_nodes_json}...")
        try:
            nodes_f = load_nodes_json(args.all_nodes_json)
            for node_obj in nodes_f:
                node_id_val = node_obj.get('node_id', node_obj.get('id'))
                sequence = node_obj.get('sequence')
                if node_id_val is None or sequence is None:
                    continue
                try:
                    nid = int(str(node_id_val))
                except Exception:
                    continue
                if nid not in node_sequences_primary:
                    node_sequences_fallback[nid] = str(sequence).upper()
            print(f"Fallback JSON provides sequences for {len(node_sequences_fallback)} node IDs.")
        except Exception as e:
            sys.exit(f"Error reading or parsing --all_nodes_json: {e}")

    # For optional view printing
    node_seq_map = {}
    node_seq_map.update(node_sequences_fallback)
    node_seq_map.update(node_sequences_primary)

    # Executor selection
    use_process = (args.executor == "process")
    Pool = ProcessPoolExecutor if use_process else ThreadPoolExecutor
    init_kw = {}
    if use_process:
        init_kw = dict(initializer=init_worker, initargs=(args.dat, args.output))

    total_tensors = 0
    processed_nodes = 0
    start_time = time.time()

    print(f"\nSubmitting tasks in waves of {args.wave_size} to {args.num_workers} {args.executor}(s)...")
    idx_stream = iter_idx_entries(args.idx)

    def make_task(node_id, offset, n_records):
        if node_id in node_sequences_primary:
            seq = node_sequences_primary[node_id]
            af_u8 = node_af_u8_primary.get(node_id, np.zeros(len(seq), dtype=np.uint8))
        else:
            seq = node_sequences_fallback.get(node_id)
            if seq is None:
                return None
            af_u8 = np.zeros(len(seq), dtype=np.uint8)
        return (node_id, offset, n_records, seq, af_u8,
                args.min_af, args.min_variants, args.min_allele_bq, args.variant_type,
                (args.view is not None), use_process,
                (None if use_process else args.dat), (None if use_process else args.output))

    with Pool(max_workers=args.num_workers, **init_kw) as ex:
        futures = []
        # prime first wave
        for node_id, offset, n_records in idx_stream:
            task = make_task(node_id, offset, n_records)
            if task is None:
                continue
            futures.append(ex.submit(process_single_node_for_pileup, task))
            if len(futures) >= args.wave_size:
                break

        # process waves
        while futures:
            for fut in as_completed(futures):
                try:
                    node_id, view_data, tensor_count = fut.result()
                    total_tensors += tensor_count
                    processed_nodes += 1

                    if args.view is not None and view_data:
                        node_seq_for_view = node_seq_map.get(node_id, "")
                        display_pileup_data(
                            view_data, str(node_id), node_seq_for_view,
                            args.max_view_reads,
                            (args.view if args.view != -1 else float('inf')),
                            show_empty_info=args.show_empty_info
                        )
                except Exception as e:
                    sys.stdout.write("\n")
                    print(f"Error processing a node: {e}", file=sys.stderr)

                elapsed = max(time.time() - start_time, 1e-9)
                rate = processed_nodes / elapsed
                _print_progress(f"Processed {processed_nodes:,} nodes | {rate:,.1f} nodes/s | tensors: {total_tensors:,} | elapsed: {elapsed:,.1f}s")

                # top up the wave
                try:
                    node_id, offset, n_records = next(idx_stream)
                    task = make_task(node_id, offset, n_records)
                    if task is not None:
                        futures.append(ex.submit(process_single_node_for_pileup, task))
                except StopIteration:
                    pass

                futures.remove(fut)

    sys.stdout.write("\n")
    print(f"Processing complete. Total nodes: {processed_nodes:,} | Total tensors: {total_tensors:,} | Elapsed: {time.time()-start_time:,.1f}s")

if __name__ == '__main__':
    main()
