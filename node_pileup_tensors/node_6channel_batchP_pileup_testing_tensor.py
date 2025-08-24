#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import math
import numpy as np
from collections import defaultdict
import re
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

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

# Globals for worker process state (used only if ProcessPool is selected)
worker_dat_file = None
worker_base_output_dir = None

# CIGAR regex/cache
CIGAR_RE = re.compile(r'(\d+)([MIDNSHPX=])')
_CIGAR_CACHE = {}

# ─────────────────────────────────────────────────────────────────────────────
# Lightweight logging / progress
def _print_progress(s: str):
    sys.stdout.write("\r\033[K" + s)
    sys.stdout.flush()

def log_once(s: str):
    sys.stderr.write(s + "\n")
    sys.stderr.flush()

# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers
def iter_idx_entries(idx_path):
    """
    Stream all entries from .idx without keeping them in memory.
    Yields tuples: (node_id, offset, n_records) in file order.
    """
    with open(idx_path, 'rb') as f:
        sz = os.fstat(f.fileno()).st_size
        if sz < 4:
            raise RuntimeError(f"Index too small: {idx_path}")
        n = struct.unpack('<I', f.read(4))[0]
        for i in range(n):
            rec = f.read(22)
            if len(rec) < 22:
                break
            node_id, offset, _unused1, n_records, _unused2 = struct.unpack('<I Q I I H', rec)
            yield (node_id, offset, n_records)

# ─────────────────────────────────────────────────────────────────────────────
# AF compaction (uint8 0..127)
def _float_to_u8_af(x: float) -> int:
    # Map AF→uint8: 0→0; else val = clamp(1..127, int(127 - 10*log10(af)))
    try:
        xf = float(x)
    except Exception:
        return 0
    if xf <= 0.0:
        return 0
    val = int(127.0 - (10.0 * math.log10(xf)))
    if val < 1:   val = 1
    if val > 127: val = 127
    return val

def af_list_to_uint8(af_list, expected_len):
    """
    Convert AF list (floats) to compact np.uint8 array length==expected_len.
    Missing/short arrays are padded with 0.
    """
    out = np.zeros(expected_len, dtype=np.uint8)
    if not af_list:
        return out
    m = min(expected_len, len(af_list))
    # Vectorized-ish loop (faster than Python list of ints for big arrays)
    for i in range(m):
        out[i] = _float_to_u8_af(af_list[i])
    return out

# ─────────────────────────────────────────────────────────────────────────────
# CIGAR / variant helpers (unchanged logic)
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
        sys.stderr.write(f"Warning: Could not parse CIGAR '{cigar_string}': {e}\n")
        return []

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def calculate_window_start(variant_pos, window_size):
    return variant_pos - (window_size // 2)

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
                            # segment_quality_values may be list or np.ndarray(uint8)
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
# Worker (process version needs an initializer; thread version just uses closures)
def init_worker(dat_file_path_for_worker, base_output_dir_for_worker):
    global worker_dat_file, worker_base_output_dir
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
    except Exception as e:
        sys.stderr.write(f"Worker init error: {e}\n"); sys.exit(1)

def process_single_node_for_pileup(task):
    """
    task: (node_id, offset, n_records, node_sequence:str, af_u8:np.ndarray|None,
           min_af, min_variants, min_allele_bq, variant_type, dat_path:str|None, base_out:str|None, use_process:bool)
    If use_process==True: read DAT via global worker_dat_file.
    If use_process==False: open a shared file descriptor per thread-group by passing dat_path and base_out and reusing (we just open per call for simplicity; OS file cache will help).
    """
    (node_id, dat_file_offset, n_records, node_sequence, af_u8,
     min_af_threshold, min_variants_threshold, min_allele_bq_threshold, variant_type_to_process,
     dat_path, base_out, use_process) = task

    # Per-exec file handle
    if use_process:
        fdat = worker_dat_file
        out_base = worker_base_output_dir
    else:
        # threads: open once per call (cheap; relies on OS page cache); avoids globals
        fdat = open(dat_path, 'rb')
        out_base = base_out

    try:
        node_len = len(node_sequence)
        node_dir = os.path.join(out_base, str(node_id))
        os.makedirs(node_dir, exist_ok=True)

        # Seek & bulk read
        fdat.seek(dat_file_offset + 10)
        buf = fdat.read(RECORD_SIZE * n_records)

        aligned_read_segments = []
        for off_from_file, raw_seq, raw_qual, raw_cigar, mapq_val, strand_byte in RECORD_STRUCT.iter_unpack(buf):
            if mapq_val < 10:
                continue
            try:
                seq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                # compact qualities
                qual_values = np.frombuffer(raw_qual.rstrip(b'\0'), dtype=np.uint8)
                cigar_str = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand_char = strand_byte.decode('ascii')
            except Exception:
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
                cur_qual = cur_qual[::-1]
                cur_ops = list(reversed(cigar_ops)) if cigar_ops else []
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

        if not aligned_read_segments:
            return node_id, 0

        # Collect candidate variants from CIGARs
        candidate_variants = defaultdict(int)
        for seg in aligned_read_segments:
            for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence
            ):
                candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

        tensors_written = 0

        for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _cnt in candidate_variants.items():
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
                tmp_af = (alt_allele_count / locus_coverage) if locus_coverage > 0 else 0.0
                if tmp_af < min_af_threshold:
                    continue
                mean_bq_tmp = float(np.mean(alt_allele_bq)) if alt_allele_bq else 0.0
                if mean_bq_tmp < min_allele_bq_threshold:
                    continue

            current_alt_freq = (alt_allele_count / locus_coverage) if locus_coverage > 0 else 0.0
            mean_alt_bq = float(np.mean(alt_allele_bq)) if alt_allele_bq else 0.0

            variant_key_string = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"
            window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
            window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)

            # Build reference row
            ref_base_indices_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
            for i, pos in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
                if 0 <= pos < node_len:
                    ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[pos].upper(), BASE_TO_INDEX['N'])

            # AF row from compact uint8 (already scaled 0..127)
            af_row = [0] * TENSOR_WINDOW_SIZE
            if af_u8 is not None and af_u8.size > 0:
                for i, pos in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
                    if 0 <= pos < af_u8.size:
                        af_row[i] = int(af_u8[pos])

            # Allocate channels
            H, W = 1 + TENSOR_MAX_READ_ROWS, TENSOR_WINDOW_SIZE
            ch1 = np.full((H, W), PADDING_BASE_INDEX, dtype=np.int8)
            ch2 = np.full((H, W), DEFAULT_QUALITY_PADDING, dtype=np.int8)
            ch3 = np.full((H, W), MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
            ch4 = np.full((H, W), DEFAULT_MAPPING_QUALITY_PADDING, dtype=np.int8)
            ch5 = np.full((H, W), CIGAR_PADDING_INDEX, dtype=np.int8)
            ch6 = np.zeros((H, W), dtype=np.int8)

            ch1[0, :] = np.asarray(ref_base_indices_row, dtype=np.int8)
            ch6[0, :] = np.asarray(af_row, dtype=np.int8)

            reads_added = 0
            variant_window_index = v_pos - window_start_pos
            for seg in aligned_read_segments:
                if reads_added >= TENSOR_MAX_READ_ROWS:
                    break
                base_row, qual_row, mapq_row, cigar_row = get_read_tensor_rows_in_window(
                    seg["cigar_ops"], seg["offset_on_node"],
                    seg["read_sequence"], seg["processed_quality_values"],
                    max(0, min(int(seg["mapping_quality"]), 127)),
                    window_start_pos, TENSOR_WINDOW_SIZE, node_len)

                if any(b != PADDING_BASE_INDEX for b in base_row):
                    r = 1 + reads_added
                    ch1[r, :] = np.asarray(base_row, dtype=np.int8)
                    ch2[r, :] = np.asarray(qual_row, dtype=np.int8)

                    # mismatch flags
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

            # Save
            arr = np.empty((6, H, W), dtype=np.int8)
            arr[0] = ch1
            arr[1] = ch2
            arr[2] = ch3
            arr[3] = ch4
            arr[4] = ch5
            arr[5] = ch6

            np.save(os.path.join(node_dir, f"{variant_key_string}.npy"), arr, allow_pickle=False)
            # Minimal sidecar; comment out if unnecessary
            with open(os.path.join(node_dir, "variant_summary.json"), 'a') as fh:
                fh.write(json.dumps({
                    "variant_key": variant_key_string,
                    "alt_allele_count": int(alt_allele_count),
                    "ref_allele_count_at_locus": int(ref_allele_count),
                    "other_allele_count_at_locus": int(other_allele_count),
                    "coverage_at_locus": int(locus_coverage),
                    "alt_allele_frequency": round(float(current_alt_freq), 4),
                    "mean_alt_allele_base_quality": round(float(mean_alt_bq), 2)
                }) + "\n")
            tensors_written += 1

        return node_id, tensors_written
    finally:
        if not use_process and fdat is not None:
            try:
                fdat.close()
            except Exception:
                pass

# ─────────────────────────────────────────────────────────────────────────────
# View (unchanged API; not generated by worker in this lean version)
def display_pileup_data(*args, **kwargs):
    pass  # intentionally omitted for the lean variant (keeps memory small)

# ─────────────────────────────────────────────────────────────────────────────
# Main with WAVES + THREADS (default)
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered tensors from alignment data (memory-lean, batched).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("dat", help=".dat alignment file")
    parser.add_argument("idx", help=".idx index file")
    parser.add_argument("output", help="Base output directory")
    parser.add_argument("candidate_variants_json",
                        help="JSON with nodes (subset) including sequences and optional genomead_af.")
    parser.add_argument("--all_nodes_json", default=None,
                        help="Optional JSON with sequences for ALL nodes (fallback; ch6=0 if no AF).")
    parser.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Parallel workers")
    parser.add_argument("--executor", choices=["thread", "process"], default="thread",
                        help="thread: share memory (low RAM); process: separate procs (higher RAM)")
    parser.add_argument("--wave_size", type=int, default=100000,
                        help="How many nodes to submit per wave")
    parser.add_argument("--chunksize", type=int, default=512,
                        help="map chunksize (mainly for process executor)")
    parser.add_argument("--min_af", type=float, default=0.1, help="Min AF for SNP to pass")
    parser.add_argument("--min_variants", type=int, default=3, help="Min alt count")
    parser.add_argument("--min_allele_bq", type=float, default=10.0, help="Min mean alt BQ")
    parser.add_argument("--variant_type", type=str, default='all', choices=['snp', 'indel', 'all'])

    args = parser.parse_args()

    for p in (args.dat, args.idx, args.candidate_variants_json):
        if not os.path.isfile(p):
            sys.exit(f"Error: missing file: {p}")
    if args.all_nodes_json and not os.path.isfile(args.all_nodes_json):
        sys.exit(f"Error: --all_nodes_json not found: {args.all_nodes_json}")
    os.makedirs(args.output, exist_ok=True)

    # Load primary JSON (kept in memory; typically a subset)
    seq_primary, af_u8_primary = {}, {}
    log_once(f"Loading nodes (primary) from {args.candidate_variants_json} ...")
    with open(args.candidate_variants_json, 'r') as f:
        data = json.load(f)
    nodes = data.get('nodes', [])
    for n in nodes:
        nid_s = n.get('node_id'); seq = n.get('sequence')
        if not (nid_s and seq): continue
        try:
            nid = int(nid_s)
        except Exception:
            continue
        seq_up = seq.upper()
        # Compact AF to uint8 (0..127)
        af_list = n.get('genomead_af', None)
        af_u8_primary[nid] = af_list_to_uint8(af_list or [], len(seq_up))
        seq_primary[nid] = seq_up
    log_once(f"Primary JSON: sequences={len(seq_primary):,}")

    # Optional fallback JSON (may be large; only store sequences, no AF)
    seq_fallback = {}
    if args.all_nodes_json:
        log_once(f"Loading fallback sequences from {args.all_nodes_json} (no AF; ch6=0) ...")
        with open(args.all_nodes_json, 'r') as f:
            data2 = json.load(f)
        for n in data2.get('nodes', []):
            nid_s = n.get('node_id'); seq = n.get('sequence')
            if not (nid_s and seq): continue
            try:
                nid = int(nid_s)
            except Exception:
                continue
            if nid not in seq_primary:
                seq_fallback[nid] = seq.upper()
        log_once(f"Fallback JSON: sequences={len(seq_fallback):,}")

    # Index iterator (streaming)
    idx_iter = iter_idx_entries(args.idx)

    # Prepare executor
    use_process = (args.executor == "process")
    Executor = ProcessPoolExecutor if use_process else ThreadPoolExecutor

    total_nodes = 0
    total_tensors = 0
    t0 = time.time()

    # Optional: process pool initializer (opens DAT once per child)
    init_kw = {}
    if use_process:
        init_kw = dict(initializer=init_worker, initargs=(args.dat, args.output))

    with Executor(max_workers=args.num_workers, **init_kw) as ex:
        wave = []
        # Build waves without holding the whole index in memory
        for node_id, offset, nrecs in idx_iter:
            # Pick sequence + AF source
            if node_id in seq_primary:
                seq = seq_primary[node_id]
                af_u8 = af_u8_primary.get(node_id, np.zeros(len(seq), dtype=np.uint8))
            else:
                seq = seq_fallback.get(node_id)
                if seq is None:
                    # No sequence at all → skip (cannot build tensors)
                    continue
                af_u8 = np.zeros(len(seq), dtype=np.uint8)  # ch6 zeros

            wave.append((
                node_id, offset, nrecs, seq, af_u8,
                args.min_af, args.min_variants, args.min_allele_bq, args.variant_type,
                (None if use_process else args.dat),
                (None if use_process else args.output),
                use_process
            ))

            if len(wave) >= args.wave_size:
                # Submit this wave
                if use_process:
                    futures = [ex.submit(process_single_node_for_pileup, t) for t in wave]
                    for fut in as_completed(futures):
                        nid, tw = fut.result()
                        total_nodes += 1
                        total_tensors += tw
                        elapsed = max(time.time() - t0, 1e-9)
                        _print_progress(f"Processed {total_nodes:,} nodes | tensors {total_tensors:,} | "
                                        f"{total_nodes/elapsed:,.1f} nodes/s | wave={args.wave_size}")
                else:
                    # threads: map is fine (no pickling)
                    for nid, tw in ex.map(process_single_node_for_pileup, wave, chunksize=args.chunksize):
                        total_nodes += 1
                        total_tensors += tw
                        elapsed = max(time.time() - t0, 1e-9)
                        _print_progress(f"Processed {total_nodes:,} nodes | tensors {total_tensors:,} | "
                                        f"{total_nodes/elapsed:,.1f} nodes/s | wave={args.wave_size}")
                wave.clear()

        # Flush remainder
        if wave:
            if use_process:
                futures = [ex.submit(process_single_node_for_pileup, t) for t in wave]
                for fut in as_completed(futures):
                    nid, tw = fut.result()
                    total_nodes += 1
                    total_tensors += tw
                    elapsed = max(time.time() - t0, 1e-9)
                    _print_progress(f"Processed {total_nodes:,} nodes | tensors {total_tensors:,} | "
                                    f"{total_nodes/elapsed:,.1f} nodes/s | wave=final")
            else:
                for nid, tw in ex.map(process_single_node_for_pileup, wave, chunksize=args.chunksize):
                    total_nodes += 1
                    total_tensors += tw
                    elapsed = max(time.time() - t0, 1e-9)
                    _print_progress(f"Processed {total_nodes:,} nodes | tensors {total_tensors:,} | "
                                    f"{total_nodes/elapsed:,.1f} nodes/s | wave=final")

    sys.stdout.write("\n")
    dt = time.time() - t0
    log_once(f"Done. Nodes processed: {total_nodes:,}; tensors: {total_tensors:,}; "
             f"elapsed: {dt:,.1f}s; rate: {total_nodes/max(dt,1e-9):,.1f} nodes/s")

if __name__ == "__main__":
    main()