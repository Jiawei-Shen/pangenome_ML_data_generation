#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np
from collections import defaultdict
import re
from concurrent.futures import ProcessPoolExecutor

# ─────────────────────────────────────────────────────────────────────────────
# Constants
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

# Globals for worker process state
worker_dat_file = None
worker_base_output_dir = None
worker_need_view = False  # gate view building

# Read-only globals (inherited by workers via fork COW)
GLOBAL_NODE_SEQS = {}
GLOBAL_NODE_AF = {}

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
def calculate_window_start(variant_pos, window_size):
    return variant_pos - (window_size // 2)

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def load_full_idx_data(idx_path):
    idx_data_map = {}
    print(f"Loading full index data from {idx_path}...")
    try:
        with open(idx_path, 'rb') as f:
            if os.fstat(f.fileno()).st_size < 4:
                sys.stderr.write(f"Error: Index file {idx_path} is too small.\n")
                return None
            num_nodes_in_idx = struct.unpack('<I', f.read(4))[0]
            print(f"  Index file reports {num_nodes_in_idx} total node entries. Reading all entries...")
            for i in range(num_nodes_in_idx):
                rec = f.read(22)
                if len(rec) < 22:
                    sys.stderr.write(f"Error: Index file ended prematurely at record {i + 1}.\n")
                    break
                node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', rec)
                idx_data_map[node_id_from_idx] = (offset, n_records)
                if (i + 1) % 5_000_000 == 0:
                    print(f"    Loaded {i + 1}/{num_nodes_in_idx} index entries...")
        print(f"Successfully loaded {len(idx_data_map)} distinct node entries.")
        return idx_data_map
    except FileNotFoundError:
        sys.stderr.write(f"Error: Index file not found at {idx_path}\n")
        return None
    except Exception as e:
        sys.stderr.write(f"Error parsing full index file {idx_path}: {e}\n")
        return None

_CIGAR_RE = re.compile(r'(\d+)([MIDNSHPX=])')
def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*':
        return []
    try:
        return [(int(L), op) for L, op in _CIGAR_RE.findall(cigar_string)]
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
            current_node_pos += length; current_read_pos += length
        elif op == 'I':
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                if current_read_pos + length <= len(read_sequence):
                    quals = read_quality_values[current_read_pos: current_read_pos + length]
                    mean_bq = (sum(quals) / len(quals)) if quals else 0.0
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
            node_pos += length; read_pos += length
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
            node_pos += L; read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                win_idx = (node_pos + i) - window_start_node
                if 0 <= win_idx < window_size:
                    window_chars[win_idx] = '*'
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L
        if node_pos >= window_start_node + window_size and read_pos > 0: break
        if read_pos >= read_seq_len: break
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
                            quals[win_idx] = segment_quality_values[r_aln]
            node_pos += L; read_pos += L
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
        if read_pos >= read_seq_len: break
    return bases, quals, mapqs, cigar_ops_indices

# ─────────────────────────────────────────────────────────────────────────────
# Worker init & task
def init_worker(dat_file_path_for_worker, base_output_dir_for_worker, need_view_flag):
    global worker_dat_file, worker_base_output_dir, worker_need_view
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
        worker_need_view = bool(need_view_flag)
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()}] init: {e}\n"); sys.exit(1)

def process_single_node_for_pileup(task_args):
    (node_id, dat_file_offset, n_records,
     min_af_threshold, min_variants_threshold, min_allele_bq_threshold,
     variant_type_to_process) = task_args

    global worker_dat_file, worker_base_output_dir, worker_need_view
    global GLOBAL_NODE_SEQS, GLOBAL_NODE_AF

    tensor_files_generated_for_node = 0
    if worker_dat_file is None or worker_base_output_dir is None:
        return node_id, None, tensor_files_generated_for_node

    node_sequence = GLOBAL_NODE_SEQS.get(node_id, "")
    if not node_sequence:
        return node_id, {}, tensor_files_generated_for_node
    genomead_af_list = GLOBAL_NODE_AF.get(node_id, [])

    node_len = len(node_sequence)
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
                cigar_str_original = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                continue
            if not seq or len(seq) != len(qual_values):
                continue

            cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
            if not cigar_ops and cigar_str_original != '*':
                continue

            cur_seq, cur_qual, cur_ops, cur_off = seq, qual_values, cigar_ops, off_from_file
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
                "original_cigar_str": cigar_str_original,
                "strand": strand_char,
                "mapping_quality": mapq_val
            })
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()} Node {node_id}]: {e}\n")
        return node_id, None, tensor_files_generated_for_node

    if not aligned_read_segments:
        return node_id, {}, tensor_files_generated_for_node

    candidate_variants = defaultdict(int)
    for seg in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    view_oriented_variant_data = {} if worker_need_view else None
    variant_headers_for_summary = []

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate_variants.items():
        if variant_type_to_process == 'snp' and v_type != 'X': continue
        if variant_type_to_process == 'indel' and v_type not in ('I', 'D'): continue

        alt_allele_count = ref_allele_count = other_allele_count = locus_coverage = 0
        alt_allele_bq = []

        expected_ref = v_ref_from_cigar
        expected_alt = v_alt_from_cigar
        ref_for_indel_ctx = None

        if v_type == 'D':
            expected_alt = "*"
            if 0 <= v_pos < node_len: expected_ref = node_sequence[v_pos]
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
                    if bq is not None: alt_allele_bq.append(bq)
                elif allele == expected_ref or (v_type in ('I','D') and allele == "REF_STATE_FOR_INDEL"):
                    ref_allele_count += 1
                else:
                    other_allele_count += 1

        if alt_allele_count < min_variants_threshold: continue
        if v_type == 'X':
            tmp_af = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
            if tmp_af < min_af_threshold: continue
            mean_bq_tmp = sum(alt_allele_bq) / len(alt_allele_bq) if alt_allele_bq else 0.0
            if mean_bq_tmp < min_allele_bq_threshold: continue

        current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
        mean_alt_bq = sum(alt_allele_bq) / len(alt_allele_bq) if alt_allele_bq else 0.0

        variant_key_string = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"
        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)

        # Reference top row
        ref_base_indices_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i, node_pos_in_window in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
            if 0 <= node_pos_in_window < node_len:
                ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[node_pos_in_window].upper(), BASE_TO_INDEX['N'])

        # Channel 6 AF row (zeros if AF not provided)
        genomead_af_row = [0] * TENSOR_WINDOW_SIZE
        if genomead_af_list:
            for i, node_pos_in_window in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
                if 0 <= node_pos_in_window < node_len:
                    af_value = genomead_af_list[node_pos_in_window]
                    if af_value == 0.0: scaled = 0
                    else:
                        phred_af = -10.0 * np.log10(af_value)
                        inverted = 127.0 - phred_af
                        scaled = max(1, min(int(inverted), 127))
                    genomead_af_row[i] = scaled

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
        for seg in aligned_read_segments:
            if reads_added >= TENSOR_MAX_READ_ROWS: break
            mapq = max(0, min(int(seg["mapping_quality"]), 127))
            base_row, qual_row, mapq_row, cigar_row = get_read_tensor_rows_in_window(
                seg["cigar_ops"], seg["offset_on_node"], seg["read_sequence"],
                seg["processed_quality_values"], mapq,
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

        tensor = np.stack([ch1, ch2, ch3, ch4, ch5, ch6], axis=0)  # (6, H, W)
        node_dir = os.path.join(worker_base_output_dir, str(node_id))
        os.makedirs(node_dir, exist_ok=True)
        np.save(os.path.join(node_dir, f"{variant_key_string}.npy"), tensor)

        variant_headers_for_summary.append({
            "variant_key": variant_key_string, "tensor_file": f"{variant_key_string}.npy",
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
# View
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if max_variants_to_display == 0: return
    if not node_data_for_display_view:
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
# Main (waves of 100k)
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered tensors from alignment data (batched waves).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dat", help=".dat alignment file")
    parser.add_argument("idx", help=".idx index file")
    parser.add_argument("output", help="Base output directory")
    parser.add_argument("candidate_variants_json",
                        help="JSON file with nodes. Supports list or {'nodes': [...]}.")
    parser.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Worker processes")
    parser.add_argument("--chunksize", type=int, default=512, help="executor.map chunksize")
    parser.add_argument("--nodes_per_wave", type=int, default=100_000,
                        help="Submit this many nodes per wave to the pool.")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N',
                        help="Print pileups for top N variants per node (-1 for all). Omit for speed.")
    parser.add_argument("--max_view_reads", type=int, default=20, help="Reads per pileup in view mode")
    parser.add_argument("--min_af", type=float, default=0.1, help="Min AF for SNP variant to pass")
    parser.add_argument("--min_variants", type=int, default=3, help="Min alt count to pass")
    parser.add_argument("--min_allele_bq", type=float, default=10.0, help="Min mean base quality for alt")
    parser.add_argument("--variant_type", type=str, default='all', choices=['snp', 'indel', 'all'],
                        help="Variant types to emit")
    args = parser.parse_args()

    if not all([os.path.isfile(args.dat), os.path.isfile(args.idx), os.path.isfile(args.candidate_variants_json)]):
        sys.exit("Error: One or more input files (dat, idx, or json) were not found.")
    os.makedirs(args.output, exist_ok=True)

    need_view = args.view is not None

    # Load candidate nodes (supports list or {'nodes': [...]})
    node_sequences = {}
    node_af_data = {}
    node_ids_to_process = set()
    print(f"Loading nodes from {args.candidate_variants_json}...")
    try:
        with open(args.candidate_variants_json, 'r') as f:
            data = json.load(f)
        nodes = data.get('nodes', data) if isinstance(data, dict) else data
        if not isinstance(nodes, list):
            raise ValueError("JSON must be a list or a dict with key 'nodes'.")
        for node_obj in nodes:
            node_id_str = node_obj.get('node_id')
            sequence = node_obj.get('sequence')
            af_list = node_obj.get('genomead_af', [])
            if node_id_str and sequence:
                try:
                    nid = int(node_id_str)
                    node_sequences[nid] = sequence.upper()
                    node_af_data[nid] = af_list
                    node_ids_to_process.add(nid)
                except ValueError:
                    pass
    except Exception as e:
        sys.exit(f"Error reading or parsing JSON file: {e}")
    print(f"Found {len(node_sequences)} nodes with integer-compatible IDs to process.")

    full_idx_data = load_full_idx_data(args.idx)
    if not full_idx_data:
        sys.exit("Failed to load index data.")

    # Prepare minimal tasks (no big payloads)
    tasks = []
    missing = 0
    for node_id in node_ids_to_process:
        if node_id in full_idx_data:
            offset, n_records = full_idx_data[node_id]
            tasks.append((node_id, offset, n_records,
                          args.min_af, args.min_variants, args.min_allele_bq, args.variant_type))
        else:
            missing += 1
    if missing:
        print(f"Warning: {missing} node IDs from JSON were not found in the index file and will be skipped.")
    if not tasks:
        sys.exit("No valid tasks to run after processing JSON and index file.")

    # Sort by DAT offset for sequential I/O
    tasks.sort(key=lambda t: t[1])

    # Publish globals for workers (COW)
    global GLOBAL_NODE_SEQS, GLOBAL_NODE_AF
    GLOBAL_NODE_SEQS = node_sequences
    GLOBAL_NODE_AF = node_af_data

    total_tasks = len(tasks)
    print(f"\nSubmitting {total_tasks:,} tasks to {args.num_workers} workers in waves of {args.nodes_per_wave:,}...")

    total_tensors = 0
    processed_total = 0
    t0 = time.time()

    # Create executor once; feed waves to the same pool
    with ProcessPoolExecutor(max_workers=args.num_workers,
                             initializer=init_worker,
                             initargs=(args.dat, args.output, need_view)) as executor:
        # iterate waves of size nodes_per_wave
        for wave_start in range(0, total_tasks, args.nodes_per_wave):
            wave_tasks = tasks[wave_start: wave_start + args.nodes_per_wave]
            wave_n = len(wave_tasks)
            wave_processed = 0
            wave_tensors = 0
            wave_t0 = time.time()

            for node_id, view_data, tensor_count in executor.map(
                    process_single_node_for_pileup, wave_tasks, chunksize=max(1, args.chunksize)):
                processed_total += 1
                wave_processed += 1
                total_tensors += tensor_count
                wave_tensors += tensor_count

                if need_view and view_data:
                    node_seq_for_view = GLOBAL_NODE_SEQS.get(node_id, "")
                    display_pileup_data(
                        view_data, str(node_id), node_seq_for_view,
                        args.max_view_reads,
                        (args.view if args.view != -1 else float('inf'))
                    )

                # Optional: light progress inside wave every 1k
                if (wave_processed % 1000) == 0 or wave_processed == wave_n:
                    dtw = max(time.time() - wave_t0, 1e-9)
                    ratew = wave_processed / dtw
                    print(f"  Wave {wave_start//args.nodes_per_wave + 1}: "
                          f"{wave_processed:,}/{wave_n:,} nodes | {ratew:,.1f} nodes/s | "
                          f"wave tensors: {wave_tensors:,}")

            # wave summary
            dtw = time.time() - wave_t0
            print(f"Wave done ({wave_start}-{wave_start+wave_n-1}). "
                  f"Nodes: {wave_n:,}, tensors: {wave_tensors:,}, time: {dtw:.1f}s")

    dt_total = time.time() - t0
    print(f"\nProcessing complete. Nodes: {processed_total:,} | "
          f"Tensors: {total_tensors:,} | Elapsed: {dt_total:,.1f}s")

if __name__ == '__main__':
    main()
