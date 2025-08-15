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
from concurrent.futures import ProcessPoolExecutor, as_completed
import torch
import shutil, sys

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

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
def _print_progress(s: str):
    # Clear the current line then write new status. Works on most terminals.
    sys.stdout.write("\r\033[K" + s)
    sys.stdout.flush()

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
                sys.stderr.write(f"Error: Index file {idx_path} is too small.\n"); return None
            num_nodes_in_idx = struct.unpack('<I', f.read(4))[0]
            print(f"  Index file reports {num_nodes_in_idx} total node entries. Reading all entries...")
            for i in range(num_nodes_in_idx):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                    sys.stderr.write(f"Error: Index file ended prematurely at record {i + 1}.\n")
                    break
                node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                idx_data_map[node_id_from_idx] = (offset, n_records)
                if (i+1) % 5_000_000 == 0:
                    print(f"    Loaded {i+1}/{num_nodes_in_idx} index entries...")
        print(f"Successfully loaded {len(idx_data_map)} distinct node entries.")
        return idx_data_map
    except FileNotFoundError:
        sys.stderr.write(f"Error: Index file not found at {idx_path}\n"); return None
    except Exception as e:
        sys.stderr.write(f"Error parsing full index file {idx_path}: {e}\n"); return None

def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*': return []
    try:
        return [(int(L), op) for L, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        sys.stderr.write(f"Warning: Could not parse CIGAR string '{cigar_string}': {e}\n"); return []

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
                    mean_bq = sum(quals)/len(quals) if quals else 0.0
                    return read_sequence[current_read_pos: current_read_pos + length].upper(), mean_bq
                return None, None
            current_read_pos += length
        elif op == 'D':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I': return "OTHER_FOR_INDEL", None
                if expected_var_type == 'D':
                    if 0 <= current_node_pos < len(node_sequence) and current_node_pos + length <= len(node_sequence):
                        del_seq = node_sequence[current_node_pos: current_node_pos + length]
                        return ("*", None) if del_seq == expected_ref_allele_for_indel else ("OTHER_FOR_INDEL", None)
                    return "OTHER_FOR_INDEL", None
                return "*", None
            current_node_pos += length
        elif op == 'S': current_read_pos += length
        elif op == 'N': current_node_pos += length
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
            if del_seq: variants.append((node_pos, 'D', "*", del_seq))
            node_pos += length
        elif op == 'S': read_pos += length
        elif op == 'N': node_pos += length
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
        if node_pos >= window_start_node + tensor_win_size and read_pos > 0: break

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
                        if r_aln < qual_len: quals[win_idx] = segment_quality_values[r_aln]
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
# Worker
def init_worker(dat_file_path_for_worker, base_output_dir_for_worker):
    global worker_dat_file, worker_base_output_dir
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
    except FileNotFoundError:
        sys.stderr.write(f"Error [Worker {os.getpid()}]: DAT file not found.\n"); sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()}] opening DAT file: {e}\n"); sys.exit(1)

def process_single_node_for_pileup(task_args):
    (node_id, dat_file_offset, n_records, node_sequence, genomead_af_list,
     min_af_threshold, min_variants_threshold, min_allele_bq_threshold, variant_type_to_process) = task_args
    global worker_dat_file, worker_base_output_dir

    tensor_files_generated_for_node = 0
    if worker_dat_file is None or worker_base_output_dir is None:
        return node_id, None, tensor_files_generated_for_node
    if not node_sequence:
        return node_id, {}, tensor_files_generated_for_node

    node_specific_output_dir = os.path.join(worker_base_output_dir, str(node_id))
    os.makedirs(node_specific_output_dir, exist_ok=True)

    node_len = len(node_sequence)
    aligned_read_segments = []
    try:
        worker_dat_file.seek(dat_file_offset + 10)
        for _ in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE: break

            off_from_file, raw_seq, raw_qual, raw_cigar, mapq_val, strand_byte = RECORD_STRUCT.unpack(data)
            if mapq_val < 10: continue

            try:
                seq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                qual_values = list(raw_qual.rstrip(b'\0'))
                cigar_str_original = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                continue

            if not seq or len(seq) != len(qual_values): continue

            original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
            if not original_decoded_cigar_ops and cigar_str_original != '*': continue

            current_read_sequence = seq
            current_quality_values = qual_values
            current_decoded_cigar_ops = original_decoded_cigar_ops
            current_offset_on_node = off_from_file

            if strand_char == '-':
                current_read_sequence = reverse_complement(seq)
                current_quality_values = qual_values[::-1]
                current_decoded_cigar_ops = list(reversed(original_decoded_cigar_ops)) if original_decoded_cigar_ops else []
                alignment_span_on_node = len(current_read_sequence)
                current_offset_on_node = node_len - alignment_span_on_node - off_from_file
                if current_offset_on_node < 0: continue

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

    # Candidate variants from CIGAR
    candidate_variants = defaultdict(int)
    for seg in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    variant_headers_for_summary = []
    view_oriented_variant_data = {}

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate_variants.items():
        if variant_type_to_process == 'snp' and v_type != 'X': continue
        if variant_type_to_process == 'indel' and v_type not in ('I', 'D'): continue

        alt_allele_count = ref_allele_count = other_allele_count = locus_coverage = 0
        alt_allele_base_qualities = []

        expected_ref_for_af = v_ref_from_cigar
        expected_alt_for_af = v_alt_from_cigar
        ref_allele_for_indel_context = None

        if v_type == 'D':
            expected_alt_for_af = "*"
            if 0 <= v_pos < node_len: expected_ref_for_af = node_sequence[v_pos]
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
                    if bq is not None: alt_allele_base_qualities.append(bq)
                elif (allele_observed == expected_ref_for_af) or (v_type in ('I','D') and allele_observed == "REF_STATE_FOR_INDEL"):
                    ref_allele_count += 1
                else:
                    other_allele_count += 1

        if alt_allele_count < min_variants_threshold: continue
        if v_type == 'X':
            current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
            if current_alt_freq < min_af_threshold: continue
            mean_alt_bq = sum(alt_allele_base_qualities)/len(alt_allele_base_qualities) if alt_allele_base_qualities else 0.0
            if mean_alt_bq < min_allele_bq_threshold: continue

        current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
        mean_alt_bq = sum(alt_allele_base_qualities)/len(alt_allele_base_qualities) if alt_allele_base_qualities else 0.0

        variant_key_string = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"
        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)

        # View JSON
        pileup_data_for_view_json = []
        for read_segment_idx, seg_data in enumerate(aligned_read_segments):
            if read_segment_idx >= TENSOR_MAX_READ_ROWS + 50: break
            row_chars = get_read_representation_in_window_for_view(
                seg_data["cigar_ops"], seg_data["offset_on_node"], seg_data["read_sequence"],
                window_start_pos, TENSOR_WINDOW_SIZE, node_len)
            if any(c != ' ' for c in row_chars):
                bases_for_view = [(PADDING_BASE_INDEX if c == ' ' else BASE_TO_INDEX.get(c.upper(), BASE_TO_INDEX['N']))
                                  for c in row_chars]
                pileup_data_for_view_json.append({
                    "bases": bases_for_view,
                    "offset": seg_data["offset_on_node"],
                    "strand": seg_data["strand"],
                    "cigar": seg_data["original_cigar_str"]
                })
        view_oriented_variant_data[variant_key_string] = {
            "pileup_reads_data": pileup_data_for_view_json[:TENSOR_MAX_READ_ROWS],
            "alt_allele_count": alt_allele_count, "ref_allele_count_at_locus": ref_allele_count,
            "other_allele_count_at_locus": other_allele_count, "coverage_at_locus": locus_coverage,
            "alt_allele_frequency": round(current_alt_freq, 4),
            "mean_alt_allele_base_quality": round(mean_alt_bq, 2)
        }

        # Tensor channels
        ch1_list, ch2_list, ch3_list, ch4_list, ch5_list, ch6_list = [], [], [], [], [], []

        # ref row (ch1 top row uses ref indices, ch3 ref row uses 0s)
        ref_base_indices_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i, node_pos_in_window in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
            if 0 <= node_pos_in_window < node_len:
                ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[node_pos_in_window].upper(), BASE_TO_INDEX['N'])

        # ch6 row: AF-scaled if provided; otherwise zeros (this is the key part)
        genomead_af_row = [0] * TENSOR_WINDOW_SIZE
        if genomead_af_list:
            for i, node_pos_in_window in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
                if 0 <= node_pos_in_window < node_len:
                    af_value = genomead_af_list[node_pos_in_window]
                    if af_value == 0.0:
                        scaled = 0
                    else:
                        phred_af = -10.0 * np.log10(af_value)
                        inverted = 127.0 - phred_af
                        scaled = max(1, min(int(inverted), 127))
                    genomead_af_row[i] = scaled

        ch1_list.append(ref_base_indices_row)
        ch2_list.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
        ch3_list.append([MISMATCH_CHANNEL_REF_ROW_VALUE] * TENSOR_WINDOW_SIZE)
        ch4_list.append([DEFAULT_MAPPING_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
        ch5_list.append([CIGAR_PADDING_INDEX] * TENSOR_WINDOW_SIZE)
        ch6_list.append(genomead_af_row)  # zeros when AF not provided

        reads_added = 0
        for seg_data in aligned_read_segments:
            if reads_added >= TENSOR_MAX_READ_ROWS: break
            mapq = max(0, min(int(seg_data["mapping_quality"]), 127))
            base_idx_row, quality_row, mapq_row, cigar_row = get_read_tensor_rows_in_window(
                seg_data["cigar_ops"], seg_data["offset_on_node"], seg_data["read_sequence"],
                seg_data["processed_quality_values"], mapq,
                window_start_pos, TENSOR_WINDOW_SIZE, node_len)

            if any(b != PADDING_BASE_INDEX for b in base_idx_row):
                ch1_list.append(base_idx_row)
                ch2_list.append(quality_row)

                variant_window_index = v_pos - window_start_pos
                mismatch_flags_row = []
                for i in range(TENSOR_WINDOW_SIZE):
                    rb, refb = base_idx_row[i], ref_base_indices_row[i]
                    if rb == PADDING_BASE_INDEX or refb == PADDING_BASE_INDEX:
                        mismatch_flags_row.append(MISMATCH_COMPARISON_PADDING_VALUE)
                    elif rb == refb:
                        mismatch_flags_row.append(0)
                    else:
                        mismatch_flags_row.append(5 if i == variant_window_index else 1)

                ch3_list.append(mismatch_flags_row)
                ch4_list.append(mapq_row)
                ch5_list.append(cigar_row)
                ch6_list.append(genomead_af_row)  # broadcast same AF row; zeros if none
                reads_added += 1

        for _ in range(TENSOR_MAX_READ_ROWS - reads_added):
            ch1_list.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            ch2_list.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            ch3_list.append([MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE)
            ch4_list.append([DEFAULT_MAPPING_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            ch5_list.append([CIGAR_PADDING_INDEX] * TENSOR_WINDOW_SIZE)
            ch6_list.append([0] * TENSOR_WINDOW_SIZE)

        try:
            tensor_chw = torch.tensor([ch1_list, ch2_list, ch3_list, ch4_list, ch5_list, ch6_list], dtype=torch.int8)
            np.save(os.path.join(node_specific_output_dir, f"{variant_key_string}.npy"), tensor_chw.numpy())
            variant_headers_for_summary.append({
                "variant_key": variant_key_string, "tensor_file": f"{variant_key_string}.npy",
                "alt_allele_count": alt_allele_count, "ref_allele_count_at_locus": ref_allele_count,
                "other_allele_count_at_locus": other_allele_count, "coverage_at_locus": locus_coverage,
                "alt_allele_frequency": round(current_alt_freq, 4),
                "mean_alt_allele_base_quality": round(mean_alt_bq, 2)
            })
            tensor_files_generated_for_node += 1
        except Exception as e:
            sys.stderr.write(f"Error creating/saving tensor for {variant_key_string}: {e}\n")

    if variant_headers_for_summary:
        with open(os.path.join(node_specific_output_dir, "variant_summary.json"), 'w') as f:
            json.dump({"node_id": node_id, "node_length": node_len,
                       "variants_passing_af_filter": variant_headers_for_summary}, f, indent=2)

    return node_id, view_oriented_variant_data, tensor_files_generated_for_node

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
        v_pos, v_type = int(variant_key.split('_')[0]), variant_key.split('_')[1]
        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)
        print(f"\n--- Variant: {variant_key} ---")
        ref_chars = [full_node_sequence[j] if 0 <= j < len(full_node_sequence) else '0'
                     for j in range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)]
        print(f"  Node Ref: {''.join(ref_chars)}")
        marker = [' '] * TENSOR_WINDOW_SIZE
        idx = v_pos - window_start_pos
        if 0 <= idx < TENSOR_WINDOW_SIZE: marker[idx] = '^'
        print(f"  Marker  : {''.join(marker)}")
        reads = node_data_for_display_view[variant_key].get("pileup_reads_data", [])
        for j, r in enumerate(reads[:max_reads_to_display_per_variant]):
            bases_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(x, '?') for x in r["bases"]])
            print(f"  Read {j + 1:3d}: {bases_str} (CIGAR:{r['cigar']})")
        meta = node_data_for_display_view[variant_key]
        print(f"  Alt Count: {meta.get('alt_allele_count','N/A')}, Ref Count: {meta.get('ref_allele_count_at_locus','N/A')}, Coverage: {meta.get('coverage_at_locus','N/A')}")
        print(f"  Alt Freq: {meta.get('alt_allele_frequency',0.0):.4f}, Mean Alt BQ: {meta.get('mean_alt_allele_base_quality',0.0):.2f}")

# ─────────────────────────────────────────────────────────────────────────────
# Main
def main():
    parser = argparse.ArgumentParser(description="Generate variant-centered tensors from alignment data.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dat", help=".dat alignment file")
    parser.add_argument("idx", help=".idx index file")
    parser.add_argument("output", help="Base output directory")
    parser.add_argument("candidate_variants_json",
                        help="JSON with nodes (subset) including sequences and optional genomead_af.")
    parser.add_argument("--all_nodes_json", default=None,
                        help="Optional JSON with sequences for ALL nodes (node_id, sequence). Used to process nodes not present in candidate_variants_json; ch6 will be zeros for those.")
    parser.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Number of worker processes")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N',
                        help="Print pileups for top N variants per node (-1 for all)")
    parser.add_argument("--max_view_reads", type=int, default=20, help="Max reads to show per pileup in view mode")
    parser.add_argument("--min_af", type=float, default=0.1, help="Minimum allele frequency for SNPs")
    parser.add_argument("--min_variants", type=int, default=3, help="Alt allele count must be >= this")
    parser.add_argument("--min_allele_bq", type=float, default=10.0, help="Minimum mean base quality for alt alleles")
    parser.add_argument("--variant_type", type=str, default='all', choices=['snp', 'indel', 'all'],
                        help="Variant types to output tensors for.")
    parser.add_argument("--show_empty_info", action="store_true",
                        help="If set, print 'No pileup data' messages for nodes with no passing variants.")

    args = parser.parse_args()

    if not (os.path.isfile(args.dat) and os.path.isfile(args.idx) and os.path.isfile(args.candidate_variants_json)):
        sys.exit("Error: One or more input files (dat, idx, or json) were not found.")
    if args.all_nodes_json and not os.path.isfile(args.all_nodes_json):
        sys.exit(f"Error: --all_nodes_json file not found: {args.all_nodes_json}")
    os.makedirs(args.output, exist_ok=True)

    # Load primary JSON (subset; has sequences and maybe AF)
    node_sequences_primary, node_af_primary = {}, {}
    print(f"Loading nodes from {args.candidate_variants_json}...")
    try:
        with open(args.candidate_variants_json, 'r') as f:
            data = json.load(f)
        for node_obj in data.get('nodes', []):
            node_id_str = node_obj.get('node_id')
            sequence = node_obj.get('sequence')
            af_list = node_obj.get('genomead_af', [])
            if node_id_str and sequence:
                try:
                    nid = int(node_id_str)
                    node_sequences_primary[nid] = sequence.upper()
                    node_af_primary[nid] = af_list
                except ValueError:
                    pass
    except Exception as e:
        sys.exit(f"Error reading or parsing JSON file: {e}")
    print(f"Primary JSON provides sequences for {len(node_sequences_primary)} node IDs.")

    # Optional fallback JSON with sequences for all nodes (no AF)
    node_sequences_fallback = {}
    if args.all_nodes_json:
        print(f"Loading fallback sequences from {args.all_nodes_json}...")
        try:
            with open(args.all_nodes_json, 'r') as f:
                data_f = json.load(f)
            for node_obj in data_f.get('nodes', []):
                node_id_str = node_obj.get('node_id')
                sequence = node_obj.get('sequence')
                if node_id_str and sequence:
                    try:
                        nid = int(node_id_str)
                        if nid not in node_sequences_primary:  # don't overwrite primary
                            node_sequences_fallback[nid] = sequence.upper()
                    except ValueError:
                        pass
            print(f"Fallback JSON provides sequences for {len(node_sequences_fallback)} node IDs.")
        except Exception as e:
            sys.exit(f"Error reading or parsing --all_nodes_json: {e}")

    # Load entire index and iterate ALL nodes
    full_idx_data = load_full_idx_data(args.idx)
    if full_idx_data is None: sys.exit("Failed to load index data.")

    tasks = []
    matched_primary = matched_fallback = skipped_no_seq = 0
    for node_id, (offset, n_records) in full_idx_data.items():
        if node_id in node_sequences_primary:
            seq = node_sequences_primary[node_id]
            af_list = node_af_primary.get(node_id, [])
            matched_primary += 1
        else:
            seq = node_sequences_fallback.get(node_id)
            if seq is None:
                skipped_no_seq += 1
                continue
            af_list = []  # ← channel 6 zeros for non-primary nodes
            matched_fallback += 1
        tasks.append((node_id, offset, n_records, seq, af_list,
                      args.min_af, args.min_variants, args.min_allele_bq, args.variant_type))

    if not tasks:
        sys.exit("No valid tasks to run: could not find sequences for any index nodes.")

    print(f"Matched {matched_primary} nodes with AF (from primary JSON) and {matched_fallback} nodes without AF (ch6=0).")
    if skipped_no_seq:
        print(f"Warning: {skipped_no_seq} index nodes were skipped due to missing sequence in both JSONs.", file=sys.stderr)

    print(f"\nSubmitting {len(tasks)} tasks to {args.num_workers} workers...")
    total_tensors = 0
    milestone_step = 10_000  # change if you want a different cadence
    next_milestone = milestone_step
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=args.num_workers, initializer=init_worker,
                             initargs=(args.dat, args.output)) as executor:
        future_to_node = {executor.submit(process_single_node_for_pileup, task): task[0] for task in tasks}

        for i, future in enumerate(as_completed(future_to_node), start=1):
            node_id = future_to_node[future]
            try:
                _, view_data, tensor_count = future.result()
                total_tensors += tensor_count

                if args.view is not None and view_data:
                    node_seq_for_view = (node_sequences.get(node_id) or
                                         node_sequences_fallback.get(node_id,
                                                                     "")) if 'node_sequences_fallback' in locals() else node_sequences.get(
                        node_id, "")
                    display_pileup_data(
                        view_data, str(node_id), node_seq_for_view,
                        args.max_view_reads,
                        (args.view if args.view != -1 else float('inf')),
                        show_empty_info=args.show_empty_info
                    )
            except Exception as e:
                # break the flashing line cleanly, log error on its own line
                sys.stdout.write("\n")
                print(f"Error processing node {node_id}: {e}", file=sys.stderr)

            elapsed = max(time.time() - start_time, 1e-9)
            rate = i / elapsed
            status = (f"Processed {i:,}/{len(tasks):,} nodes  |  {rate:,.1f} nodes/s  "
                      f"|  tensors: {total_tensors:,}  |  elapsed: {elapsed:,.1f}s")

            # <<< flashing single line
            _print_progress(status)

            # milestone (persistent line), then resume flashing
            if i >= next_milestone:
                sys.stdout.write("\n")  # finalize the flashing line
                print(f"[Milestone] {i:,} nodes processed  ({rate:,.1f} nodes/s).  "
                      f"Total tensors: {total_tensors:,}.  Elapsed: {elapsed:,.1f}s.")
                next_milestone += milestone_step
                # immediately re-render the flashing line (so it keeps updating)
                _print_progress(status)

    # end-of-run: finalize with a newline so the shell prompt isn’t stuck on the status line
    sys.stdout.write("\n")
    print(f"Processing complete. Total tensors generated: {total_tensors:,}.")


if __name__ == '__main__':
    main()
