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
import torch

# ─────────────────────────────────────────────────────────────────────────────
# Constants
RECORD_STRUCT = struct.Struct("<h150s150s20shc")  # Read offset, sequence, RAW QUALITIES, CIGAR, MAPQ, strand
RECORD_SIZE = RECORD_STRUCT.size

BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 5, ' ': 6, '-': 6}
INDEX_TO_BASE_FOR_VIEW = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N', 5: '*', 6: ' '}

TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200
PADDING_BASE_INDEX = BASE_TO_INDEX[' ']
DEFAULT_QUALITY_PADDING = 0
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = 0

# Globals for worker process state
worker_dat_file = None
worker_base_output_dir = None


# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions
# (reverse_complement, parse_idx_file_for_single_node, load_node_sequence_from_gfa,
#  decode_cigar_to_int_ops, get_allele_from_read_at_node_pos, detect_variants_from_cigar,
#  get_read_representation_in_window_for_view, get_read_tensor_rows_in_window remain the same)
# ─────────────────────────────────────────────────────────────────────────────

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]


def parse_idx_file_for_single_node(idx_path, target_node_id):
    node_info = None
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4:
                print(f"❌ Error: Index file {idx_path} is too small.", file=sys.stderr)
                return None
            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4:
                print(f"❌ Error: Could not read number of nodes from {idx_path}.", file=sys.stderr)
                return None
            num_nodes = struct.unpack('<I', num_nodes_bytes)[0]
            print(f"🔹 Index file contains {num_nodes} nodes. Searching for node {target_node_id}...")
            found = False
            for i in range(num_nodes):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                    print(f"❌ Error: Index file ended prematurely while reading record {i + 1}/{num_nodes}.",
                          file=sys.stderr)
                    break
                node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                if node_id_from_idx == target_node_id:
                    node_info = (offset, n_records)
                    found = True
                    print(f"✔ Found node {target_node_id} in index: Offset={offset}, N_Records={n_records}")
                    break
            if not found:
                print(f"❌ Error: Node ID {target_node_id} not found in the index file {idx_path}.", file=sys.stderr)
                return None
        return node_info
    except FileNotFoundError:
        print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr)
        # Changed from sys.exit(1) to allow continuation for multiple nodes
        return None
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
        # Changed from sys.exit(1)
        return None


def load_node_sequence_from_gfa(gfa_path, target_node_id):
    node_sequence = None
    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file to find sequence for node {target_node_id}: {gfa_path}")
            line_counter = 0
            for line in f:
                line_counter += 1
                if line_counter % 5_000_000 == 0:
                    print(f"  Checked {line_counter:,} lines in GFA file...")
                if not line.startswith('S\t'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                try:
                    nid = int(parts[1])
                except ValueError:
                    continue
                if nid == target_node_id:
                    node_sequence = parts[2]
                    print(f"✔ Found sequence for node {target_node_id} in GFA.")
                    break
            if line_counter > 0: print(f"✔ Finished GFA scan after {line_counter:,} lines.")
            if node_sequence is None:
                print(f"❌ Error: Sequence for node ID {target_node_id} not found in GFA file {gfa_path}.",
                      file=sys.stderr)
    except FileNotFoundError:
        print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr)
        # Changed from sys.exit(1)
        return None
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr)
        # Changed from sys.exit(1)
        return None
    return node_sequence


def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*':
        return []
    ops = []
    try:
        for length_str, op_char in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string):
            ops.append((int(length_str), op_char))
        return ops
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr)
        return []


def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence, read_cigar_ops_decoded,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None, expected_ref_allele_for_indel=None):
    current_node_pos = read_offset_on_node
    current_read_pos = 0
    for length, op in read_cigar_ops_decoded:
        if op == 'M' or op == '=' or op == 'X':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I': return "REF_STATE_FOR_INDEL"
                if expected_var_type == 'D': return "REF_STATE_FOR_INDEL"
                offset_in_block = target_node_pos - current_node_pos
                if current_read_pos + offset_in_block < len(read_sequence):
                    return read_sequence[current_read_pos + offset_in_block].upper()
                return None
            current_node_pos += length
            current_read_pos += length
        elif op == 'I':
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                return read_sequence[current_read_pos: current_read_pos + length].upper()
            current_read_pos += length
        elif op == 'D':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I': return "OTHER_FOR_INDEL"
                if expected_var_type == 'D':
                    deleted_seq_in_read = node_sequence[current_node_pos: current_node_pos + length]
                    if deleted_seq_in_read == expected_ref_allele_for_indel:
                        return "*"
                    else:
                        return "OTHER_FOR_INDEL"
                return "*"
            current_node_pos += length
        elif op == 'S':
            current_read_pos += length
        elif op == 'N':
            current_node_pos += length
        if current_node_pos > target_node_pos + 1 and op in ('M', '=', 'X', 'D', 'N'):
            if not (expected_var_type == 'I' and (current_node_pos - 1) <= target_node_pos):
                break
    return None


def detect_variants_from_cigar(offset_on_node, cigar_ops_decoded, read_sequence, node_sequence):
    variants = []
    node_pos = offset_on_node
    read_pos = 0
    for length, op in cigar_ops_decoded:
        if op == 'M' or op == '=' or op == 'X':
            for i in range(length):
                current_node_pos = node_pos + i
                current_read_pos = read_pos + i
                if current_node_pos < len(node_sequence) and current_read_pos < len(read_sequence):
                    node_base = node_sequence[current_node_pos].upper()
                    read_base = read_sequence[current_read_pos].upper()
                    if node_base != read_base and op != '=':  # Ensure mismatch op (X) or actual base diff for M
                        variants.append((current_node_pos, 'X', read_base, node_base))
            node_pos += length
            read_pos += length
        elif op == 'I':
            inserted_sequence = read_sequence[read_pos: read_pos + length].upper()
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0
            ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if 0 <= ref_anchor_pos < len(
                node_sequence) else ""
            variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor if ref_base_at_anchor else "*"))
            read_pos += length
        elif op == 'D':
            deleted_sequence_from_ref = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= len(
                node_sequence) else ""
            if deleted_sequence_from_ref:
                variants.append((node_pos, 'D', '*', deleted_sequence_from_ref))
            node_pos += length
        elif op == 'S':
            read_pos += length
        elif op == 'N':
            node_pos += length
    return variants


def get_read_representation_in_window_for_view(segment_cigar_ops, segment_offset_on_node, segment_read_sequence,
                                               window_start_node, window_size, node_len):
    window_char_representation = [' '] * window_size
    current_node_pos_in_read = segment_offset_on_node
    current_read_pos_in_read = 0
    for cigar_len, cigar_op in segment_cigar_ops:
        if cigar_op in ('M', '=', 'X'):
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                read_aln_pos = current_read_pos_in_read + i
                if window_start_node <= node_aln_pos < window_start_node + window_size:
                    window_idx = node_aln_pos - window_start_node
                    if read_aln_pos < len(segment_read_sequence):
                        window_char_representation[window_idx] = segment_read_sequence[read_aln_pos].upper()
            current_node_pos_in_read += cigar_len
            current_read_pos_in_read += cigar_len
        elif cigar_op == 'D' or cigar_op == 'N':
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                if window_start_node <= node_aln_pos < window_start_node + window_size:
                    window_idx = node_aln_pos - window_start_node
                    window_char_representation[window_idx] = '*'
            current_node_pos_in_read += cigar_len
        elif cigar_op == 'I' or cigar_op == 'S':
            current_read_pos_in_read += cigar_len
        if current_node_pos_in_read >= window_start_node + window_size and cigar_op in (
        'M', '=', 'X', 'D', 'N'):  # Check if we've passed the window
            break
    return window_char_representation


def get_read_tensor_rows_in_window(segment_cigar_ops, segment_offset_on_node,
                                   segment_read_sequence, segment_quality_str,
                                   window_start_node, current_tensor_window_size, node_len):
    base_indices_row = [PADDING_BASE_INDEX] * current_tensor_window_size
    quality_scores_row = [DEFAULT_QUALITY_PADDING] * current_tensor_window_size
    current_node_pos_in_read = segment_offset_on_node
    current_read_pos_in_read = 0
    for cigar_len, cigar_op in segment_cigar_ops:
        if current_node_pos_in_read >= window_start_node + current_tensor_window_size and cigar_op in (
                'M', 'D', 'N', '=', 'X'):
            break
        # Removed check for current_read_pos_in_read >= len(segment_read_sequence) here,
        # as it's handled per-base inside M/=/X and per-op for I/S
        if cigar_op in ('M', '=', 'X'):
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                read_aln_pos = current_read_pos_in_read + i
                if read_aln_pos >= len(segment_read_sequence): break  # Ensure we don't go past read sequence
                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_char = segment_read_sequence[read_aln_pos].upper()
                    base_indices_row[window_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])
                    if read_aln_pos < len(segment_quality_str):  # Ensure we don't go past quality string
                        try:
                            quality_scores_row[window_idx] = ord(segment_quality_str[read_aln_pos]) - 33
                        except (TypeError, IndexError):  # TypeError for None, IndexError if qual_str too short
                            quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
                    else:
                        quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
            current_node_pos_in_read += cigar_len
            current_read_pos_in_read += cigar_len
        elif cigar_op == 'D' or cigar_op == 'N':
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_indices_row[window_idx] = BASE_TO_INDEX['*']
                    quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING  # Deletions have no quality in read
            current_node_pos_in_read += cigar_len
        elif cigar_op == 'I' or cigar_op == 'S':
            # Check if current_read_pos_in_read will exceed sequence/quality length
            if current_read_pos_in_read + cigar_len > len(segment_read_sequence):
                # This indicates a CIGAR/read length mismatch, a soft clip or insertion
                # might be described longer than the available sequence.
                # We should probably log this or handle more gracefully if it's an issue.
                # For now, just advance by the CIGAR length. Tensor generation relies on aligned positions.
                pass  # current_read_pos_in_read will be advanced below.
            current_read_pos_in_read += cigar_len

        # Break if read is exhausted, relevant after advancing current_read_pos_in_read
        if current_read_pos_in_read >= len(segment_read_sequence) and cigar_op in ('M', '=', 'X', 'I', 'S'):
            break
    return base_indices_row, quality_scores_row


## ─────────────────────────────────────────────────────────────────────────────
## Worker Process Initialization and Target Function
## ─────────────────────────────────────────────────────────────────────────────

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker):
    global worker_dat_file, worker_base_output_dir
    # This function is called once per worker process when the pool starts.
    try:
        # It's crucial that each worker has its own file handle if multiple workers are
        # operating on the same file, even for read-only. Opening here ensures this.
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
        # print(f"INFO [Worker {os.getpid()}]: Initialized. DAT: {dat_file_path_for_worker}, OutDir: {base_output_dir_for_worker}")
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path_for_worker}\n")
        # Exiting worker process if essential resource like DAT file is missing.
        # The ProcessPoolExecutor will raise an exception in the main process if a worker dies.
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path_for_worker}: {e}\n")
        sys.exit(1)


def process_single_node_for_pileup(task_args_with_af_thresh):
    node_id, dat_file_offset, n_records, node_sequence, min_af_threshold = task_args_with_af_thresh
    global worker_dat_file, worker_base_output_dir

    if worker_dat_file is None or worker_base_output_dir is None:
        # This case should ideally not be reached if init_worker is successful.
        # If it is, it indicates a problem with worker initialization or state.
        sys.stderr.write(
            f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Worker not initialized properly (dat_file or output_dir missing).\n")
        return node_id, None  # Return None for view_data to indicate a critical worker problem

    if not node_sequence:
        sys.stderr.write(
            f"ℹ️ [Worker {os.getpid()} for Node {node_id}]: No sequence provided. Skipping tensor/summary generation.\n")
        return node_id, {}  # Return empty dict for view_data, not a critical error

    node_specific_output_dir = os.path.join(worker_base_output_dir, str(node_id))
    try:
        os.makedirs(node_specific_output_dir, exist_ok=True)
    except OSError as e:
        sys.stderr.write(
            f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Could not create directory {node_specific_output_dir}: {e}\n")
        return node_id, None  # Cannot proceed without output directory

    node_len = len(node_sequence)
    variant_headers_for_node = []
    view_oriented_variant_data = {}  # Initialize as empty dict

    aligned_read_segments = []
    try:
        worker_dat_file.seek(dat_file_offset + 10)  # Assuming the +10 is a known constant skip
        for record_idx in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE:
                # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()} for Node {node_id}]: Expected {RECORD_SIZE} bytes, got {len(data)}. End of records for node?\n")
                break
            off_from_file, raw_seq, raw_qual, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)

            if mapq < 10: continue  # Filter by MAPQ

            try:
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')
                qual_str = raw_qual.rstrip(b'\x00').decode('ascii', errors='replace')
                cigar_str_original = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()} for Node {node_id}]: Unicode decode error for a read. Skipping.\n")
                continue

            if len(seq) == 0 or len(seq) != len(qual_str):
                # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()} for Node {node_id}]: Mismatch length seq vs qual or empty seq. Skipping read.\n")
                continue

            original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
            if not original_decoded_cigar_ops and cigar_str_original != '*':  # Skip if CIGAR is bad but not explicitly '*'
                # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()} for Node {node_id}]: Bad CIGAR '{cigar_str_original}'. Skipping read.\n")
                continue

            current_read_sequence = seq
            current_quality_str = qual_str
            current_decoded_cigar_ops = list(original_decoded_cigar_ops)  # Make a mutable copy
            current_offset_on_node = off_from_file

            if strand_char == '-':
                current_read_sequence = reverse_complement(seq)
                current_quality_str = qual_str[::-1]

                # Reverse CIGAR operations for negative strand processing
                temp_rev_cigar_ops = []
                for length, op_code in original_decoded_cigar_ops:  # Iterate original CIGAR
                    temp_rev_cigar_ops.insert(0, (length, op_code))  # Prepend to reverse order
                current_decoded_cigar_ops = temp_rev_cigar_ops

                # Recalculate offset for negative strand based on alignment span on reference from CIGAR
                ref_span_from_cigar = sum(l for l, op in current_decoded_cigar_ops if op in ('M', 'D', 'N', '=', 'X'))
                current_offset_on_node = node_len - ref_span_from_cigar - off_from_file
                if current_offset_on_node < 0:
                    # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()} Node {node_id}]: Negative offset for strand '-', off {off_from_file}, ref_span {ref_span_from_cigar}, node_len {node_len}. Read: {seq[:20]}... CIGAR: {cigar_str_original}. Skipping.\n")
                    continue

            aligned_read_segments.append({
                "offset_on_node": current_offset_on_node,
                "read_sequence": current_read_sequence,
                "processed_quality_str": current_quality_str,
                "cigar_ops": current_decoded_cigar_ops,
                "original_cigar_str": cigar_str_original,  # For display/debug
                "strand": strand_char
            })
    except Exception as e:
        # This error is critical for the processing of this node
        sys.stderr.write(f"❌ Error [Worker {os.getpid()} for Node {node_id}] reading records: {e}\n")
        return node_id, None  # Indicate failure for this node

    candidate_variants = defaultdict(int)
    for segment in aligned_read_segments:
        variants_in_read = detect_variants_from_cigar(
            segment["offset_on_node"], segment["cigar_ops"],
            segment["read_sequence"], node_sequence
        )
        for v_pos, v_type, v_alt, v_ref in variants_in_read:
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    variant_processing_window_size = TENSOR_WINDOW_SIZE  # Same window for view and tensor data gathering
    half_window = variant_processing_window_size // 2

    for (v_pos, v_type, v_ref_defined, v_alt_defined), read_support_count in candidate_variants.items():
        # AF calculation
        af_alt_count, af_ref_count, af_other_count, af_locus_coverage = 0, 0, 0, 0
        # Determine the reference base at the variant position for AF calculation, esp. for indels
        ref_allele_for_indel_af_check = None
        if v_type == 'D':  # Deletion: ref is the deleted sequence from node_sequence
            ref_allele_for_indel_af_check = v_ref_defined  # This is node_sequence[v_pos : v_pos+len(del)]
        elif v_type == 'I':  # Insertion: ref is the base on node_sequence *before* insertion
            if 0 <= v_pos < node_len:  # v_pos is the anchor base on ref
                ref_allele_for_indel_af_check = node_sequence[v_pos]
            # else: ref_allele_for_indel_af_check remains None if v_pos is out of bounds

        for segment in aligned_read_segments:
            allele = get_allele_from_read_at_node_pos(
                segment["offset_on_node"], segment["read_sequence"], segment["cigar_ops"],
                v_pos, node_sequence,
                expected_var_type=v_type, expected_ref_allele_for_indel=ref_allele_for_indel_af_check
            )
            if allele is not None:
                af_locus_coverage += 1
                if v_type == 'X':  # SNP/Mismatch
                    if allele == v_alt_defined:
                        af_alt_count += 1
                    elif allele == v_ref_defined:
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'I':  # Insertion
                    if allele == v_alt_defined:
                        af_alt_count += 1
                    elif allele == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1  # Read matches ref at site of insertion
                    else:
                        af_other_count += 1  # Other insertion or SNP at this site
                elif v_type == 'D':  # Deletion
                    if allele == "*":
                        af_alt_count += 1  # Read shows deletion
                    elif allele == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1  # Read does not show deletion
                    else:
                        af_other_count += 1  # Other variant at this site (e.g. SNP instead of Del)

        alt_freq = af_alt_count / af_locus_coverage if af_locus_coverage > 0 else 0.0
        if alt_freq < min_af_threshold:
            continue  # Skip variant if AF is below threshold

        variant_key_str = f"{v_pos}_{v_type}_{v_ref_defined}_{v_alt_defined}"

        # Determine window start for both view and tensor based on variant type and position
        if v_type == 'I':  # For insertions, center window *after* the anchor base
            window_center_on_node = v_pos + 1
        else:  # For SNPs and Deletions, center window *at* the variant position
            window_center_on_node = v_pos
        current_window_start_on_node = window_center_on_node - half_window

        # --- Data for Pileup View ---
        pileup_reads_data_for_view = []
        for segment in aligned_read_segments:  # Limit reads for viewing if necessary later
            row_chars_for_view = get_read_representation_in_window_for_view(
                segment["cigar_ops"], segment["offset_on_node"], segment["read_sequence"],
                current_window_start_on_node, variant_processing_window_size, node_len
            )
            if any(c != ' ' for c in row_chars_for_view):  # Only add if read has content in window
                row_indices_for_view = [BASE_TO_INDEX.get(char.upper(), BASE_TO_INDEX['N']) for char in
                                        row_chars_for_view]
                pileup_reads_data_for_view.append({
                    "bases": row_indices_for_view, "offset": segment["offset_on_node"],
                    "strand": segment["strand"], "cigar": segment["original_cigar_str"]  # Store original CIGAR
                })

        view_oriented_variant_data[variant_key_str] = {
            "pileup_reads_data": pileup_reads_data_for_view[:TENSOR_MAX_READ_ROWS],  # Limit reads in json too
            "alt_allele_count": af_alt_count, "ref_allele_count_at_locus": af_ref_count,
            "other_allele_count_at_locus": af_other_count, "coverage_at_locus": af_locus_coverage,
            "alt_allele_frequency": round(alt_freq, 4)
        }

        # --- Data for Tensor ---
        tensor_ch1_bases_list = []
        tensor_ch2_qualities_list = []
        tensor_ch3_mismatches_list = []

        # Reference Row for Tensor
        ref_base_indices_row_tensor = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        ref_qual_scores_row_tensor = [DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE  # CH2 Ref Quals (can be all 0)
        ref_mismatch_row_tensor = [MISMATCH_CHANNEL_REF_ROW_VALUE] * TENSOR_WINDOW_SIZE  # CH3 Ref Mismatches (all 0)
        for i in range(TENSOR_WINDOW_SIZE):
            actual_node_pos = current_window_start_on_node + i  # Use same window_start as view
            if 0 <= actual_node_pos < node_len:
                ref_base_indices_row_tensor[i] = BASE_TO_INDEX.get(node_sequence[actual_node_pos].upper(),
                                                                   BASE_TO_INDEX['N'])
        tensor_ch1_bases_list.append(ref_base_indices_row_tensor)
        tensor_ch2_qualities_list.append(ref_qual_scores_row_tensor)
        tensor_ch3_mismatches_list.append(ref_mismatch_row_tensor)

        reads_added_to_tensor = 0
        for segment in aligned_read_segments:  # Iterate through all relevant reads
            if reads_added_to_tensor >= TENSOR_MAX_READ_ROWS:
                break

            base_indices_row_tensor, quality_scores_row_tensor = get_read_tensor_rows_in_window(
                segment["cigar_ops"], segment["offset_on_node"],
                segment["read_sequence"], segment["processed_quality_str"],
                current_window_start_on_node, TENSOR_WINDOW_SIZE, node_len  # Use same window_start
            )

            # Only add read to tensor if it has some coverage in the window
            if any(b != PADDING_BASE_INDEX for b in base_indices_row_tensor):
                tensor_ch1_bases_list.append(base_indices_row_tensor)
                tensor_ch2_qualities_list.append(quality_scores_row_tensor)

                # Mismatch Channel for this read
                mismatch_row_for_read_tensor = [MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE
                for i in range(TENSOR_WINDOW_SIZE):
                    read_base_idx = base_indices_row_tensor[i]
                    ref_base_idx_for_comp = ref_base_indices_row_tensor[i]  # From the ref row calculated above

                    if read_base_idx == PADDING_BASE_INDEX or ref_base_idx_for_comp == PADDING_BASE_INDEX:
                        mismatch_row_for_read_tensor[
                            i] = MISMATCH_COMPARISON_PADDING_VALUE  # Padding vs padding/anything
                    elif read_base_idx == ref_base_idx_for_comp:
                        mismatch_row_for_read_tensor[i] = 0  # Match
                    else:
                        mismatch_row_for_read_tensor[i] = 1  # Mismatch
                tensor_ch3_mismatches_list.append(mismatch_row_for_read_tensor)
                reads_added_to_tensor += 1

        # Pad tensor if fewer than TENSOR_MAX_READ_ROWS reads were added
        for _ in range(
                TENSOR_MAX_READ_ROWS - reads_added_to_tensor):  # Number of reads to pad = Max - (actual reads + ref_row)
            tensor_ch1_bases_list.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            tensor_ch2_qualities_list.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            tensor_ch3_mismatches_list.append([MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE)

        try:
            # Ensure all sublists for each channel have the same length (TENSOR_MAX_READ_ROWS + 1 for ref)
            # And each inner list has TENSOR_WINDOW_SIZE
            # This should be guaranteed by construction.
            final_tensor = torch.tensor([
                tensor_ch1_bases_list,
                tensor_ch2_qualities_list,
                tensor_ch3_mismatches_list
            ], dtype=torch.int8)  # Or float if qualities/mismatches are to be floats

            tensor_filename = f"{variant_key_str}.pth"
            tensor_filepath = os.path.join(node_specific_output_dir, tensor_filename)
            torch.save(final_tensor, tensor_filepath)

            variant_headers_for_node.append({
                "variant_key": variant_key_str,
                "tensor_file": tensor_filename,
                "alt_allele_count": af_alt_count,
                "ref_allele_count_at_locus": af_ref_count,
                "other_allele_count_at_locus": af_other_count,
                "coverage_at_locus": af_locus_coverage,
                "alt_allele_frequency": round(alt_freq, 4)
            })
        except Exception as e:
            sys.stderr.write(
                f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Failed to create or save tensor for {variant_key_str}: {e}\n")
            # Continue to next variant if one tensor fails

    # After processing all variants for the node
    if variant_headers_for_node:  # Only write summary if some variants were processed and tensorized
        summary_json_path = os.path.join(node_specific_output_dir, "variant_summary.json")
        try:
            with open(summary_json_path, 'w') as sjf:
                json.dump({
                    "node_id": node_id,
                    "node_length": node_len,
                    "node_sequence_preview": node_sequence[:100] + "..." if node_len > 100 else node_sequence,
                    # Preview
                    "variants": variant_headers_for_node  # List of dicts
                }, sjf, indent=2)
        except Exception as e:
            sys.stderr.write(f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Failed to write summary JSON: {e}\n")
            # view_oriented_variant_data is still returned

    return node_id, view_oriented_variant_data


## ─────────────────────────────────────────────────────────────────────────────
## Pileup Viewing Function (display_pileup_data)
## ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if not node_data_for_display_view or not isinstance(node_data_for_display_view, dict):
        print(f"ℹ️ No valid pileup data to display for node {node_id_str_for_display}.", file=sys.stderr)
        return

    print(
        f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")

    if not node_data_for_display_view:  # Checks if the dict is empty
        print(f"ℹ️ No variants found or pileups generated for this node (or all filtered by AF).")
        return

    variants_displayed_count = 0
    # Sort variants by position then type (e.g., 10_X, 10_I, 12_X)
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))

    display_window_size = TENSOR_WINDOW_SIZE  # Use consistent window size
    half_display_window = display_window_size // 2

    for variant_key in sorted_variant_keys:
        if variants_displayed_count >= max_variants_to_display:
            print(
                f"\n  ... (and {len(node_data_for_display_view) - variants_displayed_count} more variants not shown due to limit)")
            break

        variant_data = node_data_for_display_view[variant_key]
        pileup_reads_display_data = variant_data.get("pileup_reads_data", [])

        # Parse variant key
        parts = variant_key.split('_')
        v_pos = int(parts[0])
        v_type = parts[1]
        # v_ref = parts[2]
        # v_alt = parts[3]

        # Determine window center for display
        if v_type == 'I':
            window_center_on_node_display = v_pos + 1  # Center after anchor base for insertions
        else:
            window_center_on_node_display = v_pos  # Center on the variant for SNP/Del
        current_window_start_on_node_display = window_center_on_node_display - half_display_window

        print(
            f"\n--- Variant: {variant_key} (Node Pos: {v_pos}, Display Window on Node: {current_window_start_on_node_display}-{current_window_start_on_node_display + display_window_size - 1}) ---")
        print(f"  Alt Count: {variant_data.get('alt_allele_count', 'N/A')}")
        print(f"  Ref Count: {variant_data.get('ref_allele_count_at_locus', 'N/A')}")
        print(f"  Other Count: {variant_data.get('other_allele_count_at_locus', 'N/A')}")
        print(f"  Coverage: {variant_data.get('coverage_at_locus', 'N/A')}")
        alt_freq_val = variant_data.get('alt_allele_frequency', 'N/A')
        if isinstance(alt_freq_val, float):
            print(f"  Alt Freq: {alt_freq_val:.4f}")
        else:
            print(f"  Alt Freq: {alt_freq_val}")

        # Display reference sequence and marker
        ref_display_parts = [' '] * display_window_size
        marker_line_parts = [' '] * display_window_size
        variant_display_idx_in_window = -1

        # Determine the variant's visual position in the window
        if v_type != 'I':  # SNP or Deletion
            if current_window_start_on_node_display <= v_pos < current_window_start_on_node_display + display_window_size:
                variant_display_idx_in_window = v_pos - current_window_start_on_node_display
        elif v_type == 'I':  # Insertion (marker between anchor and inserted base)
            # Insertion is at v_pos, meaning it's *after* node_sequence[v_pos].
            # The marker points between v_pos and v_pos+1 conceptually.
            # If v_pos is the last base in window, marker might be tricky.
            # Let's mark at v_pos to indicate the anchor.
            if current_window_start_on_node_display <= v_pos < current_window_start_on_node_display + display_window_size:
                variant_display_idx_in_window = v_pos - current_window_start_on_node_display

        for i in range(display_window_size):
            actual_node_pos_in_window = current_window_start_on_node_display + i
            if 0 <= actual_node_pos_in_window < len(full_node_sequence):
                ref_display_parts[i] = full_node_sequence[actual_node_pos_in_window]

            if i == variant_display_idx_in_window:
                if v_type == 'I':
                    marker_line_parts[i] = "I"  # Mark the anchor base
                    if i + 1 < display_window_size:
                        marker_line_parts[i + 1] = "^"  # Point between anchor and next base
                    elif i == display_window_size - 1:
                        marker_line_parts[i] = ">"  # if at the end
                else:  # SNP/Deletion
                    marker_line_parts[i] = "^"  # Mark the affected base(s)

        print(f"  Node Ref: {''.join(ref_display_parts)}")
        print(f"  Marker  : {''.join(marker_line_parts)}")

        if not pileup_reads_display_data:
            print("  (No reads data in pileup for this variant's window or to display)")
        else:
            displayed_reads_count = 0
            for i, read_info in enumerate(
                    pileup_reads_display_data):  # pileup_reads_display_data already limited by TENSOR_MAX_READ_ROWS
                if displayed_reads_count >= max_reads_to_display_per_variant:
                    print(
                        f"  ... (and {len(pileup_reads_display_data) - displayed_reads_count} more reads not shown due to view limit)")
                    break
                base_indices = read_info["bases"]
                read_offset = read_info["offset"]
                read_strand = read_info["strand"]
                read_cigar = read_info.get("cigar", "N/A")  # Get CIGAR if available
                pileup_row_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in base_indices])
                print(
                    f"  Read {i + 1:3d}: {pileup_row_str}  (Offset: {read_offset}, Strand: {read_strand}, CIGAR: {read_cigar})")
                displayed_reads_count += 1
        variants_displayed_count += 1
    print("\n")


## ─────────────────────────────────────────────────────────────────────────────
## Main function
## ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered pileups for specified node(s) and optionally view them.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("dat", help=".dat file path (read alignment data)")
    parser.add_argument("idx", help=".idx file path (index for .dat file)")
    parser.add_argument("output", help="Base output directory for node-specific folders, tensors, and summaries.")

    # Mutually exclusive group for node_id or node_id_file
    input_node_group = parser.add_mutually_exclusive_group(required=True)
    input_node_group.add_argument("--node_id", type=int, help="The specific node ID to process.")
    input_node_group.add_argument("--node_id_file",
                                  help="Path to a text file containing node IDs (one per line) to process.")

    parser.add_argument("--gfa",
                        help="GFA graph file path (required if node sequence cache is not used/built for any specified node).")
    parser.add_argument("--load-cache", help="Load node sequences from this JSON cache file.")
    parser.add_argument("--save-cache",
                        help="Save/update node sequences to this JSON cache file (used if --gfa is provided to fetch new sequences).")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N',
                        help="Print generated pileups to console. Optionally specify N for first N variants per node. -1 for all.")
    parser.add_argument("--max_view_reads", type=int, default=20,
                        help="Max reads to display per pileup in console view.")
    parser.add_argument("--min_af", type=float, default=0.1,
                        help="Min allele frequency for a variant to be processed and included in output.")
    args = parser.parse_args()

    # --- Argument Validation ---
    if not os.path.isfile(args.dat): sys.exit(f"❌ Error: DAT file not found: {args.dat}")
    if not os.path.isfile(args.idx): sys.exit(f"❌ Error: Index file not found: {args.idx}")
    if not args.load_cache and not args.gfa:
        # If neither cache nor GFA is provided, sequences cannot be obtained.
        sys.exit("❌ Error: Must provide --gfa or --load-cache to obtain node sequences.")
    if args.load_cache and not os.path.isfile(args.load_cache):
        # If --load-cache is specified but file doesn't exist, it's an issue unless --gfa can cover all.
        # For simplicity, error out if specified load-cache is missing. User can omit --load-cache if it's optional.
        sys.exit(f"❌ Error: Specified --load-cache file not found: {args.load_cache}")
    if args.gfa and not os.path.isfile(args.gfa):
        sys.exit(f"❌ Error: GFA file not found: {args.gfa}")
    if not (0.0 <= args.min_af <= 1.0):  # Pythonic check for range
        sys.exit("❌ Error: --min_af must be between 0.0 and 1.0.")

    try:
        os.makedirs(args.output, exist_ok=True)
        print(f"🔹 Base output directory: {args.output}")
    except OSError as e:
        sys.exit(f"❌ Error: Could not create base output directory {args.output}: {e}")

    # --- Determine Target Node IDs ---
    target_node_ids = []
    if args.node_id_file:
        try:
            with open(args.node_id_file, 'r') as f_nodes:
                for line_num, line in enumerate(f_nodes, 1):
                    line = line.strip()
                    if line and not line.startswith("#"):  # Ignore empty lines and comments
                        try:
                            target_node_ids.append(int(line))
                        except ValueError:
                            print(
                                f"⚠️ Warning: Invalid non-integer node ID '{line}' in {args.node_id_file} at line {line_num}. Skipping.",
                                file=sys.stderr)
            if not target_node_ids:
                sys.exit(f"❌ Error: No valid node IDs found or specified in {args.node_id_file}.")
            print(f"🔹 Will process {len(target_node_ids)} node ID(s) from file: {args.node_id_file}")
        except FileNotFoundError:
            sys.exit(f"❌ Error: Node ID file not found: {args.node_id_file}")
    elif args.node_id is not None:
        target_node_ids.append(args.node_id)
        print(f"🔹 Will process single target node ID: {args.node_id}")

    # --- Initialize In-Memory Sequence Cache ---
    node_sequences_map = {}  # Stores node_id_str -> sequence
    if args.load_cache:  # This arg implies the file *should* exist (checked above)
        start_time_cache_load = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                node_sequences_map = json.load(cf)  # Loads as str keys
            print(
                f"✔ Loaded {len(node_sequences_map)} sequences from cache '{args.load_cache}' in {time.time() - start_time_cache_load:.2f}s.")
        except json.JSONDecodeError as e:
            print(
                f"⚠️ Warning: Error decoding JSON from cache file {args.load_cache}: {e}. Starting with an empty sequence map for loaded cache.",
                file=sys.stderr)
        except Exception as e:  # Other file reading errors
            print(
                f"⚠️ Warning: Error loading cache {args.load_cache}: {e}. Starting with an empty sequence map for loaded cache.",
                file=sys.stderr)

    # --- Main Loop for Processing Each Node ---
    processed_nodes_count = 0
    successful_nodes_with_output = 0
    overall_start_time = time.time()

    for target_node_id in target_node_ids:
        print(f"\n═════════════════════════════ PROCESSING NODE ID: {target_node_id} ═════════════════════════════")
        node_specific_start_time = time.time()
        node_sequence = None  # Reset for current node

        # 1. Get Node Data Offset from Index
        idx_parse_start_time = time.time()
        node_dat_info = parse_idx_file_for_single_node(args.idx, target_node_id)
        if not node_dat_info:
            print(f"❌ Skipping node {target_node_id} due to index parsing failure.", file=sys.stderr)
            continue  # Skip to the next node_id
        dat_offset, n_records = node_dat_info
        print(f"✔ Index parsing for node {target_node_id} took {time.time() - idx_parse_start_time:.2f}s.")

        # 2. Obtain Node Sequence (Cache -> GFA)
        node_id_str = str(target_node_id)  # Use string for dict keys
        if node_id_str in node_sequences_map:
            node_sequence = node_sequences_map[node_id_str]
            print(f"✔ Sequence for node {target_node_id} found in pre-loaded/memory cache.")
        elif args.gfa:
            print(f"ℹ️ Sequence for node {target_node_id} not in cache, attempting to load from GFA: {args.gfa}")
            gfa_load_start_time = time.time()
            node_sequence = load_node_sequence_from_gfa(args.gfa, target_node_id)
            if node_sequence:
                print(
                    f"✔ GFA sequence loading for node {target_node_id} took {time.time() - gfa_load_start_time:.2f}s.")
                node_sequences_map[node_id_str] = node_sequence  # Add to in-memory map for potential save later
            else:
                print(f"❌ Skipping node {target_node_id}: Failed to load sequence from GFA.", file=sys.stderr)
                continue
        else:  # Not in cache and no GFA to fall back on
            print(f"❌ Skipping node {target_node_id}: Sequence not in cache and --gfa not provided.", file=sys.stderr)
            continue

        if not node_sequence:  # Should be caught above, but as a safeguard
            print(f"❌ Skipping node {target_node_id}: Sequence could not be obtained.", file=sys.stderr)
            continue

        # 3. Prepare and Execute Task
        task_payload = (target_node_id, dat_offset, n_records, node_sequence, args.min_af)
        print(
            f"🔹 Preparing task for node {target_node_id} (len: {len(node_sequence)}bp, {n_records} records, min_AF: {args.min_af}).")

        node_processing_start_time = time.time()
        processed_node_id_result = None
        view_data_for_node = None  # This will hold dict from worker, or None on severe worker error

        try:
            # Using a ProcessPoolExecutor for each node to keep changes minimal from original structure.
            # For many small nodes, a single persistent pool might be more performant.
            with ProcessPoolExecutor(max_workers=1, initializer=init_worker,
                                     initargs=(args.dat, args.output)) as executor:
                future = executor.submit(process_single_node_for_pileup, task_payload)
                processed_node_id_result, view_data_for_node = future.result()  # Wait for result

            if processed_node_id_result is None:  # Implies critical error in worker before returning node_id
                print(
                    f"❌ Critical Error: Worker failed to return results for node {target_node_id}. Might be an issue in init_worker or early in process_single_node_for_pileup.",
                    file=sys.stderr)
                continue  # Skip this node

            # view_data_for_node being an empty dict {} is normal if no variants passed AF
            # view_data_for_node being None means a problem inside the worker task (e.g. couldn't create output dir)
            if view_data_for_node is None:
                print(
                    f"⚠️ Warning: Processing for node {target_node_id} returned no view data, suggesting an error within the worker task. Check logs.",
                    file=sys.stderr)
                # Still, some partial output like .pth might exist if error was late.

        except Exception as e:
            print(f"❌ An error occurred launching or during processing for node {target_node_id}: {e}", file=sys.stderr)
            # This could be an error in ProcessPoolExecutor itself or from worker dying (e.g. from sys.exit(1) in init_worker)
            continue  # Skip to next node

        print(f"✔ Worker task for node {target_node_id} finished in {time.time() - node_processing_start_time:.2f}s.")
        processed_nodes_count += 1

        # 4. Display Pileup (if requested)
        if args.view is not None:
            if view_data_for_node and isinstance(view_data_for_node,
                                                 dict) and view_data_for_node:  # If dict is not empty
                max_v_to_show = float('inf') if args.view == -1 else args.view
                view_limit_msg = f"first {int(max_v_to_show)}" if max_v_to_show != float('inf') else "all"
                print(
                    f"🔹 Displaying {view_limit_msg} pileups for node {target_node_id} (max {args.max_view_reads} reads/variant)...")
                display_pileup_data(view_data_for_node, node_id_str, node_sequence, args.max_view_reads, max_v_to_show)
            elif view_data_for_node is not None:  # view_data_for_node is {} (empty dict)
                print(
                    f"ℹ️ --view specified for node {target_node_id}, but no variants met AF threshold or were found to display.")
            # If view_data_for_node is None, warning already printed.

        # 5. Check for Output Summary
        node_summary_file_path = os.path.join(args.output, node_id_str, "variant_summary.json")
        if os.path.exists(node_summary_file_path):
            print(f"✅ Output generated for node {target_node_id} in: {os.path.dirname(node_summary_file_path)}")
            successful_nodes_with_output += 1
        elif view_data_for_node is not None and not view_data_for_node:  # Processed, no errors, but no variants passed AF
            print(f"ℹ️ No variants met AF threshold for node {target_node_id}. No tensor or summary files generated.")
        else:  # No summary, and view_data_for_node might be None (error) or {} (no variants but also no summary written for some reason)
            print(
                f"ℹ️ No output summary file found for node {target_node_id} at {node_summary_file_path}. Check logs if output was expected.")

        print(f"✔ Completed processing for node {target_node_id} in {time.time() - node_specific_start_time:.2f}s.")

    # --- End of Main Loop ---
    print("\n═════════════════════════════ PROCESSING COMPLETE ═════════════════════════════")

    # --- Save Updated Sequence Cache (if specified) ---
    if args.save_cache:
        if node_sequences_map:  # Only save if there's something in our map
            print(f"\n🔹 Saving/Updating {len(node_sequences_map)} sequence(s) to cache file: {args.save_cache}...")
            try:
                with open(args.save_cache, 'w') as wcf:
                    json.dump(node_sequences_map, wcf, indent=2)
                print(f"✔ Sequences saved to cache '{args.save_cache}'.")
            except Exception as e:
                print(f"❌ Error saving sequences to cache {args.save_cache}: {e}", file=sys.stderr)
        elif os.path.exists(
                args.save_cache):  # Cache file exists but map is empty (e.g. no nodes processed, or no sequences loaded/found)
            print(
                f"ℹ️ --save-cache specified ({args.save_cache}), but no sequences were loaded or fetched into memory. The cache file will not be modified if it exists and was not loaded from, or will be overwritten if it was the source of an empty load.")
            # To ensure consistency if the cache was loaded and then emptied by processing:
            if args.load_cache == args.save_cache and not node_sequences_map:
                try:  # Overwrite with empty if it was the loaded cache and now map is empty
                    with open(args.save_cache, 'w') as wcf:
                        json.dump({}, wcf, indent=2)
                    print(f"✔ Cache file {args.save_cache} effectively emptied as no sequences remain in map.")
                except Exception as e:
                    print(f"❌ Error writing empty cache {args.save_cache}: {e}", file=sys.stderr)

        else:  # No sequences in map and save_cache file does not exist
            print(f"ℹ️ --save-cache specified, but no sequences were loaded or fetched. No cache file created/updated.")

    total_script_time = time.time() - overall_start_time
    print(f"\nProcessed {processed_nodes_count}/{len(target_node_ids)} nodes.")
    if successful_nodes_with_output > 0:
        print(f"✅ Output (summary/tensors) successfully generated for {successful_nodes_with_output} node(s).")
    if processed_nodes_count < len(target_node_ids):
        print(
            f"⚠️ {len(target_node_ids) - processed_nodes_count} node(s) were skipped due to errors before task submission.")

    print(f"🏁 Script finished in {total_script_time:.2f} seconds.")


if __name__ == '__main__':
    main()