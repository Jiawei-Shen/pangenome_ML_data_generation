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
                sys.stderr.write(f"❌ Error: Index file {idx_path} is too small for node {target_node_id}.\n")
                return None
            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4:
                sys.stderr.write(
                    f"❌ Error: Could not read number of nodes from {idx_path} for node {target_node_id}.\n")
                return None
            num_nodes = struct.unpack('<I', num_nodes_bytes)[0]
            # Removed verbose print for every node lookup, can be re-added if debugging specific index issues
            # print(f"🔹 Index file contains {num_nodes} nodes. Searching for node {target_node_id}...")
            found = False
            for i in range(num_nodes):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                    sys.stderr.write(
                        f"❌ Error: Index file ended prematurely while reading record {i + 1}/{num_nodes} for node {target_node_id}.\n")
                    break
                node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                if node_id_from_idx == target_node_id:
                    node_info = (offset, n_records)
                    found = True
                    # print(f"✔ Found node {target_node_id} in index: Offset={offset}, N_Records={n_records}")
                    break
            if not found:
                sys.stderr.write(f"❌ Error: Node ID {target_node_id} not found in the index file {idx_path}.\n")
                return None
        return node_info
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: Index file not found at {idx_path}\n")
        return None
    except Exception as e:
        sys.stderr.write(f"❌ Error parsing index file {idx_path} for node {target_node_id}: {e}\n")
        return None


def load_multiple_node_sequences_from_gfa(gfa_path, target_node_ids_set):
    """
    Reads a GFA file once and extracts sequences for a specific set of target node IDs.
    Args:
        gfa_path (str): Path to the GFA file.
        target_node_ids_set (set): A set of integer node IDs to extract sequences for.
    Returns:
        dict: A dictionary mapping node_id (str) to sequence (str) for found nodes.
    """
    node_sequences = {}
    if not target_node_ids_set:
        return node_sequences

    nodes_to_find = target_node_ids_set.copy()  # Work on a copy

    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file to find sequences for {len(nodes_to_find)} nodes: {gfa_path}")
            line_counter = 0
            for line in f:
                line_counter += 1
                if line_counter % 10_000_000 == 0:  # Adjusted print frequency
                    print(
                        f"  Checked {line_counter:,} lines in GFA file... {len(nodes_to_find)} nodes remaining to find.")

                if not line.startswith('S\t'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue

                try:
                    nid_str_from_gfa = parts[1]
                    nid_int_from_gfa = int(nid_str_from_gfa)
                except ValueError:
                    # Potentially noisy if GFA has non-integer IDs not in target set
                    # print(f"⚠️ Warning: Could not parse node ID '{parts[1]}' as int in GFA line: {line.strip()}", file=sys.stderr)
                    continue

                if nid_int_from_gfa in nodes_to_find:
                    node_sequences[str(nid_int_from_gfa)] = parts[2]  # Store with str key
                    nodes_to_find.remove(nid_int_from_gfa)
                    if not nodes_to_find:
                        print(
                            f"✔ Found all {len(target_node_ids_set)} requested node sequences in GFA after checking {line_counter:,} lines.")
                        break

            if line_counter > 0:
                found_count = len(node_sequences)
                requested_count = len(target_node_ids_set)
                if found_count == requested_count:
                    if not nodes_to_find:  # Already printed success
                        pass
                    else:  # Should not happen if counts match and nodes_to_find is empty
                        print(
                            f"✔ Finished GFA scan after {line_counter:,} lines. Found all {found_count} requested sequences.")
                else:
                    print(
                        f"✔ Finished GFA scan after {line_counter:,} lines. Found {found_count}/{requested_count} requested sequences.")

            if nodes_to_find:  # Some requested nodes were not found
                print(
                    f"⚠️ Warning: Could not find GFA sequences for {len(nodes_to_find)} node ID(s). Examples: {list(nodes_to_find)[:5]}")

    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: GFA file not found at {gfa_path}\n")
        return {}
    except Exception as e:
        sys.stderr.write(f"❌ Error reading GFA file {gfa_path}: {e}\n")
        return node_sequences  # Return what was found so far

    return node_sequences


def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*':
        return []
    ops = []
    try:
        for length_str, op_char in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string):
            ops.append((int(length_str), op_char))
        return ops
    except Exception as e:
        sys.stderr.write(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}\n")
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
                return None  # Read ends before covering the target position within the match block
            current_node_pos += length
            current_read_pos += length
        elif op == 'I':
            # For an insertion, the variant position is the base *before* the insertion.
            # The allele is the inserted sequence itself.
            if expected_var_type == 'I' and (
                    current_node_pos - 1) == target_node_pos:  # target_node_pos is the anchor base
                return read_sequence[current_read_pos: current_read_pos + length].upper()
            current_read_pos += length
        elif op == 'D':
            # For a deletion, the variant position is the first deleted base on the reference.
            # The allele is '*' if the read has this deletion.
            if current_node_pos <= target_node_pos < current_node_pos + length:  # target_node_pos is within the deleted region
                if expected_var_type == 'I': return "OTHER_FOR_INDEL"  # Read has deletion where we expect insertion
                if expected_var_type == 'D':
                    # Check if the deletion in the read matches the expected deleted sequence
                    deleted_seq_in_read_context = node_sequence[current_node_pos: current_node_pos + length]
                    if deleted_seq_in_read_context == expected_ref_allele_for_indel:
                        return "*"  # Represents the deletion variant
                    else:
                        return "OTHER_FOR_INDEL"  # Deletion is present, but of different bases than variant
                return "*"  # Default for SNV check, a deletion is like a '*' allele at that ref pos
            current_node_pos += length
        elif op == 'S':  # Soft clip
            current_read_pos += length
        elif op == 'N':  # Skipped region from reference
            current_node_pos += length

        # Optimization: if current_node_pos has moved well past target_node_pos
        # (and it's not an insertion being checked at target_node_pos)
        if current_node_pos > target_node_pos + 1 and op in ('M', '=', 'X', 'D', 'N'):
            if not (expected_var_type == 'I' and (current_node_pos - 1) <= target_node_pos):  # allow I to check anchor
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
                    if node_base != read_base and op != '=':  # op can be 'M' or 'X'
                        variants.append((current_node_pos, 'X', read_base, node_base))
            node_pos += length
            read_pos += length
        elif op == 'I':  # Insertion to the reference
            inserted_sequence = read_sequence[read_pos: read_pos + length].upper()
            # Anchor position for insertion is the base *before* the insertion on the reference
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0
            ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if 0 <= ref_anchor_pos < len(
                node_sequence) else ""
            variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor if ref_base_at_anchor else "*"))
            read_pos += length
        elif op == 'D':  # Deletion from the reference
            deleted_sequence_from_ref = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= len(
                node_sequence) else ""
            if deleted_sequence_from_ref:  # Only add if the deletion is within node bounds
                # For deletion, ref allele is deleted sequence, alt is '*' (or sometimes first base of deleted seq)
                variants.append((node_pos, 'D', '*', deleted_sequence_from_ref))
            node_pos += length
        elif op == 'S':  # Soft clip, consumes read bases
            read_pos += length
        elif op == 'N':  # Skipped region, consumes reference
            node_pos += length
        # H and P do not consume read or reference bases in this context
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

        # Optimization: if the current reference position is already past the window, stop.
        if current_node_pos_in_read >= window_start_node + window_size and cigar_op in ('M', '=', 'X', 'D', 'N'):
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
        # Early exit if current reference position is beyond the window
        if current_node_pos_in_read >= window_start_node + current_tensor_window_size and cigar_op in ('M', 'D', 'N',
                                                                                                       '=', 'X'):
            break

        if cigar_op in ('M', '=', 'X'):  # Match, mismatch, or exact match
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                read_aln_pos = current_read_pos_in_read + i

                if read_aln_pos >= len(segment_read_sequence): break  # protect against CIGAR/read length mismatch

                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_char = segment_read_sequence[read_aln_pos].upper()
                    base_indices_row[window_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])

                    if read_aln_pos < len(segment_quality_str):  # Ensure qual string is long enough
                        try:
                            quality_scores_row[window_idx] = ord(segment_quality_str[read_aln_pos]) - 33
                        except (TypeError,
                                IndexError):  # TypeError if qual_str[read_aln_pos] is bad, IndexError if short
                            quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
                    else:
                        quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
            current_node_pos_in_read += cigar_len
            current_read_pos_in_read += cigar_len
        elif cigar_op == 'D' or cigar_op == 'N':  # Deletion or Skip in reference
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_indices_row[window_idx] = BASE_TO_INDEX['*']  # Deletion represented by '*'
                    quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING  # No quality for deletion
            current_node_pos_in_read += cigar_len
        elif cigar_op == 'I' or cigar_op == 'S':  # Insertion or Soft clip (consumes read bases only)
            current_read_pos_in_read += cigar_len

        # Early exit if read is consumed
        if current_read_pos_in_read >= len(segment_read_sequence) and cigar_op in ('M', '=', 'X', 'I', 'S'):
            break

    return base_indices_row, quality_scores_row


## ─────────────────────────────────────────────────────────────────────────────
## Worker Process Initialization and Target Function
## ─────────────────────────────────────────────────────────────────────────────

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker):
    global worker_dat_file, worker_base_output_dir
    # Each worker process will have its own file handle for the DAT file
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path_for_worker}\n")
        # Instead of sys.exit, let the error propagate or return a status
        # so the main process can handle it. However, ProcessPoolExecutor
        # might just hang or crash the worker. Best to ensure files exist before starting.
        # For now, keep exit, but main process will see this as a failed future.
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path_for_worker}: {e}\n")
        sys.exit(1)


def process_single_node_for_pileup(task_args_with_af_thresh):
    node_id, dat_file_offset, n_records, node_sequence, min_af_threshold = task_args_with_af_thresh
    global worker_dat_file, worker_base_output_dir  # These are set by init_worker

    pth_files_generated_for_node = 0

    if worker_dat_file is None or worker_base_output_dir is None:
        # This case should ideally not happen if init_worker is successful.
        sys.stderr.write(
            f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Worker not initialized (DAT file or output_dir missing).\n")
        return node_id, None, pth_files_generated_for_node  # Signal error with None for view_data

    if not node_sequence:  # Should be caught before submitting task, but double check.
        sys.stderr.write(f"ℹ️ [Worker {os.getpid()} for Node {node_id}]: No sequence provided. Skipping.\n")
        return node_id, {}, pth_files_generated_for_node  # Return empty dict for view_data

    node_specific_output_dir = os.path.join(worker_base_output_dir, str(node_id))
    try:
        os.makedirs(node_specific_output_dir, exist_ok=True)
    except OSError as e:
        sys.stderr.write(
            f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Could not create dir {node_specific_output_dir}: {e}\n")
        return node_id, None, pth_files_generated_for_node

    node_len = len(node_sequence)
    variant_headers_for_node = []
    view_oriented_variant_data = {}  # For the --view option

    aligned_read_segments = []
    try:
        # Each worker has its own file handle, so seek is safe.
        worker_dat_file.seek(dat_file_offset + 10)  # Assuming +10 is for a header in the .dat node block
        for record_idx in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE:
                sys.stderr.write(
                    f"⚠️ Warning [Node {node_id}]: Truncated record {record_idx + 1}/{n_records}. Read {len(data)}/{RECORD_SIZE} bytes.\n")
                break

            off_from_file, raw_seq, raw_qual, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)

            if mapq < 10: continue  # Filter by MAPQ

            try:
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')
                qual_str = raw_qual.rstrip(b'\x00').decode('ascii', errors='replace')
                cigar_str_original = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                # sys.stderr.write(f"⚠️ Warning [Node {node_id}]: Unicode decode error in record {record_idx}. Skipping read.\n")
                continue

            if len(seq) == 0 or len(seq) != len(qual_str):  # Basic sanity check
                # sys.stderr.write(f"⚠️ Warning [Node {node_id}]: Empty seq or seq/qual length mismatch. Seq: {len(seq)}, Qual: {len(qual_str)}. Skipping read.\n")
                continue

            original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
            if not original_decoded_cigar_ops and cigar_str_original != '*':  # '*' CIGAR is valid (unmapped or no ops)
                # sys.stderr.write(f"⚠️ Warning [Node {node_id}]: Could not decode CIGAR '{cigar_str_original}'. Skipping read.\n")
                continue  # Skip if CIGAR is invalid but not '*'

            current_read_sequence = seq
            current_quality_str = qual_str
            current_decoded_cigar_ops = list(original_decoded_cigar_ops)  # Make a mutable copy
            current_offset_on_node = off_from_file

            if strand_char == '-':
                current_read_sequence = reverse_complement(seq)
                current_quality_str = qual_str[::-1]

                # Reverse CIGAR operations
                temp_rev_cigar_ops = []
                if original_decoded_cigar_ops:  # Ensure it's not empty
                    for length, op_code in original_decoded_cigar_ops:
                        temp_rev_cigar_ops.insert(0, (length, op_code))
                current_decoded_cigar_ops = temp_rev_cigar_ops

                # Adjust offset for reverse strand based on alignment span on reference
                # The original script had a "User-confirmed fix" for alignment_span_on_node.
                # We should use CIGAR to determine reference span for robust offset calculation.
                alignment_span_on_ref = sum(l for l, op in current_decoded_cigar_ops if op in ('M', 'D', 'N', '=', 'X'))
                if alignment_span_on_ref > 0:
                    # This is a common way to adjust start for reverse strand alignments
                    # if off_from_file is original start on forward strand of node
                    current_offset_on_node = node_len - (off_from_file + alignment_span_on_ref)
                    # If original script's logic "node_len - len(current_read_sequence) - off_from_file" was correct
                    # it implies 'off_from_file' for '-' strand has a different meaning or
                    # 'len(current_read_sequence)' was a proxy for span.
                    # Reverting to user's confirmed logic for safety, but with a note:
                    # current_offset_on_node = node_len - len(current_read_sequence) - off_from_file
                    # This assumes off_from_file is from the "left" end even on reverse strand mapping.
                    # A more standard SAM-like approach: if POS is always 5'-most on reference strand sense:
                    # current_offset_on_node = off_from_file # POS is POS. CIGAR handles strand.
                    # Given the original code, will stick to its reverse strand offset logic.
                    # Original: alignment_span_on_node = len(current_read_sequence)
                    # current_offset_on_node = node_len - alignment_span_on_node - off_from_file
                    # This means off_from_file for a reverse strand read is its distance from the *end* of the node segment it aligns to.
                    # Let's use the user's logic:
                    alignment_span_on_node_user_logic = len(current_read_sequence)  # As per original user fix comment
                    current_offset_on_node = node_len - alignment_span_on_node_user_logic - off_from_file

                    if current_offset_on_node < 0:
                        # This can happen if node_len is too small or off_from_file is large.
                        # sys.stderr.write(f"⚠️ Node {node_id}, read (rev strand): calculated negative offset {current_offset_on_node}. Skipping.\n")
                        continue

            aligned_read_segments.append({
                "offset_on_node": current_offset_on_node,
                "read_sequence": current_read_sequence,
                "processed_quality_str": current_quality_str,
                "cigar_ops": current_decoded_cigar_ops,
                "original_cigar_str": cigar_str_original,
                "strand": strand_char
            })
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()} for Node {node_id}] reading records: {e}\n")
        # import traceback
        # traceback.print_exc()
        return node_id, None, pth_files_generated_for_node  # Error state

    if not aligned_read_segments:
        # print(f"ℹ️ [Node {node_id}] No aligned read segments after processing DAT records (e.g. all low MAPQ or errors).")
        return node_id, {}, pth_files_generated_for_node

    # --- Variant Detection and Processing ---
    candidate_variants = defaultdict(int)
    for segment in aligned_read_segments:
        variants_in_read = detect_variants_from_cigar(
            segment["offset_on_node"], segment["cigar_ops"],
            segment["read_sequence"], node_sequence
        )
        for v_pos, v_type, v_alt, v_ref in variants_in_read:
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1  # v_ref from CIGAR, v_alt from CIGAR

    variant_processing_window_size = TENSOR_WINDOW_SIZE  # Same window for tensor and view gen
    half_window = variant_processing_window_size // 2

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), read_support_count in candidate_variants.items():
        # AF calculation loop
        af_alt_count, af_ref_count, af_other_count, af_locus_coverage = 0, 0, 0, 0

        # Define expected ref/alt based on variant type for AF counting
        expected_ref_for_af = ""
        expected_alt_for_af = ""

        if v_type == 'X':  # SNP or MNP
            expected_ref_for_af = v_ref_from_cigar  # The actual reference base(s) at this position
            expected_alt_for_af = v_alt_from_cigar  # The variant base(s) seen in reads
        elif v_type == 'D':  # Deletion
            # For deletion, ref is the deleted sequence, alt is often represented as '*' or first base of ref.
            # get_allele_from_read_at_node_pos returns '*' for deletion.
            expected_ref_for_af = node_sequence[v_pos] if v_pos < node_len else ""  # anchor base on ref
            expected_alt_for_af = "*"
        elif v_type == 'I':  # Insertion
            # For insertion, ref is the base before insertion, alt is the inserted sequence.
            expected_ref_for_af = v_ref_from_cigar  # anchor base on ref (could be empty if at start)
            expected_alt_for_af = v_alt_from_cigar  # the inserted sequence

        for segment in aligned_read_segments:
            allele = get_allele_from_read_at_node_pos(
                segment["offset_on_node"], segment["read_sequence"], segment["cigar_ops"],
                v_pos, node_sequence,
                expected_var_type=v_type, expected_ref_allele_for_indel=v_ref_from_cigar if v_type == 'D' else (
                    node_sequence[v_pos] if v_type == 'I' and v_pos < node_len else None)
            )

            if allele is not None:
                af_locus_coverage += 1
                if v_type == 'X':
                    if allele == expected_alt_for_af:
                        af_alt_count += 1
                    elif allele == expected_ref_for_af:
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'I':  # Allele is inserted sequence or "REF_STATE_FOR_INDEL"
                    if allele == expected_alt_for_af:
                        af_alt_count += 1
                    elif allele == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1
                    else:
                        af_other_count += 1  # e.g. different insertion
                elif v_type == 'D':  # Allele is '*' for deletion or "REF_STATE_FOR_INDEL"
                    if allele == expected_alt_for_af:
                        af_alt_count += 1  # allele is '*'
                    elif allele == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1
                    else:
                        af_other_count += 1  # e.g. SNV at site of deletion call.

        alt_freq = af_alt_count / af_locus_coverage if af_locus_coverage > 0 else 0.0
        if alt_freq < min_af_threshold:
            continue

        variant_key_str = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"  # Use CIGAR-derived alleles for key

        # Windowing logic for tensor/view
        if v_type == 'I':  # Insertion happens AFTER v_pos
            window_center_on_node = v_pos + 1
        else:  # SNP/Deletion are AT v_pos
            window_center_on_node = v_pos
        current_window_start_on_node = window_center_on_node - half_window
        # Ensure window start is not negative
        current_window_start_on_node = max(0, current_window_start_on_node)

        # --- Prepare data for --view option ---
        pileup_reads_data_for_view = []
        if aligned_read_segments:  # Check if there are any reads to process for view
            for segment in aligned_read_segments:  # Iterate up to TENSOR_MAX_READ_ROWS for view consistency?
                row_chars_for_view = get_read_representation_in_window_for_view(
                    segment["cigar_ops"], segment["offset_on_node"], segment["read_sequence"],
                    current_window_start_on_node, variant_processing_window_size, node_len
                )
                if any(c != ' ' for c in row_chars_for_view):  # If read overlaps window
                    row_indices_for_view = [BASE_TO_INDEX.get(char.upper(), BASE_TO_INDEX['N']) for char in
                                            row_chars_for_view]
                    pileup_reads_data_for_view.append({
                        "bases": row_indices_for_view, "offset": segment["offset_on_node"],
                        "strand": segment["strand"], "cigar": segment["original_cigar_str"]
                    })

        view_oriented_variant_data[variant_key_str] = {
            "pileup_reads_data": pileup_reads_data_for_view[:TENSOR_MAX_READ_ROWS],  # Match tensor depth for view
            "alt_allele_count": af_alt_count,
            "ref_allele_count_at_locus": af_ref_count,
            "other_allele_count_at_locus": af_other_count,
            "coverage_at_locus": af_locus_coverage,
            "alt_allele_frequency": round(alt_freq, 4)
        }

        # --- Prepare Tensor Data ---
        tensor_ch1_bases_list = []
        tensor_ch2_qualities_list = []
        tensor_ch3_mismatches_list = []

        # Reference Row for Tensor
        ref_base_indices_row_tensor = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        ref_qual_scores_row_tensor = [DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE  # Qualities for ref are padding
        ref_mismatch_row_tensor = [MISMATCH_CHANNEL_REF_ROW_VALUE] * TENSOR_WINDOW_SIZE  # Ref matches itself perfectly

        for i in range(TENSOR_WINDOW_SIZE):
            actual_node_pos = current_window_start_on_node + i
            if 0 <= actual_node_pos < node_len:
                ref_base_indices_row_tensor[i] = BASE_TO_INDEX.get(node_sequence[actual_node_pos].upper(),
                                                                   BASE_TO_INDEX['N'])

        tensor_ch1_bases_list.append(ref_base_indices_row_tensor)
        tensor_ch2_qualities_list.append(ref_qual_scores_row_tensor)
        tensor_ch3_mismatches_list.append(ref_mismatch_row_tensor)

        reads_added_to_tensor = 0
        if aligned_read_segments:  # Check if there are any reads to process for tensor
            for segment in aligned_read_segments:
                if reads_added_to_tensor >= TENSOR_MAX_READ_ROWS:
                    break

                base_indices_row_tensor, quality_scores_row_tensor = get_read_tensor_rows_in_window(
                    segment["cigar_ops"], segment["offset_on_node"],
                    segment["read_sequence"], segment["processed_quality_str"],
                    current_window_start_on_node, TENSOR_WINDOW_SIZE, node_len
                )

                if any(b != PADDING_BASE_INDEX for b in
                       base_indices_row_tensor):  # If read overlaps window meaningfully
                    tensor_ch1_bases_list.append(base_indices_row_tensor)
                    tensor_ch2_qualities_list.append(quality_scores_row_tensor)

                    # Mismatch Channel calculation for this read against reference row
                    mismatch_row_for_read_tensor = [MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE
                    for i in range(TENSOR_WINDOW_SIZE):
                        read_base_idx = base_indices_row_tensor[i]
                        ref_base_idx_for_comp = ref_base_indices_row_tensor[i]  # The ref row already prepared

                        if read_base_idx == PADDING_BASE_INDEX or ref_base_idx_for_comp == PADDING_BASE_INDEX:
                            mismatch_row_for_read_tensor[i] = MISMATCH_COMPARISON_PADDING_VALUE
                        elif read_base_idx == ref_base_idx_for_comp:  # Match
                            mismatch_row_for_read_tensor[i] = 0
                        else:  # Mismatch (includes actual mismatches, deletions in read vs ref, insertions etc.)
                            mismatch_row_for_read_tensor[i] = 1
                    tensor_ch3_mismatches_list.append(mismatch_row_for_read_tensor)
                    reads_added_to_tensor += 1

        # Pad tensor if fewer than TENSOR_MAX_READ_ROWS were added
        for _ in range(TENSOR_MAX_READ_ROWS - reads_added_to_tensor):  # (MAX_READ_ROWS for reads + 1 for ref)
            tensor_ch1_bases_list.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            tensor_ch2_qualities_list.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            tensor_ch3_mismatches_list.append([MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE)

        try:
            # Final tensor should have shape (Channels, Reads+Ref, WindowSize)
            # Channels = 3 (bases, quals, mismatches)
            # Reads+Ref = TENSOR_MAX_READ_ROWS + 1
            final_tensor = torch.tensor([
                tensor_ch1_bases_list,  # Shape: (MAX_READS+1, WINDOW_SIZE)
                tensor_ch2_qualities_list,  # Shape: (MAX_READS+1, WINDOW_SIZE)
                tensor_ch3_mismatches_list  # Shape: (MAX_READS+1, WINDOW_SIZE)
            ], dtype=torch.int8)  # Use int8 for memory efficiency

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
            pth_files_generated_for_node += 1
        except Exception as e:
            sys.stderr.write(
                f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Failed to create/save tensor for {variant_key_str}: {e}\n")
            # import traceback
            # traceback.print_exc()

    # After processing all variants for the node:
    if variant_headers_for_node:  # Only write summary if variants were found and processed
        summary_json_path = os.path.join(node_specific_output_dir, "variant_summary.json")
        try:
            with open(summary_json_path, 'w') as sjf:
                json.dump({
                    "node_id": node_id,
                    "node_length": node_len,
                    "node_sequence_preview": node_sequence[:100] + ("..." if node_len > 100 else ""),
                    "variants_passing_af_filter": variant_headers_for_node  # List of variant info dicts
                }, sjf, indent=2)
        except Exception as e:
            sys.stderr.write(f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Failed to write summary JSON: {e}\n")

    return node_id, view_oriented_variant_data, pth_files_generated_for_node


## ─────────────────────────────────────────────────────────────────────────────
## Pileup Viewing Function (display_pileup_data)
## ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if not node_data_for_display_view or not isinstance(node_data_for_display_view, dict):
        print(f"ℹ️ No valid pileup data to display for node {node_id_str_for_display}.", file=sys.stderr)
        return

    print(f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")

    if not node_data_for_display_view:  # Empty dict
        print(f"ℹ️ No variants met AF threshold or were found to display for node {node_id_str_for_display}.")
        return

    variants_displayed_count = 0
    # Sort variant keys: primarily by position (integer part), then by type (string part)
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))

    display_window_size = TENSOR_WINDOW_SIZE  # Match tensor window for consistency
    half_display_window = display_window_size // 2

    for variant_key in sorted_variant_keys:
        if variants_displayed_count >= max_variants_to_display:
            print(
                f"\n  ... (and {len(node_data_for_display_view) - variants_displayed_count} more variants for node {node_id_str_for_display} not shown due to display limit)")
            break

        variant_data = node_data_for_display_view[variant_key]
        pileup_reads_display_data = variant_data.get("pileup_reads_data", [])
        parts = variant_key.split('_')
        v_pos = int(parts[0])
        v_type = parts[1]
        # v_ref = parts[2]
        # v_alt = parts[3]

        # Determine window center based on variant type for display
        if v_type == 'I':  # Insertion is after v_pos
            window_center_on_node_display = v_pos + 1
        else:  # SNV/Deletion is at v_pos
            window_center_on_node_display = v_pos

        current_window_start_on_node_display = window_center_on_node_display - half_display_window
        current_window_start_on_node_display = max(0, current_window_start_on_node_display)  # Don't go < 0

        print(f"\n--- Variant: {variant_key} (Node Pos: {v_pos}, Type: {v_type}) ---")
        print(
            f"  Display Window on Node: {current_window_start_on_node_display} - {current_window_start_on_node_display + display_window_size - 1}")
        print(f"  Alt Count: {variant_data.get('alt_allele_count', 'N/A')}")
        print(f"  Ref Count: {variant_data.get('ref_allele_count_at_locus', 'N/A')}")
        print(f"  Other Count: {variant_data.get('other_allele_count_at_locus', 'N/A')}")
        print(f"  Coverage: {variant_data.get('coverage_at_locus', 'N/A')}")
        alt_freq_val = variant_data.get('alt_allele_frequency', 'N/A')
        if isinstance(alt_freq_val, float):
            print(f"  Alt Freq: {alt_freq_val:.4f}")
        else:
            print(f"  Alt Freq: {alt_freq_val}")

        # Reference sequence display
        ref_display_parts = [' '] * display_window_size
        marker_line_parts = [' '] * display_window_size  # For '^' and 'I' markers

        # Determine index of the variant within the current display window
        variant_display_idx_in_window = -1
        if v_type != 'I':  # SNV or Deletion
            if current_window_start_on_node_display <= v_pos < current_window_start_on_node_display + display_window_size:
                variant_display_idx_in_window = v_pos - current_window_start_on_node_display
        elif v_type == 'I':  # Insertion (anchored at v_pos, occurs after it)
            # The marker for insertion is typically at the base *before* the insertion.
            if current_window_start_on_node_display <= v_pos < current_window_start_on_node_display + display_window_size:
                variant_display_idx_in_window = v_pos - current_window_start_on_node_display

        for i in range(display_window_size):
            actual_node_pos_in_window = current_window_start_on_node_display + i
            if 0 <= actual_node_pos_in_window < len(full_node_sequence):
                ref_display_parts[i] = full_node_sequence[actual_node_pos_in_window]

            if i == variant_display_idx_in_window:
                if v_type == 'I':
                    marker_line_parts[i] = "I"  # Mark position before insertion
                    # Optionally point between i and i+1 for insertion
                    # if i + 1 < display_window_size: marker_line_parts[i+1] = "^"
                else:  # SNV, Deletion
                    marker_line_parts[i] = "^"  # Mark the affected base

        print(f"  Node Ref: {''.join(ref_display_parts)}")
        print(f"  Marker  : {''.join(marker_line_parts)}")

        if not pileup_reads_display_data:
            print("  (No reads data available in this window for display)")
        else:
            displayed_reads_count_for_variant = 0
            for i, read_info in enumerate(pileup_reads_display_data):
                if displayed_reads_count_for_variant >= max_reads_to_display_per_variant:
                    print(
                        f"  ... (and {len(pileup_reads_display_data) - displayed_reads_count_for_variant} more reads for this variant not shown due to view limit)")
                    break
                base_indices = read_info["bases"]
                read_offset = read_info["offset"]
                read_strand = read_info["strand"]
                read_cigar = read_info.get("cigar", "N/A")  # Get CIGAR if available

                pileup_row_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in base_indices])
                print(
                    f"  Read {i + 1:3d}: {pileup_row_str}  (Offset: {read_offset}, Strand: {read_strand}, CIGAR: {read_cigar})")
                displayed_reads_count_for_variant += 1
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

    input_node_group = parser.add_mutually_exclusive_group(required=True)
    input_node_group.add_argument("--node_id", type=int, help="The specific node ID to process.")
    input_node_group.add_argument("--node_id_file",
                                  help="Path to a text file containing node IDs (one per line) to process.")

    parser.add_argument("--gfa",
                        help="GFA graph file path (required if node sequence cache is not used/built for any specified node).")
    parser.add_argument("--load-cache", help="Load node sequences from this JSON cache file.")
    parser.add_argument("--save-cache",
                        help="Save/update node sequences to this JSON cache file (used if --gfa is provided to fetch new sequences).")

    parser.add_argument("--num_workers", type=int, default=None,
                        help="Number of worker processes for parallel node processing. Defaults to number of CPU cores.")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N_VARIANTS',
                        help="Print generated pileups to console. If used without value, shows all variants. Optionally specify N for first N variants per node. -1 for all.")
    parser.add_argument("--max_view_reads", type=int, default=20,
                        help="Max reads to display per pileup in console view.")
    parser.add_argument("--min_af", type=float, default=0.1,
                        help="Min allele frequency for a variant to be processed and included in output.")
    args = parser.parse_args()

    if not os.path.isfile(args.dat): sys.exit(f"❌ Error: DAT file not found: {args.dat}")
    if not os.path.isfile(args.idx): sys.exit(f"❌ Error: Index file not found: {args.idx}")
    if not args.load_cache and not args.gfa:
        sys.exit("❌ Error: Must provide --gfa (to fetch sequences) or --load-cache (with all needed sequences).")
    if args.load_cache and not os.path.isfile(args.load_cache) and os.path.exists(
            args.load_cache):  # exists but not a file
        sys.exit(f"❌ Error: Specified --load-cache path '{args.load_cache}' is not a file.")
    if args.gfa and not os.path.isfile(args.gfa):
        sys.exit(f"❌ Error: GFA file not found: {args.gfa}")
    if not (0.0 <= args.min_af <= 1.0):
        sys.exit("❌ Error: --min_af must be between 0.0 and 1.0.")

    num_workers = args.num_workers if args.num_workers and args.num_workers > 0 else os.cpu_count()
    if num_workers is None: num_workers = 1  # Fallback if os.cpu_count() returns None
    print(f"🔹 Using {num_workers} worker process(es) for parallel node processing.")

    try:
        os.makedirs(args.output, exist_ok=True)
        print(f"🔹 Base output directory: {args.output}")
    except OSError as e:
        sys.exit(f"❌ Error: Could not create base output directory {args.output}: {e}")

    target_node_ids_set = set()
    if args.node_id_file:
        try:
            with open(args.node_id_file, 'r') as f_nodes:
                for line_num, line in enumerate(f_nodes, 1):
                    line = line.strip()
                    if line and not line.startswith("#"):
                        try:
                            target_node_ids_set.add(int(line))
                        except ValueError:
                            sys.stderr.write(
                                f"⚠️ Warning: Invalid non-integer node ID '{line}' in {args.node_id_file} at line {line_num}. Skipping.\n")
            if not target_node_ids_set:
                sys.exit(f"❌ Error: No valid node IDs found or specified in {args.node_id_file}.")
            print(f"🔹 Will process {len(target_node_ids_set)} unique node ID(s) from file: {args.node_id_file}")
        except FileNotFoundError:
            sys.exit(f"❌ Error: Node ID file not found: {args.node_id_file}")
    elif args.node_id is not None:
        target_node_ids_set.add(args.node_id)
        print(f"🔹 Will process single target node ID: {args.node_id}")

    if not target_node_ids_set:
        print("ℹ️ No target node IDs specified. Exiting.")
        sys.exit(0)

    # --- Sequence Loading ---
    node_sequences_map = {}  # Stores node_id_str -> sequence
    if args.load_cache and os.path.isfile(args.load_cache):  # Check if it's a file before trying to open
        start_time_cache_load = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                node_sequences_map = json.load(cf)  # Assumes keys are strings
            print(
                f"✔ Loaded {len(node_sequences_map)} sequences from cache '{args.load_cache}' in {time.time() - start_time_cache_load:.2f}s.")
        except json.JSONDecodeError as e:
            sys.stderr.write(
                f"⚠️ Warning: Error decoding JSON from cache file {args.load_cache}: {e}. Cache not loaded.\n")
        except Exception as e:
            sys.stderr.write(f"⚠️ Warning: Error loading cache {args.load_cache}: {e}. Cache not loaded.\n")

    nodes_needing_gfa_sequence = set()
    if args.gfa:  # Only try to fetch from GFA if GFA is provided
        for node_id_int in target_node_ids_set:
            if str(node_id_int) not in node_sequences_map:
                nodes_needing_gfa_sequence.add(node_id_int)

    if nodes_needing_gfa_sequence and args.gfa:
        print(f"🔹 {len(nodes_needing_gfa_sequence)} node(s) require sequence fetching from GFA.")
        gfa_fetch_start_time = time.time()
        fetched_sequences = load_multiple_node_sequences_from_gfa(args.gfa, nodes_needing_gfa_sequence)
        if fetched_sequences:
            node_sequences_map.update(fetched_sequences)  # Add newly fetched sequences to the map
            print(
                f"✔ Fetched {len(fetched_sequences)} new sequences from GFA in {time.time() - gfa_fetch_start_time:.2f}s. Total sequences in map: {len(node_sequences_map)}.")
        else:
            print(
                f"ℹ️ No new sequences were fetched from GFA for the {len(nodes_needing_gfa_sequence)} nodes needing them.")
    elif nodes_needing_gfa_sequence and not args.gfa:
        sys.stderr.write(
            f"⚠️ Warning: {len(nodes_needing_gfa_sequence)} node(s) need sequences, but --gfa was not provided. These nodes will be skipped if not in cache.\n")

    # --- Prepare Tasks for Workers ---
    tasks_to_submit = []
    skipped_nodes_pre_submit = set()
    for node_id_int in target_node_ids_set:  # Iterate over the original set of targets
        node_id_str = str(node_id_int)

        node_sequence = node_sequences_map.get(node_id_str)
        if not node_sequence:
            if node_id_int not in skipped_nodes_pre_submit:  # Avoid duplicate messages
                # sys.stderr.write(f"ℹ️ Skipping node {node_id_int}: Sequence not available (not in cache or GFA).\n")
                skipped_nodes_pre_submit.add(node_id_int)
            continue

        idx_parse_start_time = time.time()
        node_dat_info = parse_idx_file_for_single_node(args.idx, node_id_int)
        # print(f"Index parsing for node {node_id_int} took {time.time() - idx_parse_start_time:.3f}s")
        if not node_dat_info:
            if node_id_int not in skipped_nodes_pre_submit:
                # sys.stderr.write(f"ℹ️ Skipping node {node_id_int}: Failed to parse index information.\n")
                skipped_nodes_pre_submit.add(node_id_int)
            continue

        dat_offset, n_records = node_dat_info
        tasks_to_submit.append((node_id_int, dat_offset, n_records, node_sequence, args.min_af))

    if skipped_nodes_pre_submit:
        print(
            f"⚠️ Skipped {len(skipped_nodes_pre_submit)} node(s) before submission due to missing sequences or index errors. Examples: {list(skipped_nodes_pre_submit)[:5]}")

    if not tasks_to_submit:
        print("ℹ️ No valid tasks to process after sequence and index checks. Exiting.")
        if args.save_cache:  # Still save cache even if no tasks
            # ... (cache saving logic as at the end of script) ...
            pass
        sys.exit(0)

    # --- Parallel Processing ---
    overall_start_time = time.time()
    total_pth_files_generated = 0
    processed_nodes_count = 0  # Nodes for which worker returned a result (not necessarily successful with output)
    successful_nodes_with_output = 0  # Nodes that generated .pth or summary.json

    results_for_viewing = {}  # Store node_id -> (view_data, node_sequence_for_view)

    print(f"\n🔹 Submitting {len(tasks_to_submit)} node tasks to {num_workers} worker(s)...")
    with ProcessPoolExecutor(max_workers=num_workers, initializer=init_worker,
                             initargs=(args.dat, args.output)) as executor:
        # Map future to node_id for easier error reporting and result association
        future_to_node_id = {executor.submit(process_single_node_for_pileup, task_payload): task_payload[0] for
                             task_payload in tasks_to_submit}

        for i, future in enumerate(as_completed(future_to_node_id)):
            original_node_id = future_to_node_id[future]
            try:
                node_id_result, view_data_for_node, pth_files_for_this_node = future.result()

                processed_nodes_count += 1  # A worker completed for this node

                if node_id_result is None:  # Indicates a severe error within the worker before returning structured data
                    sys.stderr.write(
                        f"❌ Critical Error: Worker task for node {original_node_id} failed to return valid ID. Check worker logs.\n")
                    continue

                total_pth_files_generated += pth_files_for_this_node

                node_output_dir = os.path.join(args.output, str(node_id_result))
                summary_file = os.path.join(node_output_dir, "variant_summary.json")
                if pth_files_for_this_node > 0 or (
                os.path.exists(summary_file)):  # consider successful if pth or summary exists
                    successful_nodes_with_output += 1

                print(
                    f"✔ [{i + 1}/{len(tasks_to_submit)}] Completed node {node_id_result}. Generated {pth_files_for_this_node} .pth file(s).")

                if args.view is not None:
                    if view_data_for_node and isinstance(view_data_for_node, dict) and view_data_for_node:
                        # Retrieve the correct node_sequence for display (it was sent to worker, not returned)
                        current_node_sequence_for_view = node_sequences_map.get(str(node_id_result))
                        if current_node_sequence_for_view:
                            results_for_viewing[node_id_result] = (view_data_for_node, current_node_sequence_for_view)
                        else:
                            # Should not happen if task was submitted with sequence
                            sys.stderr.write(
                                f"⚠️ Warning: Could not retrieve sequence for node {node_id_result} for viewing, though task completed.\n")
                    elif view_data_for_node is not None:  # view_data_for_node is {}
                        # print(f"ℹ️ --view specified for node {node_id_result}, but worker found no variants meeting AF or to display.")
                        pass


            except Exception as exc:
                processed_nodes_count += 1  # Still count as processed attempt that failed
                sys.stderr.write(f"❌ Error processing node {original_node_id} in worker: {exc}\n")
                # import traceback
                # traceback.print_exc()

    # --- Display Pileups (after all processing done, if --view) ---
    if args.view is not None and results_for_viewing:
        print("\n═════════════════════════════ VIEWING PILEUPS ═════════════════════════════")
        # Sort by node_id for consistent view order if multiple nodes were processed
        for node_id_to_view in sorted(results_for_viewing.keys()):
            view_data, node_seq_for_view = results_for_viewing[node_id_to_view]

            max_v_to_show = float('inf') if args.view == -1 else args.view  # Max variants per node
            view_limit_msg = f"first {int(max_v_to_show)}" if max_v_to_show != float('inf') else "all"

            # This print is outside display_pileup_data to show context before potential long output
            # print(f"🔹 Displaying {view_limit_msg} variant pileups for node {node_id_to_view} (max {args.max_view_reads} reads/variant)...")
            display_pileup_data(view_data, str(node_id_to_view), node_seq_for_view, args.max_view_reads, max_v_to_show)
    elif args.view is not None and not results_for_viewing:
        print(
            f"ℹ️ --view specified, but no pileup data was generated or collected for display from any processed nodes.")

    # --- Final Summary and Cache Saving ---
    print("\n═════════════════════════════ PROCESSING COMPLETE ═════════════════════════════")
    total_script_time = time.time() - overall_start_time

    if args.save_cache:
        if node_sequences_map:  # Only save if there's something to save
            print(f"\n🔹 Saving/Updating {len(node_sequences_map)} sequence(s) to cache file: {args.save_cache}...")
            try:
                with open(args.save_cache, 'w') as wcf:
                    json.dump(node_sequences_map, wcf, indent=2)  # Keys are already strings
                print(f"✔ Sequences saved to cache '{args.save_cache}'.")
            except Exception as e:
                sys.stderr.write(f"❌ Error saving sequences to cache {args.save_cache}: {e}\n")
        elif os.path.exists(args.save_cache) and args.load_cache == args.save_cache and not node_sequences_map:
            # If loaded an empty/invalid cache and no new sequences, effectively empty it
            print(
                f"ℹ️ --save-cache specified ({args.save_cache}), but no sequences were in map (e.g., loaded empty/failed or GFA fetch failed).")
            try:
                with open(args.save_cache, 'w') as wcf:
                    json.dump({}, wcf)
                print(f"✔ Cache file {args.save_cache} written as empty (as no sequences were in memory).")
            except Exception as e:
                sys.stderr.write(f"❌ Error writing empty cache to {args.save_cache}: {e}\n")
        else:  # No map and not an overwrite of loaded empty cache
            print(
                f"ℹ️ --save-cache specified, but no sequences were in memory to save. Cache file '{args.save_cache}' not created/updated.")

    print(f"\nSummary:")
    print(f"  Targeted {len(target_node_ids_set)} unique node ID(s).")
    if skipped_nodes_pre_submit:
        print(f"  Skipped {len(skipped_nodes_pre_submit)} node(s) before submission (missing sequence/index).")
    print(f"  Submitted {len(tasks_to_submit)} node(s) to workers.")
    print(f"  Worker tasks completed for {processed_nodes_count}/{len(tasks_to_submit)} submitted node(s).")
    if successful_nodes_with_output > 0:
        print(f"✅ Output (summary/tensor files) generated for {successful_nodes_with_output} node(s).")
    else:
        print(f"ℹ️ No output files (summary/tensor) were generated for any node (check AF threshold or worker logs).")

    print(f"ℹ️ A grand total of {total_pth_files_generated} .pth tensor files were generated across all nodes.")
    print(f"🏁 Script finished in {total_script_time:.2f} seconds.")


if __name__ == '__main__':
    main()