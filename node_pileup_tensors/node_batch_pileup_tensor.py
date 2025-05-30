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
            # print(f"🔹 Index file contains {num_nodes} nodes. Searching for node {target_node_id}...") # Less verbose
            found = False
            for i in range(num_nodes):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                    # print(f"❌ Error: Index file ended prematurely while reading record {i + 1}/{num_nodes}.", # Less verbose
                    #       file=sys.stderr)
                    break
                node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                if node_id_from_idx == target_node_id:
                    node_info = (offset, n_records)
                    found = True
                    # print(f"✔ Found node {target_node_id} in index: Offset={offset}, N_Records={n_records}") # Less verbose
                    break
            if not found:
                print(f"❌ Error: Node ID {target_node_id} not found in the index file {idx_path}.", file=sys.stderr)
                return None
        return node_info
    except FileNotFoundError:
        print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
        return None


def load_node_sequence_from_gfa(gfa_path, target_node_id):
    node_sequence = None
    try:
        with open(gfa_path, 'r') as f:
            # print(f"🔹 Reading GFA file to find sequence for node {target_node_id}: {gfa_path}") # Less verbose
            line_counter = 0
            for line in f:
                line_counter += 1
                # if line_counter % 5_000_000 == 0: # Less verbose
                #     print(f"  Checked {line_counter:,} lines in GFA file...")
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
                    # print(f"✔ Found sequence for node {target_node_id} in GFA.") # Less verbose
                    break
            # if line_counter > 0: print(f"✔ Finished GFA scan after {line_counter:,} lines.") # Less verbose
            if node_sequence is None:
                print(f"❌ Error: Sequence for node ID {target_node_id} not found in GFA file {gfa_path}.",
                      file=sys.stderr)
    except FileNotFoundError:
        print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr)
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
        # print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr) # Can be noisy
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
                return None  # Read sequence exhausted
            current_node_pos += length
            current_read_pos += length
        elif op == 'I':
            if expected_var_type == 'I' and (
                    current_node_pos - 1) == target_node_pos:  # Insertion occurs after target_node_pos
                return read_sequence[current_read_pos: current_read_pos + length].upper()
            current_read_pos += length
        elif op == 'D':
            if current_node_pos <= target_node_pos < current_node_pos + length:  # Deletion spans target_node_pos
                if expected_var_type == 'I': return "OTHER_FOR_INDEL"  # Read has deletion where an Insertion is queried
                if expected_var_type == 'D':  # Querying a deletion
                    # The deletion in the read must match the expected deletion from reference
                    # This implies node_sequence[current_node_pos : current_node_pos + length]
                    # is the sequence that is deleted from the reference *according to this read*.
                    # And expected_ref_allele_for_indel is the specific deleted sequence we are checking for.
                    if node_sequence[current_node_pos: current_node_pos + length] == expected_ref_allele_for_indel:
                        return "*"  # Represents the deletion of expected_ref_allele_for_indel
                    else:
                        return "OTHER_FOR_INDEL"  # Read has a deletion, but not the one we are testing
                return "*"  # If not querying an INDEL specifically, a D at this position means the base is deleted
            current_node_pos += length
        elif op == 'S':  # Soft clip
            current_read_pos += length
        elif op == 'N':  # Skipped region from reference
            current_node_pos += length
        # Optimization: if we've passed the target position and it's not an insertion check
        if current_node_pos > target_node_pos + 1 and op in ('M', '=', 'X', 'D', 'N'):  # +1 for insertions adjacent
            if not (expected_var_type == 'I' and (
                    current_node_pos - 1) <= target_node_pos):  # allow one more step for I
                break
    return None  # Target position not covered or specific INDEL not found


def detect_variants_from_cigar(offset_on_node, cigar_ops_decoded, read_sequence, node_sequence):
    variants = []
    node_pos = offset_on_node  # Current position on the node sequence
    read_pos = 0  # Current position on the read sequence

    for length, op in cigar_ops_decoded:
        if op == 'M' or op == '=' or op == 'X':  # Match, Alignment match, or Mismatch
            for i in range(length):
                current_node_pos = node_pos + i
                current_read_pos = read_pos + i
                if current_node_pos < len(node_sequence) and current_read_pos < len(read_sequence):
                    node_base = node_sequence[current_node_pos].upper()
                    read_base = read_sequence[current_read_pos].upper()
                    if node_base != read_base and read_base != 'N' and node_base != 'N':  # Don't call N vs base as mismatch variant
                        variants.append((current_node_pos, 'X', read_base, node_base))
            node_pos += length
            read_pos += length
        elif op == 'I':  # Insertion to the reference
            inserted_sequence = read_sequence[read_pos: read_pos + length].upper()
            # Anchor position for insertion is the base *before* the insertion on the reference
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0  # Ensure not negative
            ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if 0 <= ref_anchor_pos < len(
                node_sequence) else "*"
            if inserted_sequence and 'N' not in inserted_sequence:  # Basic filter for N-only insertions
                variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor))
            read_pos += length
            # Node position does not advance for an insertion
        elif op == 'D':  # Deletion from the reference
            # Anchor position for deletion is the base *before* the deletion on the reference (same as SNP/Mismatch)
            # Or simply node_pos if we define deletion by its first deleted base
            deleted_sequence_from_ref = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= len(
                node_sequence) else ""
            if deleted_sequence_from_ref and 'N' not in deleted_sequence_from_ref:  # Basic filter
                variants.append((node_pos, 'D', "*", deleted_sequence_from_ref))  # ALT is '*', REF is deleted seq
            node_pos += length
            # Read position does not advance for a deletion
        elif op == 'S':  # Soft clip, consumes read bases
            read_pos += length
        elif op == 'N':  # Skipped region from reference, consumes reference bases
            node_pos += length
        # H (Hard clip) and P (Padding) do not consume bases in read or ref in terms of alignment coordinate progression
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
        # Optimization: if node position has passed window and it's a ref-consuming op
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
        # Condition to break early if the current ref position is beyond the window
        if current_node_pos_in_read >= window_start_node + current_tensor_window_size and cigar_op in ('M', 'D', 'N',
                                                                                                       '=', 'X'):
            break
        # Condition to break if read sequence is exhausted (can happen with trailing S or I ops beyond read length)
        # This check is now more targeted for read-consuming ops
        # if current_read_pos_in_read >= len(segment_read_sequence) and cigar_op in ('M','=','X','I','S'):
        #     break

        if cigar_op in ('M', '=', 'X'):
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                read_aln_pos = current_read_pos_in_read + i
                if read_aln_pos >= len(segment_read_sequence): break  # Exhausted read
                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_char = segment_read_sequence[read_aln_pos].upper()
                    base_indices_row[window_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])
                    if read_aln_pos < len(segment_quality_str):  # Qual string might be shorter if '*'
                        try:
                            quality_scores_row[window_idx] = ord(segment_quality_str[read_aln_pos]) - 33
                        except (TypeError, IndexError):  # Handles if qual is not string or index out of bounds
                            quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
                    else:
                        quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
            current_node_pos_in_read += cigar_len
            current_read_pos_in_read += cigar_len
        elif cigar_op == 'D' or cigar_op == 'N':  # Consumes reference
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_indices_row[window_idx] = BASE_TO_INDEX['*']  # Deletion marker
                    quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING  # No quality for deletion
            current_node_pos_in_read += cigar_len
        elif cigar_op == 'I' or cigar_op == 'S':  # Consumes read
            # Ensure we don't try to read past the end of the read sequence for quality, though it's not used for I/S in tensor directly.
            # This loop primarily advances current_read_pos_in_read.
            # The bases/qualities for insertions are not directly placed in these rows here.
            current_read_pos_in_read += cigar_len

        # Break if read is consumed, especially after I or S ops
        if current_read_pos_in_read >= len(segment_read_sequence) and cigar_op in ('M', '=', 'X', 'I', 'S'):
            break

    return base_indices_row, quality_scores_row


## ─────────────────────────────────────────────────────────────────────────────
## Worker Process Initialization and Target Function
## ─────────────────────────────────────────────────────────────────────────────

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker):
    global worker_dat_file, worker_base_output_dir
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path_for_worker}\n")
        sys.exit(1)  # Critical error for worker
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path_for_worker}: {e}\n")
        sys.exit(1)  # Critical error for worker


def process_single_node_for_pileup(task_args_with_af_thresh):
    node_id, dat_file_offset, n_records, node_sequence, min_af_threshold = task_args_with_af_thresh
    global worker_dat_file, worker_base_output_dir

    npy_files_generated_for_node = 0  # Initialize counter for .npy files

    if worker_dat_file is None or worker_base_output_dir is None:
        sys.stderr.write(
            f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Worker not initialized properly (dat_file or output_dir missing).\n")
        return node_id, None, npy_files_generated_for_node  # Return node_id, None for view_data, and 0 for npy_files

    if not node_sequence:  # Should be caught by main, but as a safeguard
        sys.stderr.write(
            f"ℹ️ [Worker {os.getpid()} for Node {node_id}]: No sequence provided. Skipping tensor/summary generation.\n")
        return node_id, {}, npy_files_generated_for_node  # Return empty dict for view_data, 0 for npy_files

    node_specific_output_dir = os.path.join(worker_base_output_dir, str(node_id))
    try:
        os.makedirs(node_specific_output_dir, exist_ok=True)
    except OSError as e:
        sys.stderr.write(
            f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Could not create directory {node_specific_output_dir}: {e}\n")
        return node_id, None, npy_files_generated_for_node

    node_len = len(node_sequence)
    variant_headers_for_node = []
    view_oriented_variant_data = {}  # For console display

    aligned_read_segments = []
    try:
        worker_dat_file.seek(dat_file_offset + 10)  # Skip the 10-byte header per node block in .dat
        for record_idx in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE:  # Should not happen if n_records is correct
                # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()} for Node {node_id}]: Unexpected EOF after {record_idx+1}/{n_records} records.\n")
                break
            off_from_file, raw_seq, raw_qual, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)

            if mapq < 10: continue  # Filter by MAPQ

            try:
                # Ensure null bytes are stripped BEFORE decoding
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')
                qual_str = raw_qual.rstrip(b'\x00').decode('ascii', errors='replace')  # Qualities are ASCII
                cigar_str_original = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()} for Node {node_id}]: Unicode decode error in record. Skipping.\n")
                continue  # Skip this problematic read

            if len(seq) == 0 or len(seq) != len(qual_str):  # Basic sanity check
                # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()} for Node {node_id}]: Empty seq or mismatched seq/qual lengths. Skipping.\n")
                continue

            original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
            if not original_decoded_cigar_ops and cigar_str_original != '*':  # If CIGAR is not '*' but fails to parse
                # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()} for Node {node_id}]: Failed to parse CIGAR '{cigar_str_original}'. Skipping read.\n")
                continue

            current_read_sequence = seq
            current_quality_str = qual_str
            current_decoded_cigar_ops = list(original_decoded_cigar_ops)  # Make a mutable copy
            current_offset_on_node = off_from_file

            if strand_char == '-':
                current_read_sequence = reverse_complement(seq)  # Reverse complement the sequence
                current_quality_str = qual_str[::-1]  # Reverse the quality string

                # Reverse CIGAR operations for negative strand
                temp_rev_cigar_ops = []
                for length, op_code in original_decoded_cigar_ops:
                    temp_rev_cigar_ops.insert(0, (length, op_code))  # Prepend to reverse
                current_decoded_cigar_ops = temp_rev_cigar_ops

                # Adjust offset for negative strand based on node length and alignment span on node
                # This logic was previously discussed and confirmed.
                alignment_span_on_node = len(
                    current_read_sequence)  # This might need adjustment based on CIGAR (M,D,N sum)
                # For simplicity, using read length for now as per original.
                # User confirmed fix: len(current_read_sequence)
                if alignment_span_on_node > 0:  # Ensure no div by zero if seq len is 0 (already checked)
                    current_offset_on_node = node_len - alignment_span_on_node - off_from_file
                    if current_offset_on_node < 0:
                        # This can happen if node_len is too small or off_from_file is large.
                        # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()} for Node {node_id}]: Negative offset calculated for '-' strand. Read offset: {off_from_file}, Read len: {len(seq)}, Node len: {node_len}. Calculated: {current_offset_on_node}. Skipping read.\n")
                        continue

            aligned_read_segments.append({
                "offset_on_node": current_offset_on_node,
                "read_sequence": current_read_sequence,
                "processed_quality_str": current_quality_str,  # Already reversed if needed
                "cigar_ops": current_decoded_cigar_ops,  # Already reversed if needed
                "original_cigar_str": cigar_str_original,  # Keep original for view
                "strand": strand_char
            })
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()} for Node {node_id}] reading records: {e}\n")
        return node_id, None, npy_files_generated_for_node

    candidate_variants = defaultdict(int)
    for segment in aligned_read_segments:
        variants_in_read = detect_variants_from_cigar(
            segment["offset_on_node"], segment["cigar_ops"],
            segment["read_sequence"], node_sequence
        )
        for v_pos, v_type, v_alt, v_ref in variants_in_read:
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1  # Note: v_ref, v_alt order for key

    variant_processing_window_size = TENSOR_WINDOW_SIZE
    half_window = variant_processing_window_size // 2

    for (v_pos, v_type, v_ref_defined, v_alt_defined), read_support_count in candidate_variants.items():
        af_alt_count, af_ref_count, af_other_count, af_locus_coverage = 0, 0, 0, 0

        # Define ref_allele_for_indel_af_check based on variant type being assessed
        ref_allele_for_indel_af_check = None
        if v_type == 'D':  # For Deletion, the ref allele is the sequence deleted from reference
            ref_allele_for_indel_af_check = v_ref_defined  # v_ref_defined is the deleted sequence
        elif v_type == 'I':  # For Insertion, the ref allele is the base *at the insertion site* on the reference
            # v_ref_defined is the base *before* insertion (anchor base)
            # For AF calculation, we need the base *at* v_pos (if it's 'X'-like comparison)
            # or simply that the read *doesn't* have the insertion.
            if 0 <= v_pos < node_len:  # v_pos is anchor for 'I'
                ref_allele_for_indel_af_check = node_sequence[v_pos]  # Base on ref at anchor site for 'I'
                # Not strictly 'ref allele' for insertion, but used by get_allele...

        for segment in aligned_read_segments:
            allele = get_allele_from_read_at_node_pos(
                segment["offset_on_node"], segment["read_sequence"], segment["cigar_ops"],
                v_pos, node_sequence,
                expected_var_type=v_type,
                expected_ref_allele_for_indel=ref_allele_for_indel_af_check if v_type in ['I', 'D'] else v_ref_defined
            )

            if allele is not None:
                af_locus_coverage += 1
                if v_type == 'X':
                    if allele == v_alt_defined:
                        af_alt_count += 1
                    elif allele == v_ref_defined:
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'I':  # v_alt_defined is the inserted sequence
                    if allele == v_alt_defined:
                        af_alt_count += 1
                    elif allele == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1  # Read doesn't have this specific insertion
                    else:
                        af_other_count += 1  # Read has a different insertion or other variant
                elif v_type == 'D':  # v_alt_defined is '*' meaning deletion
                    if allele == "*":
                        af_alt_count += 1  # Read confirms the deletion
                    elif allele == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1  # Read has the reference bases (no deletion)
                    else:
                        af_other_count += 1  # Read has something else (e.g. SNP at that site)

        alt_freq = af_alt_count / af_locus_coverage if af_locus_coverage > 0 else 0.0
        if alt_freq < min_af_threshold:
            continue  # Skip this variant if AF is too low

        # Construct a unique key string for the variant
        variant_key_str = f"{v_pos}_{v_type}_{v_ref_defined}_{v_alt_defined}"

        # Determine window center for view and tensor (they are aligned)
        if v_type == 'I':  # For insertions, center window *after* the anchor base
            window_center_on_node = v_pos + 1
        else:  # For SNPs and Deletions, center window *at* the variant position
            window_center_on_node = v_pos
        current_window_start_on_node = window_center_on_node - half_window

        # --- Prepare data for console view (pileup_reads_data_for_view) ---
        pileup_reads_data_for_view = []  # Limited to TENSOR_MAX_READ_ROWS for consistency if displayed
        reads_for_view_count = 0
        for segment in aligned_read_segments:
            if reads_for_view_count >= TENSOR_MAX_READ_ROWS: break  # Limit reads for view consistency
            row_chars_for_view = get_read_representation_in_window_for_view(
                segment["cigar_ops"], segment["offset_on_node"], segment["read_sequence"],
                current_window_start_on_node, variant_processing_window_size, node_len  # node_len for bounds
            )
            # Only add if there's some content in the window for this read
            if any(c != ' ' for c in row_chars_for_view):
                row_indices_for_view = [BASE_TO_INDEX.get(char.upper(), BASE_TO_INDEX['N']) for char in
                                        row_chars_for_view]
                pileup_reads_data_for_view.append({
                    "bases": row_indices_for_view, "offset": segment["offset_on_node"],
                    "strand": segment["strand"], "cigar": segment["original_cigar_str"]
                })
                reads_for_view_count += 1

        view_oriented_variant_data[variant_key_str] = {
            "pileup_reads_data": pileup_reads_data_for_view,  # Already limited
            "alt_allele_count": af_alt_count,
            "ref_allele_count_at_locus": af_ref_count,
            "other_allele_count_at_locus": af_other_count,
            "coverage_at_locus": af_locus_coverage,
            "alt_allele_frequency": round(alt_freq, 4)
        }

        # --- Prepare data for tensor ---
        tensor_ch1_bases_list = []  # Channel 1: Base indices
        tensor_ch2_qualities_list = []  # Channel 2: Quality scores
        tensor_ch3_mismatches_list = []  # Channel 3: Mismatch to reference indicator

        # Reference Row for Tensor (Row 0 in each channel)
        ref_base_indices_row_tensor = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        ref_qual_scores_row_tensor = [DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE  # Quality for ref row is padding
        ref_mismatch_row_tensor = [MISMATCH_CHANNEL_REF_ROW_VALUE] * TENSOR_WINDOW_SIZE  # Ref row is 0 (no mismatch with itself)

        for i in range(TENSOR_WINDOW_SIZE):
            actual_node_pos = current_window_start_on_node + i
            if 0 <= actual_node_pos < node_len:
                ref_base_indices_row_tensor[i] = BASE_TO_INDEX.get(node_sequence[actual_node_pos].upper(),
                                                                   BASE_TO_INDEX['N'])
        tensor_ch1_bases_list.append(ref_base_indices_row_tensor)
        tensor_ch2_qualities_list.append(ref_qual_scores_row_tensor)
        tensor_ch3_mismatches_list.append(ref_mismatch_row_tensor)

        reads_added_to_tensor = 0
        for segment in aligned_read_segments:  # Iterate through all aligned reads again for tensor
            if reads_added_to_tensor >= TENSOR_MAX_READ_ROWS:
                break  # Stop if max read rows for tensor is reached

            base_indices_row_tensor, quality_scores_row_tensor = get_read_tensor_rows_in_window(
                segment["cigar_ops"], segment["offset_on_node"],
                segment["read_sequence"], segment["processed_quality_str"],
                current_window_start_on_node, TENSOR_WINDOW_SIZE, node_len
            )

            # Only add if there's some actual base information (not all padding)
            # This check might be too strict if a read legitimately aligns outside and only shows padding
            # For now, keeping it to ensure some content, or else it's just a padding row.
            if any(b != PADDING_BASE_INDEX for b in base_indices_row_tensor):
                tensor_ch1_bases_list.append(base_indices_row_tensor)
                tensor_ch2_qualities_list.append(quality_scores_row_tensor)

                # Mismatch Channel (Channel 3) for this read
                mismatch_row_for_read_tensor = [MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE
                for i in range(TENSOR_WINDOW_SIZE):
                    read_base_idx = base_indices_row_tensor[i]
                    ref_base_idx_for_comp = ref_base_indices_row_tensor[i]  # From ref row calculated above

                    if read_base_idx == PADDING_BASE_INDEX or ref_base_idx_for_comp == PADDING_BASE_INDEX:
                        mismatch_row_for_read_tensor[i] = MISMATCH_COMPARISON_PADDING_VALUE
                    elif read_base_idx == ref_base_idx_for_comp:
                        mismatch_row_for_read_tensor[i] = 0  # Match
                    else:
                        mismatch_row_for_read_tensor[i] = 1  # Mismatch
                tensor_ch3_mismatches_list.append(mismatch_row_for_read_tensor)
                reads_added_to_tensor += 1

        # Pad tensor if fewer than TENSOR_MAX_READ_ROWS were processed
        for _ in range(TENSOR_MAX_READ_ROWS - reads_added_to_tensor):
            tensor_ch1_bases_list.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            tensor_ch2_qualities_list.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            tensor_ch3_mismatches_list.append([MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE)

        try:
            # Create PyTorch tensor
            final_tensor = torch.tensor([
                tensor_ch1_bases_list,
                tensor_ch2_qualities_list,
                tensor_ch3_mismatches_list
            ], dtype=torch.int8)

            # Convert to NumPy array
            numpy_array = final_tensor.numpy()

            # Save as .npy file
            tensor_filename = f"{variant_key_str}.npy"  # Changed extension
            tensor_filepath = os.path.join(node_specific_output_dir, tensor_filename)
            np.save(tensor_filepath, numpy_array)  # Use np.save

            variant_headers_for_node.append({
                "variant_key": variant_key_str,
                "tensor_file": tensor_filename,  # Filename now has .npy
                "alt_allele_count": af_alt_count,
                "ref_allele_count_at_locus": af_ref_count,
                "other_allele_count_at_locus": af_other_count,
                "coverage_at_locus": af_locus_coverage,
                "alt_allele_frequency": round(alt_freq, 4)
            })
        except Exception as e:
            sys.stderr.write(
                f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Failed to create or save tensor for {variant_key_str}: {e}\n")

    npy_files_generated_for_node = len(variant_headers_for_node)  # Count of .npy files
    if variant_headers_for_node:  # Only create summary if there are variants passing AF
        summary_json_path = os.path.join(node_specific_output_dir, "variant_summary.json")
        try:
            with open(summary_json_path, 'w') as sjf:
                json.dump({
                    "node_id": node_id,
                    "node_length": node_len,
                    "node_sequence_preview": node_sequence[:100] + "..." if node_len > 100 else node_sequence,
                    # Keep preview
                    "variants": variant_headers_for_node  # List of dicts, each with "tensor_file" key
                }, sjf, indent=2)
        except Exception as e:
            sys.stderr.write(f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Failed to write summary JSON: {e}\n")
            # If summary fails, .npy files might still exist.

    return node_id, view_oriented_variant_data, npy_files_generated_for_node


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

    if not node_data_for_display_view:  # Empty dictionary
        print(f"ℹ️ No variants found or pileups generated for this node (or all filtered by AF).")
        return

    variants_displayed_count = 0
    # Sort by position then by type (e.g., 10_D_..., 10_I_..., 10_X_...)
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
        pileup_reads_display_data = variant_data.get("pileup_reads_data",
                                                     [])  # This is already limited by TENSOR_MAX_READ_ROWS

        parts = variant_key.split('_')  # pos_type_ref_alt
        v_pos = int(parts[0])
        v_type = parts[1]
        # v_ref = parts[2]
        # v_alt = parts[3]

        # Determine window center for display (consistent with tensor generation)
        if v_type == 'I':
            window_center_on_node_display = v_pos + 1
        else:
            window_center_on_node_display = v_pos
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

        # Display Reference sequence and Marker line
        ref_display_parts = [' '] * display_window_size
        marker_line_parts = [' '] * display_window_size
        variant_display_idx_in_window = -1  # Relative to window start

        # For 'I', the v_pos is the anchor base *before* the insertion.
        # The marker points between v_pos and v_pos+1 (or at v_pos if it's the first base of insertion symbol like 'I^')
        # For 'X' and 'D', v_pos is the affected base.
        if v_type != 'I':  # SNP or Deletion
            if current_window_start_on_node_display <= v_pos < current_window_start_on_node_display + display_window_size:
                variant_display_idx_in_window = v_pos - current_window_start_on_node_display
        elif v_type == 'I':  # Insertion, v_pos is anchor
            # The insertion happens *after* v_pos. So we mark between v_pos and v_pos+1
            # If v_pos is in window, mark it.
            if current_window_start_on_node_display <= v_pos < current_window_start_on_node_display + display_window_size:
                variant_display_idx_in_window = v_pos - current_window_start_on_node_display

        for i in range(display_window_size):
            actual_node_pos_in_window = current_window_start_on_node_display + i
            if 0 <= actual_node_pos_in_window < len(full_node_sequence):
                ref_display_parts[i] = full_node_sequence[actual_node_pos_in_window]

            if i == variant_display_idx_in_window:
                if v_type == 'I':
                    marker_line_parts[i] = "I"  # Mark at the anchor base
                    if i + 1 < display_window_size:  # If there's space for '^'
                        marker_line_parts[i + 1] = "^"
                    elif i == display_window_size - 1:  # Insertion at the very end of window view
                        marker_line_parts[i] = ">"  # Indicate insertion continues
                else:  # SNP or Deletion
                    marker_line_parts[i] = "^"

        print(f"  Node Ref: {''.join(ref_display_parts)}")
        print(f"  Marker  : {''.join(marker_line_parts)}")

        if not pileup_reads_display_data:
            print("  (No reads data in pileup for this variant's window or to display)")
        else:
            displayed_reads_count = 0
            for i, read_info in enumerate(pileup_reads_display_data):  # pileup_reads_display_data is already limited
                if displayed_reads_count >= max_reads_to_display_per_variant:  # Apply further limit if arg is smaller
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
        description="Generate variant-centered pileups (.npy tensors and JSON summaries) for specified node(s) and optionally view them.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter  # Helps with showing defaults in --help
    )
    parser.add_argument("dat", help=".dat file path (read alignment data)")
    parser.add_argument("idx", help=".idx file path (index for .dat file)")
    parser.add_argument("output", help="Base output directory for node-specific folders, tensors, and summaries.")

    # Group for mutually exclusive node ID input
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
                        # -1 for all, N for N variants
                        help="Print generated pileups to console. Optionally specify N for first N variants per node. Use -1 for all variants. Default: No view.")
    parser.add_argument("--max_view_reads", type=int, default=20,
                        help="Max reads to display per pileup in console view.")
    parser.add_argument("--min_af", type=float, default=0.1,
                        help="Min allele frequency for a variant to be processed and included in output.")
    # Removed --num_workers as it's fixed to 1 per node in this script structure
    args = parser.parse_args()

    # --- Argument Validation ---
    if not os.path.isfile(args.dat): sys.exit(f"❌ Error: DAT file not found: {args.dat}")
    if not os.path.isfile(args.idx): sys.exit(f"❌ Error: Index file not found: {args.idx}")
    if not args.load_cache and not args.gfa:
        sys.exit("❌ Error: Must provide --gfa or --load-cache to obtain node sequences.")
    if args.load_cache and not os.path.isfile(args.load_cache):  # if specified, it must exist
        sys.exit(f"❌ Error: Specified --load-cache file not found: {args.load_cache}")
    if args.gfa and not os.path.isfile(args.gfa):  # if specified, it must exist
        sys.exit(f"❌ Error: GFA file not found: {args.gfa}")
    if not (0.0 <= args.min_af <= 1.0):
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
                    if line and not line.startswith("#"):  # Skip empty lines and comments
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
    elif args.node_id is not None:  # This elif is correct due to mutually_exclusive_group
        target_node_ids.append(args.node_id)
        print(f"🔹 Will process single target node ID: {args.node_id}")
    # No else needed as one of them is required

    # --- Load Initial Sequence Cache ---
    node_sequences_map = {}  # Holds node_id_str -> sequence
    if args.load_cache:  # args.load_cache is path to existing cache
        start_time_cache_load = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                node_sequences_map = json.load(cf)  # Loads into {node_id_str: sequence}
            print(
                f"✔ Loaded {len(node_sequences_map)} sequences from cache '{args.load_cache}' in {time.time() - start_time_cache_load:.2f}s.")
        except json.JSONDecodeError as e:
            print(
                f"⚠️ Warning: Error decoding JSON from cache file {args.load_cache}: {e}. Starting with an empty sequence map for loaded cache.",
                file=sys.stderr)
            node_sequences_map = {}  # Ensure it's an empty dict on error
        except Exception as e:  # Catch other potential errors like permission issues
            print(
                f"⚠️ Warning: Error loading cache {args.load_cache}: {e}. Starting with an empty sequence map for loaded cache.",
                file=sys.stderr)
            node_sequences_map = {}

    # --- Main Processing Loop for Each Node ID ---
    processed_nodes_count = 0  # Nodes that worker was submitted for and returned (even if error in worker)
    successful_nodes_with_output = 0  # Nodes that had a summary.json written
    overall_start_time = time.time()
    total_npy_files_generated = 0  # Grand total of .npy files
    nodes_iterated_count = 0  # Nodes attempted from the input list

    for target_node_id in target_node_ids:
        nodes_iterated_count += 1
        print(
            f"\n═════════════════════════════ PROCESSING NODE ID: {target_node_id} ({nodes_iterated_count}/{len(target_node_ids)}) ═════════════════════════════")
        node_specific_start_time = time.time()
        node_sequence = None

        # 1. Get DAT file offset and record count for the node
        idx_parse_start_time = time.time()
        node_dat_info = parse_idx_file_for_single_node(args.idx, target_node_id)
        if not node_dat_info:
            print(f"❌ Skipping node {target_node_id} due to index parsing failure.", file=sys.stderr)
            if args.node_id_file and nodes_iterated_count > 0 and nodes_iterated_count % 100 == 0:
                print(
                    f"ℹ️ Progress after iterating {nodes_iterated_count} node(s) from file: {total_npy_files_generated} .npy files generated so far.")
            continue  # Move to the next node_id in the list
        dat_offset, n_records = node_dat_info
        print(
            f"✔ Index parsing for node {target_node_id} took {time.time() - idx_parse_start_time:.2f}s. (Offset: {dat_offset}, Records: {n_records})")

        # 2. Obtain Node Sequence (Cache -> GFA)
        node_id_str = str(target_node_id)
        if node_id_str in node_sequences_map:
            node_sequence = node_sequences_map[node_id_str]
            print(f"✔ Sequence for node {target_node_id} found in pre-loaded/memory cache.")
        elif args.gfa:  # If not in cache, try GFA
            print(f"ℹ️ Sequence for node {target_node_id} not in cache, attempting to load from GFA: {args.gfa}")
            gfa_load_start_time = time.time()
            node_sequence = load_node_sequence_from_gfa(args.gfa, target_node_id)
            if node_sequence:
                print(
                    f"✔ GFA sequence loading for node {target_node_id} took {time.time() - gfa_load_start_time:.2f}s.")
                node_sequences_map[node_id_str] = node_sequence  # Add to in-memory map for potential save-cache later
            else:
                print(f"❌ Skipping node {target_node_id}: Failed to load sequence from GFA.", file=sys.stderr)
                if args.node_id_file and nodes_iterated_count > 0 and nodes_iterated_count % 100 == 0:
                    print(
                        f"ℹ️ Progress after iterating {nodes_iterated_count} node(s) from file: {total_npy_files_generated} .npy files generated so far.")
                continue
        else:  # Not in cache AND no GFA provided
            print(f"❌ Skipping node {target_node_id}: Sequence not in cache and --gfa not provided.", file=sys.stderr)
            if args.node_id_file and nodes_iterated_count > 0 and nodes_iterated_count % 100 == 0:
                print(
                    f"ℹ️ Progress after iterating {nodes_iterated_count} node(s) from file: {total_npy_files_generated} .npy files generated so far.")
            continue

        if not node_sequence:  # Should be redundant due to continues above, but as a final check
            print(f"❌ Skipping node {target_node_id}: Sequence could not be obtained.", file=sys.stderr)
            if args.node_id_file and nodes_iterated_count > 0 and nodes_iterated_count % 100 == 0:
                print(
                    f"ℹ️ Progress after iterating {nodes_iterated_count} node(s) from file: {total_npy_files_generated} .npy files generated so far.")
            continue

        # 3. Prepare and Submit Task to Worker
        task_payload = (target_node_id, dat_offset, n_records, node_sequence, args.min_af)
        print(
            f"🔹 Preparing task for node {target_node_id} (len: {len(node_sequence)}bp, {n_records} records, min_AF: {args.min_af}).")

        node_processing_start_time = time.time()
        processed_node_id_result = None  # To check if worker returned node_id
        view_data_for_node = None  # To hold dict for viewing
        npy_files_for_this_node = 0  # Files generated by this worker call

        try:
            # Using a new ProcessPoolExecutor for each node ensures clean worker state if a worker crashes,
            # but adds overhead. For single worker use, it's okay.
            # If you were to use multiple workers for multiple nodes simultaneously, a single pool would be better.
            with ProcessPoolExecutor(max_workers=1, initializer=init_worker,
                                     initargs=(args.dat, args.output)) as executor:
                future = executor.submit(process_single_node_for_pileup, task_payload)
                # Unpack results: node_id, view_data_dict, count_of_npy_files
                processed_node_id_result, view_data_for_node, npy_files_for_this_node = future.result()

            total_npy_files_generated += npy_files_for_this_node  # Add to grand total

            if processed_node_id_result is None:  # Should ideally not happen if future.result() succeeded without exception
                print(
                    f"❌ Critical Error: Worker failed to return results properly for node {target_node_id}. Check worker logs.",
                    file=sys.stderr)
                if args.node_id_file and nodes_iterated_count > 0 and nodes_iterated_count % 100 == 0:
                    print(
                        f"ℹ️ Progress after iterating {nodes_iterated_count} node(s) from file: {total_npy_files_generated} .npy files generated so far.")
                continue  # Skip this node for display and counting as "processed"

            # view_data_for_node can be {} if no variants passed AF, or None if worker had internal error before populating it
            if view_data_for_node is None:
                print(
                    f"⚠️ Warning: Processing for node {target_node_id} returned no view data (may indicate a worker-side issue post-start). .npy files for this node: {npy_files_for_this_node}",
                    file=sys.stderr)

        except Exception as e:  # Covers errors in ProcessPoolExecutor, task submission, or future.result()
            print(f"❌ An error occurred launching or during processing for node {target_node_id}: {e}", file=sys.stderr)
            if args.node_id_file and nodes_iterated_count > 0 and nodes_iterated_count % 100 == 0:
                print(
                    f"ℹ️ Progress after iterating {nodes_iterated_count} node(s) from file: {total_npy_files_generated} .npy files generated so far.")
            continue  # Skip to next node

        print(
            f"✔ Worker task for node {target_node_id} finished in {time.time() - node_processing_start_time:.2f}s. Generated {npy_files_for_this_node} .npy file(s).")
        processed_nodes_count += 1  # Count nodes for which worker processing was attempted and returned a result (even if problematic)

        # 4. Display Pileup Data (Optional)
        if args.view is not None:
            if view_data_for_node and isinstance(view_data_for_node,
                                                 dict) and view_data_for_node:  # If dict is not empty
                max_v_to_show = float('inf') if args.view == -1 else args.view  # -1 means all
                if args.view < -1:  # Invalid negative number other than -1
                    print(
                        f"⚠️ Warning: Invalid number for --view ({args.view}). Showing all variants for node {target_node_id}.",
                        file=sys.stderr)
                    max_v_to_show = float('inf')

                view_limit_msg = f"first {int(max_v_to_show)}" if max_v_to_show != float('inf') else "all"
                print(
                    f"🔹 Displaying {view_limit_msg} pileups for node {target_node_id} (max {args.max_view_reads} reads/variant)...")
                display_pileup_data(view_data_for_node, node_id_str, node_sequence, args.max_view_reads, max_v_to_show)
            elif view_data_for_node is not None:  # view_data_for_node is {} (empty dict means no variants passed AF)
                print(
                    f"ℹ️ --view specified for node {target_node_id}, but no variants met AF threshold or were found to display.")
            # If view_data_for_node is None, a warning was already printed.

        # 5. Report Output Generation
        node_summary_file_path = os.path.join(args.output, node_id_str, "variant_summary.json")
        if os.path.exists(node_summary_file_path):
            print(f"✅ Output generated for node {target_node_id} in: {os.path.dirname(node_summary_file_path)}")
            successful_nodes_with_output += 1
        elif view_data_for_node is not None and not view_data_for_node and npy_files_for_this_node == 0:  # No variants passed AF
            print(f"ℹ️ No variants met AF threshold for node {target_node_id}. No tensor or summary files generated.")
        else:  # Summary doesn't exist, but npy_files might, or view_data was None
            print(
                f"ℹ️ No output summary file found for node {target_node_id} at {node_summary_file_path}. (Generated {npy_files_for_this_node} .npy files for this node).")

        print(f"✔ Completed processing for node {target_node_id} in {time.time() - node_specific_start_time:.2f}s.")

        # Progress report for long lists of nodes
        if args.node_id_file and nodes_iterated_count > 0 and nodes_iterated_count % 100 == 0:
            print(
                f"\nℹ️ Progress after iterating {nodes_iterated_count} node(s) from file: {total_npy_files_generated} .npy files generated so far.\n")

    print("\n═════════════════════════════ PROCESSING COMPLETE ═════════════════════════════")

    # --- Save/Update Sequence Cache (if specified) ---
    if args.save_cache:  # args.save_cache is path to save to
        if node_sequences_map:  # Only save if there's something in the map (newly fetched or originally loaded)
            print(f"\n🔹 Saving/Updating {len(node_sequences_map)} sequence(s) to cache file: {args.save_cache}...")
            try:
                # Prune empty sequences before saving if any crept in (should not happen with current logic)
                # final_sequences_to_save = {k: v for k, v in node_sequences_map.items() if v}
                with open(args.save_cache, 'w') as wcf:
                    json.dump(node_sequences_map, wcf, indent=2)  # Save the potentially updated map
                print(f"✔ Sequences saved to cache '{args.save_cache}'.")
            except Exception as e:
                print(f"❌ Error saving sequences to cache {args.save_cache}: {e}", file=sys.stderr)
        elif os.path.exists(args.save_cache):  # If map is empty but file exists (e.g. loaded empty, fetched none)
            print(
                f"ℹ️ --save-cache specified ({args.save_cache}), but no sequences were loaded or fetched into memory to update it.")
            # Optionally, clear the cache file if it was meant to be a fresh build and nothing was found
            if args.load_cache == args.save_cache and not node_sequences_map:  # If loading and saving to same file and map ended up empty
                try:
                    with open(args.save_cache, 'w') as wcf:
                        json.dump({}, wcf, indent=2)
                    print(f"✔ Cache file {args.save_cache} effectively emptied as no sequences remain in map.")
                except Exception as e:
                    print(f"❌ Error writing empty cache {args.save_cache}: {e}", file=sys.stderr)
        else:  # Map is empty and no pre-existing save_cache file to worry about overwriting implicitly
            print(
                f"ℹ️ --save-cache specified, but no sequences were loaded or fetched. No cache file created/updated as map is empty.")

    # --- Final Summary ---
    total_script_time = time.time() - overall_start_time
    print(
        f"\nProcessed {processed_nodes_count}/{len(target_node_ids)} node IDs from input list that were submitted to worker.")
    if successful_nodes_with_output > 0:
        print(f"✅ Output summary files successfully generated for {successful_nodes_with_output} node(s).")
    # If nodes_iterated_count is greater than processed_nodes_count, it implies some nodes were skipped before worker submission
    if processed_nodes_count < nodes_iterated_count:
        print(
            f"⚠️ {nodes_iterated_count - processed_nodes_count} node(s) were skipped before worker submission (e.g., due to index parsing or sequence loading errors).")

    print(f"ℹ️ A grand total of {total_npy_files_generated} .npy files were generated.")
    print(f"🏁 Script finished in {total_script_time:.2f} seconds.")


if __name__ == '__main__':
    main()