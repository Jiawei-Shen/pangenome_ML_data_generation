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
# Removed ProcessPoolExecutor as per request for no parallel run in this revision
import torch  # Still needed for tensor creation before converting to NumPy

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
MISMATCH_COMPARISON_PADDING_VALUE = -1  # Using -1 for padding in mismatch channel for clarity


# Globals for worker-like function state (if any were needed, but DAT file will be opened per node)
# For serial processing, these might not be 'globals' in the same way.
# worker_dat_file = None # Will be opened locally in process_node_serially
# worker_base_output_dir = None # Will be passed as arg

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions
# ─────────────────────────────────────────────────────────────────────────────

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]


def load_full_idx_data(idx_path):
    """Loads all node offset information from the .idx file into a dictionary."""
    idx_data_map = {}  # node_id (int) -> (offset, n_records)
    print(f"🔹 Loading full index data from {idx_path}...")
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4:
                sys.stderr.write(f"❌ Error: Index file {idx_path} is too small (size: {file_size} bytes).\n")
                return None
            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4:
                sys.stderr.write(f"❌ Error: Could not read number of nodes from {idx_path}.\n")
                return None
            num_nodes_in_idx = struct.unpack('<I', num_nodes_bytes)[0]
            print(f"  Index file reports {num_nodes_in_idx} total node entries. Reading all entries...")

            if num_nodes_in_idx == 0:
                print(f"  Index file contains 0 node entries according to its header.")
                return idx_data_map

            expected_min_size = 4 + (num_nodes_in_idx * 22)
            if file_size < expected_min_size:
                sys.stderr.write(
                    f"⚠️ Warning: Index file size ({file_size} bytes) is smaller than expected ({expected_min_size} bytes) for {num_nodes_in_idx} records. File may be truncated.\n")

            processed_entries = 0
            for i in range(num_nodes_in_idx):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                    sys.stderr.write(
                        f"❌ Error: Index file ended prematurely while reading record {i + 1}/{num_nodes_in_idx}. Loaded {processed_entries} entries.\n")
                    break
                node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                idx_data_map[node_id_from_idx] = (offset, n_records)
                processed_entries += 1
                if processed_entries > 0 and processed_entries % 2000000 == 0:
                    print(f"    Loaded {processed_entries}/{num_nodes_in_idx} index entries...")

            print(f"✔ Successfully loaded {len(idx_data_map)} distinct node entries from index file {idx_path}.")
        return idx_data_map
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: Index file not found at {idx_path}\n")
        return None
    except Exception as e:
        sys.stderr.write(f"❌ Error parsing full index file {idx_path}: {e}\n")
        return None


def load_multiple_node_sequences_from_gfa(gfa_path, target_node_ids_set):
    """Reads GFA once for sequences of target_node_ids_set (integers). Returns dict: node_id_str -> sequence."""
    node_sequences = {}  # node_id (str) -> sequence
    if not target_node_ids_set:
        return node_sequences

    nodes_to_find = target_node_ids_set.copy()  # Work on a copy to remove found IDs

    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file to find sequences for {len(nodes_to_find)} target node(s): {gfa_path}")
            line_counter = 0
            for line in f:
                line_counter += 1
                if line_counter % 10_000_000 == 0:  # Progress for very large GFA files
                    print(
                        f"  Checked {line_counter:,} lines in GFA file... {len(nodes_to_find)} nodes remaining to find.")

                if not line.startswith('S\t'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 3:  # Need at least S, name, sequence
                    continue

                try:
                    nid_int_from_gfa = int(parts[1])  # GFA node IDs are often strings, but script uses ints
                except ValueError:
                    continue  # Skip if node ID in GFA is not purely numeric if our targets are

                if nid_int_from_gfa in nodes_to_find:
                    node_sequences[str(nid_int_from_gfa)] = parts[2]  # Store with string key
                    nodes_to_find.remove(nid_int_from_gfa)
                    if not nodes_to_find:  # All target nodes found
                        print(
                            f"✔ Found all {len(target_node_ids_set)} requested node sequences in GFA after checking {line_counter:,} lines.")
                        break  # Stop reading GFA

            found_count = len(node_sequences)
            requested_count = len(target_node_ids_set)
            if not nodes_to_find:  # All found, message already printed
                pass
            else:  # Some not found
                print(
                    f"✔ Finished GFA scan ({line_counter:,} lines). Found {found_count}/{requested_count} of the targeted sequences.")
                if nodes_to_find:  # Log which ones were missing if any
                    print(
                        f"⚠️ Warning: Could not find GFA sequences for {len(nodes_to_find)} node ID(s). Examples: {list(nodes_to_find)[:5]}")

    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: GFA file not found at {gfa_path}\n")
        return {}  # Return empty dict if critical error
    except Exception as e:
        sys.stderr.write(f"❌ Error reading GFA file {gfa_path}: {e}\n")
        return node_sequences  # Return whatever was found so far

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
        # sys.stderr.write(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}\n") # Can be noisy
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
                    # For querying a deletion, the 'ref' is what's deleted from the reference genome.
                    # The read shows a deletion ('*') if its CIGAR aligns it over these ref bases with a 'D' op.
                    # The expected_ref_allele_for_indel should be the actual sequence deleted from the reference.
                    ref_deleted_segment = node_sequence[current_node_pos: current_node_pos + length]
                    if ref_deleted_segment == expected_ref_allele_for_indel:
                        return "*"  # This read confirms the specific deletion we're looking for
                    else:
                        return "OTHER_FOR_INDEL"  # Read has a deletion, but of different bases
                return "*"  # If just checking a position that happens to be deleted in this read
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
    node_seq_len = len(node_sequence)
    read_seq_len = len(read_sequence)

    for length, op in cigar_ops_decoded:
        if op == 'M' or op == '=' or op == 'X':
            for i in range(length):
                current_node_pos = node_pos + i
                current_read_pos = read_pos + i
                if current_node_pos < node_seq_len and current_read_pos < read_seq_len:
                    node_base = node_sequence[current_node_pos].upper()
                    read_base = read_sequence[current_read_pos].upper()
                    # Consider Ns as non-variant inducing or handle as per specific needs
                    if node_base != read_base and read_base != 'N' and node_base != 'N':
                        variants.append((current_node_pos, 'X', read_base, node_base))
                else:
                    break
            node_pos += length
            read_pos += length
        elif op == 'I':
            inserted_sequence = read_sequence[read_pos: read_pos + length].upper()
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0
            ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if 0 <= ref_anchor_pos < node_seq_len else "*"
            if inserted_sequence and 'N' not in inserted_sequence:  # Basic filter for N-only insertions
                variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor))
            read_pos += length
        elif op == 'D':
            deleted_sequence_from_ref = node_sequence[
                                        node_pos: node_pos + length].upper() if node_pos + length <= node_seq_len else ""
            if deleted_sequence_from_ref and 'N' not in deleted_sequence_from_ref:  # Basic filter for N-only deletions in ref
                variants.append((node_pos, 'D', "*", deleted_sequence_from_ref))
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
        if current_node_pos_in_read >= window_start_node + current_tensor_window_size and cigar_op in ('M', 'D', 'N',
                                                                                                       '=', 'X'):
            break

        if cigar_op in ('M', '=', 'X'):
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                read_aln_pos = current_read_pos_in_read + i
                if read_aln_pos >= len(segment_read_sequence): break
                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_char = segment_read_sequence[read_aln_pos].upper()
                    base_indices_row[window_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])
                    if read_aln_pos < len(segment_quality_str):
                        try:
                            quality_scores_row[window_idx] = ord(segment_quality_str[read_aln_pos]) - 33
                        except:
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
                    quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
            current_node_pos_in_read += cigar_len
        elif cigar_op == 'I' or cigar_op == 'S':
            current_read_pos_in_read += cigar_len

        if current_read_pos_in_read >= len(segment_read_sequence) and cigar_op in ('M', '=', 'X', 'I', 'S'):
            break
    return base_indices_row, quality_scores_row


## ─────────────────────────────────────────────────────────────────────────────
## Core Processing Logic (formerly worker function)
## ─────────────────────────────────────────────────────────────────────────────
def process_node_serially(dat_file_path, base_output_dir,
                          node_id, dat_file_offset, n_records, node_sequence, min_af_threshold):
    """
    Processes a single node: reads DAT, calls variants, generates tensors (.npy) and summary.
    This function is called directly in a serial loop.
    """
    npy_files_generated_for_node = 0
    view_oriented_variant_data = {}  # For console display

    node_specific_output_dir = os.path.join(base_output_dir, str(node_id))
    try:
        os.makedirs(node_specific_output_dir, exist_ok=True)
    except OSError as e:
        sys.stderr.write(f"❌ Error creating directory {node_specific_output_dir} for Node {node_id}: {e}\n")
        return node_id, None, npy_files_generated_for_node  # Return None for view_data to signal error

    if not node_sequence:
        sys.stderr.write(f"ℹ️ Node {node_id}: No sequence provided. Skipping tensor/summary generation.\n")
        return node_id, {}, npy_files_generated_for_node

    node_len = len(node_sequence)
    aligned_read_segments = []

    try:
        with open(dat_file_path, 'rb') as dat_f:  # Open DAT file for this node's processing
            dat_f.seek(dat_file_offset + 10)  # Skip 10-byte header for this node's block
            for record_idx in range(n_records):
                data = dat_f.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE:
                    sys.stderr.write(f"⚠️ Node {node_id}: Unexpected EOF after {record_idx + 1}/{n_records} records.\n")
                    break
                off_from_file, raw_seq, raw_qual, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)

                if mapq < 10: continue

                try:
                    seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')
                    qual_str = raw_qual.rstrip(b'\x00').decode('ascii', errors='replace')
                    cigar_str_original = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                    strand_char = strand_byte.decode('ascii')
                except UnicodeDecodeError:
                    continue

                if len(seq) == 0 or len(seq) != len(qual_str): continue

                original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
                if not original_decoded_cigar_ops and cigar_str_original != '*': continue

                current_read_sequence, current_quality_str = seq, qual_str
                current_decoded_cigar_ops = list(original_decoded_cigar_ops)
                current_offset_on_node = off_from_file

                if strand_char == '-':
                    current_read_sequence = reverse_complement(seq)
                    current_quality_str = qual_str[::-1]
                    current_decoded_cigar_ops = [op for op in reversed(
                        original_decoded_cigar_ops)] if original_decoded_cigar_ops else []

                    alignment_span_on_node = len(current_read_sequence)  # User-confirmed logic for offset
                    if alignment_span_on_node > 0:
                        current_offset_on_node = node_len - alignment_span_on_node - off_from_file
                        if current_offset_on_node < 0: continue

                aligned_read_segments.append({
                    "offset_on_node": current_offset_on_node, "read_sequence": current_read_sequence,
                    "processed_quality_str": current_quality_str, "cigar_ops": current_decoded_cigar_ops,
                    "original_cigar_str": cigar_str_original, "strand": strand_char
                })
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: DAT file {dat_file_path} not found during processing for node {node_id}.\n")
        return node_id, None, npy_files_generated_for_node
    except Exception as e:
        sys.stderr.write(f"❌ Error reading records for node {node_id} from {dat_file_path}: {e}\n")
        return node_id, None, npy_files_generated_for_node

    # --- Variant Calling and AF Calculation ---
    candidate_variants = defaultdict(int)
    for segment in aligned_read_segments:
        variants_in_read = detect_variants_from_cigar(
            segment["offset_on_node"], segment["cigar_ops"],
            segment["read_sequence"], node_sequence)
        for v_pos, v_type, v_alt, v_ref in variants_in_read:
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    variant_headers_for_node = []  # To store info for summary.json
    variant_processing_window_size = TENSOR_WINDOW_SIZE
    half_window = variant_processing_window_size // 2

    for (v_pos, v_type, v_ref_defined, v_alt_defined), read_support_count in candidate_variants.items():
        af_alt_count, af_ref_count, af_other_count, af_locus_coverage = 0, 0, 0, 0
        ref_allele_for_indel_af_check = None
        if v_type == 'D':
            ref_allele_for_indel_af_check = v_ref_defined
        elif v_type == 'I':
            if 0 <= v_pos < node_len: ref_allele_for_indel_af_check = node_sequence[v_pos]

        for segment in aligned_read_segments:
            allele = get_allele_from_read_at_node_pos(
                segment["offset_on_node"], segment["read_sequence"], segment["cigar_ops"],
                v_pos, node_sequence, expected_var_type=v_type,
                expected_ref_allele_for_indel=ref_allele_for_indel_af_check if v_type in ['I', 'D'] else v_ref_defined)
            if allele is not None:
                af_locus_coverage += 1
                if v_type == 'X':
                    if allele == v_alt_defined:
                        af_alt_count += 1
                    elif allele == v_ref_defined:
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'I':
                    if allele == v_alt_defined:
                        af_alt_count += 1
                    elif allele == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'D':
                    if allele == "*":
                        af_alt_count += 1
                    elif allele == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1
                    else:
                        af_other_count += 1

        alt_freq = af_alt_count / af_locus_coverage if af_locus_coverage > 0 else 0.0
        if alt_freq < min_af_threshold: continue

        variant_key_str = f"{v_pos}_{v_type}_{v_ref_defined}_{v_alt_defined}"
        window_center_on_node = v_pos + 1 if v_type == 'I' else v_pos
        current_window_start_on_node = max(0,
                                           window_center_on_node - half_window)  # Ensure window start is not negative

        # --- Prepare data for console view (pileup_reads_data_for_view) ---
        pileup_reads_data_for_view = []
        reads_for_view_count = 0
        for segment in aligned_read_segments:
            if reads_for_view_count >= TENSOR_MAX_READ_ROWS: break
            row_chars = get_read_representation_in_window_for_view(
                segment["cigar_ops"], segment["offset_on_node"], segment["read_sequence"],
                current_window_start_on_node, variant_processing_window_size, node_len)
            if any(c != ' ' for c in row_chars):
                pileup_reads_data_for_view.append({
                    "bases": [BASE_TO_INDEX.get(char.upper(), BASE_TO_INDEX['N']) for char in row_chars],
                    "offset": segment["offset_on_node"], "strand": segment["strand"],
                    "cigar": segment["original_cigar_str"]})
                reads_for_view_count += 1

        view_oriented_variant_data[variant_key_str] = {
            "pileup_reads_data": pileup_reads_data_for_view,
            "alt_allele_count": af_alt_count, "ref_allele_count_at_locus": af_ref_count,
            "other_allele_count_at_locus": af_other_count, "coverage_at_locus": af_locus_coverage,
            "alt_allele_frequency": round(alt_freq, 4)}

        # --- Prepare data for tensor ---
        tensor_ch1_bases_list, tensor_ch2_qualities_list, tensor_ch3_mismatches_list = [], [], []
        ref_base_indices_row_tensor = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i in range(TENSOR_WINDOW_SIZE):
            actual_node_pos = current_window_start_on_node + i
            if 0 <= actual_node_pos < node_len:
                ref_base_indices_row_tensor[i] = BASE_TO_INDEX.get(node_sequence[actual_node_pos].upper(),
                                                                   BASE_TO_INDEX['N'])

        tensor_ch1_bases_list.append(ref_base_indices_row_tensor)
        tensor_ch2_qualities_list.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
        tensor_ch3_mismatches_list.append([MISMATCH_CHANNEL_REF_ROW_VALUE] * TENSOR_WINDOW_SIZE)

        reads_added_to_tensor = 0
        for segment in aligned_read_segments:
            if reads_added_to_tensor >= TENSOR_MAX_READ_ROWS: break
            base_indices_row, quality_scores_row = get_read_tensor_rows_in_window(
                segment["cigar_ops"], segment["offset_on_node"],
                segment["read_sequence"], segment["processed_quality_str"],
                current_window_start_on_node, TENSOR_WINDOW_SIZE, node_len)
            if any(b != PADDING_BASE_INDEX for b in base_indices_row):
                tensor_ch1_bases_list.append(base_indices_row)
                tensor_ch2_qualities_list.append(quality_scores_row)
                mismatch_row = [MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE
                for i in range(TENSOR_WINDOW_SIZE):
                    read_base_idx, ref_base_idx = base_indices_row[i], ref_base_indices_row_tensor[i]
                    if read_base_idx == PADDING_BASE_INDEX or ref_base_idx == PADDING_BASE_INDEX:
                        mismatch_row[i] = MISMATCH_COMPARISON_PADDING_VALUE
                    else:
                        mismatch_row[i] = 0 if read_base_idx == ref_base_idx else 1
                tensor_ch3_mismatches_list.append(mismatch_row)
                reads_added_to_tensor += 1

        for _ in range(TENSOR_MAX_READ_ROWS - reads_added_to_tensor):  # Pad remaining rows
            tensor_ch1_bases_list.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            tensor_ch2_qualities_list.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            tensor_ch3_mismatches_list.append([MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE)

        try:
            final_tensor = torch.tensor([tensor_ch1_bases_list, tensor_ch2_qualities_list, tensor_ch3_mismatches_list],
                                        dtype=torch.int8)
            numpy_array = final_tensor.numpy()  # Convert to NumPy array
            tensor_filename = f"{variant_key_str}.npy"  # Use .npy extension
            tensor_filepath = os.path.join(node_specific_output_dir, tensor_filename)
            np.save(tensor_filepath, numpy_array)  # Save as .npy

            variant_headers_for_node.append({
                "variant_key": variant_key_str, "tensor_file": tensor_filename,
                "alt_allele_count": af_alt_count, "ref_allele_count_at_locus": af_ref_count,
                "other_allele_count_at_locus": af_other_count, "coverage_at_locus": af_locus_coverage,
                "alt_allele_frequency": round(alt_freq, 4)})
        except Exception as e:
            sys.stderr.write(f"❌ Node {node_id}: Error creating/saving tensor for {variant_key_str}: {e}\n")

    npy_files_generated_for_node = len(variant_headers_for_node)
    if variant_headers_for_node:
        summary_json_path = os.path.join(node_specific_output_dir, "variant_summary.json")
        try:
            with open(summary_json_path, 'w') as sjf:
                json.dump({
                    "node_id": node_id, "node_length": node_len,
                    "node_sequence_preview": node_sequence[:100] + "..." if node_len > 100 else node_sequence,
                    "variants": variant_headers_for_node}, sjf, indent=2)
        except Exception as e:
            sys.stderr.write(f"❌ Node {node_id}: Error writing summary JSON: {e}\n")

    return node_id, view_oriented_variant_data, npy_files_generated_for_node


## ─────────────────────────────────────────────────────────────────────────────
## Pileup Viewing Function (display_pileup_data)
## ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if not node_data_for_display_view or not isinstance(node_data_for_display_view, dict):
        print(f"ℹ️ No valid pileup data to display for node {node_id_str_for_display}.", file=sys.stderr)
        return
    print(f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")
    if not node_data_for_display_view:
        print(f"ℹ️ No variants met AF threshold or found for node {node_id_str_for_display}.")
        return

    variants_displayed_count = 0
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))
    display_window_size, half_display_window = TENSOR_WINDOW_SIZE, TENSOR_WINDOW_SIZE // 2

    for variant_key in sorted_variant_keys:
        if variants_displayed_count >= max_variants_to_display:
            print(
                f"\n  ... (and {len(node_data_for_display_view) - variants_displayed_count} more variants for node {node_id_str_for_display} not shown due to display limit)")
            break

        variant_data = node_data_for_display_view[variant_key]
        v_pos, v_type = int(variant_key.split('_')[0]), variant_key.split('_')[1]
        window_center = v_pos + 1 if v_type == 'I' else v_pos
        window_start = max(0, window_center - half_display_window)  # Ensure window start is not negative

        print(f"\n--- Variant: {variant_key} (Node Pos: {v_pos}, Type: {v_type}) ---")
        print(f"  Display Window on Node: {window_start} - {window_start + display_window_size - 1}")
        for stat_key, stat_name in [
            ("alt_allele_count", "Alt Count"), ("ref_allele_count_at_locus", "Ref Count"),
            ("other_allele_count_at_locus", "Other Count"), ("coverage_at_locus", "Coverage")]:
            print(f"  {stat_name}: {variant_data.get(stat_key, 'N/A')}")
        alt_freq = variant_data.get('alt_allele_frequency', 'N/A')
        print(f"  Alt Freq: {alt_freq:.4f}" if isinstance(alt_freq, float) else f"  Alt Freq: {alt_freq}")

        ref_display, marker_line = [' '] * display_window_size, [' '] * display_window_size
        var_idx_in_window = (v_pos - window_start) if window_start <= v_pos < window_start + display_window_size else -1

        for i in range(display_window_size):
            actual_pos = window_start + i
            if 0 <= actual_pos < len(full_node_sequence): ref_display[i] = full_node_sequence[actual_pos]
            if i == var_idx_in_window:
                marker_line[i] = "I" if v_type == 'I' else "^"
                if v_type == 'I' and i + 1 < display_window_size:
                    marker_line[i + 1] = "^"
                elif v_type == 'I' and i == display_window_size - 1:
                    marker_line[i] = ">"

        print(f"  Node Ref: {''.join(ref_display)}")
        print(f"  Marker  : {''.join(marker_line)}")

        pileup_reads = variant_data.get("pileup_reads_data", [])
        if not pileup_reads:
            print("  (No reads data in window for display)")
        else:
            for i, read_info in enumerate(pileup_reads):  # Already limited by TENSOR_MAX_READ_ROWS in worker
                if i >= max_reads_to_display_per_variant:  # Apply user's view limit if smaller
                    print(f"  ... ({len(pileup_reads) - i} more reads not shown due to view limit)")
                    break
                pileup_row_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in read_info["bases"]])
                print(
                    f"  Read {i + 1:3d}: {pileup_row_str}  (Off: {read_info['offset']}, Str: {read_info['strand']}, CIG: {read_info.get('cigar', 'N/A')})")
        variants_displayed_count += 1
    print("\n")


## ─────────────────────────────────────────────────────────────────────────────
## Main function
## ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered pileups (.npy tensors and JSON summaries) for specified node(s) serially. Reads GFA and IDX once.",
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
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N_VARIANTS',
                        help="Print generated pileups to console. No value or -1 for all variants. Optionally specify N for first N variants per node. Default: No view.")
    parser.add_argument("--max_view_reads", type=int, default=20,
                        help="Max reads to display per pileup in console view.")
    parser.add_argument("--min_af", type=float, default=0.1,
                        help="Min allele frequency for a variant to be processed and included in output.")
    args = parser.parse_args()

    # --- Argument Validation ---
    if not os.path.isfile(args.dat): sys.exit(f"❌ Error: DAT file not found: {args.dat}")
    if not os.path.isfile(args.idx): sys.exit(f"❌ Error: Index file not found: {args.idx}")
    if not args.load_cache and not args.gfa:
        sys.exit("❌ Error: Must provide --gfa or --load-cache to obtain node sequences.")
    if args.load_cache and not os.path.isfile(args.load_cache) and os.path.exists(args.load_cache):
        sys.exit(f"❌ Error: Specified --load-cache path '{args.load_cache}' is not a file.")
    if args.gfa and not os.path.isfile(args.gfa):
        sys.exit(f"❌ Error: GFA file not found: {args.gfa}")
    if not (0.0 <= args.min_af <= 1.0):
        sys.exit("❌ Error: --min_af must be between 0.0 and 1.0.")

    try:
        os.makedirs(args.output, exist_ok=True)
        print(f"🔹 Base output directory: {args.output}")
    except OSError as e:
        sys.exit(f"❌ Error: Could not create base output directory {args.output}: {e}")

    # --- Determine Target Node IDs ---
    target_node_ids_input = set()  # Use a set to store unique integer node IDs
    if args.node_id_file:
        try:
            with open(args.node_id_file, 'r') as f_nodes:
                for line_num, line in enumerate(f_nodes, 1):
                    line = line.strip()
                    if line and not line.startswith("#"):
                        try:
                            target_node_ids_input.add(int(line))
                        except ValueError:
                            sys.stderr.write(
                                f"⚠️ Warning: Invalid non-integer node ID '{line}' in {args.node_id_file} at line {line_num}. Skipping.\n")
            if not target_node_ids_input:
                sys.exit(f"❌ Error: No valid node IDs found or specified in {args.node_id_file}.")
            print(f"🔹 Will process {len(target_node_ids_input)} unique node ID(s) from file: {args.node_id_file}")
        except FileNotFoundError:
            sys.exit(f"❌ Error: Node ID file not found: {args.node_id_file}")
    elif args.node_id is not None:
        target_node_ids_input.add(args.node_id)
        print(f"🔹 Will process single target node ID: {args.node_id}")

    if not target_node_ids_input:  # Should not happen due to required group, but good check
        print("ℹ️ No target node IDs specified. Exiting.")
        sys.exit(0)

    # --- Load Initial Data (IDX, Sequence Cache, GFA) ---
    overall_start_time = time.time()

    idx_data_map = load_full_idx_data(args.idx)
    if idx_data_map is None:  # Error occurred during IDX loading
        sys.exit("❌ Critical error loading IDX data. Exiting.")

    node_sequences_map = {}  # Holds node_id_str -> sequence
    if args.load_cache and os.path.isfile(args.load_cache):
        s_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                node_sequences_map = json.load(cf)
            print(
                f"✔ Loaded {len(node_sequences_map)} sequences from cache '{args.load_cache}' in {time.time() - s_time:.2f}s.")
        except Exception as e:
            sys.stderr.write(f"⚠️ Error loading cache {args.load_cache}: {e}. Will proceed without cached sequences.\n")
            node_sequences_map = {}  # Ensure it's empty on error

    # Identify nodes needing sequences from GFA
    nodes_needing_gfa_sequence = {nid_int for nid_int in target_node_ids_input if
                                  str(nid_int) not in node_sequences_map}
    if nodes_needing_gfa_sequence and args.gfa:
        s_time = time.time()
        fetched_sequences = load_multiple_node_sequences_from_gfa(args.gfa, nodes_needing_gfa_sequence)
        if fetched_sequences:
            node_sequences_map.update(fetched_sequences)  # Add newly fetched to our main map
            print(
                f"✔ Fetched {len(fetched_sequences)} new sequences from GFA in {time.time() - s_time:.2f}s. Total sequences in map: {len(node_sequences_map)}.")
    elif nodes_needing_gfa_sequence and not args.gfa:
        sys.stderr.write(
            f"⚠️ {len(nodes_needing_gfa_sequence)} node(s) need sequences, but --gfa not provided. These will be skipped if not cached.\n")

    # --- Main Serial Processing Loop ---
    processed_nodes_count = 0
    successful_nodes_with_output = 0
    total_npy_files_generated = 0

    # Sort target_node_ids_input for deterministic processing order (optional but good practice)
    sorted_target_node_ids = sorted(list(target_node_ids_input))

    for i, target_node_id_int in enumerate(sorted_target_node_ids):
        node_id_str = str(target_node_id_int)
        print(f"\n══════ PROCESSING NODE: {target_node_id_int} ({i + 1}/{len(sorted_target_node_ids)}) ══════")
        node_specific_start_time = time.time()

        # Get DAT info from pre-loaded map
        node_dat_info = idx_data_map.get(target_node_id_int)
        if not node_dat_info:
            sys.stderr.write(f"❌ Node {target_node_id_int}: Not found in loaded IDX data. Skipping.\n")
            continue
        dat_offset, n_records = node_dat_info

        # Get sequence from pre-loaded map
        node_sequence = node_sequences_map.get(node_id_str)
        if not node_sequence:
            sys.stderr.write(f"❌ Node {target_node_id_int}: Sequence not available (not in cache or GFA). Skipping.\n")
            continue

        print(f"  Node {target_node_id_int}: Len={len(node_sequence)}bp, Records={n_records}, MinAF={args.min_af}")

        # Call the processing function directly (serially)
        try:
            _, view_data_for_node, npy_files_for_this_node = process_node_serially(
                args.dat, args.output, target_node_id_int, dat_offset, n_records, node_sequence, args.min_af
            )
            processed_nodes_count += 1
            total_npy_files_generated += npy_files_for_this_node

            if npy_files_for_this_node > 0:  # or if summary exists
                node_summary_file_path = os.path.join(args.output, node_id_str, "variant_summary.json")
                if os.path.exists(node_summary_file_path):
                    print(f"  ✅ Output generated for node {target_node_id_int}. Summary: {node_summary_file_path}")
                    successful_nodes_with_output += 1
                else:  # npy files might exist even if summary doesn't (e.g. write error for summary)
                    print(
                        f"  Generated {npy_files_for_this_node} .npy file(s) for node {target_node_id_int}, but summary JSON might be missing.")
                    # still count as successful if npy files were made
                    if npy_files_for_this_node > 0: successful_nodes_with_output += 1


            elif view_data_for_node is not None and not view_data_for_node:  # No variants passed AF
                print(f"  ℹ️ No variants met AF threshold for node {target_node_id_int}.")
            elif view_data_for_node is None:  # Error in processing
                print(f"  ⚠️ Processing error for node {target_node_id_int}, check logs above.")

            # Display Pileup Data (Optional)
            if args.view is not None:
                if view_data_for_node and isinstance(view_data_for_node, dict) and view_data_for_node:
                    max_v_to_show = float('inf') if args.view == -1 else args.view
                    if args.view < -1: max_v_to_show = float('inf')
                    view_limit_msg = f"first {int(max_v_to_show)}" if max_v_to_show != float('inf') else "all"
                    print(
                        f"  🔹 Displaying {view_limit_msg} pileups for node {target_node_id_int} (max {args.max_view_reads} reads/variant)...")
                    display_pileup_data(view_data_for_node, node_id_str, node_sequence, args.max_view_reads,
                                        max_v_to_show)
                elif view_data_for_node is not None:  # Empty dict {}
                    print(f"  ℹ️ --view: No variants met AF for node {target_node_id_int} to display.")

            print(f"  ✔ Node {target_node_id_int} processing took {time.time() - node_specific_start_time:.2f}s.")

        except Exception as e:
            sys.stderr.write(f"❌ CRITICAL ERROR processing node {target_node_id_int} in main loop: {e}\n")
            # import traceback; traceback.print_exc(file=sys.stderr) # For detailed debugging
            # Decide if you want to continue with other nodes or exit
            # For now, it will continue to the next node.

    print("\n═════════════════════════════ OVERALL PROCESSING COMPLETE ═════════════════════════════")

    # --- Save/Update Sequence Cache (if specified) ---
    if args.save_cache:
        if node_sequences_map:
            print(f"\n🔹 Saving/Updating {len(node_sequences_map)} sequence(s) to cache file: {args.save_cache}...")
            try:
                with open(args.save_cache, 'w') as wcf:
                    json.dump(node_sequences_map, wcf, indent=2)
                print(f"✔ Sequences saved to cache '{args.save_cache}'.")
            except Exception as e:
                sys.stderr.write(f"❌ Error saving sequences to cache {args.save_cache}: {e}\n")
        elif os.path.exists(args.save_cache) and args.load_cache == args.save_cache and not node_sequences_map:
            try:  # Empty the cache if it was loaded as empty and nothing new was added
                with open(args.save_cache, 'w') as wcf:
                    json.dump({}, wcf, indent=2)
                print(f"✔ Cache file {args.save_cache} effectively emptied as no sequences remain in map.")
            except Exception as e:
                sys.stderr.write(f"❌ Error writing empty cache {args.save_cache}: {e}\n")
        else:
            print(
                f"ℹ️ --save-cache specified, but no sequences were in memory to save. Cache file '{args.save_cache}' not created/updated.")

    # --- Final Summary ---
    total_script_time = time.time() - overall_start_time
    print(f"\nSummary:")
    print(f"  Attempted to process: {len(sorted_target_node_ids)} node(s).")
    print(f"  Successfully processed (data operations completed): {processed_nodes_count} node(s).")
    skipped_count = len(sorted_target_node_ids) - processed_nodes_count
    if skipped_count > 0:
        print(f"  Skipped due to missing data or critical errors before/during processing: {skipped_count} node(s).")
    print(f"  Output files (summary.json and/or .npy) generated for: {successful_nodes_with_output} node(s).")
    print(f"  Total .npy tensor files generated: {total_npy_files_generated}.")
    print(f"🏁 Total script execution time: {total_script_time:.2f} seconds.")


if __name__ == '__main__':
    main()