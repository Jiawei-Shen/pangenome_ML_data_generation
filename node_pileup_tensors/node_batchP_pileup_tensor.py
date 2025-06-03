#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np  # Not explicitly used in latest, but often kept for bioinformatics
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
TENSOR_MAX_READ_ROWS = 200  # Max reads per pileup (excluding reference)
PADDING_BASE_INDEX = BASE_TO_INDEX[' ']
DEFAULT_QUALITY_PADDING = 0
MISMATCH_CHANNEL_REF_ROW_VALUE = 0  # Mismatch value for ref base compared to itself (0=match)
MISMATCH_COMPARISON_PADDING_VALUE = -1  # Value for padding in mismatch channel where no comparison happens

# Globals for worker process state (set by initializer)
worker_dat_file = None
worker_base_output_dir = None


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
            if file_size < 4:  # Must have at least num_nodes
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
                return idx_data_map  # Return empty map if header says 0 nodes

            # Check if file size is consistent with num_nodes_in_idx
            expected_min_size = 4 + (num_nodes_in_idx * 22)  # 4 for header, 22 per record
            if file_size < expected_min_size:
                sys.stderr.write(
                    f"⚠️ Warning: Index file size ({file_size} bytes) is smaller than expected ({expected_min_size} bytes) for {num_nodes_in_idx} records. File may be truncated.\n")
                # We can still try to read what's there.

            processed_entries = 0
            for i in range(num_nodes_in_idx):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                    sys.stderr.write(
                        f"❌ Error: Index file ended prematurely while reading record {i + 1}/{num_nodes_in_idx}. Loaded {processed_entries} entries.\n")
                    break
                    # node_id_from_idx (I), offset (Q), num_uniq_kmers (I), n_records (I), avg_records_per_kmer (H)
                node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                idx_data_map[node_id_from_idx] = (offset, n_records)
                processed_entries += 1
                if processed_entries > 0 and processed_entries % 2000000 == 0:  # Progress for very large index files
                    print(f"    Loaded {processed_entries}/{num_nodes_in_idx} index entries...")

            if processed_entries != num_nodes_in_idx and len(
                    record_bytes) == 22:  # If loop finished but not all read (should not happen if f.read worked)
                sys.stderr.write(
                    f"⚠️ Warning: Read {processed_entries} entries, but index header indicated {num_nodes_in_idx}.\n")
            print(f"✔ Successfully loaded {len(idx_data_map)} distinct node entries from index file {idx_path}.")
        return idx_data_map
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: Index file not found at {idx_path}\n")
        return None
    except Exception as e:
        sys.stderr.write(f"❌ Error parsing full index file {idx_path}: {e}\n")
        # import traceback; traceback.print_exc() # Uncomment for detailed debug
        return None


def load_multiple_node_sequences_from_gfa(gfa_path, target_node_ids_set):
    """Reads GFA once for sequences of target_node_ids_set (integers). Returns dict: node_id_str -> sequence."""
    node_sequences = {}
    if not target_node_ids_set:
        return node_sequences
    nodes_to_find = target_node_ids_set.copy()

    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file to find sequences for {len(nodes_to_find)} nodes: {gfa_path}")
            line_counter = 0
            for line in f:
                line_counter += 1
                if line_counter % 10_000_000 == 0:
                    print(
                        f"  Checked {line_counter:,} lines in GFA file... {len(nodes_to_find)} nodes remaining to find.")

                if not line.startswith('S\t'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                try:
                    nid_int_from_gfa = int(parts[1])
                except ValueError:
                    continue
                if nid_int_from_gfa in nodes_to_find:
                    node_sequences[str(nid_int_from_gfa)] = parts[2]
                    nodes_to_find.remove(nid_int_from_gfa)
                    if not nodes_to_find:
                        print(
                            f"✔ Found all {len(target_node_ids_set)} requested node sequences in GFA after checking {line_counter:,} lines.")
                        break

            found_count = len(node_sequences)
            requested_count = len(target_node_ids_set)
            if not nodes_to_find:  # All found
                if found_count != requested_count:  # Should not happen if nodes_to_find is empty
                    sys.stderr.write(
                        f"⚠️ GFA Load: Mismatch - found_count {found_count}, requested {requested_count}, but all marked found.\n")
            else:  # Some not found
                print(f"✔ Finished GFA scan ({line_counter:,} lines). Found {found_count}/{requested_count} sequences.")
                print(
                    f"⚠️ Warning: Could not find GFA sequences for {len(nodes_to_find)} node ID(s). Examples: {list(nodes_to_find)[:5]}")
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: GFA file not found at {gfa_path}\n")
        return {}
    except Exception as e:
        sys.stderr.write(f"❌ Error reading GFA file {gfa_path}: {e}\n")
        return node_sequences
    return node_sequences


def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*': return []
    ops = []
    try:
        for length_str, op_char in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string):
            ops.append((int(length_str), op_char))
    except Exception as e:
        sys.stderr.write(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}\n")
        return []
    return ops


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
                    deleted_seq_in_read_context = node_sequence[current_node_pos: current_node_pos + length]
                    if deleted_seq_in_read_context == expected_ref_allele_for_indel:
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
                    if node_base != read_base and op != '=':
                        variants.append((current_node_pos, 'X', read_base, node_base))
                else:
                    break  # Alignment goes off end of sequence
            node_pos += length
            read_pos += length
        elif op == 'I':
            inserted_sequence = read_sequence[read_pos: read_pos + length].upper()
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0
            ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if 0 <= ref_anchor_pos < node_seq_len else "*"
            variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor))
            read_pos += length
        elif op == 'D':
            deleted_sequence_from_ref = node_sequence[
                                        node_pos: node_pos + length].upper() if node_pos + length <= node_seq_len else ""
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


# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function
# ─────────────────────────────────────────────────────────────────────────────

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker):
    global worker_dat_file, worker_base_output_dir
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path_for_worker}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path_for_worker}: {e}\n")
        sys.exit(1)


def process_single_node_for_pileup(task_args_with_af_thresh):
    node_id, dat_file_offset, n_records, node_sequence, min_af_threshold = task_args_with_af_thresh
    global worker_dat_file, worker_base_output_dir

    pth_files_generated_for_node = 0
    if worker_dat_file is None or worker_base_output_dir is None:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Worker not initialized.\n")
        return node_id, None, pth_files_generated_for_node
    if not node_sequence:
        sys.stderr.write(f"ℹ️ [Worker {os.getpid()} for Node {node_id}]: No sequence. Skipping.\n")
        return node_id, {}, pth_files_generated_for_node

    node_specific_output_dir = os.path.join(worker_base_output_dir, str(node_id))
    try:
        os.makedirs(node_specific_output_dir, exist_ok=True)
    except OSError as e:
        sys.stderr.write(
            f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Dir creation {node_specific_output_dir}: {e}\n")
        return node_id, None, pth_files_generated_for_node

    node_len = len(node_sequence)
    view_oriented_variant_data = {}
    aligned_read_segments = []
    try:
        worker_dat_file.seek(dat_file_offset + 10)
        for record_idx in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE: break
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
                current_decoded_cigar_ops = [op for op in
                                             reversed(original_decoded_cigar_ops)] if original_decoded_cigar_ops else []

                # Using user's original logic for reverse strand offset calculation
                alignment_span_on_node_user_logic = len(current_read_sequence)
                current_offset_on_node = node_len - alignment_span_on_node_user_logic - off_from_file
                if current_offset_on_node < 0: continue

            aligned_read_segments.append({
                "offset_on_node": current_offset_on_node, "read_sequence": current_read_sequence,
                "processed_quality_str": current_quality_str, "cigar_ops": current_decoded_cigar_ops,
                "original_cigar_str": cigar_str_original, "strand": strand_char
            })
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()} for Node {node_id}] reading DAT: {e}\n")
        return node_id, None, pth_files_generated_for_node

    if not aligned_read_segments:
        return node_id, {}, pth_files_generated_for_node

    candidate_variants = defaultdict(int)
    for segment in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                segment["offset_on_node"], segment["cigar_ops"], segment["read_sequence"], node_sequence):
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    variant_headers_for_node = []
    variant_processing_window_size = TENSOR_WINDOW_SIZE
    half_window = variant_processing_window_size // 2

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate_variants.items():
        af_alt_count, af_ref_count, af_other_count, af_locus_coverage = 0, 0, 0, 0
        expected_ref_for_af, expected_alt_for_af = "", ""
        ref_allele_for_indel_af_check = None

        if v_type == 'X':
            expected_ref_for_af, expected_alt_for_af = v_ref_from_cigar, v_alt_from_cigar
        elif v_type == 'D':
            expected_ref_for_af = node_sequence[v_pos] if v_pos < node_len else ""
            expected_alt_for_af = "*"
            ref_allele_for_indel_af_check = v_ref_from_cigar
        elif v_type == 'I':
            expected_ref_for_af = v_ref_from_cigar
            expected_alt_for_af = v_alt_from_cigar
            ref_allele_for_indel_af_check = node_sequence[v_pos] if v_pos < node_len else None

        for segment in aligned_read_segments:
            allele = get_allele_from_read_at_node_pos(
                segment["offset_on_node"], segment["read_sequence"], segment["cigar_ops"],
                v_pos, node_sequence, expected_var_type=v_type,
                expected_ref_allele_for_indel=ref_allele_for_indel_af_check
            )
            if allele is not None:
                af_locus_coverage += 1
                if allele == expected_alt_for_af:
                    af_alt_count += 1
                elif allele == expected_ref_for_af or (v_type in 'ID' and allele == "REF_STATE_FOR_INDEL"):
                    af_ref_count += 1
                else:
                    af_other_count += 1

        alt_freq = af_alt_count / af_locus_coverage if af_locus_coverage > 0 else 0.0
        if alt_freq < min_af_threshold: continue

        variant_key_str = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"
        window_center_on_node = v_pos + 1 if v_type == 'I' else v_pos
        current_window_start_on_node = max(0, window_center_on_node - half_window)

        pileup_reads_data_for_view = []
        for seg_idx, segment in enumerate(aligned_read_segments):
            if seg_idx >= TENSOR_MAX_READ_ROWS * 2: break  # Limit reads for view data generation too
            row_chars = get_read_representation_in_window_for_view(
                segment["cigar_ops"], segment["offset_on_node"], segment["read_sequence"],
                current_window_start_on_node, variant_processing_window_size, node_len)
            if any(c != ' ' for c in row_chars):
                pileup_reads_data_for_view.append({
                    "bases": [BASE_TO_INDEX.get(c.upper(), BASE_TO_INDEX['N']) for c in row_chars],
                    "offset": segment["offset_on_node"], "strand": segment["strand"],
                    "cigar": segment["original_cigar_str"]})

        view_oriented_variant_data[variant_key_str] = {
            "pileup_reads_data": pileup_reads_data_for_view[:TENSOR_MAX_READ_ROWS],
            "alt_allele_count": af_alt_count, "ref_allele_count_at_locus": af_ref_count,
            "other_allele_count_at_locus": af_other_count, "coverage_at_locus": af_locus_coverage,
            "alt_allele_frequency": round(alt_freq, 4)}

        tensor_ch1_bases, tensor_ch2_qualities, tensor_ch3_mismatches = [], [], []
        ref_bases_tensor = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i in range(TENSOR_WINDOW_SIZE):
            actual_node_pos = current_window_start_on_node + i
            if 0 <= actual_node_pos < node_len:
                ref_bases_tensor[i] = BASE_TO_INDEX.get(node_sequence[actual_node_pos].upper(), BASE_TO_INDEX['N'])

        tensor_ch1_bases.append(ref_bases_tensor)
        tensor_ch2_qualities.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
        tensor_ch3_mismatches.append([MISMATCH_CHANNEL_REF_ROW_VALUE] * TENSOR_WINDOW_SIZE)

        reads_added_to_tensor = 0
        for segment in aligned_read_segments:
            if reads_added_to_tensor >= TENSOR_MAX_READ_ROWS: break
            base_idx_row, qual_row = get_read_tensor_rows_in_window(
                segment["cigar_ops"], segment["offset_on_node"], segment["read_sequence"],
                segment["processed_quality_str"], current_window_start_on_node, TENSOR_WINDOW_SIZE, node_len)
            if any(b != PADDING_BASE_INDEX for b in base_idx_row):
                tensor_ch1_bases.append(base_idx_row)
                tensor_ch2_qualities.append(qual_row)
                mismatch_row = [MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE
                for i in range(TENSOR_WINDOW_SIZE):
                    if base_idx_row[i] == PADDING_BASE_INDEX or ref_bases_tensor[i] == PADDING_BASE_INDEX: continue
                    mismatch_row[i] = 0 if base_idx_row[i] == ref_bases_tensor[i] else 1
                tensor_ch3_mismatches.append(mismatch_row)
                reads_added_to_tensor += 1

        for _ in range(TENSOR_MAX_READ_ROWS - reads_added_to_tensor):
            tensor_ch1_bases.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            tensor_ch2_qualities.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            tensor_ch3_mismatches.append([MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE)

        try:
            final_tensor = torch.tensor([tensor_ch1_bases, tensor_ch2_qualities, tensor_ch3_mismatches],
                                        dtype=torch.int8)
            tensor_filename = f"{variant_key_str}.pth"
            torch.save(final_tensor, os.path.join(node_specific_output_dir, tensor_filename))
            variant_headers_for_node.append({
                "variant_key": variant_key_str, "tensor_file": tensor_filename,
                "alt_allele_count": af_alt_count, "ref_allele_count_at_locus": af_ref_count,
                "other_allele_count_at_locus": af_other_count, "coverage_at_locus": af_locus_coverage,
                "alt_allele_frequency": round(alt_freq, 4)})
            pth_files_generated_for_node += 1
        except Exception as e:
            sys.stderr.write(f"❌ Error [Worker Node {node_id}]: Tensor save for {variant_key_str}: {e}\n")

    if variant_headers_for_node:
        summary_path = os.path.join(node_specific_output_dir, "variant_summary.json")
        try:
            with open(summary_path, 'w') as sjf:
                json.dump({"node_id": node_id, "node_length": node_len,
                           "node_sequence_preview": node_sequence[:100] + ("..." if node_len > 100 else ""),
                           "variants_passing_af_filter": variant_headers_for_node}, sjf, indent=2)
        except Exception as e:
            sys.stderr.write(f"❌ Error [Worker Node {node_id}]: Summary JSON write: {e}\n")

    return node_id, view_oriented_variant_data, pth_files_generated_for_node


# ─────────────────────────────────────────────────────────────────────────────
# Pileup Viewing Function
# ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if not node_data_for_display_view:
        print(f"ℹ️ No pileup data to display for node {node_id_str_for_display}.", file=sys.stderr)
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
        window_start = max(0, window_center - half_display_window)

        print(f"\n--- Variant: {variant_key} (Node Pos: {v_pos}, Type: {v_type}) ---")
        print(f"  Display Window on Node: {window_start} - {window_start + display_window_size - 1}")
        for stat in ["alt_allele_count", "ref_allele_count_at_locus", "other_allele_count_at_locus",
                     "coverage_at_locus"]:
            print(f"  {stat.replace('_', ' ').capitalize()}: {variant_data.get(stat, 'N/A')}")
        alt_freq = variant_data.get('alt_allele_frequency', 'N/A')
        print(f"  Alt allele frequency: {alt_freq:.4f}" if isinstance(alt_freq, float) else f"  Alt Freq: {alt_freq}")

        ref_display, marker_line = [' '] * display_window_size, [' '] * display_window_size
        var_idx_in_window = (v_pos - window_start) if window_start <= v_pos < window_start + display_window_size else -1

        for i in range(display_window_size):
            actual_pos = window_start + i
            if 0 <= actual_pos < len(full_node_sequence): ref_display[i] = full_node_sequence[actual_pos]
            if i == var_idx_in_window: marker_line[i] = "I^"[v_type == 'I']  # "I" or "^"

        print(f"  Node Ref: {''.join(ref_display)}")
        print(f"  Marker  : {''.join(marker_line)}")

        pileup_reads = variant_data.get("pileup_reads_data", [])
        if not pileup_reads:
            print("  (No reads data in window for display)")
        else:
            for i, read_info in enumerate(pileup_reads):
                if i >= max_reads_to_display_per_variant:
                    print(f"  ... ({len(pileup_reads) - i} more reads not shown due to view limit)")
                    break
                pileup_row_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in read_info["bases"]])
                print(
                    f"  Read {i + 1:3d}: {pileup_row_str}  (Off: {read_info['offset']}, Str: {read_info['strand']}, CIG: {read_info.get('cigar', 'N/A')})")
        variants_displayed_count += 1
    print("\n")


# ─────────────────────────────────────────────────────────────────────────────
# Main function
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered pileups for specified node(s) and optionally view them.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dat", help=".dat file path")
    parser.add_argument("idx", help=".idx file path")
    parser.add_argument("output", help="Base output directory")
    input_node_group = parser.add_mutually_exclusive_group(required=True)
    input_node_group.add_argument("--node_id", type=int, help="Specific node ID to process.")
    input_node_group.add_argument("--node_id_file", help="File with node IDs (one per line).")
    parser.add_argument("--gfa", help="GFA graph file path.")
    parser.add_argument("--load-cache", help="Load node sequences from JSON cache.")
    parser.add_argument("--save-cache", help="Save/update node sequences to JSON cache.")
    parser.add_argument("--num_workers", type=int, default=None,
                        help="Number of worker processes. Defaults to CPU cores.")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N_VARIANTS',
                        help="Print pileups. No value or -1 for all variants, N for first N variants per node.")
    parser.add_argument("--max_view_reads", type=int, default=20, help="Max reads per pileup in console view.")
    parser.add_argument("--min_af", type=float, default=0.1, help="Min allele frequency for variant processing.")
    args = parser.parse_args()

    # --- Initial Validations ---
    if not os.path.isfile(args.dat): sys.exit(f"❌ DAT file not found: {args.dat}")
    if not os.path.isfile(args.idx): sys.exit(f"❌ Index file not found: {args.idx}")
    if not args.load_cache and not args.gfa: sys.exit("❌ Must provide --gfa or --load-cache.")
    if args.load_cache and not os.path.isfile(args.load_cache) and os.path.exists(args.load_cache):
        sys.exit(f"❌ Cache path '{args.load_cache}' is not a file.")
    if args.gfa and not os.path.isfile(args.gfa): sys.exit(f"❌ GFA file not found: {args.gfa}")
    if not (0.0 <= args.min_af <= 1.0): sys.exit("❌ --min_af must be between 0.0 and 1.0.")

    num_workers = args.num_workers if args.num_workers and args.num_workers > 0 else (os.cpu_count() or 1)
    print(f"🔹 Using {num_workers} worker process(es).")
    os.makedirs(args.output, exist_ok=True)
    print(f"🔹 Base output directory: {args.output}")

    # --- Target Node IDs ---
    target_node_ids_set = set()
    if args.node_id_file:
        try:
            with open(args.node_id_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        try:
                            target_node_ids_set.add(int(line))
                        except ValueError:
                            sys.stderr.write(f"⚠️ Invalid node ID '{line}' in file. Skipping.\n")
            if not target_node_ids_set: sys.exit(f"❌ No valid node IDs in {args.node_id_file}.")
            print(f"🔹 Will process {len(target_node_ids_set)} unique node ID(s) from {args.node_id_file}")
        except FileNotFoundError:
            sys.exit(f"❌ Node ID file not found: {args.node_id_file}")
    elif args.node_id is not None:
        target_node_ids_set.add(args.node_id)
        print(f"🔹 Will process single target node ID: {args.node_id}")
    if not target_node_ids_set: sys.exit("ℹ️ No target node IDs. Exiting.")

    # --- Sequence Loading (Cache & GFA) ---
    node_sequences_map = {}  # node_id_str -> sequence
    if args.load_cache and os.path.isfile(args.load_cache):
        s_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                node_sequences_map = json.load(cf)
            print(
                f"✔ Loaded {len(node_sequences_map)} sequences from cache '{args.load_cache}' in {time.time() - s_time:.2f}s.")
        except Exception as e:
            sys.stderr.write(f"⚠️ Error loading cache {args.load_cache}: {e}. Cache not loaded.\n")

    nodes_needing_gfa_sequence = {nid for nid in target_node_ids_set if str(nid) not in node_sequences_map}
    if nodes_needing_gfa_sequence and args.gfa:
        print(f"🔹 {len(nodes_needing_gfa_sequence)} node(s) require sequence fetching from GFA.")
        s_time = time.time()
        fetched = load_multiple_node_sequences_from_gfa(args.gfa, nodes_needing_gfa_sequence)
        node_sequences_map.update(fetched)
        print(
            f"✔ Fetched {len(fetched)} new sequences from GFA in {time.time() - s_time:.2f}s. Total sequences in map: {len(node_sequences_map)}.")
    elif nodes_needing_gfa_sequence and not args.gfa:
        sys.stderr.write(
            f"⚠️ {len(nodes_needing_gfa_sequence)} nodes need sequences, but --gfa not provided. These will be skipped if not cached.\n")

    # --- Load Full Index Data ---
    idx_load_start_time = time.time()
    full_idx_map = load_full_idx_data(args.idx)  # int_node_id -> (offset, n_records)
    if full_idx_map is None: sys.exit(f"❌ Failed to load index data from {args.idx}. Exiting.")
    print(f"✔ Index data ({len(full_idx_map)} entries) loaded in {time.time() - idx_load_start_time:.2f}s.")

    # --- Prepare Tasks for Workers ---
    tasks_to_submit, skipped_nodes_pre_submit = [], set()
    print(f"🔹 Preparing tasks for {len(target_node_ids_set)} target nodes...")
    task_prep_start_time = time.time()
    for i, node_id_int in enumerate(target_node_ids_set):
        if (i + 1) % 50000 == 0: print(f"  Prepared tasks for {i + 1}/{len(target_node_ids_set)} nodes...")
        node_seq = node_sequences_map.get(str(node_id_int))
        node_dat_info = full_idx_map.get(node_id_int)
        if not node_seq or not node_dat_info:
            skipped_nodes_pre_submit.add(node_id_int)
            continue
        tasks_to_submit.append((node_id_int, node_dat_info[0], node_dat_info[1], node_seq, args.min_af))
    print(f"✔ Task preparation completed in {time.time() - task_prep_start_time:.2f}s.")
    if skipped_nodes_pre_submit:
        print(
            f"⚠️ Skipped {len(skipped_nodes_pre_submit)} nodes (missing sequence/index data). Examples: {list(skipped_nodes_pre_submit)[:3]}")
    if not tasks_to_submit: sys.exit("ℹ️ No valid tasks to process. Exiting.")

    # --- Parallel Processing ---
    overall_start_time = time.time()  # Marks the start of the parallel execution phase
    total_pth_files_generated = 0  # Cumulative total .pth files
    processed_nodes_count = 0  # Cumulative count of futures that completed (success or failure)
    successful_nodes_output = 0  # Cumulative count of nodes for which output files were made

    results_for_viewing = {}  # node_id_int -> (view_data, node_sequence)

    # For batch reporting
    batch_start_time = time.time()
    batch_nodes_processed_count = 0
    batch_pth_files_generated = 0

    print(f"\n🔹 Submitting {len(tasks_to_submit)} node tasks to {num_workers} worker(s)...")
    with ProcessPoolExecutor(max_workers=num_workers, initializer=init_worker,
                             initargs=(args.dat, args.output)) as executor:
        future_to_node_id = {executor.submit(process_single_node_for_pileup, task): task[0] for task in tasks_to_submit}

        for future_idx, future in enumerate(as_completed(future_to_node_id)):
            current_completed_total_count = future_idx + 1  # 1-indexed count of all completed futures
            original_node_id_for_future = future_to_node_id[future]

            try:
                # Attempt to get the result from the worker
                res_node_id, view_data, pth_files_for_this_node = future.result()

                processed_nodes_count += 1  # A future completed and returned a result (even if res_node_id is None)

                if res_node_id is None:
                    # Worker signaled an internal issue by returning None for node_id
                    sys.stderr.write(
                        f"❌ Worker task for node {original_node_id_for_future} indicated failure (returned None ID).\n")
                else:
                    # Worker processed the node and returned data
                    total_pth_files_generated += pth_files_for_this_node
                    batch_pth_files_generated += pth_files_for_this_node  # Accumulate for current batch report

                    # Check if output was actually generated for this node
                    summary_file_path = os.path.join(args.output, str(res_node_id), "variant_summary.json")
                    if pth_files_for_this_node > 0 or os.path.exists(summary_file_path):
                        successful_nodes_output += 1

                    # Store data for --view option if enabled
                    if args.view is not None and view_data:  # view_data might be {}
                        sequence_for_viewing = node_sequences_map.get(str(res_node_id))
                        if sequence_for_viewing:
                            results_for_viewing[res_node_id] = (view_data, sequence_for_viewing)
                        # else: # Should be rare if node_sequences_map is correctly populated
                        #     sys.stderr.write(f"⚠️ Warning: Sequence for node {res_node_id} not found in map for viewing.\n")

            except Exception as exc:
                processed_nodes_count += 1  # A future completed by raising an exception
                sys.stderr.write(
                    f"❌ Error processing node {original_node_id_for_future} (future raised exception): {exc}\n")
                # For detailed debugging, you might want to uncomment the next two lines:
                # import traceback
                # traceback.print_exc(file=sys.stderr)

            batch_nodes_processed_count += 1

            # --- Batch Reporting Logic ---
            # Report every 1000 nodes processed in the batch, or if it's the very last task
            if batch_nodes_processed_count == 1000 or current_completed_total_count == len(tasks_to_submit):
                batch_end_time = time.time()
                batch_duration_seconds = batch_end_time - batch_start_time
                nodes_per_second_in_batch = batch_nodes_processed_count / batch_duration_seconds if batch_duration_seconds > 0 else 0

                print(f"  Processed batch of {batch_nodes_processed_count} nodes "
                      f"({current_completed_total_count}/{len(tasks_to_submit)} total completed) "
                      f"in {batch_duration_seconds:.2f}s ({nodes_per_second_in_batch:.2f} nodes/sec). "
                      f"Generated {batch_pth_files_generated} .pth files in this batch.")

                # Reset for the next batch
                batch_start_time = time.time()  # Mark start of a new batch
                batch_nodes_processed_count = 0
                batch_pth_files_generated = 0

    # --- Display Pileups (if --view) ---
    if args.view is not None and results_for_viewing:
        print("\n══════════ VIEWING PILEUPS ══════════")
        for node_id_view in sorted(results_for_viewing.keys()):
            view_data, node_seq = results_for_viewing[node_id_view]
            max_v = float('inf') if args.view == -1 else args.view
            display_pileup_data(view_data, str(node_id_view), node_seq, args.max_view_reads, max_v)
    elif args.view is not None:
        print(f"ℹ️ --view specified, but no pileup data to display.")

    # --- Final Summary and Cache Saving ---
    print("\n══════════ PROCESSING COMPLETE ══════════")
    if args.save_cache and node_sequences_map:
        print(f"\n🔹 Saving {len(node_sequences_map)} sequences to cache: {args.save_cache}...")
        try:
            with open(args.save_cache, 'w') as wcf:
                json.dump(node_sequences_map, wcf, indent=2)
            print(f"✔ Sequences saved to cache.")
        except Exception as e:
            sys.stderr.write(f"❌ Error saving to cache {args.save_cache}: {e}\n")
    elif args.save_cache:
        print(f"ℹ️ --save-cache: No sequences in map to save.")

    print(f"\nSummary:")
    print(f"  Targeted: {len(target_node_ids_set)} unique node IDs.")
    if skipped_nodes_pre_submit: print(f"  Skipped pre-submission: {len(skipped_nodes_pre_submit)} node(s).")
    print(f"  Submitted to workers: {len(tasks_to_submit)} node(s).")
    print(f"  Worker tasks completed: {processed_nodes}/{len(tasks_to_submit)}.")
    print(f"  Output files generated for: {successful_nodes_output} node(s).")
    print(f"  Total .pth tensor files: {total_pth}.")
    print(f"🏁 Script finished in {time.time() - overall_start_time:.2f} seconds (parallel processing phase).")


if __name__ == '__main__':
    main()