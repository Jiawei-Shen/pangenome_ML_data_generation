#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np
from collections import defaultdict
import re  # Import re at the top level
from concurrent.futures import ProcessPoolExecutor

# from os import cpu_count # For defaulting workers, can be set explicitly

# ─────────────────────────────────────────────────────────────────────────────
# Constants
RECORD_STRUCT = struct.Struct("<h150s150s20shc")  # Read offset, sequence, base qualities, CIGAR, MAPQ, strand
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 4,
                 ' ': 4}  # '*' for deletions, ' ' for no coverage in window
INDEX_TO_BASE_FOR_VIEW = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N', '-': '*'}  # '-' for viewing deletions

# Global for worker process state (file handle)
worker_dat_file = None


# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions

def reverse_complement(sequence):
    """Returns the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]


def parse_idx_file_for_single_node(idx_path, target_node_id):
    """
    Parses the .idx file to find the offset and record count for a specific target_node_id.
    """
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
                record_bytes = f.read(22)  # node_id (I), offset (Q), block_size (I), n_records (I), padding (H)
                if len(record_bytes) < 22:
                    print(f"❌ Error: Index file ended prematurely while reading record {i + 1}/{num_nodes}.",
                          file=sys.stderr)
                    break
                node_id, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                if node_id == target_node_id:
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
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
        sys.exit(1)


def load_node_sequence_from_gfa(gfa_path, target_node_id):
    """
    Loads the sequence for a single target_node_id from a GFA file.
    """
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
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr)
        sys.exit(1)
    return node_sequence


def decode_cigar_to_int_ops(cigar_string):
    """Decodes a CIGAR string into a list of (int_length, operation_str) tuples."""
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
    """
    Determines the allele presented by a read at a specific target_node_pos.
    read_cigar_ops_decoded is a list of (int_length, str_op) tuples.
    For indels, expected_ref_allele_for_indel is the sequence deleted from ref, or anchor base for ins.
    """
    current_node_pos = read_offset_on_node
    current_read_pos = 0

    for length, op in read_cigar_ops_decoded:
        if op == 'M' or op == '=' or op == 'X':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I': return "REF_STATE_FOR_INDEL"  # Read matches over anchor, no insertion here
                if expected_var_type == 'D': return "REF_STATE_FOR_INDEL"  # Read matches over deleted region

                offset_in_block = target_node_pos - current_node_pos
                if current_read_pos + offset_in_block < len(read_sequence):
                    return read_sequence[current_read_pos + offset_in_block].upper()
                return None  # Should not happen
            current_node_pos += length
            current_read_pos += length
        elif op == 'I':
            # Insertion occurs after (current_node_pos - 1) on reference
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                return read_sequence[current_read_pos: current_read_pos + length].upper()
            current_read_pos += length
        elif op == 'D':
            # Deletion spans from current_node_pos to current_node_pos + length - 1 on reference
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I': return "OTHER_FOR_INDEL"  # Deletion where insertion expected
                if expected_var_type == 'D':
                    # The read shows a deletion. Does it match the *specific* deletion we're testing?
                    # The specific deletion is defined by v_ref_defined (expected_ref_allele_for_indel)
                    # and its length.
                    deleted_seq_in_read = node_sequence[current_node_pos: current_node_pos + length]
                    if deleted_seq_in_read == expected_ref_allele_for_indel:
                        return "*"  # Represents this specific deletion
                    else:
                        return "OTHER_FOR_INDEL"  # Different deletion
                return "*"  # General case: read shows some deletion here
            current_node_pos += length
        elif op == 'S':
            current_read_pos += length
        elif op == 'N':
            current_node_pos += length

        # Optimization to stop early if past the target
        if current_node_pos > target_node_pos + 1 and op in ('M', '=', 'X', 'D', 'N'):  # Small buffer for indel anchors
            # For insertions, target_node_pos is the anchor *before* the insertion.
            # If current_node_pos (start of next ref-consuming block) is beyond target_node_pos,
            # and we expected an insertion at target_node_pos, we've passed the point of check.
            if not (expected_var_type == 'I' and (current_node_pos - 1) <= target_node_pos):
                break
    return None


def detect_variants_from_cigar(offset_on_node, cigar_ops_decoded, read_sequence, node_sequence):
    """
    Detects variants from a single read.
    Returns list of (pos, type, alt, ref)
    """
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
                    if node_base != read_base and op != '=':
                        variants.append((current_node_pos, 'X', read_base, node_base))
            node_pos += length
            read_pos += length
        elif op == 'I':
            inserted_sequence = read_sequence[read_pos: read_pos + length].upper()
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0
            ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if ref_anchor_pos < len(
                node_sequence) and ref_anchor_pos >= 0 else ""
            variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor if ref_base_at_anchor else "*"))
            read_pos += length
        elif op == 'D':
            deleted_sequence_from_ref = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= len(
                node_sequence) else ""
            if deleted_sequence_from_ref:  # Only report if valid deletion from ref
                variants.append((node_pos, 'D', '*', deleted_sequence_from_ref))
            node_pos += length
        elif op == 'S':
            read_pos += length
        elif op == 'N':
            node_pos += length
    return variants


def get_read_representation_in_window(segment_cigar_ops, segment_offset_on_node, segment_read_sequence,
                                      window_start_node, window_size, node_len):
    """
    For a single read, generates its character representation across a defined window on the node.
    Returns a list of characters of length window_size.
    ' ' for no coverage, '*' for deletion, base for match/mismatch.
    Insertions are not explicitly expanded in this node-coordinate based window.
    """
    window_char_representation = [' '] * window_size  # Initialize with space for no coverage

    current_node_pos_in_read = segment_offset_on_node
    current_read_pos_in_read = 0

    for cigar_len, cigar_op in segment_cigar_ops:
        if cigar_op in ('M', '=', 'X'):  # Match, Equal, Mismatch
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                read_aln_pos = current_read_pos_in_read + i

                if node_aln_pos >= window_start_node and node_aln_pos < window_start_node + window_size:
                    window_idx = node_aln_pos - window_start_node
                    if read_aln_pos < len(segment_read_sequence):
                        window_char_representation[window_idx] = segment_read_sequence[read_aln_pos].upper()
            current_node_pos_in_read += cigar_len
            current_read_pos_in_read += cigar_len
        elif cigar_op == 'D' or cigar_op == 'N':  # Deletion from ref or Skip in ref
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                if node_aln_pos >= window_start_node and node_aln_pos < window_start_node + window_size:
                    window_idx = node_aln_pos - window_start_node
                    window_char_representation[window_idx] = '*'
            current_node_pos_in_read += cigar_len
        elif cigar_op == 'I':  # Insertion to ref
            # Insertions occur between node bases. They don't occupy a position *in* this node-based window.
            # The base *before* the insertion will be shown.
            current_read_pos_in_read += cigar_len
        elif cigar_op == 'S':  # Soft clip
            current_read_pos_in_read += cigar_len

        # Optimization: if we've processed beyond the window on the node
        if current_node_pos_in_read >= window_start_node + window_size:
            break

    return window_char_representation


# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function

def init_worker(dat_file_path_for_worker):
    global worker_dat_file
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
    except FileNotFoundError:
        print(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path_for_worker}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path_for_worker}: {e}", file=sys.stderr)
        sys.exit(1)


def process_single_node_for_pileup(task_args):
    node_id, dat_file_offset, n_records, node_sequence = task_args
    global worker_dat_file

    if worker_dat_file is None: return node_id, {}
    if not node_sequence: return node_id, {}

    node_len = len(node_sequence)
    final_variant_output = {}  # Store final data including AF and new pileup matrix
    aligned_read_segments = []  # Store all reads aligned to this node

    try:
        worker_dat_file.seek(dat_file_offset + 10)
        for record_idx in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE: break
            off, raw_seq, _, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)
            if mapq < 10: continue
            try:
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')
                cigar_str = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                continue

            current_read_sequence = reverse_complement(seq) if strand_char == '-' else seq
            decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str)  # Decode to (int_len, str_op)

            aligned_read_segments.append({
                "offset_on_node": off, "read_sequence": current_read_sequence,
                "cigar_ops": decoded_cigar_ops,
            })
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] reading records for node {node_id}: {e}", file=sys.stderr)
        return node_id, {}

        # 1. Identify all unique candidate variants from all reads
    candidate_variants = defaultdict(int)  # Key: (pos, type, ref, alt), Value: count of reads initially defining it
    for segment in aligned_read_segments:
        variants_in_read = detect_variants_from_cigar(
            segment["offset_on_node"], segment["cigar_ops"],
            segment["read_sequence"], node_sequence
        )
        for v_pos, v_type, v_alt, v_ref in variants_in_read:
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    window_size = 50  # As per user request
    half_window = window_size // 2

    # 2. For each candidate variant, calculate AF and build the new pileup matrix
    for (v_pos, v_type, v_ref_defined, v_alt_defined), _ in candidate_variants.items():

        # --- Allele Frequency Calculation ---
        af_alt_count = 0
        af_ref_count = 0
        af_other_count = 0
        af_locus_coverage = 0

        # For deletions, v_ref_defined is the deleted sequence. For AF, we need its length.
        # For insertions, v_ref_defined is the base before the insertion.
        ref_allele_for_indel_af = v_ref_defined if v_type == 'D' else (
            node_sequence[v_pos] if v_type == 'I' and 0 <= v_pos < node_len else None)

        for segment in aligned_read_segments:  # Iterate ALL reads
            allele_in_segment_at_v_pos = get_allele_from_read_at_node_pos(
                segment["offset_on_node"], segment["read_sequence"], segment["cigar_ops"],
                v_pos, node_sequence, v_type, ref_allele_for_indel_af
            )

            if allele_in_segment_at_v_pos is not None:  # Read informs the locus
                af_locus_coverage += 1
                if v_type == 'X':
                    if allele_in_segment_at_v_pos == v_alt_defined:
                        af_alt_count += 1
                    elif allele_in_segment_at_v_pos == v_ref_defined:
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'I':
                    if allele_in_segment_at_v_pos == v_alt_defined:
                        af_alt_count += 1
                    elif allele_in_segment_at_v_pos == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'D':
                    if allele_in_segment_at_v_pos == "*":  # Read shows the defined deletion
                        af_alt_count += 1
                    elif allele_in_segment_at_v_pos == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1
                    else:
                        af_other_count += 1

        alt_freq = af_alt_count / af_locus_coverage if af_locus_coverage > 0 else 0.0

        # --- Pileup Matrix Generation (from ALL covering reads for the window) ---
        pileup_matrix_for_json = []
        # Define window around v_pos. If v_pos is variant start.
        # Window: [v_pos - half_window + 1, v_pos + half_window] if v_pos is center of interest
        # Let's make window start at v_pos - (window_size/2 -1) and end at v_pos + window_size/2
        # This makes v_pos roughly in the middle.
        # If window_size = 50, half_window = 25.
        # Window: [v_pos - 24, v_pos + 25] (inclusive, 50bp)
        window_start_on_node = v_pos - (half_window - 1) if v_type != 'I' else (v_pos + 1) - (
                    half_window - 1)  # Adjust for insertion anchor
        # For insertion at v_pos (after base v_pos), center window around v_pos+1

        for segment in aligned_read_segments:
            row_chars = get_read_representation_in_window(
                segment["cigar_ops"], segment["offset_on_node"], segment["read_sequence"],
                window_start_on_node, window_size, node_len
            )
            if row_chars:  # If the read overlaps and provides content for the window
                row_indices = [BASE_TO_INDEX.get(char.upper(), BASE_TO_INDEX['N']) for char in row_chars]
                pileup_matrix_for_json.append(row_indices)

        variant_key_str = f"{v_pos}_{v_type}_{v_ref_defined}_{v_alt_defined}"
        final_variant_output[variant_key_str] = {
            "pileup_matrix": pileup_matrix_for_json if pileup_matrix_for_json else [],  # Ensure it's a list
            "alt_allele_count": af_alt_count,
            "ref_allele_count_at_locus": af_ref_count,
            "other_allele_count_at_locus": af_other_count,
            "coverage_at_locus": af_locus_coverage,
            "alt_allele_frequency": round(alt_freq, 4)
        }

    return node_id, final_variant_output


# ─────────────────────────────────────────────────────────────────────────────
# Pileup Viewing Function (Integrated)
# ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(pileup_data_for_node, node_id_str_for_display, max_reads_to_display_per_variant,
                        max_variants_to_display=float('inf')):
    if not pileup_data_for_node or not isinstance(pileup_data_for_node, dict):
        print(f"ℹ️ No valid pileup data to display for node {node_id_str_for_display}.", file=sys.stderr)
        return

    node_length = pileup_data_for_node.get("node_length")
    variants_dict = pileup_data_for_node.get("variants", {})

    print(
        f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {node_length if node_length is not None else 'N/A'}) ===")

    if not variants_dict:
        print(f"ℹ️ No variants found or pileups generated for this node.")
        return

    variants_displayed_count = 0
    sorted_variant_keys = sorted(variants_dict.keys(), key=lambda x: int(x.split('_')[0]))

    for variant_key in sorted_variant_keys:
        if variants_displayed_count >= max_variants_to_display:
            print(f"\n  ... (and {len(variants_dict) - variants_displayed_count} more variants not shown due to limit)")
            break

        variant_data = variants_dict[variant_key]
        pileup_matrix_indices = variant_data.get("pileup_matrix", [])

        print(f"\n--- Variant: {variant_key} ---")
        print(f"  Alt Count: {variant_data.get('alt_allele_count', 'N/A')}")
        print(f"  Ref Count: {variant_data.get('ref_allele_count_at_locus', 'N/A')}")
        print(f"  Other Count: {variant_data.get('other_allele_count_at_locus', 'N/A')}")
        print(f"  Coverage: {variant_data.get('coverage_at_locus', 'N/A')}")
        alt_freq_val = variant_data.get('alt_allele_frequency', 'N/A')
        if isinstance(alt_freq_val, float):
            print(f"  Alt Freq: {alt_freq_val:.4f}")
        else:
            print(f"  Alt Freq: {alt_freq_val}")

        if not pileup_matrix_indices:
            print("  (No reads in pileup matrix for this variant's window)")
        else:
            displayed_reads_count = 0
            for i, row_indices in enumerate(pileup_matrix_indices):
                if displayed_reads_count >= max_reads_to_display_per_variant:
                    print(
                        f"  ... (and {len(pileup_matrix_indices) - max_reads_to_display_per_variant} more reads for this variant's pileup window)")
                    break
                # Use INDEX_TO_BASE_FOR_VIEW which maps deletion index to '*'
                pileup_row_str = "".join(
                    [INDEX_TO_BASE_FOR_VIEW.get(idx, '?') if idx != BASE_TO_INDEX[' '] else ' ' for idx in row_indices])
                print(f"  Read {i + 1:3d}: {pileup_row_str}")
                displayed_reads_count += 1
        variants_displayed_count += 1
    print("\n")


# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered pileups for a single specified node and optionally view them.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("dat", help=".dat file path (read alignment data)")
    parser.add_argument("idx", help=".idx file path (index for .dat file)")
    parser.add_argument("output", help="JSON output file path for pileups")
    parser.add_argument("--node_id", type=int, required=True, help="The specific node ID to process.")
    parser.add_argument("--gfa", help="GFA graph file path (required if node sequence cache is not used/built).")
    parser.add_argument("--load-cache", help="Load node sequence from this JSON cache file.")
    parser.add_argument("--save-cache", help="Save node sequence to this JSON cache file (used if --gfa is provided).")

    parser.add_argument(
        "--view",
        nargs='?', const=-1, default=None, type=int, metavar='N',
        help="Print generated pileups to console. Optionally specify N to view the first N variants. If no N, all variants are shown."
    )
    parser.add_argument(
        "--max_view_reads", type=int, default=20,
        help="Maximum number of reads to display per pileup matrix in console view (if --view is used)."
    )

    args = parser.parse_args()

    if not os.path.isfile(args.dat): sys.exit(f"❌ Error: DAT file not found: {args.dat}")
    if not os.path.isfile(args.idx): sys.exit(f"❌ Error: Index file not found: {args.idx}")
    if not args.load_cache and not args.gfa: sys.exit("❌ Error: Must provide --gfa or --load-cache.")
    if args.load_cache and not os.path.isfile(args.load_cache): sys.exit(
        f"❌ Error: Cache file not found: {args.load_cache}")
    if args.gfa and not os.path.isfile(args.gfa): sys.exit(f"❌ Error: GFA file not found: {args.gfa}")
    if args.gfa and args.load_cache: print("🔹 Info: Both --gfa and --load-cache provided. Cache preferred.")

    target_node_id = args.node_id
    print(f"🔹 Processing single target node ID: {target_node_id}")

    start_time = time.time()
    node_dat_info = parse_idx_file_for_single_node(args.idx, target_node_id)
    if not node_dat_info: sys.exit(f"❌ Error: Failed to get index info for node {target_node_id}.")
    dat_offset, n_records = node_dat_info
    print(f"✔ Index parsing for node {target_node_id} took {time.time() - start_time:.2f}s.")

    node_sequence = None
    if args.load_cache and os.path.isfile(args.load_cache):
        start_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                node_sequence = json.load(cf).get(str(target_node_id))
            if node_sequence:
                print(f"✔ Loaded sequence for node {target_node_id} from cache in {time.time() - start_time:.2f}s.")
            else:
                print(f"⚠️ Warning: Node {target_node_id} not in cache {args.load_cache}.")
                if not args.gfa: sys.exit(f"❌ Error: Node {target_node_id} not in cache and no GFA. Exiting.")
        except Exception as e:
            print(f"❌ Error loading cache {args.load_cache}: {e}", file=sys.stderr)
            if not args.gfa: sys.exit(1)
            node_sequence = None

    if node_sequence is None and args.gfa:
        start_time = time.time()
        node_sequence = load_node_sequence_from_gfa(args.gfa, target_node_id)
        if node_sequence and args.save_cache:
            print(f"🔹 Saving sequence for node {target_node_id} to cache: {args.save_cache}...")
            try:
                existing_cache = {}
                if os.path.isfile(args.save_cache):
                    with open(args.save_cache, 'r') as rcf:
                        try:
                            existing_cache = json.load(rcf)
                        except json.JSONDecodeError:
                            print(f"⚠️ Corrupt cache {args.save_cache}, overwriting.", file=sys.stderr)
                existing_cache[str(target_node_id)] = node_sequence
                with open(args.save_cache, 'w') as wcf:
                    json.dump(existing_cache, wcf)
                print(f"✔ Saved sequence for node {target_node_id} to cache.")
            except Exception as e:
                print(f"❌ Error saving to cache {args.save_cache}: {e}", file=sys.stderr)
        elif node_sequence:
            print(f"✔ Sequence loading for node {target_node_id} from GFA took {time.time() - start_time:.2f}s.")

    if not node_sequence: sys.exit(f"❌ Error: Failed to obtain sequence for node {target_node_id}. Exiting.")

    node_len_for_output = len(node_sequence)

    task = (target_node_id, dat_offset, n_records, node_sequence)
    print(f"🔹 Prepared task for node {target_node_id}.")

    output_data_for_json = {}

    print(f"🔹 Processing node {target_node_id} using 1 worker...")
    start_proc_time = time.time()

    try:
        with ProcessPoolExecutor(max_workers=1, initializer=init_worker, initargs=(args.dat,)) as executor:
            future = executor.submit(process_single_node_for_pileup, task)
            processed_node_id, variants_dict_from_worker = future.result()

            if variants_dict_from_worker is not None:
                output_data_for_json[str(processed_node_id)] = {
                    "node_length": node_len_for_output,
                    "variants": variants_dict_from_worker  # This is the dict of variants
                }
            else:
                print(f"⚠️ Warning: Processing for node {processed_node_id} did not yield a variants dictionary.",
                      file=sys.stderr)
                output_data_for_json[str(processed_node_id)] = {
                    "node_length": node_len_for_output,
                    "variants": {}
                }

    except Exception as pool_exc:
        sys.exit(f"\n❌ An error occurred during processing node {target_node_id}: {pool_exc}")

    total_elapsed_time = time.time() - start_proc_time
    print(f"✔ Node {target_node_id} processing finished in {total_elapsed_time:.2f}s.")

    # Check if the node specific data was actually populated
    node_result_data = output_data_for_json.get(str(target_node_id))

    if node_result_data and node_result_data.get("variants"):
        print(f"🔹 Writing pileup results for node {target_node_id} to JSON output: {args.output}")
        start_write_time = time.time()
        try:
            with open(args.output, 'w') as out_f:
                json.dump(output_data_for_json, out_f, indent=2)
            print(f"✔ Output written in {time.time() - start_write_time:.2f}s. ✅ Pileup JSON saved to {args.output}")

            if args.view is not None:
                max_v_show = float('inf') if args.view == -1 else (args.view if args.view >= 0 else float('inf'))
                if args.view < -1: print("⚠️ Warning: Invalid number for --view. Showing all.", file=sys.stderr)

                view_msg = f"all variants (max {args.max_view_reads} reads/variant)..."
                if max_v_show != float(
                    'inf'): view_msg = f"first {int(max_v_show)} variants (max {args.max_view_reads} reads/variant)..."
                print(f"\n🔹 Displaying pileups: {view_msg}")
                # Pass the specific node's data to the display function
                display_pileup_data({str(target_node_id): node_result_data}, str(target_node_id), args.max_view_reads,
                                    max_v_show)

        except Exception as e:
            sys.exit(f"❌ Error writing/viewing output: {e}")
    else:
        print(
            f"ℹ️ No variant pileup data generated or processed for node {target_node_id}. Output file not written or is empty of variants.")

    print("✅ Script finished.")


if __name__ == '__main__':
    main()
