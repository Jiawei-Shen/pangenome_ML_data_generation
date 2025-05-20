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
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 4}  # '*' for deletions in pileup vis
INDEX_TO_BASE_FOR_VIEW = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N'}  # For viewing, '*' will show as 'N'

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
        sys.exit(1)  # Critical error
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
        sys.exit(1)  # Critical error


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
                if line_counter % 5_000_000 == 0:  # Progress update
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
                    break  # Found the target node, no need to parse further

            if line_counter > 0: print(f"✔ Finished GFA scan after {line_counter:,} lines.")
            if node_sequence is None:
                print(f"❌ Error: Sequence for node ID {target_node_id} not found in GFA file {gfa_path}.",
                      file=sys.stderr)

    except FileNotFoundError:
        print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr)
        sys.exit(1)  # Critical error
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr)
        sys.exit(1)  # Critical error
    return node_sequence


def decode_cigar(cigar_string):
    """Decodes a CIGAR string into a list of (length, operation) tuples."""
    if not cigar_string or cigar_string == '*':
        return []
    try:
        return re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)  # Common CIGAR ops
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr)
        return []


def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence, read_cigar_ops,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None, expected_ref_len_for_del=0):
    """
    Determines the allele presented by a read at a specific target_node_pos.
    Returns:
        - The base character if it's a match/mismatch at target_node_pos.
        - The inserted sequence string if an insertion occurs *after* target_node_pos.
        - A '*' character if a deletion covering target_node_pos occurs.
        - None if the read doesn't informatively cover the target_node_pos for the expected_var_type.
    """
    current_node_pos = read_offset_on_node
    current_read_pos = 0

    for length, op in read_cigar_ops:
        if op == 'M' or op == '=' or op == 'X':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                # Target position is within this aligned block
                if expected_var_type == 'I':  # If we expect an INS, a match here means REF state for INS
                    return "REF_STATE_FOR_INDEL"
                if expected_var_type == 'D':  # If we expect a DEL, a match here means REF state for DEL
                    return "REF_STATE_FOR_INDEL"
                # For SNP or general query
                offset_in_block = target_node_pos - current_node_pos
                return read_sequence[current_read_pos + offset_in_block]
            current_node_pos += length
            current_read_pos += length
        elif op == 'I':
            # Insertion in read occurs after current_node_pos-1 on reference
            # If target_node_pos is the anchor base *before* the insertion
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                return read_sequence[current_read_pos: current_read_pos + length]
            current_read_pos += length
        elif op == 'D':
            # Deletion from reference starts at current_node_pos and spans 'length' bases
            if current_node_pos <= target_node_pos < current_node_pos + length:
                # The target_node_pos is part of this deletion in the read
                if expected_var_type == 'I':  # If we expect an INS, a DEL here means OTHER state for INS
                    return "OTHER_FOR_INDEL"  # Or handle as non-informative for specific INS
                if expected_var_type == 'D':
                    # Check if this deletion matches the one we are assessing
                    # This requires comparing the length of deletion
                    if length == expected_ref_len_for_del:  # Simple check, could be more specific
                        return "*"  # Represents the deletion allele
                    else:
                        return "OTHER_FOR_INDEL"  # Different deletion
                return "*"  # General case: read shows a deletion at target_node_pos
            current_node_pos += length
        elif op == 'S':
            current_read_pos += length
        elif op == 'N':
            current_node_pos += length

        # If we've passed the target position without an informative operation
        if current_node_pos > target_node_pos + 1 and (op in ['M', '=', 'X', 'D', 'N']):  # Add some buffer for indels
            # For indels, target_node_pos can be an anchor, so check current_node_pos > target_node_pos
            if not (expected_var_type == 'I' and (current_node_pos - 1) <= target_node_pos):
                break  # Optimization: stop if we've passed the relevant part of the read

    return None  # Read does not informatively cover the position for the given type or is complex


def detect_variants_from_cigar(offset_on_node, cigar_string, read_sequence, node_sequence):
    """
    Detects variants (mismatches, insertions, deletions) based on CIGAR and sequence alignment.
    Returns a list of tuples: (position_on_node, variant_type, alt_allele, ref_allele)
    """
    variants = []
    node_pos = offset_on_node
    read_pos = 0
    cigar_ops = decode_cigar(cigar_string)

    for length_str, op in cigar_ops:
        try:
            length = int(length_str)
        except ValueError:
            print(f"⚠️ Warning: Invalid length in CIGAR operation '{length_str}{op}' in string '{cigar_string}'",
                  file=sys.stderr)
            continue

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
            inserted_sequence = read_sequence[read_pos: read_pos + length]
            # Anchor position for insertion is the base *before* it on the reference
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0  # Handle insertion at start of node alignment
            # Ref allele for an insertion is often represented by the base at anchor_pos or empty/placeholder
            # For simplicity, using '*' as ref for indel events.
            ref_base_at_anchor = node_sequence[ref_anchor_pos] if ref_anchor_pos < len(
                node_sequence) and ref_anchor_pos >= 0 else ""
            variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor if ref_base_at_anchor else "*"))
            read_pos += length
        elif op == 'D':
            deleted_sequence_from_ref = node_sequence[node_pos: node_pos + length]
            variants.append((node_pos, 'D', '*', deleted_sequence_from_ref))  # Alt is '*', Ref is deleted part
            node_pos += length
        elif op == 'S':
            read_pos += length
        elif op == 'N':
            node_pos += length
    return variants


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
    final_variant_data = {}  # Store final data including AF
    aligned_read_segments = []

    try:  # Read all segments for the node
        worker_dat_file.seek(dat_file_offset + 10)
        for record_idx in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE: break
            off, raw_seq, _, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)
            if mapq < 10: continue
            try:
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')
                cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                continue
            current_read_sequence = reverse_complement(seq) if strand_char == '-' else seq
            aligned_read_segments.append({
                "offset_on_node": off, "read_sequence": current_read_sequence,
                "cigar_string": cigar, "cigar_ops": decode_cigar(cigar),  # Pre-decode CIGAR
                "original_strand": strand_char
            })
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] reading records for node {node_id}: {e}", file=sys.stderr)
        return node_id, {}

    # 1. Detect all unique variants defined by reads
    # Key: (pos, type, ref, alt), Value: list of read_info for pileup matrix
    raw_variant_occurrences = defaultdict(list)
    for segment in aligned_read_segments:
        variants_in_read = detect_variants_from_cigar(
            segment["offset_on_node"], segment["cigar_string"],
            segment["read_sequence"], node_sequence
        )
        for v_pos, v_type, v_alt, v_ref in variants_in_read:
            raw_variant_occurrences[(v_pos, v_type, v_ref, v_alt)].append({
                "read_offset_on_node": segment["offset_on_node"],
                "read_sequence": segment["read_sequence"]
            })

    # 2. For each unique variant, calculate AF and build pileup
    window_size = 60
    half_window = window_size // 2

    for (v_pos, v_type, v_ref_defined, v_alt_defined), reads_defining_alt in raw_variant_occurrences.items():
        alt_allele_count = len(reads_defining_alt)
        ref_allele_count = 0
        other_allele_count = 0
        coverage_at_locus = 0

        expected_ref_len_for_del = len(v_ref_defined) if v_type == 'D' else 0

        for segment in aligned_read_segments:  # Iterate all reads covering the locus
            # Check if this read informs the locus v_pos for v_type
            # This requires checking if the read spans v_pos and what allele it presents

            # Simplified check: does read span v_pos?
            # A more precise check would use CIGAR to see if v_pos is matched/mismatched or part of an indel
            read_start_on_node = segment["offset_on_node"]
            read_end_on_node = read_start_on_node  # Needs CIGAR to determine actual end

            # Quick CIGAR parse to find read span on reference
            temp_node_pos_for_span = segment["offset_on_node"]
            for l_span, op_span in segment["cigar_ops"]:
                if op_span in ('M', '=', 'X', 'D', 'N'):
                    temp_node_pos_for_span += l_span
            read_end_on_node = temp_node_pos_for_span

            # Check if v_pos (or anchor for insertion) is covered
            locus_covered = False
            if v_type == 'I':  # v_pos is anchor, insertion is after v_pos
                # Read needs to span the base before v_pos and the base after v_pos (or end of read)
                # to inform about presence/absence of insertion
                if read_start_on_node <= v_pos < read_end_on_node:  # Simplified: covers anchor
                    locus_covered = True
            else:  # SNP or Deletion, v_pos is the start of the event
                if read_start_on_node <= v_pos < read_end_on_node:
                    locus_covered = True

            if not locus_covered:
                continue

            # If covered, determine what this read supports at v_pos for v_type
            coverage_at_locus += 1

            # Does this read define the exact (v_pos, v_type, v_ref_defined, v_alt_defined)?
            # We know it does if it's in reads_defining_alt, but those are already counted in alt_allele_count.
            # We need to check reads NOT in reads_defining_alt for their support of REF or OTHER.

            # Check if this read *also* generated the current unique_variant
            # This is a bit tricky because a read can generate multiple variants.
            # We are interested if *this specific variant being iterated* is present in the read.

            is_this_read_defining_current_alt = False
            # This check is not perfect because reads_defining_alt contains dicts, not segments.
            # For simplicity, assume if a read is NOT among those that defined the ALT, it's either REF or OTHER.
            # A more robust way is to re-evaluate the read against this specific variant.

            # Re-evaluate what this read presents for the *specific variant* (v_pos, v_type, v_ref_defined, v_alt_defined)
            read_allele_type = "OTHER"  # Default

            # Simplified logic for now: if it's not one of the reads that *defined* this alt,
            # assume it's ref or other. This is an approximation.
            # A full re-evaluation of each read against each variant type is needed for perfect AF.

            # For now, let's use a proxy:
            # If this read is not in `reads_defining_alt` (hard to check directly with current structure)
            # we assume it's either REF or OTHER.
            # This part needs significant refinement for accurate ref/other counts for indels.

            # Let's count ref based on reads that *don't* have this specific variant
            # but *do* cover the site and match the reference state.
            # This is the most complex part.

            # Simplified AF: alt_count is known. Coverage is reads spanning the site.
            # Ref_count = coverage - alt_count - other_count (where other is not easy yet)

            # Let's use the get_allele_from_read_at_node_pos for a more direct check
            allele_in_read = get_allele_from_read_at_node_pos(
                segment["offset_on_node"], segment["read_sequence"], segment["cigar_ops"],
                v_pos, node_sequence, v_type, expected_ref_len_for_del
            )

            if allele_in_read is None:  # Read doesn't inform this specific locus/type well
                coverage_at_locus -= 1  # Don't count it in coverage for AF
                continue

            # Check if this read supports the v_alt_defined
            supported_alt_in_this_read = False
            if v_type == 'X' and allele_in_read == v_alt_defined:
                supported_alt_in_this_read = True
            elif v_type == 'I' and allele_in_read == v_alt_defined:  # allele_in_read is inserted seq
                supported_alt_in_this_read = True
            elif v_type == 'D' and allele_in_read == "*":  # and v_ref_defined matches deletion from read
                # This check needs to be more robust for specific deleted sequence if needed
                supported_alt_in_this_read = True  # Assuming '*' means it supports *a* deletion here

            if supported_alt_in_this_read:
                # This read contributes to alt_allele_count, which is already len(reads_defining_alt)
                # If this read was NOT in reads_defining_alt but still supports it, this logic is flawed.
                # Assuming reads_defining_alt correctly captures all reads for this exact alt.
                pass
            elif allele_in_read == "REF_STATE_FOR_INDEL" or \
                    (v_type == 'X' and allele_in_read == v_ref_defined):
                ref_allele_count += 1
            else:  # Supports neither the defined alt nor the defined ref
                other_allele_count += 1

        # Adjust coverage: it's sum of alt (already counted), ref, and other
        # The initial coverage_at_locus was just reads spanning.
        # More accurate coverage for AF is alt_allele_count + ref_allele_count + other_allele_count
        effective_coverage = alt_allele_count + ref_allele_count + other_allele_count
        alt_freq = alt_allele_count / effective_coverage if effective_coverage > 0 else 0.0

        # Build pileup matrix from reads that defined this alternate variant
        pileup_matrix = np.full((len(reads_defining_alt), window_size), BASE_TO_INDEX['N'], dtype=np.uint8)
        for i, read_info in enumerate(reads_defining_alt):
            read_offset_on_node = read_info["read_offset_on_node"]
            read_seq = read_info["read_sequence"]
            read_len = len(read_seq)
            for j in range(window_size):
                window_pos_on_node = (v_pos - half_window) + j
                read_char_idx = window_pos_on_node - read_offset_on_node
                if 0 <= read_char_idx < read_len:
                    base = read_seq[read_char_idx].upper()
                    pileup_matrix[i, j] = BASE_TO_INDEX.get(base, BASE_TO_INDEX['N'])

        variant_key_str = f"{v_pos}_{v_type}_{v_ref_defined}_{v_alt_defined}"
        final_variant_data[variant_key_str] = {
            "pileup_matrix": pileup_matrix.tolist(),
            "alt_allele_count": alt_allele_count,
            "ref_allele_count_at_locus": ref_allele_count,
            "other_allele_count_at_locus": other_allele_count,
            "coverage_at_locus": effective_coverage,
            "alt_allele_frequency": round(alt_freq, 4)  # Round for display
        }

    return node_id, final_variant_data


# ─────────────────────────────────────────────────────────────────────────────
# Pileup Viewing Function (Integrated)
# ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(pileup_data_dict, max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if not pileup_data_dict:
        print("ℹ️ No pileup data to display.", file=sys.stderr)
        return

    for node_id_str, variants_dict in pileup_data_dict.items():
        if not isinstance(variants_dict, dict):
            print(f"⚠️ Warning: Data for node '{node_id_str}' is not in the expected format. Skipping display.",
                  file=sys.stderr)
            continue
        if not variants_dict:
            print(f"ℹ️ No variants found or pileups generated for node ID: {node_id_str}")
            continue

        print(f"\n=== Displaying Pileups for Node ID: {node_id_str} ===")
        variants_displayed_count = 0
        # Sort variants by position (first part of the key)
        sorted_variant_keys = sorted(variants_dict.keys(), key=lambda x: int(x.split('_')[0]))

        for variant_key in sorted_variant_keys:
            if variants_displayed_count >= max_variants_to_display:
                print(
                    f"\n  ... (and {len(variants_dict) - variants_displayed_count} more variants not shown due to limit)")
                break

            variant_data = variants_dict[variant_key]
            pileup_matrix_indices = variant_data.get("pileup_matrix", [])

            print(f"\n--- Variant: {variant_key} ---")
            print(f"  Alt Count: {variant_data.get('alt_allele_count', 'N/A')}")
            print(f"  Ref Count: {variant_data.get('ref_allele_count_at_locus', 'N/A')}")
            print(f"  Other Count: {variant_data.get('other_allele_count_at_locus', 'N/A')}")
            print(f"  Coverage: {variant_data.get('coverage_at_locus', 'N/A')}")
            print(f"  Alt Freq: {variant_data.get('alt_allele_frequency', 'N/A'):.4f}")

            if not pileup_matrix_indices:
                print("  (No reads in pileup matrix for this variant)")
            else:
                displayed_reads_count = 0
                for i, row_indices in enumerate(pileup_matrix_indices):
                    if displayed_reads_count >= max_reads_to_display_per_variant:
                        print(
                            f"  ... (and {len(pileup_matrix_indices) - max_reads_to_display_per_variant} more reads for this variant's pileup matrix)")
                        break
                    pileup_row_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in row_indices])
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

    task = (target_node_id, dat_offset, n_records, node_sequence)
    print(f"🔹 Prepared task for node {target_node_id}.")
    results = {}
    print(f"🔹 Processing node {target_node_id} using 1 worker...")
    start_proc_time = time.time()

    try:
        with ProcessPoolExecutor(max_workers=1, initializer=init_worker, initargs=(args.dat,)) as executor:
            future = executor.submit(process_single_node_for_pileup, task)
            processed_node_id, pileup_dict = future.result()
            if pileup_dict is not None:
                results[str(processed_node_id)] = pileup_dict
            else:
                print(f"⚠️ Warning: Processing for node {processed_node_id} did not yield results.", file=sys.stderr)
    except Exception as pool_exc:
        sys.exit(f"\n❌ An error occurred during processing node {target_node_id}: {pool_exc}")

    total_elapsed_time = time.time() - start_proc_time
    print(f"✔ Node {target_node_id} processing finished in {total_elapsed_time:.2f}s.")

    if results:
        print(f"🔹 Writing pileup results for node {target_node_id} to JSON output: {args.output}")
        start_write_time = time.time()
        try:
            with open(args.output, 'w') as out_f:
                json.dump(results, out_f, indent=2)
            print(f"✔ Output written in {time.time() - start_write_time:.2f}s. ✅ Pileup JSON saved to {args.output}")
            if args.view is not None:
                max_v_show = float('inf') if args.view == -1 else (args.view if args.view >= 0 else float('inf'))
                if args.view < -1: print("⚠️ Warning: Invalid number for --view. Showing all.", file=sys.stderr)

                view_msg = f"all variants (max {args.max_view_reads} reads/variant)..."
                if max_v_show != float(
                    'inf'): view_msg = f"first {int(max_v_show)} variants (max {args.max_view_reads} reads/variant)..."
                print(f"\n🔹 Displaying pileups: {view_msg}")
                display_pileup_data(results, args.max_view_reads, max_v_show)
        except Exception as e:
            sys.exit(f"❌ Error writing/viewing output: {e}")
    else:
        print(f"ℹ️ No pileup data generated for node {target_node_id}. Output file not written.")
    print("✅ Script finished.")


if __name__ == '__main__':
    main()
