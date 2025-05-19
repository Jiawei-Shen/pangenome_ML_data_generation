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
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 4}  # Added '*' for deletions in pileup vis

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


def detect_variants_from_cigar(offset_on_node, cigar_string, read_sequence, node_sequence):
    """
    Detects variants (mismatches, insertions, deletions) based on CIGAR and sequence alignment.
    Returns a list of tuples: (position_on_node, variant_type, alt_allele, ref_allele)
    For insertions, alt_allele is the inserted sequence, ref_allele is '*'.
    For deletions, alt_allele is '*', ref_allele is the deleted sequence from node.
    For mismatches, alt_allele is read base, ref_allele is node base.
    """
    variants = []
    node_pos = offset_on_node  # Current 0-based position on the node sequence
    read_pos = 0  # Current 0-based position on the read sequence
    cigar_ops = decode_cigar(cigar_string)

    for length_str, op in cigar_ops:
        try:
            length = int(length_str)
        except ValueError:
            print(f"⚠️ Warning: Invalid length in CIGAR operation '{length_str}{op}' in string '{cigar_string}'",
                  file=sys.stderr)
            continue

        if op == 'M' or op == '=' or op == 'X':  # Match, Sequence Match, or Mismatch
            for i in range(length):
                current_node_pos = node_pos + i
                current_read_pos = read_pos + i
                if current_node_pos < len(node_sequence) and current_read_pos < len(read_sequence):
                    node_base = node_sequence[current_node_pos].upper()
                    read_base = read_sequence[current_read_pos].upper()
                    if node_base != read_base and op != '=':  # Mismatch (M or X)
                        variants.append((current_node_pos, 'X', read_base, node_base))
                # else: alignment extends beyond sequence, ignore
            node_pos += length
            read_pos += length
        elif op == 'I':  # Insertion to the reference (node)
            # Insertion occurs *after* node_pos-1 on the node.
            # The inserted sequence is from read_pos to read_pos + length - 1 in the read.
            inserted_sequence = read_sequence[read_pos: read_pos + length]
            # Position of insertion is typically the base *before* the insertion on the reference.
            # If node_pos is 0, it's an insertion at the beginning.
            ref_anchor_pos = node_pos - 1  # Base before insertion
            variants.append((ref_anchor_pos, 'I', inserted_sequence, '*'))
            read_pos += length  # Consumes read bases
            # Does not consume node bases
        elif op == 'D':  # Deletion from the reference (node)
            # Deletion starts *at* node_pos on the node.
            # The deleted sequence is from node_pos to node_pos + length - 1 on the node.
            deleted_sequence = node_sequence[node_pos: node_pos + length]
            variants.append((node_pos, 'D', '*', deleted_sequence))
            node_pos += length  # Consumes node bases
            # Does not consume read bases
        elif op == 'S':  # Soft clipping
            read_pos += length  # Consumes read bases
        # Other CIGAR ops like H (hard clipping), P (padding), N (skipped region) are ignored for basic pileup
        # H and P do not consume read or reference for alignment purposes here. N consumes reference.
        elif op == 'N':  # Skipped region from reference
            node_pos += length

    return variants


# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function

def init_worker(dat_file_path_for_worker):
    """Initializer for each worker process: opens the .dat file."""
    global worker_dat_file
    # print(f"[Worker {os.getpid()}] Initializing and opening {dat_file_path_for_worker}")
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
    except FileNotFoundError:
        print(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path_for_worker}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path_for_worker}: {e}", file=sys.stderr)
        sys.exit(1)


def process_single_node_for_pileup(task_args):
    """
    Function executed by the worker process for the single target node.
    Processes the node to find variants and generate pileups.
    """
    node_id, dat_file_offset, n_records, node_sequence = task_args
    global worker_dat_file

    if worker_dat_file is None:
        print(f"❌ Error [Worker {os.getpid()}]: Worker DAT file handle not initialized for node {node_id}.",
              file=sys.stderr)
        return node_id, {}

    if not node_sequence:
        print(f"⚠️ Warning [Worker {os.getpid()}]: Empty sequence provided for node {node_id}. Cannot generate pileup.",
              file=sys.stderr)
        return node_id, {}

    node_len = len(node_sequence)
    # segments = [] # Store (offset_on_node, read_sequence, cigar_string)

    # Pileup matrix: rows = positions on node, columns = A, C, G, T, N, Ins, Del
    # For simplicity, let's create a dictionary of pileups for each variant position.
    # pileups_at_variant_sites = defaultdict(lambda: np.zeros(5, dtype=np.uint32)) # A,C,G,T,N counts

    # More detailed pileup: for each variant, store a window of aligned read bases
    variant_pileups = {}  # Key: "pos_type", Value: list of lists (pileup matrix)

    # Store all read segments that align to this node
    aligned_read_segments = []

    try:
        # The first 10 bytes of a block in .dat are metadata (node_id_block, offset_block, block_size_block)
        # The actual alignment records start after this.
        # The offset from .idx file is the start of the block for this node.
        worker_dat_file.seek(dat_file_offset + 10)  # Skip block header

        for record_idx in range(n_records):
            try:
                data = worker_dat_file.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE:
                    print(
                        f"⚠️ Warning [Worker {os.getpid()}]: Short read ({len(data)} bytes) for node {node_id}, record {record_idx + 1}/{n_records}. Stopping.",
                        file=sys.stderr)
                    break

                off, raw_seq, _, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)

                if mapq < 10:  # Basic MAPQ filter
                    continue

                try:
                    seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')  # Replace invalid bytes
                    cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                    strand_char = strand_byte.decode('ascii')
                except UnicodeDecodeError as ude:
                    print(
                        f"⚠️ Warning [Worker {os.getpid()}]: Unicode decode error in record {record_idx + 1} for node {node_id}: {ude}. Skipping.",
                        file=sys.stderr)
                    continue

                read_len_from_seq = len(seq)  # Actual length of decoded sequence

                # Adjust offset and reverse complement sequence for reverse strand reads
                current_offset_on_node = off
                current_read_sequence = seq

                if strand_char == '-':
                    # The offset 'off' from the .dat file for a reverse strand read
                    # is typically the start position of the alignment on the FORWARD strand
                    # representation of the node.
                    # To map this to the reverse complemented node sequence (if we were to use one),
                    # or to correctly place the reverse complemented read on the forward node sequence,
                    # we need to adjust.
                    # If the node_sequence is always forward, and we reverse_complement the read,
                    # the 'off' still refers to the start on the forward node.
                    # The original script had: adj_off = node_len - off - read_len_from_seq
                    # This adj_off would be the start of the revcomp read on a revcomp node.
                    # Let's keep node_sequence as forward. If read is '-', revcomp it.
                    # The 'off' from .dat should be the 0-based start on the forward node.
                    current_read_sequence = reverse_complement(seq)
                    # 'off' remains the start position on the forward reference node.

                # Store the processed read segment for later pileup generation
                aligned_read_segments.append({
                    "offset_on_node": current_offset_on_node,
                    "read_sequence": current_read_sequence,  # Already reverse_complemented if needed
                    "cigar_string": cigar,
                    "mapq": mapq,
                    "strand": strand_char  # Original strand
                })

            except struct.error as se:
                print(
                    f"❌ Error [Worker {os.getpid()}]: Failed to unpack record {record_idx + 1} for node {node_id}: {se}. Stopping.",
                    file=sys.stderr)
                break
            except Exception as e_inner:
                print(
                    f"❌ Error [Worker {os.getpid()}]: Unexpected error processing record {record_idx + 1} for node {node_id}: {e_inner}",
                    file=sys.stderr)
                continue

    except IOError as ioe:
        print(
            f"❌ Error [Worker {os.getpid()}]: I/O error reading data for node {node_id} at offset {dat_file_offset}: {ioe}",
            file=sys.stderr)
        return node_id, {}
    except Exception as e_outer:
        print(f"❌ Error [Worker {os.getpid()}]: Unexpected error seeking/reading node {node_id}: {e_outer}",
              file=sys.stderr)
        return node_id, {}

    # --- Variant Detection and Pileup Generation from collected segments ---
    # First, identify all unique variant sites from all reads aligned to this node
    all_variant_sites = defaultdict(lambda: {"count": 0, "reads_info": []})  # Store reads supporting each variant

    for segment in aligned_read_segments:
        # Variants are (pos_on_node, type, alt, ref)
        variants_in_read = detect_variants_from_cigar(
            segment["offset_on_node"],
            segment["cigar_string"],
            segment["read_sequence"],  # Use the (potentially revcomped) read sequence
            node_sequence  # Use the forward node sequence
        )
        for v_pos, v_type, v_alt, v_ref in variants_in_read:
            variant_key = f"{v_pos}_{v_type}_{v_ref}_{v_alt}"  # More specific key
            all_variant_sites[variant_key]["count"] += 1
            # Store read info for pileup generation around this variant
            all_variant_sites[variant_key]["reads_info"].append({
                "read_offset_on_node": segment["offset_on_node"],
                "read_sequence": segment["read_sequence"]
                # "original_strand": segment["strand"] # If needed for strand bias, etc.
            })

    # Now, generate pileup windows around each identified variant site
    window_size = 60  # Characters in the pileup window (e.g., 30 left, variant, 30 right)
    half_window = window_size // 2

    for variant_key, data in all_variant_sites.items():
        parts = variant_key.split('_')
        v_pos = int(parts[0])  # 0-based position of the variant on the node
        # v_type = parts[1]
        # v_ref = parts[2]
        # v_alt = parts[3]

        # Pileup matrix for this variant: (num_supporting_reads, window_size)
        # Filled with index of 'N' initially
        pileup_matrix = np.full((len(data["reads_info"]), window_size), BASE_TO_INDEX['N'], dtype=np.uint8)

        for i, read_info in enumerate(data["reads_info"]):
            read_offset_on_node = read_info["read_offset_on_node"]
            read_seq = read_info["read_sequence"]
            read_len = len(read_seq)

            # Determine the portion of the read that aligns to the pileup window
            # The pileup window is centered around v_pos on the node sequence
            # Window start on node: v_pos - half_window + 1 (if v_pos is the variant itself)
            # Or, if v_pos is the first base of a multi-base variant:
            # Let's center the window around v_pos.

            # For each position in the window (0 to window_size - 1)
            for j in range(window_size):
                # Position in the window relative to the center (v_pos)
                # window_pos_on_node is the 0-based coordinate on the node sequence
                # for the j-th base of the pileup window.
                window_pos_on_node = (v_pos - half_window) + j

                # Corresponding position in the current read
                # read_char_idx = (position on node) - (read's start offset on node)
                read_char_idx = window_pos_on_node - read_offset_on_node

                if 0 <= read_char_idx < read_len:
                    base = read_seq[read_char_idx].upper()
                    pileup_matrix[i, j] = BASE_TO_INDEX.get(base, BASE_TO_INDEX['N'])
                # else: it's outside the read, remains 'N'

        variant_pileups[variant_key] = pileup_matrix.tolist()  # Convert numpy array to list for JSON

    # print(f"[Worker {os.getpid()}] Finished node {node_id}. Found {len(variant_pileups)} variant sites with pileups.")
    return node_id, variant_pileups


# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered pileups for a single specified node.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("dat", help=".dat file path (read alignment data)")
    parser.add_argument("idx", help=".idx file path (index for .dat file)")
    parser.add_argument("output", help="JSON output file path for pileups")
    parser.add_argument("--node_id", type=int, required=True, help="The specific node ID to process.")
    parser.add_argument("--gfa", help="GFA graph file path (required if node sequence cache is not used/built).")
    parser.add_argument("--load-cache", help="Load node sequence from this JSON cache file.")
    parser.add_argument("--save-cache", help="Save node sequence to this JSON cache file (used if --gfa is provided).")
    # Workers argument is kept for consistency but will effectively be 1 for a single node.
    # parser.add_argument("-w", "--workers", type=int, default=1, help="Number of worker processes (should be 1 for single node).")

    args = parser.parse_args()

    # --- Input Validation ---
    if not os.path.isfile(args.dat):
        print(f"❌ Error: DAT file not found: {args.dat}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.idx):
        print(f"❌ Error: Index file not found: {args.idx}", file=sys.stderr)
        sys.exit(1)

    if not args.load_cache and not args.gfa:
        print(
            "❌ Error: You must provide either a GFA file (`--gfa`) or a sequence cache (`--load-cache`) for the target node.",
            file=sys.stderr)
        sys.exit(1)
    if args.load_cache and not os.path.isfile(args.load_cache):
        print(f"❌ Error: Specified cache file to load does not exist: {args.load_cache}", file=sys.stderr)
        sys.exit(1)
    if args.gfa and not os.path.isfile(args.gfa):
        print(f"❌ Error: GFA file not found: {args.gfa}", file=sys.stderr)
        sys.exit(1)
    if args.gfa and args.load_cache:
        print("🔹 Info: Both --gfa and --load-cache provided. Cache will be preferred if node sequence is found there.")

    target_node_id = args.node_id
    print(f"🔹 Processing single target node ID: {target_node_id}")

    # --- Load Index for the Single Node ---
    print(f"🔹 Parsing index file for node {target_node_id}...")
    start_time = time.time()
    node_dat_info = parse_idx_file_for_single_node(args.idx, target_node_id)
    if not node_dat_info:
        print(f"❌ Error: Failed to get index information for node {target_node_id}. Exiting.", file=sys.stderr)
        sys.exit(1)
    dat_offset, n_records = node_dat_info
    print(f"✔ Index parsing for node {target_node_id} took {time.time() - start_time:.2f} seconds.")

    # --- Load Sequence for the Single Node ---
    node_sequence = None
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"🔹 Loading node sequence for {target_node_id} from cache: {args.load_cache}...")
        start_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                cached_sequences = json.load(cf)
                # Keys in JSON cache are strings
                node_sequence = cached_sequences.get(str(target_node_id))
            if node_sequence:
                print(
                    f"✔ Loaded sequence for node {target_node_id} from cache in {time.time() - start_time:.2f} seconds.")
            else:
                print(f"⚠️ Warning: Node {target_node_id} not found in cache {args.load_cache}.")
                if not args.gfa:  # If GFA is not provided as fallback
                    print(f"❌ Error: Node {target_node_id} not in cache and no GFA provided. Exiting.", file=sys.stderr)
                    sys.exit(1)
        except Exception as e:
            print(f"❌ Error loading cache file {args.load_cache}: {e}", file=sys.stderr)
            # Potentially fallback to GFA if provided
            if not args.gfa: sys.exit(1)
            node_sequence = None  # Ensure it's reset if cache load failed

    if node_sequence is None and args.gfa:  # Fallback to GFA or primary source
        print(f"🔹 Loading node sequence for {target_node_id} from GFA: {args.gfa}...")
        start_time = time.time()
        node_sequence = load_node_sequence_from_gfa(args.gfa, target_node_id)
        if node_sequence and args.save_cache:
            print(f"🔹 Saving node sequence for {target_node_id} to cache: {args.save_cache}...")
            try:
                # Load existing cache if it exists to update it, otherwise create new
                existing_cache_data = {}
                if os.path.isfile(args.save_cache):
                    with open(args.save_cache, 'r') as scf_read:
                        try:
                            existing_cache_data = json.load(scf_read)
                        except json.JSONDecodeError:
                            print(f"⚠️ Warning: Existing cache file {args.save_cache} is corrupted. Overwriting.",
                                  file=sys.stderr)

                existing_cache_data[str(target_node_id)] = node_sequence  # Add/update current node

                with open(args.save_cache, 'w') as scf_write:
                    json.dump(existing_cache_data, scf_write)
                print(f"✔ Saved sequence for node {target_node_id} to cache.")
            except Exception as e:
                print(f"❌ Error saving sequence for node {target_node_id} to cache {args.save_cache}: {e}",
                      file=sys.stderr)
        elif node_sequence:
            print(f"✔ Sequence loading for node {target_node_id} from GFA took {time.time() - start_time:.2f} seconds.")

    if not node_sequence:
        print(f"❌ Error: Failed to obtain sequence for target node {target_node_id}. Exiting.", file=sys.stderr)
        sys.exit(1)

    # --- Prepare Single Task ---
    task = (target_node_id, dat_offset, n_records, node_sequence)
    print(f"🔹 Prepared task for node {target_node_id}.")

    # --- Execute Task (using ProcessPoolExecutor with 1 worker for consistency with init_worker) ---
    results = {}
    print(f"🔹 Processing node {target_node_id} using 1 worker...")
    start_proc_time = time.time()

    try:
        # Using max_workers=1 as we are processing a single node.
        # The initializer will open the .dat file for this single worker.
        with ProcessPoolExecutor(max_workers=1,
                                 initializer=init_worker,
                                 initargs=(args.dat,)) as executor:
            # Submit the single task
            future = executor.submit(process_single_node_for_pileup, task)
            # Get the result (this will block until the task is done)
            processed_node_id, pileup_dict = future.result()

            if pileup_dict is not None:
                results[str(processed_node_id)] = pileup_dict  # Store with string key for JSON
            else:
                print(f"⚠️ Warning: Processing for node {processed_node_id} did not yield results.", file=sys.stderr)

    except Exception as pool_exc:
        print(f"\n❌ An error occurred during processing node {target_node_id}: {pool_exc}", file=sys.stderr)
        sys.exit(1)  # Exit if processing fails for the single node

    total_elapsed_time = time.time() - start_proc_time
    print(f"✔ Node {target_node_id} processing finished in {total_elapsed_time:.2f} seconds.")

    # --- Write Output ---
    if results:
        print(f"🔹 Writing pileup results for node {target_node_id} to JSON output: {args.output}")
        start_write_time = time.time()
        try:
            with open(args.output, 'w') as out_f:
                json.dump(results, out_f, indent=2)  # results will contain one key: target_node_id (as string)
            write_elapsed_time = time.time() - start_write_time
            print(f"✔ Output written in {write_elapsed_time:.2f} seconds.")
            print(f"✅ Done. Output saved to {args.output}")
        except IOError as ioe:
            print(f"❌ Error writing output JSON to {args.output}: {ioe}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"❌ Unexpected error writing output JSON: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print(
            f"ℹ️ No pileup data generated for node {target_node_id} (this might be expected if no variants or reads). Output file not written if empty.")


if __name__ == '__main__':
    main()
