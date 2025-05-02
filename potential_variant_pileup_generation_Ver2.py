#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np
from collections import defaultdict
import re # Import re at the top level
from concurrent.futures import ProcessPoolExecutor
# from os import cpu_count # Uncomment if using cpu_count for default workers

# ─────────────────────────────────────────────────────────────────────────────
# Constants
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4} # N used for padding/unknown

# Global for worker process state (file handle)
worker_dat_file = None

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions

def reverse_complement(sequence):
    """Computes the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def parse_idx_file(idx_path):
    """Parses the .idx file to get node offsets and record counts."""
    node_index = {}
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4:
                 print(f"❌ Error: Index file {idx_path} is too small.", file=sys.stderr)
                 sys.exit(1)

            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4:
                print(f"❌ Error: Could not read number of nodes from {idx_path}.", file=sys.stderr)
                sys.exit(1)
            num_nodes = struct.unpack('<I', num_nodes_bytes)[0]

            expected_min_size = 4 + num_nodes * 22
            if file_size < expected_min_size:
                print(f"❌ Error: Index file {idx_path} appears truncated. Expected at least {expected_min_size} bytes for {num_nodes} nodes, found {file_size}.", file=sys.stderr)
                # Optionally continue if partial processing is desired
                # sys.exit(1)

            print(f"🔹 Reading index for {num_nodes} nodes...")
            for i in range(num_nodes):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                     print(f"❌ Error: Index file ended prematurely while reading record {i+1}/{num_nodes}.", file=sys.stderr)
                     break
                # <I Q I I H => node_id, offset, block_size, n_records, padding
                node_id, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                if n_records > 0: # Only store nodes that have records
                    node_index[node_id] = (offset, n_records)

        print(f"✔ Parsed {len(node_index)} nodes with records from index file.")
        if len(node_index) != num_nodes and file_size >= expected_min_size:
             print(f"⚠️ Warning: Index header indicated {num_nodes} nodes, but {num_nodes - len(node_index)} had 0 records or were missing. Processing {len(node_index)} nodes.", file=sys.stderr)
        return node_index
    except FileNotFoundError:
        print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
        sys.exit(1)


def load_node_sequences_from_gfa(gfa_path, target_node_ids):
    """Loads sequences for specified node IDs from a GFA file."""
    node_sequences = {}
    target_node_set = set(target_node_ids)
    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file: {gfa_path}")
            line_counter = 0
            parsed_count = 0
            start_time = time.time()
            for line in f:
                line_counter += 1
                if line_counter % 20_000_000 == 0: # Print progress less often
                    elapsed = time.time() - start_time
                    print(f"  Checked {line_counter:,} lines in GFA file... ({elapsed:.1f}s)", end='\r')

                if not line.startswith('S\t'):
                    continue

                parts = line.strip().split('\t')
                # GFA S line: S <Name> <Sequence> [Optional tags]
                if len(parts) < 3:
                    continue

                try:
                    # Node ID might not always be numeric in GFA, but assume integer based on .idx
                    nid = int(parts[1])
                except ValueError:
                    continue # Skip nodes with non-integer IDs if index uses integers

                if nid in target_node_set:
                    node_sequences[nid] = parts[2]
                    parsed_count += 1
                    # Optimization: Remove found node ID from set
                    target_node_set.remove(nid)
                    if not target_node_set:
                       print(f"\n✔ Found all {len(target_node_ids)} target sequences after checking {line_counter:,} lines.")
                       break # Stop reading GFA if all targets found

            print(f"\n✔ Checked {line_counter:,} total lines in GFA.")
            print(f"✔ Loaded {len(node_sequences)} target sequences from GFA.")
            if target_node_set: # Check if any targets were *not* found
                 print(f"⚠️ Warning: Could not find sequences for {len(target_node_set)} target nodes in the GFA file (e.g., node {next(iter(target_node_set))}).", file=sys.stderr)

    except FileNotFoundError:
        print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr)
        sys.exit(1)

    return node_sequences

def decode_cigar(cigar_string):
    """Decodes a CIGAR string into a list of (length, operation) tuples."""
    if not cigar_string or cigar_string == '*':
        return []
    try:
        # Allow M, I, D, X, = (standard operations involved in alignment differences)
        # Others like S, H, P, N are ignored by findall
        return re.findall(r'(\d+)([MIDX=])', cigar_string)
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr)
        return []


def detect_variants_from_cigar(offset, cigar_string, node_len):
    """
    Identifies potential variant positions relative to the node sequence based on CIGAR.
    Returns a list of tuples: (position_on_node, variant_type).
    Variant Types: 'X' (mismatch), 'I' (insertion *after* position), 'D' (deletion *at* position).
    """
    variants = []
    ref_pos = offset # Current position on the reference node sequence (0-based)
    cigar_ops = decode_cigar(cigar_string)

    for length_str, op in cigar_ops:
        try:
            length = int(length_str)
        except ValueError:
            print(f"⚠️ Warning: Invalid length in CIGAR operation '{length_str}{op}' in string '{cigar_string}'", file=sys.stderr)
            continue

        if op == 'M' or op == '=': # Match or sequence match
            ref_pos += length
        elif op == 'X': # Sequence mismatch
            # Report variant site for each mismatch position
            for i in range(length):
                current_pos = ref_pos + i
                if 0 <= current_pos < node_len: # Check boundary
                   variants.append((current_pos, 'X'))
            ref_pos += length
        elif op == 'I': # Insertion to the reference
            # Insertion occurs *after* the reference position `ref_pos - 1`
            # We need to associate it with a valid reference position.
            # Conventionally, it's often linked to the base *before* the insertion.
            current_pos = ref_pos - 1
            if current_pos >= -1 and current_pos < node_len: # Allow pos -1 if insertion is at the beginning
                 variants.append((current_pos, 'I'))
            # Insertions do not consume reference bases
        elif op == 'D': # Deletion from the reference
            # Deletion starts *at* the reference position `ref_pos`
            # Report variant site for the start of the deletion
            if 0 <= ref_pos < node_len: # Check boundary for start pos
                 variants.append((ref_pos, 'D'))
            ref_pos += length # Deletions consume reference bases

        # Check if current ref_pos has exceeded node length prematurely (e.g., bad CIGAR/offset)
        # This might indicate an issue upstream, but helps prevent invalid variant positions later.
        # if ref_pos > node_len:
        #    print(f"⚠️ Warning: CIGAR processing resulted in ref_pos ({ref_pos}) exceeding node_len ({node_len}). CIGAR: '{cigar_string}', offset: {offset}", file=sys.stderr)
        #    break # Stop processing CIGAR for this read

    return variants

# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function

def init_worker(dat_file_path):
    """Initializer for each worker process: opens the .dat file."""
    global worker_dat_file
    try:
        worker_dat_file = open(dat_file_path, 'rb')
        # Optional: print(f"[Worker {os.getpid()}] Initialized with {dat_file_path}")
    except FileNotFoundError:
         print(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path}", file=sys.stderr)
         sys.exit(1) # Worker cannot proceed
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path}: {e}", file=sys.stderr)
        sys.exit(1) # Worker cannot proceed


def process_node_parallel(task_args):
    """
    Function executed by each worker process.
    Processes reads for a single node, identifies variants, and generates pileups.
    """
    node_id, file_offset, n_records, sequence, window_size = task_args # Added window_size
    global worker_dat_file

    if worker_dat_file is None:
         print(f"❌ Error [Worker {os.getpid()}]: Worker DAT file handle not initialized for node {node_id}.", file=sys.stderr)
         return node_id, {}

    if not sequence:
         # This might happen if GFA loading failed for this specific node
         # print(f"ℹ️ Info [Worker {os.getpid()}]: No sequence for node {node_id}, skipping.", file=sys.stderr)
         return node_id, {} # Cannot process without sequence

    node_len = len(sequence)
    if node_len == 0:
        # print(f"ℹ️ Info [Worker {os.getpid()}]: Empty sequence for node {node_id}, skipping.", file=sys.stderr)
        return node_id, {} # Cannot process empty sequence

    segments = [] # Store (offset_on_node, read_sequence, cigar_string)

    try:
        # Seek to the correct position for this node's records
        # Assuming index offset points directly to the start of the block,
        # and the first 10 bytes are metadata (node_id, block_offset, block_size).
        # If the .dat format is different, this needs adjustment.
        # Let's assume offset is the start of the block, and we skip metadata
        block_metadata_size = 10 # Example: size of node_id (4) + offset (8?) + block_size (4?) -> Adjust if known!
        # Correction: The index file format unpacks as <I Q I I H (22 bytes total)
        # node_id (4), offset (8), block_size (4), n_records (4), padding (2)
        # The offset (Q, 8 bytes) likely points to the START of the data block in the .dat file.
        # The first N bytes of that block might be metadata before the actual records start.
        # Let's ASSUME the first record starts immediately at 'file_offset'. If there's block metadata, adjust seek.
        # Example: If there's a 4-byte node_id at the start of the block:
        # worker_dat_file.seek(file_offset + 4)
        # For now, assume no extra metadata offset:
        worker_dat_file.seek(file_offset)

        for record_idx in range(n_records):
            try:
                data = worker_dat_file.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE:
                    # print(f"⚠️ Warning [Worker {os.getpid()}]: Short read ({len(data)} bytes) for node {node_id}, record {record_idx+1}/{n_records}. Assuming end of records for node.", file=sys.stderr)
                    break

                # Unpack record data: off, seq, bq, cigar, mapq, strand
                off, raw_seq, _, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)

                # Basic filtering
                if mapq < 10: # Example MAPQ filter
                    continue

                # Decode fields
                try:
                    # Use 'replace' or 'ignore' for potential non-ASCII bytes
                    seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
                    cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
                    strand_char = strand_byte.decode('ascii')
                except UnicodeDecodeError as ude:
                    # print(f"⚠️ Warning [Worker {os.getpid()}]: Unicode decode error in record {record_idx+1} for node {node_id}: {ude}. Skipping record.", file=sys.stderr)
                    continue

                if not seq or not cigar: # Skip reads with empty sequence or CIGAR
                    continue

                read_len = len(seq)

                # *** CORRECTED Reverse Strand Handling ***
                # The offset 'off' should refer to the position on the forward strand node sequence.
                # If the read aligns to the reverse strand, we only need to reverse complement the read sequence.
                current_seq = seq
                if strand_char == '-':
                    current_seq = reverse_complement(seq)
                    # Do NOT adjust offset based on node_len here. 'off' is the coordinate.

                # Validate offset before adding segment
                # Offset should be within the node bounds considering the alignment length from CIGAR later
                # A simple check: offset should be non-negative. More complex checks could parse CIGAR length.
                if off < 0:
                    # print(f"⚠️ Warning [Worker {os.getpid()}]: Negative offset ({off}) encountered for node {node_id}, record {record_idx+1}. Skipping read.", file=sys.stderr)
                    continue
                # Potentially add check: off < node_len?

                segments.append((off, current_seq, cigar)) # Store offset, potentially RC'd sequence, CIGAR

            except struct.error as se:
                 print(f"❌ Error [Worker {os.getpid()}]: Failed to unpack record {record_idx+1} for node {node_id}: {se}. Stopping reads for this node.", file=sys.stderr)
                 break
            except Exception as e_inner:
                 print(f"❌ Error [Worker {os.getpid()}]: Unexpected error processing record {record_idx+1} for node {node_id}: {e_inner}", file=sys.stderr)
                 continue # Skip problematic record

    except IOError as ioe:
         print(f"❌ Error [Worker {os.getpid()}]: I/O error reading data for node {node_id} at offset {file_offset}: {ioe}", file=sys.stderr)
         return node_id, {}
    except Exception as e_outer:
        print(f"❌ Error [Worker {os.getpid()}]: Unexpected error seeking/reading node {node_id}: {e_outer}", file=sys.stderr)
        return node_id, {}

    # --- Variant Detection and Pileup Generation ---
    reads_by_variant = defaultdict(list)
    for read_offset, read_seq, cigar in segments:
        # Detect variants based on CIGAR string relative to this read's alignment
        variants_in_read = detect_variants_from_cigar(read_offset, cigar, node_len)

        # Store the read sequence associated with each variant it covers
        for vpos, vtype in variants_in_read:
            # Key: tuple (variant_position_on_node, variant_type)
            # Value: list of (read_alignment_start_offset, read_sequence) tuples
            reads_by_variant[(vpos, vtype)].append((read_offset, read_seq))


    # --- Generate Pileups ---
    pileups = {}
    # Use the window_size passed from args
    half_window = window_size // 2

    for (vpos, vtype), supporting_reads in reads_by_variant.items():
        if not supporting_reads:
             continue

        # Create pileup matrix for this variant
        # Shape: (number_of_supporting_reads, window_size)
        # Initialize with 'N' index
        pileup_matrix = np.full((len(supporting_reads), window_size), BASE_TO_INDEX['N'], dtype=np.uint8)

        for i, (read_offset, read_seq) in enumerate(supporting_reads):
            read_len = len(read_seq)

            # Calculate the index *in the read* corresponding to the *start* of the window on the node.
            # Window start on node = vpos - half_window
            # Position of window start relative to read start = (vpos - half_window) - read_offset
            window_start_in_read = vpos - read_offset - half_window

            # Fill the row in the pileup matrix
            for j in range(window_size): # j is the column index in the pileup matrix (0 to window_size-1)
                # Corresponding position on the node = (vpos - half_window) + j
                # Corresponding index in the read = window_start_in_read + j
                read_idx = window_start_in_read + j

                # Check if this position falls within the bounds of the current read
                if 0 <= read_idx < read_len:
                    base = read_seq[read_idx].upper()
                    pileup_matrix[i, j] = BASE_TO_INDEX.get(base, BASE_TO_INDEX['N'])
                # else: Position is outside the read, leave as 'N' (already initialized)

        # Store pileup matrix (as list of lists for JSON)
        # Key format: "position_variantType" (e.g., "123_X", "455_I", "678_D")
        pileups[f"{vpos}_{vtype}"] = pileup_matrix.tolist()

    # Optional: print(f"[Worker {os.getpid()}] Finished node {node_id}. Found {len(pileups)} variant sites.")
    return node_id, pileups

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Parallel variant-centered read segment pileup generation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("dat", help=".dat file path (read alignment data)")
    parser.add_argument("idx", help=".idx file path (index for .dat file)")
    parser.add_argument("output", help="JSON output file path for pileups")
    parser.add_argument("--gfa", help="GFA graph file path (needed if node sequence cache is not used/built)")
    parser.add_argument("--load-cache", help="Load node sequences from this JSON cache file")
    parser.add_argument("--save-cache", help="Save node sequences to this JSON cache file (used if --gfa is provided)")
    parser.add_argument("-w", "--workers", type=int, default=8, # Consider os.cpu_count()
                        help="Number of worker processes to use")
    # *** ADDED window argument ***
    parser.add_argument("--window", type=int, default=75,
                        help="Pileup window size centered on the variant position")
    parser.add_argument("-c", "--chunksize", type=int, default=100, # Adjusted default chunksize
                        help="Number of nodes processed by a worker before returning results (approx)")
    args = parser.parse_args()

    # --- Input Validation ---
    if not os.path.isfile(args.dat):
         print(f"❌ Error: DAT file not found: {args.dat}", file=sys.stderr)
         sys.exit(1)
    if not os.path.isfile(args.idx):
         print(f"❌ Error: Index file not found: {args.idx}", file=sys.stderr)
         sys.exit(1)
    # Validate window size
    if args.window <= 0 or args.window % 2 == 0:
        print(f"❌ Error: Pileup window size (--window) must be a positive odd integer. Found: {args.window}", file=sys.stderr)
        sys.exit(1)

    # Cache/GFA validation
    if not args.load_cache and not args.gfa:
        print("❌ Error: You must provide either a GFA file (`--gfa`) or a node sequence cache (`--load-cache`).", file=sys.stderr)
        sys.exit(1)
    if args.gfa and not args.save_cache and not args.load_cache:
        print("❌ Error: When providing `--gfa` without `--load-cache`, you must also specify `--save-cache`.", file=sys.stderr)
        sys.exit(1)
    if args.load_cache and not os.path.isfile(args.load_cache):
        print(f"❌ Error: Cache file to load does not exist: {args.load_cache}", file=sys.stderr)
        sys.exit(1)
    if args.gfa and not os.path.isfile(args.gfa):
         print(f"❌ Error: GFA file not found: {args.gfa}", file=sys.stderr)
         sys.exit(1)

    # --- Load Index ---
    print("🔹 Parsing index file...")
    start_time = time.time()
    node_index = parse_idx_file(args.idx)
    if not node_index:
         print("❌ Error: Failed to parse node index. Exiting.", file=sys.stderr)
         sys.exit(1)
    print(f"✔ Index parsing took {time.time() - start_time:.2f} seconds. {len(node_index)} nodes to process.")


    # --- Load or Build Node Sequences ---
    node_sequences = {}
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"🔹 Loading node sequences from cache: {args.load_cache}...")
        start_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                loaded_data = json.load(cf)
                # Ensure keys are integers after loading from JSON
                node_sequences = {int(k): v for k, v in loaded_data.items()}
            print(f"✔ Loaded {len(node_sequences)} sequences from cache in {time.time() - start_time:.2f} seconds.")
        except Exception as e:
             print(f"❌ Error loading cache file {args.load_cache}: {e}", file=sys.stderr)
             sys.exit(1) # Critical error if cache load fails
    elif args.gfa:
        print(f"🔹 Loading/Building node sequences from GFA: {args.gfa}...")
        start_time = time.time()
        # Load sequences only for nodes present in the index
        node_sequences = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        print(f"✔ Sequence loading/building took {time.time() - start_time:.2f} seconds.")
        if args.save_cache:
            print(f"🔹 Saving node sequence cache to: {args.save_cache}...")
            start_time = time.time()
            try:
                # Ensure keys are saved as strings in JSON
                save_data = {str(k): v for k, v in node_sequences.items()}
                with open(args.save_cache, 'w') as cf:
                    json.dump(save_data, cf) # No indent for smaller file size
                print(f"✔ Saved cache ({len(save_data)} nodes) in {time.time() - start_time:.2f} seconds.")
            except Exception as e:
                 print(f"❌ Error saving cache file {args.save_cache}: {e}", file=sys.stderr)
                 # Non-critical, continue execution

    # --- Prepare Tasks for Parallel Processing ---
    print("🔹 Preparing tasks for parallel processing...")
    tasks = []
    nodes_missing_sequence = 0
    nodes_with_records = 0
    for node_id, (offset, n_records) in node_index.items():
        sequence = node_sequences.get(node_id)
        if sequence:
            # *** ADDED args.window to task arguments ***
            tasks.append((node_id, offset, n_records, sequence, args.window))
            nodes_with_records += 1
        else:
            nodes_missing_sequence += 1

    if nodes_missing_sequence > 0:
        print(f"⚠️ Warning: Skipped {nodes_missing_sequence} nodes present in index but missing sequence data (from GFA/cache).")
    if not tasks:
        print("❌ Error: No tasks to process (no nodes found with both sequence data and records in index). Exiting.", file=sys.stderr)
        sys.exit(1)

    total_tasks = len(tasks)
    num_workers = min(args.workers, total_tasks, os.cpu_count() or args.workers) # Be reasonable with workers
    if num_workers <= 0: num_workers = 1

    print(f"🔹 Processing {total_tasks} nodes using {num_workers} workers (window: {args.window}, chunksize: {args.chunksize})...")

    # --- Execute in Parallel ---
    results = {} # Store {node_id: pileup_dict}
    start_proc_time = time.time()

    try:
        with ProcessPoolExecutor(max_workers=num_workers,
                                 initializer=init_worker,
                                 initargs=(args.dat,)) as executor:

            # Use executor.map for simplicity, results arrive in submission order (approx)
            future_results = executor.map(process_node_parallel, tasks, chunksize=args.chunksize)

            processed_count = 0
            variants_found_count = 0
            # Iterate through results as they become available
            for node_id, pileup_dict in future_results:
                 processed_count += 1
                 if pileup_dict: # Store result only if pileups were generated
                     results[node_id] = pileup_dict
                     variants_found_count += len(pileup_dict)

                 # Progress Update
                 if processed_count % (args.chunksize * num_workers / 2) == 0 or processed_count == total_tasks: # Update periodically
                     elapsed = time.time() - start_proc_time
                     rate = processed_count / elapsed if elapsed > 0 else 0
                     eta = (total_tasks - processed_count) / rate if rate > 0 else 0
                     print(f"  Processed: {processed_count}/{total_tasks} nodes ({rate:.1f}/s) | Variants Found: {variants_found_count} | Elapsed: {elapsed:.1f}s | ETA: {eta:.1f}s", end='\r')

    except Exception as pool_exc:
         print(f"\n❌ An error occurred during parallel processing: {pool_exc}", file=sys.stderr)
         print("⚠️ Attempting to write any partial results obtained...", file=sys.stderr)
    finally:
        # Ensure newline after progress indicator
        print() # Moves to the next line after the final progress update

    total_elapsed_time = time.time() - start_proc_time
    print(f"✔ Parallel processing finished in {total_elapsed_time:.2f} seconds.")

    # --- Write Output ---
    if not results:
        print("⚠️ No pileups were generated. Output file will be empty or not created if it doesn't exist.")
        print("✅ Done.")
        sys.exit(0) # Exit cleanly, no results to write

    print(f"🔹 Writing pileups for {len(results)} nodes (total {variants_found_count} variant sites) to JSON: {args.output}")
    start_write_time = time.time()
    try:
        with open(args.output, 'w') as out_f:
            # Save node IDs as strings for JSON keys
            json.dump({str(k): v for k, v in results.items()}, out_f) # No indent for smaller file size
        write_elapsed_time = time.time() - start_write_time
        print(f"✔ Output written in {write_elapsed_time:.2f} seconds.")
        print(f"✅ Done. Output saved to {args.output}")
    except Exception as e:
         print(f"❌ Error writing output JSON to {args.output}: {e}", file=sys.stderr)
         sys.exit(1)


if __name__ == '__main__':
    main()