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
import traceback # Import traceback for detailed error logging
from concurrent.futures import ProcessPoolExecutor, as_completed # Import as_completed

# ─────────────────────────────────────────────────────────────────────────────
# Constants
RECORD_STRUCT = struct.Struct("<h150s150s20shc") # Assumed structure (verify!)
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# Global for worker process state
worker_dat_file = None

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions (Keep functions from your script: reverse_complement, parse_idx_file, load_node_sequences_from_gfa, decode_cigar, detect_variants_from_cigar)

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def parse_idx_file(idx_path):
    # --- Using the robust binary index parsing from your script ---
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
                # sys.exit(1) # Optionally exit

            print(f"🔹 Reading index for {num_nodes} nodes...")
            for i in range(num_nodes):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                     print(f"❌ Error: Index file ended prematurely while reading record {i+1}/{num_nodes}.", file=sys.stderr)
                     break
                # --- Assuming index struct is: node_id(I), offset(Q), block_size(I), n_records(I), padding(H) ---
                # --- Verify this matches your actual index file structure ---
                try:
                     node_id, offset, block_size, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                     node_index[node_id] = (offset, n_records)
                except struct.error as se:
                     print(f"❌ Error: Failed unpacking index record {i+1}/{num_nodes}: {se}. Check index format.", file=sys.stderr)
                     # Decide whether to break or continue
                     break

        parsed_count = len(node_index)
        print(f"✔ Parsed {parsed_count} nodes from index file.")
        if parsed_count != num_nodes:
             print(f"⚠️ Warning: Expected {num_nodes} nodes based on header, but parsed {parsed_count}. File might be corrupt or truncated.", file=sys.stderr)
        return node_index
    except FileNotFoundError:
        print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr) # Print traceback for unexpected errors
        sys.exit(1)


def load_node_sequences_from_gfa(gfa_path, target_node_ids):
    # --- Using the GFA parsing logic from your script ---
    node_sequences = {}
    target_node_set = set(target_node_ids)
    try:
        with open(gfa_path, 'r', encoding='utf-8') as f: # Specify UTF-8 for GFA
            print(f"🔹 Reading GFA file: {gfa_path}")
            line_counter = 0
            parsed_count = 0
            start_gfa_time = time.time()
            for line in f:
                line_counter += 1
                if line_counter % 10_000_000 == 0:
                    elapsed = time.time() - start_gfa_time
                    rate = line_counter / elapsed if elapsed > 0 else 0
                    print(f"  Checked {line_counter:,} lines in GFA file... ({rate:.0f} lines/s)")

                if not line.startswith('S\t'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue

                try:
                    nid = int(parts[1])
                except ValueError:
                    continue

                if nid in target_node_set:
                    node_sequences[nid] = parts[2]
                    parsed_count += 1
                    # Optimization: stop early if all found (uncomment if desired)
                    # target_node_set.remove(nid)
                    # if not target_node_set:
                    #    print(f"\n✔ Found all {len(target_node_ids)} target sequences.")
                    #    break

            print(f"✔ Checked {line_counter:,} lines in GFA.")
            print(f"✔ Loaded {len(node_sequences)} target sequences from GFA.")
            if len(node_sequences) != len(target_node_ids):
                 missing_count = len(target_node_ids) - len(node_sequences)
                 print(f"⚠️ Warning: Could not find sequences for {missing_count} target nodes in the GFA.", file=sys.stderr)

    except FileNotFoundError:
        print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

    return node_sequences

def decode_cigar(cigar_string):
    # --- Using the regex CIGAR parsing from your script ---
    if not cigar_string or cigar_string == '*':
        return []
    try:
        # Include common CIGAR ops M, I, D, X, = (S, H, N, P often ignored in pileups)
        return re.findall(r'(\d+)([MIDX=])', cigar_string)
    except Exception as e:
        # Add error handling here too
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr)
        return []


def detect_variants_from_cigar(offset, cigar_string):
    # --- Using the variant detection logic from your script ---
    variants = []
    pos = offset
    cigar_ops = decode_cigar(cigar_string)

    for length_str, op in cigar_ops:
        try:
            length = int(length_str)
        except ValueError:
            print(f"⚠️ Warning: Invalid length in CIGAR operation '{length_str}{op}' in string '{cigar_string}'", file=sys.stderr)
            continue

        if op == 'M' or op == '=':
            pos += length
        elif op == 'X':
            for i in range(length):
                variants.append((pos + i, 'X'))
            pos += length
        elif op == 'I':
             # Insertion *after* pos-1 on reference. Store pos-1 as the anchor.
             variants.append((pos - 1, 'I'))
             # pos does not change
        elif op == 'D':
             # Deletion *at* pos on reference. Store pos as the anchor.
             variants.append((pos, 'D'))
             pos += length
        # Other ops ignored

    return variants

# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function

def init_worker(dat_file_path):
    """Initializer for each worker process: opens the .dat file."""
    global worker_dat_file
    pid = os.getpid()
    # print(f"[Worker {pid}] Initializing...") # Debug
    try:
        worker_dat_file = open(dat_file_path, 'rb')
        # print(f"[Worker {pid}] Opened {dat_file_path}") # Debug
    except FileNotFoundError:
         print(f"❌ Error [Worker {pid}]: DAT file not found at {dat_file_path}. Worker exiting.", file=sys.stderr)
         sys.exit(1) # Exit worker if essential file is missing
    except Exception as e:
        print(f"❌ Error [Worker {pid}] opening DAT file {dat_file_path}: {e}. Worker exiting.", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

def process_node_parallel(task_args):
    """
    Function executed by each worker process. Enhanced error handling.
    Returns (node_id, pileup_dict) on success, or (node_id, {}) on failure.
    """
    node_id, offset, n_records, sequence = task_args
    global worker_dat_file
    pid = os.getpid()

    # --- Pre-check worker state ---
    if worker_dat_file is None:
         print(f"❌ Error [Worker {pid}][Node {node_id}]: Worker DAT file handle not initialized.", file=sys.stderr)
         return node_id, {} # Indicate failure
    if not sequence:
         # print(f"⚠️ Warning [Worker {pid}][Node {node_id}]: Empty sequence provided. Skipping.", file=sys.stderr)
         return node_id, {} # Cannot process without sequence
    node_len = len(sequence)

    # --- Main processing block with enhanced error handling ---
    segments = []
    pileups = {} # Final result for this node
    try:
        # --- Seek and Read Records ---
        try:
            # Use offset from index directly, add 10 bytes for block header
            worker_dat_file.seek(offset + 10)
        except IOError as ioe:
            print(f"❌ Error [Worker {pid}][Node {node_id}]: Failed to seek to offset {offset + 10}: {ioe}", file=sys.stderr)
            return node_id, {} # Fail task

        for record_idx in range(n_records):
            try:
                data = worker_dat_file.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE:
                    print(f"⚠️ Warning [Worker {pid}][Node {node_id}]: Short read ({len(data)} bytes) for record {record_idx+1}/{n_records}. Stopping reads.", file=sys.stderr)
                    break # Stop processing this node

                # --- Unpack Record Data ---
                try:
                    # --- Use the structure defined in your script ---
                    off, raw_seq, raw_bq, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)
                except struct.error as se:
                     print(f"❌ Error [Worker {pid}][Node {node_id}]: Struct unpack error, record {record_idx+1}: {se}. Skipping record.", file=sys.stderr)
                     continue # Skip this problematic record

                # --- Decode and Filter ---
                # Basic filtering
                if mapq < 10:
                    continue

                try:
                    # Use 'replace' or 'ignore' for robustness
                    seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')
                    cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                    strand_char = strand_byte.decode('ascii')
                except Exception as decode_err: # Catch potential errors in decode too
                    print(f"⚠️ Warning [Worker {pid}][Node {node_id}]: Decode error, record {record_idx+1}: {decode_err}. Skipping.", file=sys.stderr)
                    continue

                read_len = len(seq)

                # --- Adjust for Reverse Strand ---
                # Using the logic from your script
                if strand_char == '-' and node_len > 0:
                     adj_off = node_len - off - read_len
                     if adj_off < 0:
                        # print(f"⚠️ Warning [Worker {pid}][Node {node_id}]: Negative adjusted offset ({adj_off}), record {record_idx+1}. Skipping.", file=sys.stderr)
                        continue # Skip inconsistent read
                     seq = reverse_complement(seq)
                     off = adj_off # Use adjusted offset

                # Append processed segment info
                segments.append((off, seq, cigar))

            except IOError as read_ioe: # Catch error during read itself
                 print(f"❌ Error [Worker {pid}][Node {node_id}]: I/O error reading record {record_idx+1}: {read_ioe}. Stopping reads.", file=sys.stderr)
                 break # Stop processing this node
            except Exception as record_err: # Catch unexpected errors processing a single record
                 print(f"❌ Error [Worker {pid}][Node {node_id}]: Unexpected error processing record {record_idx+1}: {record_err}", file=sys.stderr)
                 traceback.print_exc(file=sys.stderr)
                 continue # Skip problematic record

        # --- End of reading loop ---

        # --- Variant Detection and Pileup Generation ---
        reads_by_variant = defaultdict(list)
        try:
            for off, seq, cigar in segments:
                if off < 0: continue # Skip segments with negative offset

                # Detect variants based on CIGAR string
                variants = detect_variants_from_cigar(off, cigar)
                for vpos, vtype in variants:
                    # Validate variant position (using logic from your script)
                    is_valid_pos = False
                    if vtype == 'I' and vpos >= -1 and vpos < node_len:
                         is_valid_pos = True
                    elif vtype in ('D', 'X') and vpos >= 0 and vpos < node_len:
                         is_valid_pos = True

                    if is_valid_pos:
                        reads_by_variant[(vpos, vtype)].append((off, seq))

        except Exception as var_err:
             print(f"❌ Error [Worker {pid}][Node {node_id}]: Error during variant detection phase: {var_err}", file=sys.stderr)
             traceback.print_exc(file=sys.stderr)
             # Decide if we can continue to pileup with potentially partial variant info
             # For safety, let's fail the node if variant detection has issues
             return node_id, {}

        # --- Generate Pileups ---
        try:
            window = 60
            half = window // 2
            for (vpos, vtype), reads in reads_by_variant.items():
                if not reads: continue

                mat = np.full((len(reads), window), BASE_TO_INDEX['N'], dtype=np.uint8)
                for i, (read_offset, read_seq) in enumerate(reads):
                    read_len = len(read_seq)
                    start_in_read = vpos - read_offset - half

                    for j in range(window):
                        read_idx = start_in_read + j
                        if 0 <= read_idx < read_len:
                            base = read_seq[read_idx].upper()
                            mat[i, j] = BASE_TO_INDEX.get(base, BASE_TO_INDEX['N'])

                # Use f-string for key, ensure it's unique
                pileups[f"{vpos}_{vtype}"] = mat.tolist()

        except Exception as pileup_err:
            print(f"❌ Error [Worker {pid}][Node {node_id}]: Error during pileup generation phase: {pileup_err}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            # Fail the node if pileup generation has issues
            return node_id, {}

    # --- Catch-all for unexpected errors in the main worker logic ---
    except Exception as e_main:
        print(f"❌❌❌ [Worker {pid}][Node {node_id}] UNCAUGHT MAJOR error during processing: {e_main}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return node_id, {} # Indicate failure

    # print(f"[Worker {pid}][Node {node_id}] Processed successfully. Found {len(pileups)} variants.") # Debug
    return node_id, pileups # Return result on success

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Parallel variant-centered pileup generation (Robust Version).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # --- Using argument parsing from your script ---
    parser.add_argument("dat", help=".dat file path (read alignment data)")
    parser.add_argument("idx", help=".idx file path (index for .dat file)")
    parser.add_argument("output", help="JSON output file path for pileups")
    parser.add_argument("--gfa", help="GFA graph file path (to build node sequence cache)")
    parser.add_argument("--load-cache", help="Load node sequences from this JSON cache file")
    parser.add_argument("--save-cache", help="Save node sequences to this JSON cache file (if using --gfa)")
    parser.add_argument("-w", "--workers", type=int, default=8,
                        help="Number of worker processes to use. Tune based on I/O tests!")
    # Chunksize is less relevant for submit/as_completed, can be removed or ignored
    # parser.add_argument("-c", "--chunksize", type=int, default=200, help="Nodes per chunk (map only)")
    args = parser.parse_args()

    # --- Input Validation (from your script) ---
    if not os.path.isfile(args.dat): print(f"❌ Error: DAT file not found: {args.dat}", file=sys.stderr); sys.exit(1)
    if not os.path.isfile(args.idx): print(f"❌ Error: Index file not found: {args.idx}", file=sys.stderr); sys.exit(1)
    if not args.load_cache and not args.gfa: print("❌ Error: Must provide --gfa or --load-cache.", file=sys.stderr); sys.exit(1)
    if args.gfa and not args.save_cache and not args.load_cache: print("❌ Error: If using --gfa without --load-cache, must specify --save-cache.", file=sys.stderr); sys.exit(1)
    if args.load_cache and not os.path.isfile(args.load_cache): print(f"❌ Error: Cache file to load not found: {args.load_cache}", file=sys.stderr); sys.exit(1)
    if args.gfa and not os.path.isfile(args.gfa): print(f"❌ Error: GFA file not found: {args.gfa}", file=sys.stderr); sys.exit(1)

    # --- Load Index (using your function) ---
    print("🔹 Parsing index file...")
    start_time = time.time()
    node_index = parse_idx_file(args.idx)
    if not node_index: print("❌ Error: Failed to parse node index.", file=sys.stderr); sys.exit(1)
    print(f"✔ Index parsing took {time.time() - start_time:.2f} seconds.")

    # --- Load or Build Node Sequences (using your logic) ---
    node_sequences = {}
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"🔹 Loading node sequences from cache: {args.load_cache}...")
        start_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                loaded_data = json.load(cf)
                node_sequences = {int(k): v for k, v in loaded_data.items()} # Load keys as int
            print(f"✔ Loaded {len(node_sequences)} sequences from cache in {time.time() - start_time:.2f} seconds.")
        except Exception as e:
             print(f"❌ Error loading cache file {args.load_cache}: {e}", file=sys.stderr); sys.exit(1)
    elif args.gfa:
        print(f"🔹 Building node sequence cache from GFA: {args.gfa}...")
        start_time = time.time()
        node_sequences = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        print(f"✔ Sequence loading from GFA took {time.time() - start_time:.2f} seconds.")
        if args.save_cache:
            print(f"🔹 Saving node sequence cache to: {args.save_cache}...")
            start_time = time.time()
            try:
                with open(args.save_cache, 'w') as cf:
                    json.dump({str(k): v for k, v in node_sequences.items()}, cf) # Save keys as str
                print(f"✔ Saved cache in {time.time() - start_time:.2f} seconds.")
            except Exception as e:
                 print(f"❌ Error saving cache file {args.save_cache}: {e}", file=sys.stderr) # Continue even if save fails

    # --- Prepare Tasks ---
    print("🔹 Preparing tasks...")
    tasks = []
    nodes_missing_sequence = 0
    for node_id, (offset, n_records) in node_index.items():
        sequence = node_sequences.get(node_id)
        if sequence is not None and n_records > 0: # Also skip nodes with 0 records
            tasks.append((node_id, offset, n_records, sequence))
        elif sequence is None:
            nodes_missing_sequence += 1

    if nodes_missing_sequence > 0: print(f"⚠️ Warning: Skipped {nodes_missing_sequence} nodes missing sequence data.")
    if not tasks: print("❌ Error: No tasks to process.", file=sys.stderr); sys.exit(1)

    total_tasks = len(tasks)
    num_workers = min(args.workers, total_tasks)
    if num_workers <= 0: num_workers = 1

    print(f"🔹 Processing {total_tasks} nodes using {num_workers} workers...")

    # --- Execute in Parallel using submit/as_completed ---
    results = {}
    start_proc_time = time.time()
    processed_count = 0
    errors_reported = 0 # Count errors reported from workers or future retrieval

    try:
        with ProcessPoolExecutor(max_workers=num_workers,
                                 initializer=init_worker,
                                 initargs=(args.dat,)) as executor:

            # Submit tasks and map future object back to node_id
            futures_map = {executor.submit(process_node_parallel, task): task[0] for task in tasks}
            print(f"--- Submitted {len(futures_map)} tasks. Waiting for results...")

            # Process results as they complete
            for future in as_completed(futures_map):
                node_id = futures_map[future] # Get node_id for this future
                try:
                    # Get result - this re-raises exceptions from worker *startup* or *pickling*,
                    # but NOT exceptions caught within process_node_parallel's try/except blocks.
                    _returned_node_id, pileup_dict = future.result()

                    # Check if worker returned a valid (non-empty) dict
                    if pileup_dict:
                         results[node_id] = pileup_dict
                    else:
                         # Worker returned {} indicating a caught error, already logged
                         errors_reported += 1

                except Exception as exc:
                    # Handle errors during future retrieval itself (e.g., worker crashed hard)
                    print(f'❌❌❌ Main thread error getting result for node {node_id}: {exc}', file=sys.stderr)
                    traceback.print_exc(file=sys.stderr)
                    errors_reported += 1

                # --- Progress Update ---
                processed_count += 1
                if processed_count % 10000 == 0 or processed_count == total_tasks:
                     elapsed = time.time() - start_proc_time
                     rate = processed_count / elapsed if elapsed > 0 else 0
                     # Report errors seen so far
                     error_msg = f"({errors_reported} errors reported)" if errors_reported > 0 else ""
                     print(f"✔ {processed_count}/{total_tasks} nodes processed {error_msg}. Rate: {rate:.1f} nodes/s. Elapsed: {elapsed:.2f}s")

    except Exception as pool_exc:
         print(f"\n❌ An error occurred within the ProcessPoolExecutor context: {pool_exc}", file=sys.stderr)
         traceback.print_exc(file=sys.stderr)
         print("⚠️ Attempting to write any partial results obtained...", file=sys.stderr)

    total_elapsed_time = time.time() - start_proc_time
    print(f"\n✔ Parallel processing finished in {total_elapsed_time:.2f} seconds.")
    print(f"Processed {processed_count} tasks. Stored results for {len(results)} nodes.")
    if errors_reported > 0:
         print(f"⚠️ Encountered {errors_reported} errors during processing (check logs above).", file=sys.stderr)

    # --- Write Output ---
    print(f"🔹 Writing {len(results)} results to JSON output: {args.output}")
    start_write_time = time.time()
    try:
        with open(args.output, 'w') as out_f:
            # Ensure keys are strings for JSON
            json.dump({str(k): v for k, v in results.items()}, out_f, indent=2) # Add indent for readability
        write_elapsed_time = time.time() - start_write_time
        print(f"✔ Output written in {write_elapsed_time:.2f} seconds.")
        print(f"✅ Done. Output saved to {args.output}")
    except Exception as e:
         print(f"❌ Error writing output JSON to {args.output}: {e}", file=sys.stderr)
         traceback.print_exc(file=sys.stderr)
         sys.exit(1)


if __name__ == '__main__':
    main()