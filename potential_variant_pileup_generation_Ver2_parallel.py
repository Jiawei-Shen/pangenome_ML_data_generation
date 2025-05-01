#!/usr/bin/env python3

import argparse
import sys
import os
import time
import json
import struct
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import numpy as np # Assuming numpy is used for pileups

# --- Constants ---
RECORD_SIZE = 322 # Size of each record in the .dat file (adjust if needed)
CIGAR_OPS = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X'}
BASES = ['A', 'C', 'G', 'T', 'N'] # Order for pileup array

# --- Global variable for worker ---
worker_dat_file = None

# ===============================================
# Worker Initialization
# ===============================================
def init_worker(dat_file_path):
    """Initializer for each worker process: opens the .dat file."""
    global worker_dat_file
    pid = os.getpid()
    print(f"--- [Worker {pid}] Initializing...")
    try:
        worker_dat_file = open(dat_file_path, 'rb')
        print(f"--- [Worker {pid}] Successfully opened {dat_file_path}")
    except Exception as e:
        print(f"❌❌❌ [Worker {pid}] FAILED to open DAT file {dat_file_path}: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        worker_dat_file = None # Ensure it's None if open failed

# ===============================================
# Helper Functions (Example CIGAR Parser)
# ===============================================
def parse_cigar(cigar_bytes):
    """Parses CIGAR bytes into a list of (operation, length) tuples."""
    # Assuming CIGAR bytes are packed integers (adjust format if needed)
    # Example: 4 bytes per operation/length pair? Or null-terminated string?
    # This needs to match your actual CIGAR data format in the .dat file
    # --- Placeholder Implementation ---
    # Replace with your actual CIGAR parsing logic
    cigar_string = cigar_bytes.rstrip(b'\x00').decode('ascii')
    if not cigar_string or cigar_string == '*':
        return []

    parsed = []
    current_len = ""
    for char in cigar_string:
        if char.isdigit():
            current_len += char
        else:
            op = CIGAR_OPS.get(int(char) if char.isdigit() else -1) # Handle potential direct op codes? Unlikely in standard CIGAR
            if char in CIGAR_OPS.values() and current_len: # More standard check
                 parsed.append((char, int(current_len)))
                 current_len = ""
            # Handle other CIGAR format variations if needed
    return parsed # Example return format

def get_reference_position(pos, cigar_tuples, query_pos):
    """ Calculates reference position given alignment start, CIGAR, and query position. """
    # --- Placeholder Implementation ---
    # Replace with your actual logic based on CIGAR
    ref_pos = pos
    q_pos = 0
    for op, length in cigar_tuples:
        if op in ('M', '=', 'X'): # Consumes query and reference
            if q_pos + length > query_pos:
                return ref_pos + (query_pos - q_pos)
            ref_pos += length
            q_pos += length
        elif op == 'D' or op == 'N': # Consumes reference only
            ref_pos += length
        elif op == 'I' or op == 'S': # Consumes query only
            q_pos += length
        if q_pos > query_pos and op != 'D' and op != 'N': # Check if query position passed
             break # Approximate position found or exact if in M/=/X block
        # H and P don't consume sequence positions

    # Return based on last alignment block if exact position not found in M/=/X
    return ref_pos # Might need refinement based on exact CIGAR semantics needed


# ===============================================
# Core Worker Function
# ===============================================
def process_node_parallel(task_args):
    """
    Function executed by each worker process.
    Processes a single node to find variants and generate pileups.
    """
    node_id, offset, n_records, node_sequence = task_args
    global worker_dat_file
    pid = os.getpid()
    # print(f"--- [Worker {pid}] Starting task for Node ID: {node_id}") # Verbose logging

    # --- Check Worker Initialization ---
    if worker_dat_file is None:
         print(f"❌ Error [Worker {pid}][Node {node_id}]: Worker DAT file handle not initialized.", file=sys.stderr)
         return node_id, {} # Return empty dict for this node

    pileups = {} # Store pileup data: {ref_pos: [A, C, G, T, N, Ins, Del]}

    # --- Main Processing Block ---
    try:
        # --- Seek and Read Records ---
        try:
            # Seek past the 10-byte header for the block (adjust if header size differs)
            worker_dat_file.seek(offset + 10)
            read_count = 0
            for _ in range(n_records):
                record_bytes = worker_dat_file.read(RECORD_SIZE)
                if len(record_bytes) != RECORD_SIZE:
                    print(f"⚠️ Warning [Worker {pid}][Node {node_id}]: Short read at record {read_count+1}/{n_records}. Expected {RECORD_SIZE}, got {len(record_bytes)}. Offset: {offset}", file=sys.stderr)
                    break # Stop processing this node if reads are truncated

                # --- Unpack Record Data ---
                # Adjust unpack format string based on your ACTUAL record structure
                # Example format (replace with yours):
                # qname (e.g., 32s), flag (H=uint16), rname (e.g., 16s), pos (i=int32), mapq (B=uint8),
                # cigar (e.g., 64s), rnext (e.g., 16s), pnext (i=int32), tlen (i=int32),
                # seq (e.g., 151s), qual (e.g., 151s), ... other fields
                # THIS IS JUST AN EXAMPLE - VERIFY YOUR FORMAT!
                try:
                    # Example: Assuming fixed fields first, then variable CIGAR/SEQ/QUAL
                    # Recalculate offsets based on actual fixed field sizes!
                    # flag = struct.unpack_from('<H', record_bytes, 32)[0] # Example offset
                    pos = struct.unpack_from('<i', record_bytes, 34)[0]  # Example offset
                    mapq = struct.unpack_from('<B', record_bytes, 38)[0] # Example offset
                    cigar_bytes = record_bytes[40:104] # Example slice for CIGAR
                    seq_bytes = record_bytes[104:255] # Example slice for SEQ
                    # Add unpacking for other fields as needed

                except struct.error as se:
                    print(f"❌ Error [Worker {pid}][Node {node_id}]: Struct unpack error processing record {read_count+1}/{n_records}: {se}. Offset: {offset}", file=sys.stderr)
                    # Optionally continue to next record or fail node? Decide policy.
                    continue # Skip this record

                read_count += 1

                # --- Filter Reads (Example) ---
                if mapq < 20: # Example MAPQ filter
                    continue
                if pos < 0: # Unmapped or invalid position
                    continue

                # --- Process Read (Placeholder Logic) ---
                try:
                    # Decode sequence, parse CIGAR, etc.
                    read_seq = seq_bytes.rstrip(b'\x00').decode('ascii')
                    cigar_tuples = parse_cigar(cigar_bytes) # Use your actual parser

                    # --- Pileup Logic ---
                    # Iterate through read bases or CIGAR ops to build pileup
                    # This is complex and highly dependent on your specific needs
                    # Example: Simplified base counting at reference positions
                    ref_idx = pos - 1 # Assuming 1-based position in file -> 0-based index
                    read_idx = 0
                    for op, length in cigar_tuples:
                        if op == 'M' or op == '=' or op == 'X': # Match or mismatch
                            for i in range(length):
                                current_ref_idx = ref_idx + i
                                current_read_idx = read_idx + i
                                if current_ref_idx < 0 or current_ref_idx >= len(node_sequence): continue # Bounds check
                                if current_read_idx < 0 or current_read_idx >= len(read_seq): continue # Bounds check

                                ref_base = node_sequence[current_ref_idx].upper()
                                read_base = read_seq[current_read_idx].upper()

                                if current_ref_idx not in pileups:
                                     # A, C, G, T, N, Ins, Del counts
                                     pileups[current_ref_idx] = np.zeros(7, dtype=np.uint32)

                                if read_base in BASES:
                                    pileups[current_ref_idx][BASES.index(read_base)] += 1
                                else:
                                    pileups[current_ref_idx][BASES.index('N')] += 1 # Count others as N

                            ref_idx += length
                            read_idx += length
                        elif op == 'D' or op == 'N': # Deletion relative to reference
                             for i in range(length):
                                 current_ref_idx = ref_idx + i
                                 if current_ref_idx < 0 or current_ref_idx >= len(node_sequence): continue
                                 if current_ref_idx not in pileups:
                                     pileups[current_ref_idx] = np.zeros(7, dtype=np.uint32)
                                 pileups[current_ref_idx][6] += 1 # Increment Del count
                             ref_idx += length
                        elif op == 'I': # Insertion relative to reference
                            # Increment Ins count at position *before* insertion
                            if ref_idx >= 0 and ref_idx < len(node_sequence):
                                if ref_idx not in pileups:
                                    pileups[ref_idx] = np.zeros(7, dtype=np.uint32)
                                pileups[ref_idx][5] += 1 # Increment Ins count
                            read_idx += length
                        elif op == 'S': # Soft clipping - consumes query
                            read_idx += length
                        # Handle H, P if necessary (usually ignored for pileups)

                except Exception as e_proc:
                     print(f"❌ Error [Worker {pid}][Node {node_id}]: Error processing read {read_count}/{n_records} (Pos: {pos}): {e_proc}", file=sys.stderr)
                     traceback.print_exc(file=sys.stderr)
                     # Decide: continue to next read or fail the node?
                     continue # Example: Skip this problematic read

        except IOError as ioe:
             print(f"❌ Error [Worker {pid}][Node {node_id}]: I/O error reading data block: {ioe}. Offset: {offset}", file=sys.stderr)
             return node_id, {} # Fail this task

        except Exception as e_read_block:
             print(f"❌ Error [Worker {pid}][Node {node_id}]: Unexpected error during record reading loop: {e_read_block}. Offset: {offset}", file=sys.stderr)
             traceback.print_exc(file=sys.stderr)
             return node_id, {} # Fail this task

    # --- Catch-all for any other unexpected error in the worker ---
    except Exception as e_main:
        print(f"❌❌❌ [Worker {pid}][Node {node_id}] UNCAUGHT MAJOR error during processing: {e_main}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return node_id, {} # Fail this task

    # --- Convert numpy arrays to lists for JSON serialization ---
    final_pileups = {pos: arr.tolist() for pos, arr in pileups.items()}

    # print(f"--- [Worker {pid}] Finished task for Node ID: {node_id}. Found {len(final_pileups)} pileup positions.") # Verbose logging
    return node_id, final_pileups

# ===============================================
# Main Execution Logic
# ===============================================
def main():
    parser = argparse.ArgumentParser(description="Generate pileup information from alignment data.")
    parser.add_argument("dat_file", help="Path to the input binary .dat file.")
    parser.add_argument("index_file", help="Path to the index file (.idx).")
    parser.add_argument("output_json", help="Path to the output JSON file for pileups.")
    parser.add_argument("--load-cache", help="Path to node sequences JSON cache.")
    parser.add_argument("-w", "--workers", type=int, default=4,
                        help="Number of worker processes. IMPORTANT: Match this to SLURM CPUs requested & tune based on I/O stability tests. Default: 4")
    # Chunksize doesn't directly apply to submit/as_completed in the same way as map
    # parser.add_argument("--chunksize", type=int, default=50, help="Number of tasks per chunk (for map, less relevant for submit).")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    print(f"Starting pileup generation at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Parameters: DAT='{args.dat_file}', Index='{args.index_file}', Output='{args.output_json}', Cache='{args.load_cache}', Workers={args.workers}")

    # --- Load Index ---
    start_time = time.time()
    node_index = {}
    print(f"🔹 Parsing index file: {args.index_file}...")
    try:
        with open(args.index_file, 'r') as f_idx:
            for line in f_idx:
                parts = line.strip().split('\t')
                if len(parts) == 3:
                    node_id, offset, n_records = parts
                    node_index[node_id] = (int(offset), int(n_records))
                else:
                    print(f"⚠️ Warning: Skipping malformed index line: {line.strip()}", file=sys.stderr)
        index_count = len(node_index)
        print(f"✔ Parsed {index_count} nodes from index file.")
    except FileNotFoundError:
        print(f"❌ Error: Index file not found: {args.index_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error parsing index file: {e}", file=sys.stderr)
        sys.exit(1)
    print(f"✔ Index parsing took {time.time() - start_time:.2f} seconds.")


    # --- Load Node Sequences ---
    start_time = time.time()
    node_sequences = {}
    if args.load_cache:
        print(f"🔹 Loading node sequences from cache: {args.load_cache}...")
        try:
            with open(args.load_cache, 'r') as f_cache:
                node_sequences = json.load(f_cache)
            print(f"✔ Loaded {len(node_sequences)} sequences from cache.")
        except FileNotFoundError:
            print(f"❌ Error: Cache file not found: {args.load_cache}", file=sys.stderr)
            sys.exit(1)
        except json.JSONDecodeError as jde:
             print(f"❌ Error: Could not decode JSON from cache file {args.load_cache}: {jde}", file=sys.stderr)
             sys.exit(1)
        except Exception as e:
             print(f"❌ Error loading cache file: {e}", file=sys.stderr)
             sys.exit(1)
    else:
        print("❌ Error: Sequence cache file must be provided via --load-cache.", file=sys.stderr)
        sys.exit(1) # Require cache for this workflow
    print(f"✔ Sequence loading took {time.time() - start_time:.2f} seconds.")

    # --- Prepare Tasks ---
    start_time = time.time()
    print("🔹 Preparing tasks for parallel processing...")
    tasks = []
    nodes_missing_sequence = 0
    skipped_no_records = 0
    for node_id, (offset, n_records) in node_index.items():
        sequence = node_sequences.get(node_id)
        if n_records <= 0:
            skipped_no_records += 1
            continue # Skip nodes with no records listed in index
        if sequence is not None:
            tasks.append((node_id, offset, n_records, sequence))
        else:
            nodes_missing_sequence += 1

    if nodes_missing_sequence > 0:
        print(f"⚠️ Warning: Skipped {nodes_missing_sequence} nodes present in index but missing from sequence cache.", file=sys.stderr)
    if skipped_no_records > 0:
        print(f"⚠️ Warning: Skipped {skipped_no_records} nodes with 0 records according to index.", file=sys.stderr)

    total_tasks = len(tasks)
    if total_tasks == 0:
        print("❌ Error: No valid tasks to process. Check index and cache files.", file=sys.stderr)
        sys.exit(1)
    print(f"✔ Prepared {total_tasks} tasks in {time.time() - start_time:.2f} seconds.")

    # --- Execute in Parallel using submit/as_completed ---
    print(f"🔹 Processing {total_tasks} nodes using {args.workers} workers...")
    results = {}
    start_proc_time = time.time()
    processed_count = 0
    errors_count = 0

    # Use context manager for ProcessPoolExecutor
    try:
        with ProcessPoolExecutor(max_workers=args.workers,
                                 initializer=init_worker,
                                 initargs=(args.dat_file,)) as executor:

            # Submit tasks and store future keyed by node_id for error reporting
            # Note: If node_id is not unique in tasks (shouldn't happen here), this needs adjustment
            futures_map = {executor.submit(process_node_parallel, task): task[0] for task in tasks}

            print(f"--- Submitted {len(futures_map)} tasks. Waiting for results...")

            for future in as_completed(futures_map):
                node_id = futures_map[future] # Get node_id associated with this future
                try:
                    # Get the result. This will re-raise exceptions *if* they occurred
                    # *during pickling* or *executor shutdown*, but not exceptions
                    # caught *inside* process_node_parallel (which returns normally).
                    _returned_node_id, pileup_dict = future.result()

                    # Check if the worker returned a valid dictionary (not the error indicator '{}')
                    if pileup_dict: # Check if dict is not empty
                        results[node_id] = pileup_dict
                    # else: Worker returned {}, indicating caught error - already logged by worker

                except Exception as exc:
                    # Catch exceptions raised *outside* the worker's main try/except,
                    # e.g., worker process crashing, pickling errors.
                    print(f'❌❌❌ Main thread caught exception for future associated with Node ID {node_id}: {exc}', file=sys.stderr)
                    # Log traceback if desired for these rarer errors
                    # traceback.print_exc(file=sys.stderr)
                    errors_count += 1

                # --- Progress Update ---
                processed_count += 1
                if processed_count % 1000 == 0 or processed_count == total_tasks:
                    elapsed = time.time() - start_proc_time
                    rate = processed_count / elapsed if elapsed > 0 else 0
                    print(f"✔ Processed {processed_count}/{total_tasks} nodes ({errors_count} errors reported). Rate: {rate:.1f} nodes/s. Elapsed: {elapsed:.2f}s", flush=True)

    except Exception as pool_exc:
         # Errors during pool creation or shutdown itself
         print(f"\n❌ An error occurred within the ProcessPoolExecutor context: {pool_exc}", file=sys.stderr)
         print("⚠️ Attempting to write any partial results obtained...", file=sys.stderr)
         traceback.print_exc(file=sys.stderr)


    # --- Final Timing and Summary ---
    total_proc_time = time.time() - start_proc_time
    print(f"\n✔ Parallel processing finished in {total_proc_time:.2f} seconds.")
    print(f"Processed {processed_count} tasks. Found pileups for {len(results)} nodes.")
    if errors_count > 0 or processed_count < total_tasks :
         print(f"⚠️ Encountered {errors_count} errors during future processing. {total_tasks - processed_count} tasks may not have completed.", file=sys.stderr)


    # --- Save Results ---
    start_time = time.time()
    print(f"🔹 Saving {len(results)} pileup results to {args.output_json}...")
    try:
        with open(args.output_json, 'w') as f_out:
            json.dump(results, f_out) # Consider indent=2 for readability if files aren't huge
        print(f"✔ Results saved successfully.")
    except IOError as ioe:
         print(f"❌ Error: Could not write output file {args.output_json}: {ioe}", file=sys.stderr)
    except Exception as e:
         print(f"❌ Error saving results to JSON: {e}", file=sys.stderr)
    print(f"✔ Result saving took {time.time() - start_time:.2f} seconds.")

    print(f"Total execution time: {time.time() - (start_proc_time - (time.time() - start_time) )} seconds") # Approximate total
    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    main()