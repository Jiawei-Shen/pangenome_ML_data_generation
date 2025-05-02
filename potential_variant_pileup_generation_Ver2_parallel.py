#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np # Using numpy for efficient pileup matrix creation (optional)
from collections import defaultdict
import re # Import re at the top level
from concurrent.futures import ProcessPoolExecutor
# from os import cpu_count # Can uncomment to default workers

# ─────────────────────────────────────────────────────────────────────────────
# Constants
RECORD_STRUCT = struct.Struct("<h150s150s20shc") # From original code
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4,
                 'a': 0, 'c': 1, 'g': 2, 't': 3, 'n': 4}
N_INDEX = BASE_TO_INDEX['N']

# Default window size
DEFAULT_WINDOW_SIZE = 100

# Global for worker process state (file handle)
worker_dat_file = None

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions (Unchanged from previous version)

def reverse_complement(sequence):
    """Computes the reverse complement of a DNA sequence."""
    if not isinstance(sequence, str):
        return ""
    try:
        complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return sequence.translate(complement_map)[::-1]
    except Exception:
        return ""

def parse_cigar(cigar_string):
    """Parses a CIGAR string into a list of (length, operation) tuples."""
    if not cigar_string or cigar_string == '*':
        return []
    try:
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr)
        return []

def detect_variants_core(ref_seq, read_seq, start_offset, cigar_ops):
    """
    Core CIGAR-aware variant detection logic. Yields specific variants.
    Yields: tuple: (ref_pos, variant_type, ref_base, alt_base)
    """
    ref_len = len(ref_seq)
    read_len = len(read_seq)
    ref_ptr = start_offset
    read_ptr = 0
    last_valid_ref_ptr = start_offset - 1

    try:
        for length, op in cigar_ops:
            if op in ('M', '=', 'X'):
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    current_read_pos = read_ptr + i
                    if current_ref_pos >= ref_len or current_read_pos >= read_len: return
                    ref_base = ref_seq[current_ref_pos].upper()
                    read_base = read_seq[current_read_pos].upper()
                    last_valid_ref_ptr = current_ref_pos
                    if op != '=' and ref_base != 'N' and read_base != 'N' and ref_base != read_base:
                        yield (current_ref_pos, 'X', ref_base, read_base)
                ref_ptr += length
                read_ptr += length
            elif op == 'I':
                if read_ptr + length > read_len: return
                inserted_bases = read_seq[read_ptr : read_ptr + length].upper()
                yield (last_valid_ref_ptr, 'I', '-', inserted_bases)
                read_ptr += length
            elif op == 'D':
                start_del_pos = ref_ptr
                if ref_ptr + length > ref_len: return
                deleted_bases = ref_seq[ref_ptr : ref_ptr + length].upper()
                yield (start_del_pos, 'D', deleted_bases, '-')
                ref_ptr += length
                last_valid_ref_ptr = ref_ptr - 1
            elif op == 'S':
                read_ptr += length
            elif op == 'N':
                 ref_ptr += length
                 last_valid_ref_ptr = ref_ptr - 1
            elif op in ('H', 'P'):
                 pass
    except IndexError:
        print(f"⚠️ Warning: IndexError during variant detection. RefLen={ref_len}, ReadLen={read_len}, Offset={start_offset}, CIGAR='{''.join(map(str,cigar_ops))}', RefPtr={ref_ptr}, ReadPtr={read_ptr}", file=sys.stderr)
    except Exception as e:
        print(f"❌ Error: Unexpected error in detect_variants_core: {e}", file=sys.stderr)

def create_pileup_row(read_seq, read_offset, cigar_ops, window_start, window_end):
    """
    Creates a single pileup row (list of base indices) based on CIGAR alignment
    within a specified reference window.
    """
    window_len = window_end - window_start
    if window_len <= 0: return None

    pileup_row = [N_INDEX] * window_len
    read_len = len(read_seq)
    ref_ptr = read_offset
    read_ptr = 0
    row_filled = False

    try:
        for length, op in cigar_ops:
            if ref_ptr >= window_end and not row_filled: return None

            if op in ('M', '=', 'X'):
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    current_read_pos = read_ptr + i
                    if current_read_pos >= read_len: return pileup_row if row_filled else None
                    if window_start <= current_ref_pos < window_end:
                        pileup_idx = current_ref_pos - window_start
                        base = read_seq[current_read_pos]
                        pileup_row[pileup_idx] = BASE_TO_INDEX.get(base.upper(), N_INDEX)
                        row_filled = True
                ref_ptr += length
                read_ptr += length
            elif op == 'I':
                if read_ptr + length > read_len: return pileup_row if row_filled else None
                read_ptr += length
            elif op == 'D':
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    if window_start <= current_ref_pos < window_end:
                        pileup_idx = current_ref_pos - window_start
                        pileup_row[pileup_idx] = N_INDEX
                        row_filled = True
                ref_ptr += length
            elif op == 'S':
                read_ptr += length
            elif op == 'N':
                 for i in range(length):
                    current_ref_pos = ref_ptr + i
                    if window_start <= current_ref_pos < window_end:
                         pileup_idx = current_ref_pos - window_start
                         pileup_row[pileup_idx] = N_INDEX
                         row_filled = True
                 ref_ptr += length
            elif op in ('H', 'P'):
                 pass
        return pileup_row if row_filled else None
    except IndexError:
         print(f"⚠️ Warning: IndexError during pileup creation. ReadLen={read_len}, Offset={read_offset}, CIGAR='{''.join(map(str,cigar_ops))}', RefPtr={ref_ptr}, ReadPtr={read_ptr}, Win=[{window_start}:{window_end}]", file=sys.stderr)
         return pileup_row if row_filled else None
    except Exception as e:
         print(f"❌ Error: Unexpected error in create_pileup_row: {e}", file=sys.stderr)
         return None

def parse_idx_file(idx_path):
     # (Function unchanged - kept for brevity)
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
                 print(f"⚠️ Warning: Index file {idx_path} may be truncated. Expected at least {expected_min_size} bytes for {num_nodes} nodes, found {file_size}.", file=sys.stderr)
             print(f"🔹 Reading index for {num_nodes} nodes...")
             nodes_parsed_count = 0
             for i in range(num_nodes):
                 record_bytes = f.read(22)
                 if len(record_bytes) < 22:
                      print(f"❌ Error: Index file ended prematurely while reading record {i+1}/{num_nodes}.", file=sys.stderr)
                      break
                 node_id, offset, block_size, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                 if n_records > 0:
                    node_index[node_id] = (offset, n_records)
                    nodes_parsed_count += 1
             print(f"✔ Parsed {nodes_parsed_count} nodes with records from index file.")
             if nodes_parsed_count == 0 and num_nodes > 0 :
                  print(f"⚠️ Warning: Header indicated {num_nodes} nodes, but 0 nodes with records were found in the index data.", file=sys.stderr)
             elif len(node_index) != num_nodes and nodes_parsed_count > 0 :
                  print(f"⚠️ Warning: Header indicated {num_nodes} nodes, but parsed {nodes_parsed_count} nodes with records.", file=sys.stderr)
             return node_index
     except FileNotFoundError:
         print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr)
         sys.exit(1)
     except Exception as e:
         print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
         sys.exit(1)


def load_node_sequences_from_gfa(gfa_path, target_node_ids):
     # (Function unchanged - kept for brevity)
     node_sequences = {}
     target_node_set = set(target_node_ids)
     nodes_found_count = 0
     try:
         with open(gfa_path, 'r') as f:
             print(f"🔹 Reading GFA file: {gfa_path}")
             line_counter = 0
             for line in f:
                 line_counter += 1
                 if line_counter % 5_000_000 == 0:
                     print(f"  Checked {line_counter:,} lines... found {nodes_found_count}/{len(target_node_ids)} target sequences.")
                 if not line.startswith('S\t'): continue
                 parts = line.strip().split('\t')
                 if len(parts) < 3: continue
                 try: nid = int(parts[1])
                 except ValueError: continue
                 if nid in target_node_set:
                     seq = parts[2].upper()
                     if seq != '*' and re.match(r'^[ACGTN]+$', seq):
                         node_sequences[nid] = seq
                         nodes_found_count += 1
                         target_node_set.remove(nid)
                         if not target_node_set:
                            print(f"✔ Found all {len(target_node_ids)} target sequences after {line_counter:,} lines.")
                            break
             print(f"✔ Finished checking {line_counter:,} lines in GFA.")
             print(f"✔ Loaded sequences for {len(node_sequences)} target nodes.")
             if target_node_set:
                  print(f"⚠️ Warning: Did not find sequences for {len(target_node_set)} target node(s), e.g., {list(target_node_set)[:5]}.", file=sys.stderr)
     except FileNotFoundError:
         print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr)
         sys.exit(1)
     except Exception as e:
         print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr)
         sys.exit(1)
     return node_sequences

# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function

def init_worker(dat_file_path):
    """Initializer for each worker process: opens the .dat file."""
    global worker_dat_file
    try:
        worker_dat_file = open(dat_file_path, 'rb')
    except FileNotFoundError:
         print(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path}", file=sys.stderr)
         sys.exit(1)
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path}: {e}", file=sys.stderr)
        sys.exit(1)


def process_node_parallel(task_args):
    """
    Function executed by each worker process. Processes a single node.
    """
    # (Function unchanged - kept for brevity, assumes correctness from previous step)
    node_id, file_read_offset, n_records, node_sequence, window_size = task_args
    global worker_dat_file
    if worker_dat_file is None: return node_id, {}
    if not node_sequence: return node_id, {}

    node_len = len(node_sequence)
    half_window = window_size // 2
    reads_by_variant = defaultdict(list)
    processed_read_count = 0
    mapq_filter = 10

    try:
        worker_dat_file.seek(file_read_offset)
        for record_idx in range(n_records):
            try:
                data = worker_dat_file.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE: break
                read_offset_on_node, raw_seq, _, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)
                if mapq < mapq_filter: continue
                try:
                    seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
                    cigar_str = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
                    strand_char = strand_byte.decode('ascii')
                except UnicodeDecodeError: continue
                if not seq or not cigar_str or cigar_str == '*': continue

                processed_read_count += 1
                oriented_seq = seq
                current_offset = read_offset_on_node
                if strand_char == '-':
                    oriented_seq = reverse_complement(seq)
                    if not oriented_seq: continue
                    # NOTE: Keep offset/CIGAR relative to forward strand reference for simplicity

                cigar_ops = parse_cigar(cigar_str)
                if not cigar_ops: continue

                variants_found_in_read = list(detect_variants_core(node_sequence, oriented_seq, current_offset, cigar_ops))
                for vpos, vtype, ref_base, alt_base in variants_found_in_read:
                    if (vtype == 'I' and -1 <= vpos < node_len) or \
                       (vtype in ('D', 'X') and 0 <= vpos < node_len):
                        variant_key = f"{vpos}_{vtype}_{ref_base}_{alt_base}"
                        reads_by_variant[variant_key].append((current_offset, oriented_seq, cigar_ops))

            except struct.error as se: break
            except Exception as e_inner: continue
    except IOError as ioe: return node_id, {}
    except Exception as e_outer: return node_id, {}

    final_pileups = {}
    if not reads_by_variant: return node_id, {}

    for variant_key, supporting_reads in reads_by_variant.items():
        if not supporting_reads: continue
        try: vpos = int(variant_key.split('_')[0])
        except (ValueError, IndexError): continue

        if node_len <= window_size:
            window_start = 0
            window_end = node_len
            actual_window_len = node_len
        else:
            center_pos = vpos
            window_start = max(0, center_pos - half_window)
            window_end = min(node_len, window_start + window_size)
            window_start = max(0, window_end - window_size)
            actual_window_len = window_end - window_start

        if actual_window_len <= 0: continue

        pileup_rows = []
        for read_offset, read_seq, read_cigar_ops in supporting_reads: # Renamed cigar_ops
            row = create_pileup_row(read_seq, read_offset, read_cigar_ops, window_start, window_end)
            if row: pileup_rows.append(row)
        if pileup_rows:
            final_pileups[variant_key] = pileup_rows

    # print(node_id, final_pileups, '\n\n')
    return node_id, final_pileups

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Parallel variant-centered pileup generation from DAT/IDX files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("dat", help=".dat file path (read alignment data)")
    parser.add_argument("idx", help=".idx file path (index for .dat file)")
    parser.add_argument("output", help="JSON output file path for pileups")
    # --- Re-added Cache Arguments ---
    parser.add_argument("--gfa", help="GFA graph file path (needed if node sequence cache is not used/built)")
    parser.add_argument("--load-cache", help="Load node sequences from this JSON cache file")
    parser.add_argument("--save-cache", help="Save node sequences to this JSON cache file (used if --gfa is provided and cache is built)")
    # --- End Re-added Cache Arguments ---
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW_SIZE, help="Pileup window size around variant")
    parser.add_argument("-w", "--workers", type=int, default=os.cpu_count(), help="Number of worker processes to use")
    parser.add_argument("-c", "--chunksize", type=int, default=1, help="Approximate number of nodes per worker task")
    args = parser.parse_args()

    # --- Input Validation ---
    if not os.path.isfile(args.dat): print(f"❌ Error: DAT file not found: {args.dat}", file=sys.stderr); sys.exit(1)
    if not os.path.isfile(args.idx): print(f"❌ Error: Index file not found: {args.idx}", file=sys.stderr); sys.exit(1)
    # --- Cache/GFA Validation ---
    if not args.load_cache and not args.gfa:
        print("❌ Error: You must provide either a GFA file (`--gfa`) to build node sequences or a pre-built cache file (`--load-cache`).", file=sys.stderr)
        sys.exit(1)
    if args.load_cache and args.gfa:
        print("ℹ️ Info: Both --gfa and --load-cache provided. Will use the cache file if it exists.")
    if args.load_cache and not os.path.isfile(args.load_cache):
        # If cache specified but not found, *and* GFA *is* provided, proceed to build from GFA.
        if args.gfa:
            print(f"⚠️ Warning: Specified cache file '{args.load_cache}' not found. Will build sequences from GFA '{args.gfa}'.")
        else:
             print(f"❌ Error: Specified cache file to load does not exist: {args.load_cache}, and no GFA provided.", file=sys.stderr)
             sys.exit(1) # Cannot proceed without sequences
    if args.gfa and not os.path.isfile(args.gfa):
        # If GFA is needed (no cache specified or cache not found) but GFA doesn't exist.
        if not (args.load_cache and os.path.isfile(args.load_cache)):
             print(f"❌ Error: GFA file not found: {args.gfa}", file=sys.stderr)
             sys.exit(1)
    # --- End Cache/GFA Validation ---
    if args.window <= 0: print(f"❌ Error: Window size (--window={args.window}) must be positive.", file=sys.stderr); sys.exit(1)

    # --- Load Index ---
    print("🔹 Parsing index file...")
    start_time = time.time()
    node_index = parse_idx_file(args.idx)
    if not node_index: print("❌ Error: Failed to parse node index. Exiting.", file=sys.stderr); sys.exit(1)
    print(f"✔ Index parsing took {time.time() - start_time:.2f} seconds.")
    target_ids = node_index.keys() # Get target IDs needed from index

    # --- Load or Build Node Sequences (with Cache Logic) ---
    node_sequences = {}
    loaded_from_cache = False

    # Try loading from cache first if specified and exists
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"🔹 Loading node sequences from cache: {args.load_cache}...")
        start_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                loaded_data = json.load(cf)
                # Convert keys back to integers and filter for needed IDs
                node_sequences = {int(k): v for k, v in loaded_data.items() if int(k) in target_ids}
            print(f"✔ Loaded {len(node_sequences)} required sequences from cache in {time.time() - start_time:.2f} seconds.")
            loaded_from_cache = True
        except json.JSONDecodeError as jde:
             print(f"❌ Error decoding JSON from cache file {args.load_cache}: {jde}. Will attempt to build from GFA if provided.", file=sys.stderr)
        except Exception as e:
             print(f"❌ Error loading cache file {args.load_cache}: {e}. Will attempt to build from GFA if provided.", file=sys.stderr)

    # If not loaded from cache, try building from GFA if provided
    if not loaded_from_cache:
        if args.gfa: # GFA must be provided if not loading from valid cache
            print(f"🔹 Building node sequence cache from GFA: {args.gfa}...")
            start_time = time.time()
            node_sequences = load_node_sequences_from_gfa(args.gfa, target_ids)
            print(f"✔ Sequence loading from GFA took {time.time() - start_time:.2f} seconds.")

            # Save to cache if requested
            if args.save_cache:
                print(f"🔹 Saving {len(node_sequences)} node sequences to cache: {args.save_cache}...")
                start_time = time.time()
                try:
                    with open(args.save_cache, 'w') as cf:
                        # Save keys as strings in JSON
                        json.dump({str(k): v for k, v in node_sequences.items()}, cf)
                    print(f"✔ Saved cache in {time.time() - start_time:.2f} seconds.")
                except IOError as ioe:
                     print(f"❌ Error saving cache file {args.save_cache}: {ioe}", file=sys.stderr)
                except Exception as e:
                     print(f"❌ Error during cache saving to {args.save_cache}: {e}", file=sys.stderr)
        else:
            # This case should be prevented by initial validation, but is a fallback.
             print("❌ Error: Cannot load sequences. No valid cache file provided or found, and no GFA file specified.", file=sys.stderr)
             sys.exit(1)

    # Final check if sequences were obtained
    if not node_sequences:
         print("❌ Error: Failed to load or build any required node sequences. Exiting.", file=sys.stderr)
         sys.exit(1)

    # --- Prepare Tasks for Parallel Processing ---
    print("🔹 Preparing tasks for parallel processing...")
    tasks = []
    nodes_missing_sequence = 0
    skipped_nodes_no_records = 0

    for node_id, (file_read_offset, n_records) in node_index.items():
        if n_records <= 0:
            skipped_nodes_no_records += 1
            continue

        sequence = node_sequences.get(node_id) # Get sequence using integer key
        if sequence is not None:
            tasks.append((node_id, file_read_offset, n_records, sequence, args.window))
        else:
            # Node in index (with records) but sequence missing from cache/GFA
            nodes_missing_sequence += 1

    if skipped_nodes_no_records > 0: print(f"ℹ️ Note: {skipped_nodes_no_records} nodes were filtered due to 0 records in index.")
    if nodes_missing_sequence > 0: print(f"⚠️ Warning: Skipped {nodes_missing_sequence} nodes present in index (with records) but sequence data was not found/loaded.")
    if not tasks: print("❌ Error: No tasks to process (no nodes with both index records and sequence data found). Exiting.", file=sys.stderr); sys.exit(1)

    total_tasks = len(tasks)
    num_workers = min(args.workers, total_tasks, os.cpu_count() or 1)
    if num_workers <= 0: num_workers = 1
    print(f"🔹 Processing {total_tasks} nodes using {num_workers} workers (chunksize: {args.chunksize})...")

    # --- Execute in Parallel ---
    results = {}
    start_proc_time = time.time()
    nodes_processed_count = 0

    try:
        with ProcessPoolExecutor(max_workers=num_workers,
                                 initializer=init_worker,
                                 initargs=(args.dat,)) as executor:
            print("start")
            future_results = executor.map(process_node_parallel, tasks, chunksize=args.chunksize)
            print(future_results)
            for node_id, pileup_dict in future_results:
                 nodes_processed_count += 1
                 if pileup_dict:
                     results[node_id] = pileup_dict # Store using integer key
                 if nodes_processed_count % max(1, total_tasks // 20000) == 0 or nodes_processed_count == total_tasks:
                     elapsed = time.time() - start_proc_time
                     rate = nodes_processed_count / elapsed if elapsed > 0 else 0
                     print(f"  Processed {nodes_processed_count}/{total_tasks} nodes ({len(results)} with pileups). Rate: {rate:.1f} nodes/s. Elapsed: {elapsed:.2f}s")
    except Exception as pool_exc:
         print(f"\n❌ An error occurred during parallel processing: {pool_exc}", file=sys.stderr)
         print("⚠️ Attempting to write any partial results obtained...", file=sys.stderr)

    total_elapsed_time = time.time() - start_proc_time
    print(f"✔ Parallel processing finished in {total_elapsed_time:.2f} seconds.")
    print(f"✔ Found pileups for {len(results)} nodes.")

    # --- Write Output ---
    print(f"🔹 Writing {len(results)} node results to JSON: {args.output}")
    start_write_time = time.time()
    try:
        with open(args.output, 'w') as out_f:
            # Save node IDs as strings in the final JSON output
            json.dump({str(k): v for k, v in results.items()}, out_f, indent=2)
        write_elapsed_time = time.time() - start_write_time
        print(f"✔ Output written in {write_elapsed_time:.2f} seconds.")
        print(f"✅ Done. Output saved to {args.output}")
    except IOError as ioe: print(f"❌ Error writing output JSON to {args.output}: {ioe}", file=sys.stderr); sys.exit(1)
    except Exception as e: print(f"❌ Unexpected error writing output JSON: {e}", file=sys.stderr); sys.exit(1)

if __name__ == '__main__':
    main()