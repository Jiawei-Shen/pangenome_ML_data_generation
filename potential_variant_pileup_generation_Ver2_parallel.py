#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import math  # For math.ceil
# import numpy as np # Not strictly needed if outputting lists of lists
from collections import defaultdict
import re
from concurrent.futures import ProcessPoolExecutor, as_completed

# --- Constants (same as your script) ---
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4,
                 'a': 0, 'c': 1, 'g': 2, 't': 3, 'n': 4}
N_INDEX = BASE_TO_INDEX['N']
DEFAULT_WINDOW_SIZE = 100

worker_dat_file = None  # Global for worker


# --- Helper Functions (ellipsis indicates these are from your script, assumed unchanged) ---
def reverse_complement(sequence):
    # ... (implementation from your script) ...
    if not isinstance(sequence, str): return ""
    try:
        complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return sequence.translate(complement_map)[::-1]
    except Exception:
        return ""


def parse_cigar(cigar_string):
    # ... (implementation from your script) ...
    if not cigar_string or cigar_string == '*': return []
    try:
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr, flush=True)
        return []


def detect_variants_core(ref_seq, read_seq, start_offset, cigar_ops):
    # ... (implementation from your script) ...
    ref_len = len(ref_seq);
    read_len = len(read_seq)
    ref_ptr = start_offset;
    read_ptr = 0
    last_valid_ref_ptr = start_offset - 1
    try:
        for length, op in cigar_ops:
            if op in ('M', '=', 'X'):
                for i in range(length):
                    current_ref_pos = ref_ptr + i;
                    current_read_pos = read_ptr + i
                    if current_ref_pos >= ref_len or current_read_pos >= read_len: return
                    ref_base = ref_seq[current_ref_pos].upper();
                    read_base = read_seq[current_read_pos].upper()
                    last_valid_ref_ptr = current_ref_pos
                    if op != '=' and ref_base != 'N' and read_base != 'N' and ref_base != read_base:
                        yield (current_ref_pos, 'X', ref_base, read_base)
                ref_ptr += length;
                read_ptr += length
            elif op == 'I':
                if read_ptr + length > read_len: return
                inserted_bases = read_seq[read_ptr: read_ptr + length].upper()
                yield (last_valid_ref_ptr, 'I', '-', inserted_bases)
                read_ptr += length
            elif op == 'D':
                start_del_pos = ref_ptr
                if ref_ptr + length > ref_len: return
                deleted_bases = ref_seq[ref_ptr: ref_ptr + length].upper()
                yield (start_del_pos, 'D', deleted_bases, '-')
                ref_ptr += length;
                last_valid_ref_ptr = ref_ptr - 1
            elif op == 'S':
                read_ptr += length
            elif op == 'N':
                ref_ptr += length; last_valid_ref_ptr = ref_ptr - 1
            elif op in ('H', 'P'):
                pass
    except IndexError:
        print(
            f"⚠️ Warning: IndexError during variant detection. RefLen={ref_len}, ReadLen={read_len}, Offset={start_offset}, CIGAR='{''.join(map(str, cigar_ops))}', RefPtr={ref_ptr}, ReadPtr={read_ptr}",
            file=sys.stderr, flush=True)
    except Exception as e:
        print(f"❌ Error: Unexpected error in detect_variants_core: {e}", file=sys.stderr, flush=True)


def create_pileup_row(read_seq, read_offset, cigar_ops, window_start, window_end):
    # ... (implementation from your script) ...
    window_len = window_end - window_start
    if window_len <= 0: return None
    pileup_row = [N_INDEX] * window_len;
    read_len = len(read_seq)
    ref_ptr = read_offset;
    read_ptr = 0;
    row_filled = False
    try:
        for length, op in cigar_ops:
            if ref_ptr >= window_end and not row_filled: return None
            if op in ('M', '=', 'X'):
                for i in range(length):
                    current_ref_pos = ref_ptr + i;
                    current_read_pos = read_ptr + i
                    if current_read_pos >= read_len: return pileup_row if row_filled else None
                    if window_start <= current_ref_pos < window_end:
                        pileup_idx = current_ref_pos - window_start
                        pileup_row[pileup_idx] = BASE_TO_INDEX.get(read_seq[current_read_pos].upper(), N_INDEX)
                        row_filled = True
                ref_ptr += length;
                read_ptr += length
            elif op == 'I':
                if read_ptr + length > read_len: return pileup_row if row_filled else None
                read_ptr += length
            elif op == 'D':
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    if window_start <= current_ref_pos < window_end:
                        pileup_row[current_ref_pos - window_start] = N_INDEX;
                        row_filled = True
                ref_ptr += length
            elif op == 'S':
                read_ptr += length
            elif op == 'N':
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    if window_start <= current_ref_pos < window_end:
                        pileup_row[current_ref_pos - window_start] = N_INDEX;
                        row_filled = True
                ref_ptr += length
            elif op in ('H', 'P'):
                pass
        return pileup_row if row_filled else None
    except IndexError:
        print(
            f"⚠️ Warning: IndexError during pileup creation. ReadLen={read_len}, Offset={read_offset}, CIGAR='{''.join(map(str, cigar_ops))}', RefPtr={ref_ptr}, ReadPtr={read_ptr}, Win=[{window_start}:{window_end}]",
            file=sys.stderr, flush=True); return pileup_row if row_filled else None
    except Exception as e:
        print(f"❌ Error: Unexpected error in create_pileup_row: {e}", file=sys.stderr, flush=True); return None


def parse_idx_file(idx_path):
    # ... (implementation from your script) ...
    node_index = {}
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4: print(f"❌ Error: Index file {idx_path} is too small.", file=sys.stderr,
                                    flush=True); sys.exit(1)
            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4: print(f"❌ Error: Could not read number of nodes from {idx_path}.",
                                               file=sys.stderr, flush=True); sys.exit(1)
            num_nodes_header = struct.unpack('<I', num_nodes_bytes)[0]
            expected_min_size = 4 + num_nodes_header * 22
            if file_size < expected_min_size: print(
                f"⚠️ Warning: Index file {idx_path} may be truncated. Expected {expected_min_size}, found {file_size}.",
                file=sys.stderr, flush=True)
            print(f"🔹 Reading index for up to {num_nodes_header} nodes (from header)...", flush=True)
            nodes_parsed_count = 0
            for i in range(num_nodes_header):
                record_bytes = f.read(22)
                if len(record_bytes) < 22: print(
                    f"❌ Error: Index file ended prematurely reading record {i + 1}/{num_nodes_header}.",
                    file=sys.stderr, flush=True); break
                node_id, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                if n_records > 0: node_index[node_id] = (offset, n_records); nodes_parsed_count += 1
            print(f"✔ Parsed {nodes_parsed_count} nodes with records from index file.")
            if nodes_parsed_count == 0 and num_nodes_header > 0:
                print(f"⚠️ Warning: Header {num_nodes_header} nodes, but 0 with records found.", file=sys.stderr,
                      flush=True)
            elif len(
                    node_index) != num_nodes_header and nodes_parsed_count > 0 and file_size >= expected_min_size:  # Check only if file wasn't truncated
                print(f"⚠️ Warning: Header {num_nodes_header} nodes, parsed {nodes_parsed_count} with records.",
                      file=sys.stderr, flush=True)
            return node_index
    except FileNotFoundError:
        print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr, flush=True); sys.exit(1)
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr, flush=True); sys.exit(1)


def load_node_sequences_from_gfa(gfa_path, target_node_ids_set):
    # ... (implementation from your script) ...
    node_sequences = {};
    nodes_found_count = 0
    current_target_ids_set = set(target_node_ids_set)
    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file: {gfa_path} for up to {len(current_target_ids_set)} target sequences...",
                  flush=True)
            line_counter = 0
            for line in f:
                line_counter += 1
                if not current_target_ids_set: break
                if line_counter % 1_000_000 == 0: print(
                    f"  GFA Checked {line_counter:,} lines... found {nodes_found_count}/{len(target_node_ids_set)} target sequences.")
                    # end='\r', flush=True)
                if not line.startswith('S\t'): continue
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                try:
                    nid = int(parts[1])
                except ValueError:
                    continue
                if nid in current_target_ids_set:
                    seq = parts[2].upper()
                    if seq != '*' and re.match(r'^[ACGTN]+$', seq):
                        node_sequences[nid] = seq;
                        nodes_found_count += 1
                        current_target_ids_set.remove(nid)
            print(f"\n✔ Finished GFA. Checked {line_counter:,} lines.", flush=True)
            print(f"✔ Loaded sequences for {len(node_sequences)} target nodes from GFA.", flush=True)
            if current_target_ids_set: print(
                f"⚠️ Warning (GFA): Did not find sequences for {len(current_target_ids_set)} target node(s), e.g., {list(current_target_ids_set)[:5]}.",
                file=sys.stderr, flush=True)
    except FileNotFoundError:
        print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr, flush=True); sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr, flush=True); sys.exit(1)
    return node_sequences


def init_worker(dat_file_path):
    # ... (implementation from your script) ...
    global worker_dat_file
    try:
        worker_dat_file = open(dat_file_path, 'rb')
    except FileNotFoundError:
        print(f"❌ Error [Worker {os.getpid()}]: DAT file not found: {dat_file_path}", file=sys.stderr,
              flush=True); sys.exit(1)
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] opening DAT: {e}", file=sys.stderr, flush=True); sys.exit(1)


def process_node_parallel(task_args):
    # ... (implementation from your script, assumed to be the core logic for one task) ...
    node_id, file_read_offset, n_records, node_sequence, window_size = task_args
    global worker_dat_file
    if worker_dat_file is None: return node_id, {"error": "worker_dat_file not initialized"}
    if not node_sequence: return node_id, {"error": "node_sequence is empty"}

    node_len = len(node_sequence);
    half_window = window_size // 2
    reads_by_variant = defaultdict(list);
    mapq_filter = 10
    final_pileups = {}
    try:
        worker_dat_file.seek(file_read_offset)
        for _ in range(n_records):
            try:
                data = worker_dat_file.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE: break
                read_offset_on_node, raw_seq, _, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)
                if mapq < mapq_filter: continue
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
                cigar_str = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
                strand_char = strand_byte.decode('ascii', errors='ignore')
                if not seq or not cigar_str or cigar_str == '*': continue
                oriented_seq = seq;
                current_offset = read_offset_on_node
                if strand_char == '-':
                    oriented_seq = reverse_complement(seq)
                    if not oriented_seq: continue
                cigar_ops = parse_cigar(cigar_str)
                if not cigar_ops: continue
                variants_found_in_read = list(
                    detect_variants_core(node_sequence, oriented_seq, current_offset, cigar_ops))
                for vpos, vtype, ref_base, alt_base in variants_found_in_read:
                    if (vtype == 'I' and -1 <= vpos < node_len) or \
                            (vtype in ('D', 'X') and 0 <= vpos < node_len):
                        variant_key = f"{vpos}_{vtype}_{ref_base}_{alt_base}"
                        reads_by_variant[variant_key].append((current_offset, oriented_seq, cigar_ops))
            except struct.error:
                break
            except UnicodeDecodeError:
                continue
            except Exception:
                continue
    except IOError as e:
        return node_id, {"error": f"IOError: {e}"}
    except Exception as e:
        return node_id, {"error": f"Outer error: {e}"}

    if not reads_by_variant: return node_id, {}
    for variant_key, supporting_reads in reads_by_variant.items():
        if not supporting_reads: continue
        try:
            vpos = int(variant_key.split('_')[0])
        except (ValueError, IndexError):
            continue
        window_start, window_end = 0, 0
        if node_len <= window_size:
            window_start, window_end = 0, node_len
        else:
            center_pos = vpos;
            window_start = max(0, center_pos - half_window)
            window_end = min(node_len, window_start + window_size)
            window_start = max(0, window_end - window_size)
        if window_start >= window_end: continue
        pileup_rows = []
        for ro, rs, rco in supporting_reads:
            row = create_pileup_row(rs, ro, rco, window_start, window_end)
            if row: pileup_rows.append(row)
        if pileup_rows: final_pileups[variant_key] = pileup_rows
    return node_id, final_pileups


# --- NEW Generator function ---
def generate_task_args(node_info_list, node_sequences_dict, window_size_arg):
    """
    Generator that yields task arguments one by one.
    This avoids creating a massive list of all task arguments in memory.
    """
    processed_count = 0
    for node_id, file_offset, n_records in node_info_list:
        sequence = node_sequences_dict.get(node_id)  # Use integer node_id
        if sequence:  # Only yield tasks for which we have a sequence
            processed_count += 1
            yield (node_id, file_offset, n_records, sequence, window_size_arg)
    print(f"\nℹ️ Task generator finished. Yielded {processed_count} tasks.", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Pileup generation.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dat", help=".dat file")
    parser.add_argument("idx", help=".idx file")
    parser.add_argument("output", help="JSON output file")
    parser.add_argument("--gfa", help="GFA file (if not using cache)")
    parser.add_argument("--load-cache", help="JSON cache file for node sequences")
    parser.add_argument("--save-cache", help="JSON cache file to save node sequences")
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW_SIZE, help="Pileup window size")
    parser.add_argument("-w", "--workers", type=int, default=os.cpu_count(), help="Number of worker processes")
    # Renamed --chunksize to --progress-interval as it's now used for reporting frequency
    parser.add_argument("--progress-interval", type=int, default=1000,
                        help="Interval for progress updates (number of tasks).")
    parser.add_argument("--max-active-futures", type=int, default=0,
                        help="Max futures in flight for as_completed (0 for num_workers * 2, min 1)")

    args = parser.parse_args()

    # --- Basic File/Arg Validation ---
    if not os.path.isfile(args.dat): print(f"❌ DAT not found: {args.dat}", file=sys.stderr, flush=True); sys.exit(1)
    if not os.path.isfile(args.idx): print(f"❌ IDX not found: {args.idx}", file=sys.stderr, flush=True); sys.exit(1)
    if not args.load_cache and not args.gfa: print("❌ Must provide --gfa or --load-cache.", file=sys.stderr,
                                                   flush=True); sys.exit(1)
    if args.load_cache and not os.path.isfile(args.load_cache) and not args.gfa: print(
        f"❌ Cache {args.load_cache} not found and no GFA fallback.", file=sys.stderr, flush=True); sys.exit(1)
    if args.gfa and not os.path.isfile(args.gfa) and not (args.load_cache and os.path.isfile(args.load_cache)): print(
        f"❌ GFA {args.gfa} not found and no valid cache.", file=sys.stderr, flush=True); sys.exit(1)
    if args.window <= 0: print("❌ Window must be > 0.", file=sys.stderr, flush=True); sys.exit(1)

    print(f"Script started at {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    overall_start_time = time.time()

    # 1. Load full node_index (node_id -> (offset, n_records) dict)
    node_index_dict = parse_idx_file(args.idx)
    if not node_index_dict: print("❌ No nodes in index. Exiting.", file=sys.stderr, flush=True); sys.exit(1)

    # Create a list of (node_id, offset, n_records) for easier iteration by the generator
    all_indexed_nodes_info = []
    for nid_int, data_tuple in node_index_dict.items():  # node_index_dict keys are already int from parse_idx
        all_indexed_nodes_info.append((nid_int, data_tuple[0], data_tuple[1]))

    if not all_indexed_nodes_info:
        print("❌ No processable nodes found in index (after filtering for n_records > 0). Exiting.", file=sys.stderr,
              flush=True)
        sys.exit(1)

    target_ids_for_sequences = {info[0] for info in all_indexed_nodes_info}  # Set of all node IDs from index

    # 2. Load ALL necessary node_sequences ONCE into a dictionary
    # (This assumes the complete all_node_sequences dictionary can fit in main process memory)
    all_node_sequences = {}  # Integer node_id -> sequence string
    loaded_from_cache = False
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"🔹 Loading node sequences from cache: {args.load_cache}...", flush=True)
        cache_load_start = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                loaded_data = json.load(cf)
            for k_str, seq_val in loaded_data.items():
                try:
                    k_int = int(k_str)
                except ValueError:
                    continue
                if k_int in target_ids_for_sequences: all_node_sequences[k_int] = seq_val
            print(
                f"✔ Loaded {len(all_node_sequences)} relevant sequences from cache in {time.time() - cache_load_start:.2f}s.",
                flush=True)
            loaded_from_cache = True
        except Exception as e:
            print(f"❌ Error loading/parsing cache {args.load_cache}: {e}. Will try GFA if provided.", file=sys.stderr,
                  flush=True)
            all_node_sequences = {}  # Reset

    if not loaded_from_cache or len(all_node_sequences) < len(target_ids_for_sequences):
        if args.gfa:
            print(
                f"🔹 {'Attempting to load remaining' if loaded_from_cache and len(all_node_sequences) > 0 else 'Loading ALL'} sequences from GFA: {args.gfa}...",
                flush=True)
            still_needed_ids = target_ids_for_sequences - set(
                all_node_sequences.keys()) if loaded_from_cache else target_ids_for_sequences

            if still_needed_ids:
                gfa_sequences = load_node_sequences_from_gfa(args.gfa, still_needed_ids)
                all_node_sequences.update(gfa_sequences)  # Add newly loaded sequences
            elif not still_needed_ids and loaded_from_cache:  # All were already in cache
                print("ℹ️ All required sequences were already found in cache.", flush=True)

            # Save if cache was not loaded OR if GFA provided new sequences AND save_cache is specified
            if args.save_cache and (not loaded_from_cache or (still_needed_ids and gfa_sequences)):
                print(f"🔹 Saving {len(all_node_sequences)} sequences to cache: {args.save_cache}...", flush=True)
                try:
                    with open(args.save_cache, 'w') as cf:
                        json.dump({str(k): v for k, v in all_node_sequences.items()}, cf, indent=2)
                    print(f"✔ Saved cache to {args.save_cache}", flush=True)
                except Exception as e:
                    print(f"❌ Error saving cache: {e}", file=sys.stderr, flush=True)
        elif not all_node_sequences:  # No GFA and cache failed/empty and not loaded
            print("❌ No sequence source (valid cache or GFA) available for required nodes. Exiting.", file=sys.stderr,
                  flush=True)
            sys.exit(1)

    # Filter node_info_list to only those for which we actually have sequences
    # This is important because GFA loading might not find all target_ids
    process_node_info_list = [info for info in all_indexed_nodes_info if info[0] in all_node_sequences]
    total_tasks_to_process = len(process_node_info_list)

    if not process_node_info_list:
        print("❌ No nodes to process after matching index with available sequences. Exiting.", flush=True)
        sys.exit(1)

    if total_tasks_to_process < len(all_indexed_nodes_info):
        print(
            f"ℹ️ Will process {total_tasks_to_process} nodes for which sequences were loaded (out of {len(all_indexed_nodes_info)} indexed nodes).",
            flush=True)

    # 3. Create the task generator
    task_arg_generator = generate_task_args(process_node_info_list, all_node_sequences, args.window)

    # --- Execute in Parallel using bounded as_completed with task generator ---
    overall_results = {}
    future_to_node_id = {}  # Map Future objects to their original node_id for error reporting

    num_workers = min(args.workers, os.cpu_count() or 1, total_tasks_to_process)
    if num_workers <= 0: num_workers = 1
    max_active_futures = args.max_active_futures if args.max_active_futures > 0 else num_workers * 2
    max_active_futures = max(1, max_active_futures)  # Ensure at least 1

    progress_update_interval = args.progress_interval

    print(f"🔹 Concurrently submitting and processing {total_tasks_to_process} tasks using {num_workers} workers.",
          flush=True)
    print(f"   Max active futures to manage: {max_active_futures}", flush=True)

    nodes_completed_count = 0
    tasks_submitted_count = 0

    with ProcessPoolExecutor(max_workers=num_workers, initializer=init_worker, initargs=(args.dat,)) as executor:
        # Submit initial batch of tasks from the generator
        for _ in range(max_active_futures):
            try:
                task_args = next(task_arg_generator)
                future = executor.submit(process_node_parallel, task_args)
                future_to_node_id[future] = task_args[0]  # Map future to node_id
                tasks_submitted_count += 1
            except StopIteration:  # No more tasks from generator
                break  # Initial batch might be smaller than max_active_futures

        if tasks_submitted_count > 0:
            print(f"  Initial submission: {tasks_submitted_count} tasks in flight.", flush=True)

        while future_to_node_id:  # Loop as long as there are futures being processed
            # Get a snapshot of current active futures for as_completed
            active_futures_list = list(future_to_node_id.keys())
            if not active_futures_list: break  # Should be caught by `while future_to_node_id`

            for future in as_completed(active_futures_list):
                nodes_completed_count += 1
                original_node_id = future_to_node_id.pop(future)  # Process one and remove from tracking

                try:
                    _ret_node_id, pileup_dict = future.result()  # This will raise exception if worker failed

                    if _ret_node_id != original_node_id:  # Sanity check
                        print(f"⚠️ Main: Node ID mismatch! Original: {original_node_id}, Returned: {_ret_node_id}",
                              file=sys.stderr, flush=True)

                    # Store result if valid
                    if pileup_dict and not (isinstance(pileup_dict, dict) and "error" in pileup_dict):
                        overall_results[original_node_id] = pileup_dict
                    elif isinstance(pileup_dict, dict) and "error" in pileup_dict:
                        print(f"  ℹ️ Worker for node {original_node_id} reported error: {pileup_dict['error']}",
                              file=sys.stderr, flush=True)

                except Exception as exc:
                    print(f"\n❌ Task for node {original_node_id} resulted in an exception: {exc}", file=sys.stderr,
                          flush=True)
                    # import traceback; traceback.print_exc(file=sys.stderr) # Uncomment for full worker traceback

                # Try to submit a new task from the generator if we have capacity
                # and tasks are remaining.
                if len(future_to_node_id) < max_active_futures:
                    try:
                        if tasks_submitted_count < total_tasks_to_process:  # Check against total actual tasks
                            new_task_args = next(task_arg_generator)
                            new_future = executor.submit(process_node_parallel, new_task_args)
                            future_to_node_id[new_future] = new_task_args[0]
                            tasks_submitted_count += 1
                    except StopIteration:
                        # Generator is exhausted, no more tasks to submit
                        pass

                # Progress Reporting
                if nodes_completed_count % progress_update_interval == 0 or \
                        nodes_completed_count == total_tasks_to_process or \
                        not future_to_node_id:  # Update on the very last one too
                    elapsed = time.time() - overall_start_time
                    rate = nodes_completed_count / elapsed if elapsed > 0 else 0
                    # ETA calculation is more complex here as submission isn't all at once
                    print(
                        f"  Progress: {nodes_completed_count}/{total_tasks_to_process} done. {len(future_to_node_id)} active. {tasks_submitted_count} submitted. Rate: {rate:.1f}/s. Elapsed: {elapsed:.2f}s",
                        end='\r', flush=True)

                # After processing one future from as_completed, we break from this inner loop.
                # This allows the outer 'while future_to_node_id:' loop to re-evaluate
                # and ensures as_completed is called on the most current set of active futures.
                break
            else:
                # This 'else' belongs to the 'for future in as_completed(...)' loop.
                # It executes if the as_completed loop finishes without any 'break'.
                # This means all futures in the current `active_futures_list` snapshot were processed.
                # If `future_to_node_id` is now empty, the outer `while` loop will terminate.
                if not future_to_node_id:
                    break

    print()  # Newline after progress bar
    total_script_elapsed = time.time() - overall_start_time
    print(f"\n🏁 Finished processing in {total_script_elapsed:.2f} seconds.")
    print(f"🔹 Found pileups for {len(overall_results)} nodes overall.")

    print(f"🔹 Writing {len(overall_results)} node results to JSON: {args.output}")
    try:
        with open(args.output, 'w') as out_f:
            json.dump({str(k): v for k, v in overall_results.items()}, out_f, indent=2)
        print(f"✔ Output written to {args.output}")
    except Exception as e:
        print(f"❌ Error writing output JSON: {e}", file=sys.stderr, flush=True); sys.exit(1)

    print(f"✅ Done. Total time: {time.time() - overall_start_time:.2f} seconds.", flush=True)


if __name__ == '__main__':
    main()