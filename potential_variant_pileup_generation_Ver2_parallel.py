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

# --- Constants ---
RECORD_STRUCT = struct.Struct("<h150s150s20shc")  # From original code
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4,
                 'a': 0, 'c': 1, 'g': 2, 't': 3, 'n': 4}
N_INDEX = BASE_TO_INDEX['N']

# Default window size
DEFAULT_WINDOW_SIZE = 100
# Default for how many nodes' tasks to prepare and submit in one processing chunk
DEFAULT_PROCESSING_CHUNK_SIZE = 50000

# Global for worker process state (file handle)
worker_dat_file = None


# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions (largely unchanged from your provided script version)

def reverse_complement(sequence):
    """Computes the reverse complement of a DNA sequence."""
    if not isinstance(sequence, str): return ""
    try:
        complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return sequence.translate(complement_map)[::-1]
    except Exception:
        return ""


def parse_cigar(cigar_string):
    """Parses a CIGAR string into a list of (length, operation) tuples."""
    if not cigar_string or cigar_string == '*': return []
    try:
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr, flush=True)
        return []


def detect_variants_core(ref_seq, read_seq, start_offset, cigar_ops):
    """
    Core CIGAR-aware variant detection logic. Yields specific variants.
    Yields: tuple: (ref_pos, variant_type, ref_base, alt_base)
    """
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
    """
    Creates a single pileup row (list of base indices) based on CIGAR alignment
    within a specified reference window.
    """
    window_len = window_end - window_start
    if window_len <= 0: return None
    pileup_row = [N_INDEX] * window_len;
    read_len = len(read_seq)
    ref_ptr = read_offset;
    read_ptr = 0;
    row_filled = False
    try:
        for length, op in cigar_ops:
            if ref_ptr >= window_end and not row_filled: return None  # Optimization
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
    node_index = {}  # Using dict for node_id -> (offset, n_records)
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
            for i in range(num_nodes_header):  # Iterate based on header count
                record_bytes = f.read(22)
                if len(record_bytes) < 22: print(
                    f"❌ Error: Index file ended prematurely reading record {i + 1}/{num_nodes_header}.",
                    file=sys.stderr, flush=True); break
                node_id, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                if n_records > 0:  # Only store/process nodes with records
                    node_index[node_id] = (offset, n_records)
                    nodes_parsed_count += 1
            print(f"✔ Parsed {nodes_parsed_count} nodes with records from index file.")
            if nodes_parsed_count == 0 and num_nodes_header > 0: print(
                f"⚠️ Warning: Header indicated {num_nodes_header} nodes, but 0 with records were found or processed.",
                file=sys.stderr, flush=True)
            return node_index
    except FileNotFoundError:
        print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr, flush=True); sys.exit(1)
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr, flush=True); sys.exit(1)


def load_node_sequences_from_gfa(gfa_path, target_node_ids_set):  # target_node_ids_set should be a set
    node_sequences = {};
    nodes_found_count = 0
    # Make a copy to modify if needed, or ensure caller passes a set if they need original
    current_target_ids_set = set(target_node_ids_set)
    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file: {gfa_path} for up to {len(current_target_ids_set)} target sequences...",
                  flush=True)
            line_counter = 0
            for line in f:
                line_counter += 1
                if not current_target_ids_set: break  # Optimization: if all targets found
                if line_counter % 1_000_000 == 0:  # Progress update
                    print(
                        f"  GFA Checked {line_counter:,} lines... found {nodes_found_count}/{len(target_node_ids_set)} target sequences.",
                        end='\r', flush=True)

                if not line.startswith('S\t'): continue
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                try:
                    nid = int(parts[1])
                except ValueError:
                    continue

                if nid in current_target_ids_set:
                    seq = parts[2].upper()
                    if seq != '*' and re.match(r'^[ACGTN]+$', seq):  # Basic sequence validation
                        node_sequences[nid] = seq
                        nodes_found_count += 1
                        current_target_ids_set.remove(nid)  # Remove from set of needed IDs
                    # else: print(f"Debug: Node {nid} in GFA has invalid sequence: {seq[:20]}...")
            print(f"\n✔ Finished GFA. Checked {line_counter:,} lines.", flush=True)  # Ensure newline after progress
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
    global worker_dat_file
    try:
        worker_dat_file = open(dat_file_path, 'rb')
        # print(f"[Worker {os.getpid()}] Initialized with {dat_file_path}", flush=True)
    except FileNotFoundError:
        print(f"❌ Error [Worker {os.getpid()}]: DAT file not found: {dat_file_path}", file=sys.stderr,
              flush=True); sys.exit(1)
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] opening DAT: {e}", file=sys.stderr, flush=True); sys.exit(1)


def process_node_parallel(task_args):
    # (This function is largely the same as the one you provided in the prompt,
    #  assuming its internal logic is what you want for processing a single node's data.
    #  Make sure it's robust.)
    node_id, file_read_offset, n_records, node_sequence, window_size = task_args
    global worker_dat_file
    if worker_dat_file is None: return node_id, {}  # Should have been initialized
    if not node_sequence: return node_id, {}  # Should have sequence

    node_len = len(node_sequence);
    half_window = window_size // 2
    reads_by_variant = defaultdict(list);
    mapq_filter = 10  # Or pass mapq_filter in task_args
    final_pileups = {}

    try:
        worker_dat_file.seek(file_read_offset)
        for _ in range(n_records):  # Use _ if record_idx is not used
            try:
                data = worker_dat_file.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE: break
                read_offset_on_node, raw_seq, _, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)

                if mapq < mapq_filter: continue

                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
                cigar_str = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
                strand_char = strand_byte.decode('ascii', errors='ignore')  # Added errors='ignore'

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
                    # Basic boundary check for variant position
                    if (vtype == 'I' and -1 <= vpos < node_len) or \
                            (vtype in ('D', 'X') and 0 <= vpos < node_len):
                        variant_key = f"{vpos}_{vtype}_{ref_base}_{alt_base}"
                        reads_by_variant[variant_key].append(
                            (current_offset, oriented_seq, cigar_ops))  # Store parsed ops
            except struct.error as se:
                # print(f"Debug: [Worker {os.getpid()}] Node {node_id}: Struct error unpacking record: {se}", flush=True)
                break  # Problem with this node's data block
            except UnicodeDecodeError:
                # print(f"Debug: [Worker {os.getpid()}] Node {node_id}: Unicode decode error.", flush=True)
                continue  # Skip problematic record
            except Exception:  # Catch any other unexpected error in loop
                # print(f"Debug: [Worker {os.getpid()}] Node {node_id}: Generic error in record loop.", flush=True)
                continue
    except IOError:
        # print(f"Debug: [Worker {os.getpid()}] Node {node_id}: IOError reading .dat file.", flush=True)
        return node_id, {}
    except Exception:
        # print(f"Debug: [Worker {os.getpid()}] Node {node_id}: Major error processing node.", flush=True)
        return node_id, {}

    if not reads_by_variant: return node_id, {}

    for variant_key, supporting_reads in reads_by_variant.items():
        if not supporting_reads: continue
        try:
            vpos = int(variant_key.split('_')[0])
        except (ValueError, IndexError):
            continue

        window_start, window_end = 0, 0
        if node_len <= window_size:  # Node is shorter or equal to window, use the whole node
            window_start, window_end = 0, node_len
        else:  # Node is longer, calculate window around vpos
            center_pos = vpos
            window_start = max(0, center_pos - half_window)
            window_end = min(node_len, window_start + window_size)
            # Adjust start if window was clamped at the end of the node, to maintain window_size if possible
            window_start = max(0, window_end - window_size)

        if window_start >= window_end: continue  # Skip if window is invalid

        pileup_rows = []
        for ro, rs, rco in supporting_reads:  # ro=read_offset, rs=read_seq, rco=read_cigar_ops
            row = create_pileup_row(rs, ro, rco, window_start, window_end)
            if row: pileup_rows.append(row)
        if pileup_rows: final_pileups[variant_key] = pileup_rows

    return node_id, final_pileups


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
    parser.add_argument("--processing-chunk-size", type=int, default=DEFAULT_PROCESSING_CHUNK_SIZE,
                        help="Number of nodes to load/prepare tasks for at a time.")
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

    # Create a list of (node_id, offset, n_records) for easier chunking by index later
    # This also filters out nodes that might be in GFA/cache but not index
    all_indexed_nodes_info = []
    for nid_int, data_tuple in node_index_dict.items():
        all_indexed_nodes_info.append((nid_int, data_tuple[0], data_tuple[1]))

    if not all_indexed_nodes_info:
        print("❌ No processable nodes found in index (after filtering for n_records > 0). Exiting.", file=sys.stderr,
              flush=True)
        sys.exit(1)

    target_ids_for_sequences = {info[0] for info in all_indexed_nodes_info}

    # 2. Load ALL necessary node_sequences ONCE into a dictionary
    all_node_sequences = {}  # Integer node_id -> sequence string
    loaded_from_cache = False
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"🔹 Loading node sequences from cache: {args.load_cache}...", flush=True)
        cache_load_start = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                loaded_data = json.load(cf)
            # Filter for only needed sequences and ensure keys are integers
            for k_str, seq_val in loaded_data.items():
                try:
                    k_int = int(k_str)
                except ValueError:
                    continue  # Skip non-integer keys in cache
                if k_int in target_ids_for_sequences: all_node_sequences[k_int] = seq_val
            print(
                f"✔ Loaded {len(all_node_sequences)} relevant sequences from cache in {time.time() - cache_load_start:.2f}s.",
                flush=True)
            loaded_from_cache = True
        except Exception as e:
            print(f"❌ Error loading/parsing cache {args.load_cache}: {e}. Will try GFA if provided.", file=sys.stderr,
                  flush=True)
            all_node_sequences = {}  # Reset on error

    if not loaded_from_cache or len(all_node_sequences) < len(target_ids_for_sequences):
        # If cache failed, or didn't contain all needed sequences, or no cache specified, try GFA
        if args.gfa:
            print(
                f"🔹 {'Attempting to load remaining' if loaded_from_cache else 'Loading ALL'} sequences from GFA: {args.gfa}...",
                flush=True)
            # Determine which sequences are still needed if some were loaded from cache
            still_needed_ids = target_ids_for_sequences - set(
                all_node_sequences.keys()) if loaded_from_cache else target_ids_for_sequences

            if still_needed_ids:
                gfa_sequences = load_node_sequences_from_gfa(args.gfa, still_needed_ids)
                all_node_sequences.update(gfa_sequences)  # Add newly loaded sequences
            else:
                print("ℹ️ All required sequences already found in cache.", flush=True)

            if args.save_cache and (
                    not loaded_from_cache or gfa_sequences):  # Save if built from GFA or updated from GFA
                print(f"🔹 Saving {len(all_node_sequences)} sequences to cache: {args.save_cache}...", flush=True)
                try:
                    with open(args.save_cache, 'w') as cf:
                        json.dump({str(k): v for k, v in all_node_sequences.items()}, cf,
                                  indent=2)  # Indent for readability
                    print(f"✔ Saved cache to {args.save_cache}", flush=True)
                except Exception as e:
                    print(f"❌ Error saving cache: {e}", file=sys.stderr, flush=True)
        elif not all_node_sequences:  # No GFA and cache failed/empty
            print("❌ No sequence source (valid cache or GFA) available for required nodes. Exiting.", file=sys.stderr,
                  flush=True)
            sys.exit(1)

    if len(all_node_sequences) < len(target_ids_for_sequences):
        missing_count = len(target_ids_for_sequences) - len(all_node_sequences)
        print(
            f"⚠️ Warning: Could not load sequences for {missing_count} nodes listed in index. Processing will continue for nodes with sequences.",
            file=sys.stderr, flush=True)
    if not all_node_sequences:
        print("❌ Failed to load any sequences for indexed nodes. Exiting.", file=sys.stderr, flush=True)
        sys.exit(1)

    # --- Main processing loop over chunks of nodes ---
    overall_results = {}
    num_total_nodes_to_process = len(all_indexed_nodes_info)  # Based on index

    # Filter all_indexed_nodes_info to only include nodes for which we have sequences
    process_node_info_list = [info for info in all_indexed_nodes_info if info[0] in all_node_sequences]
    actual_nodes_to_process_count = len(process_node_info_list)
    if actual_nodes_to_process_count < num_total_nodes_to_process:
        print(
            f"ℹ️ Will process {actual_nodes_to_process_count} nodes for which sequences were found (out of {num_total_nodes_to_process} in index).",
            flush=True)

    if not process_node_info_list:
        print("❌ No nodes to process after matching index with available sequences. Exiting.", flush=True)
        sys.exit(1)

    num_chunks = math.ceil(actual_nodes_to_process_count / args.processing_chunk_size)

    num_workers = min(args.workers, os.cpu_count() or 1,
                      actual_nodes_to_process_count)  # Don't use more workers than tasks
    if num_workers <= 0: num_workers = 1

    max_active_futures = args.max_active_futures if args.max_active_futures > 0 else num_workers * 2
    max_active_futures = max(1, max_active_futures)  # Ensure at least 1

    print(f"🔹 Will process {actual_nodes_to_process_count} nodes in {num_chunks} chunks using {num_workers} workers.",
          flush=True)
    print(f"   Max active futures to manage: {max_active_futures}", flush=True)

    progress_update_interval = max(1, int(args.processing_chunk_size / (
                num_workers * 2)))  # Heuristic, update a few times per chunk

    with ProcessPoolExecutor(max_workers=num_workers, initializer=init_worker, initargs=(args.dat,)) as executor:
        for i_chunk in range(num_chunks):
            chunk_start_node_idx = i_chunk * args.processing_chunk_size
            chunk_end_node_idx = min((i_chunk + 1) * args.processing_chunk_size, actual_nodes_to_process_count)
            current_node_info_chunk = process_node_info_list[chunk_start_node_idx:chunk_end_node_idx]

            if not current_node_info_chunk: continue  # Should not happen if num_chunks is correct

            print(
                f"\n--- Processing Chunk {i_chunk + 1}/{num_chunks} (Nodes {chunk_start_node_idx + 1}-{chunk_end_node_idx} of processable list) ---",
                flush=True)

            tasks_for_this_chunk = []
            for node_id, file_offset, n_records in current_node_info_chunk:
                sequence = all_node_sequences.get(node_id)  # Already filtered to have sequences
                if sequence:  # Should always be true now
                    tasks_for_this_chunk.append((node_id, file_offset, n_records, sequence, args.window))

            if not tasks_for_this_chunk:
                print(f"  No tasks in chunk {i_chunk + 1}. Skipping.", flush=True)  # Should not happen
                continue

            print(f"  Submitting {len(tasks_for_this_chunk)} tasks for chunk {i_chunk + 1}...", flush=True)

            future_to_node_id = {}
            task_iterator_chunk = iter(tasks_for_this_chunk)
            tasks_submitted_this_chunk = 0
            nodes_completed_this_chunk = 0

            # Submit initial batch for this chunk
            for _ in range(min(max_active_futures, len(tasks_for_this_chunk))):
                try:
                    task_args = next(task_iterator_chunk)
                    future = executor.submit(process_node_parallel, task_args)
                    future_to_node_id[future] = task_args[0]  # Map future to node_id
                    tasks_submitted_this_chunk += 1
                except StopIteration:
                    break

            if tasks_submitted_this_chunk > 0:
                print(f"  Initial submission for chunk {i_chunk + 1}: {tasks_submitted_this_chunk} tasks in flight.",
                      flush=True)

            while future_to_node_id:  # Loop as long as there are active futures for this chunk
                active_futures_list = list(future_to_node_id.keys())
                if not active_futures_list: break

                for future in as_completed(active_futures_list):
                    nodes_completed_this_chunk += 1
                    original_node_id = future_to_node_id.pop(future)  # Get node_id and remove from active
                    try:
                        _ret_node_id, pileup_dict = future.result()  # Get result or exception
                        if _ret_node_id != original_node_id:
                            print(
                                f"⚠️ Chunk {i_chunk + 1}: Node ID mismatch! Original: {original_node_id}, Returned: {_ret_node_id}",
                                file=sys.stderr, flush=True)
                        if pileup_dict:
                            if isinstance(pileup_dict, dict) and "error" in pileup_dict:
                                print(
                                    f"  ⚠️ Worker for node {original_node_id} (chunk {i_chunk + 1}) reported: {pileup_dict['error']}",
                                    file=sys.stderr, flush=True)
                            elif pileup_dict:  # Ensure it's not empty {}
                                overall_results[original_node_id] = pileup_dict
                    except Exception as exc:
                        print(f"\n❌ Task for node {original_node_id} (chunk {i_chunk + 1}) generated exception: {exc}",
                              file=sys.stderr, flush=True)
                        # import traceback; traceback.print_exc(file=sys.stderr) # For more detail

                    # Submit new task if available FOR THIS CHUNK
                    try:
                        if tasks_submitted_this_chunk < len(tasks_for_this_chunk):
                            new_task_args = next(task_iterator_chunk)
                            new_future = executor.submit(process_node_parallel, new_task_args)
                            future_to_node_id[new_future] = new_task_args[0]  # Add new future
                            tasks_submitted_this_chunk += 1
                    except StopIteration:
                        pass  # No more tasks in this chunk's iterator

                    if nodes_completed_this_chunk % progress_update_interval == 0 or \
                            nodes_completed_this_chunk == len(tasks_for_this_chunk) or \
                            not future_to_node_id:
                        print(
                            f"  Chunk {i_chunk + 1} Progress: {nodes_completed_this_chunk}/{len(tasks_for_this_chunk)} done. {len(future_to_node_id)} active.",
                            end='\r', flush=True)
                    break  # Re-evaluate active futures in the while loop
                else:  # Inner for loop completed without break
                    if not future_to_node_id: break

            print(
                f"\n  ✔ Chunk {i_chunk + 1} complete. Processed {nodes_completed_this_chunk} tasks. Total results collected so far: {len(overall_results)}",
                flush=True)
            del tasks_for_this_chunk  # Help GC

    # --- End of processing all chunks ---
    total_script_elapsed = time.time() - overall_start_time
    print(f"\n🏁 Finished all chunks in {total_script_elapsed:.2f} seconds.")
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