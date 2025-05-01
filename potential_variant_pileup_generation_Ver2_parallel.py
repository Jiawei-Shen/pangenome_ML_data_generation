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
# !!! IMPORTANT: Verify this struct matches your ACTUAL .dat record format !!!
# Example: <h (offset?) 150s (seq?) 150s (qual?) 20s (cigar?) h (mapq?) c (strand?)
# Adjust sizes (150, 150, 20) and types (h, s, c) precisely.
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# Global for worker process state
worker_dat_file = None

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions

def reverse_complement(sequence):
    """Computes the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def parse_idx_file(idx_path):
    """Parses the binary index file."""
    node_index = {}
    print(f"🔹 Parsing index file: {idx_path}")
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4:
                 print(f"❌ Error: Index file {idx_path} is too small.", file=sys.stderr)
                 sys.exit(1)

            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4:
                print(f"❌ Error: Could not read number of nodes header from {idx_path}.", file=sys.stderr)
                sys.exit(1)
            num_nodes = struct.unpack('<I', num_nodes_bytes)[0]
            print(f"🔹 Reading index header indicates {num_nodes:,} nodes...")


            # Verify expected file size based on header
            # Structure: node_id(I=4), offset(Q=8), block_size(I=4), n_records(I=4), padding(H=2) = 22 bytes/record
            INDEX_RECORD_SIZE = 22
            expected_min_size = 4 + num_nodes * INDEX_RECORD_SIZE
            if file_size < expected_min_size:
                print(f"⚠️ Warning: Index file {idx_path} size ({file_size:,} bytes) is smaller than expected ({expected_min_size:,} bytes) for {num_nodes:,} nodes. File might be truncated.", file=sys.stderr)
                # Continue cautiously, reading only available records

            parsed_count = 0
            while True:
                record_bytes = f.read(INDEX_RECORD_SIZE)
                if not record_bytes: # End of file
                    break
                if len(record_bytes) < INDEX_RECORD_SIZE:
                     print(f"❌ Error: Index file ended prematurely while reading record {parsed_count+1}/{num_nodes}. Read {len(record_bytes)} bytes, expected {INDEX_RECORD_SIZE}.", file=sys.stderr)
                     break # Stop processing further records

                try:
                     # --- Verify this unpack format matches your ACTUAL index file structure ---
                     node_id, offset, block_size, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                     node_index[node_id] = (offset, n_records)
                     parsed_count += 1
                except struct.error as se:
                     print(f"❌ Error: Failed unpacking index record {parsed_count+1}: {se}. Check index format.", file=sys.stderr)
                     break # Stop if index format is wrong

                if parsed_count % 1_000_000 == 0:
                    print(f"    Parsed {parsed_count:,} index records...")


        print(f"✔ Parsed {parsed_count:,} nodes from index file.")
        if parsed_count != num_nodes:
             print(f"⚠️ Warning: Header indicated {num_nodes:,} nodes, but parsed {parsed_count:,}. Index file might be corrupt or truncated.", file=sys.stderr)
        return node_index
    except FileNotFoundError:
        print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr) # Print traceback for unexpected errors
        sys.exit(1)


def load_node_sequences_from_gfa(gfa_path, target_node_ids):
    """Loads required node sequences from a GFA file."""
    node_sequences = {}
    target_node_set = set(target_node_ids)
    print(f"🔹 Reading GFA file: {gfa_path}")
    try:
        # Try to detect encoding, but default to utf-8 or latin-1 if issues arise
        try:
            f = open(gfa_path, 'r', encoding='utf-8')
            f.readline() # Try reading a line to detect encoding issues early
            f.seek(0)    # Reset cursor
        except UnicodeDecodeError:
            print("⚠️ Warning: UTF-8 decode failed for GFA, trying Latin-1...", file=sys.stderr)
            f = open(gfa_path, 'r', encoding='latin-1')

        with f:
            line_counter = 0
            parsed_count = 0
            found_count = 0
            start_gfa_time = time.time()

            for line in f:
                line_counter += 1
                if line_counter % 10_000_000 == 0:
                    elapsed = time.time() - start_gfa_time
                    rate = line_counter / elapsed if elapsed > 0 else 0
                    print(f"    Checked {line_counter:,} lines in GFA... ({found_count:,} target seqs found, {rate:.0f} lines/s)")

                # Optimization: Skip parsing if line doesn't start with 'S'
                if not line.startswith('S\t'):
                    continue

                parts = line.strip().split('\t')
                # Check for expected format S <node_id> <sequence> [*optional tags]
                if len(parts) < 3:
                    continue # Skip malformed S lines

                try:
                    # Assume node IDs are integers based on index format
                    nid = int(parts[1])
                except ValueError:
                    # print(f"⚠️ Warning: Skipping S line {line_counter} with non-integer ID: {parts[1]}", file=sys.stderr)
                    continue # Skip nodes with non-integer IDs

                if nid in target_node_set:
                    node_sequences[nid] = parts[2] # Store the sequence
                    found_count += 1
                    # Optimization: remove found ID from set to speed up future checks
                    target_node_set.remove(nid)
                    # Optimization: stop early if all target nodes are found
                    if not target_node_set:
                       print(f"\n✔ Found all {len(target_node_ids):,} target sequences after checking {line_counter:,} lines.")
                       break

            print(f"✔ Finished checking {line_counter:,} lines in GFA.")
            print(f"✔ Loaded sequences for {found_count:,} target nodes.")
            if target_node_set: # Check if any nodes were not found
                 missing_count = len(target_node_set)
                 print(f"⚠️ Warning: Could not find sequences for {missing_count:,} target nodes present in the index file.", file=sys.stderr)
                 # Optionally list missing nodes if count is small
                 # if missing_count < 20:
                 #    print(f"   Missing node IDs: {list(target_node_set)}")

    except FileNotFoundError:
        print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

    return node_sequences

def decode_cigar(cigar_string):
    """Parses a CIGAR string into operations and lengths."""
    if not cigar_string or cigar_string == '*':
        return []
    try:
        # Common CIGAR ops for alignment analysis: M, I, D, X, =
        # Others like S, H, N, P might indicate clipping or introns
        return re.findall(r'(\d+)([MIDX=])', cigar_string)
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr)
        return []


def detect_variants_from_cigar(offset, cigar_string):
    """Identifies simple variant events (Ins/Del/Mismatch) from CIGAR."""
    variants = []
    pos = offset # 0-based position on the reference node sequence
    cigar_ops = decode_cigar(cigar_string)

    for length_str, op in cigar_ops:
        try:
            length = int(length_str)
        except ValueError:
            # print(f"⚠️ Warning: Invalid length in CIGAR op '{length_str}{op}' in '{cigar_string}'", file=sys.stderr)
            continue # Skip this malformed operation

        if op == 'M' or op == '=': # Match or sequence match
            pos += length
        elif op == 'X': # Sequence mismatch
            # Report one 'X' event per base in the mismatch block
            for i in range(length):
                variants.append((pos + i, 'X'))
            pos += length
        elif op == 'I': # Insertion to the reference
            # Insertion occurs *after* the reference position `pos - 1`. Anchor to that position.
            variants.append((pos - 1, 'I'))
            # Insertions do not consume reference bases, so pos does not change here.
        elif op == 'D': # Deletion from the reference
            # Deletion starts *at* the reference position `pos`. Anchor to that position.
            variants.append((pos, 'D'))
            pos += length # Deletions consume reference bases
        # Other CIGAR ops (like S, H, N, P) are ignored for this variant detection

    return variants

# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function

def init_worker(dat_file_path):
    """Initializer for each worker process: opens the .dat file."""
    global worker_dat_file
    pid = os.getpid()
    # print(f"[Worker {pid}] Initializing...") # Verbose debug
    try:
        worker_dat_file = open(dat_file_path, 'rb')
        # print(f"[Worker {pid}] Opened {dat_file_path}") # Verbose debug
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
    Processes reads for a single node, finds variants, generates pileups.
    Returns (node_id, pileup_dict) on success, or (node_id, {}) on failure.
    """
    node_id, offset, n_records, sequence = task_args
    global worker_dat_file
    pid = os.getpid()
    # print(f"[Worker {pid}] Processing Node {node_id}...") # Verbose debug

    # --- Pre-check worker state ---
    if worker_dat_file is None:
         print(f"❌ Error [Worker {pid}][Node {node_id}]: Worker DAT file handle not initialized.", file=sys.stderr)
         return node_id, {} # Indicate failure
    if not sequence:
         # This check might be redundant if generate_tasks filters Nones, but good safeguard
         # print(f"⚠️ Warning [Worker {pid}][Node {node_id}]: Empty sequence provided. Skipping.", file=sys.stderr)
         return node_id, {} # Cannot process without sequence
    node_len = len(sequence)

    # --- Main processing block with enhanced error handling ---
    segments = []
    pileups = {} # Final result for this node
    try:
        # --- Seek and Read Records ---
        try:
            # Offset is the start of the block containing node info + records
            # Assume first 10 bytes are node_id, offset, block_size header (VERIFY THIS)
            record_data_offset = offset + 10
            worker_dat_file.seek(record_data_offset)
        except IOError as ioe:
            print(f"❌ Error [Worker {pid}][Node {node_id}]: Failed to seek to offset {record_data_offset}: {ioe}", file=sys.stderr)
            return node_id, {} # Fail task if seek fails

        for record_idx in range(n_records):
            # --- Read Record ---
            try:
                data = worker_dat_file.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE:
                    # Check if we simply reached the end of the *file* unexpectedly
                    current_pos = worker_dat_file.tell()
                    print(f"⚠️ Warning [Worker {pid}][Node {node_id}]: Short read ({len(data)} bytes, expected {RECORD_SIZE}) for record {record_idx+1}/{n_records}. Current file pos: {current_pos}. Stopping reads for node.", file=sys.stderr)
                    break # Stop processing this node if data seems truncated
            except IOError as read_ioe: # Catch error during read itself
                 print(f"❌ Error [Worker {pid}][Node {node_id}]: I/O error reading record {record_idx+1}: {read_ioe}. Stopping reads.", file=sys.stderr)
                 break # Stop processing this node

            # --- Unpack Record Data ---
            try:
                # !!! Use the EXACT structure defined globally and verified against your data !!!
                read_off, raw_seq, raw_bq, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)
            except struct.error as se:
                 # This often indicates RECORD_STRUCT or RECORD_SIZE is wrong, or file corruption
                 print(f"❌ Error [Worker {pid}][Node {node_id}]: Struct unpack error, record {record_idx+1}. Expected size {RECORD_SIZE}. Error: {se}. Check RECORD_STRUCT/SIZE or data corruption. Skipping record.", file=sys.stderr)
                 continue # Skip this problematic record

            # --- Decode and Filter ---
            try:
                # Basic filtering example (adjust as needed)
                if mapq < 10: # Example MAPQ filter
                    continue

                # Decode fields, handling potential non-ASCII bytes robustly
                # Use 'replace' or 'ignore' to avoid crashing on unexpected bytes
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')
                cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                strand_char = strand_byte.decode('ascii', errors='replace') # Should be '+' or '-'
            except Exception as decode_err: # Catch potential errors in decode/filter too
                print(f"⚠️ Warning [Worker {pid}][Node {node_id}]: Decode/Filter error, record {record_idx+1}: {decode_err}. Skipping.", file=sys.stderr)
                continue

            read_len = len(seq) # Get length after decoding

            # --- Adjust for Reverse Strand ---
            # Use the logic from your script (ensure it matches coordinate system)
            current_off = read_off # Use the offset unpacked from the record
            if strand_char == '-' and node_len > 0:
                 # This calculation assumes 'read_off' is the start position on the forward strand
                 # And the alignment covers 'read_len' bases on the forward strand.
                 # The start on the reverse complement would be node_len - (forward_start + forward_len)
                 adj_off = node_len - current_off - read_len
                 if adj_off < 0:
                    # print(f"⚠️ Warning [Worker {pid}][Node {node_id}]: Negative adjusted offset ({adj_off}), record {record_idx+1}. Orig off={current_off}, node_len={node_len}, read_len={read_len}. Skipping.", file=sys.stderr)
                    continue # Skip reads with inconsistent coordinates
                 seq = reverse_complement(seq)
                 current_off = adj_off # Use adjusted offset for subsequent logic

            # Append processed segment info for variant detection/pileup
            segments.append((current_off, seq, cigar))

            # --- Catch unexpected errors processing a single record ---
            # This except block is now less likely to be hit due to specific catches above,
            # but remains as a safeguard.
            # except Exception as record_err:
            #      print(f"❌ Error [Worker {pid}][Node {node_id}]: Unexpected error processing record {record_idx+1}: {record_err}", file=sys.stderr)
            #      traceback.print_exc(file=sys.stderr)
            #      continue # Skip problematic record

        # --- End of reading loop for this node ---

        # --- Variant Detection ---
        reads_by_variant = defaultdict(list)
        try:
            for read_offset, read_seq, read_cigar in segments:
                if read_offset < 0: continue # Skip segments with invalid offset

                # Detect variants based on CIGAR string relative to the node sequence
                variants = detect_variants_from_cigar(read_offset, read_cigar)
                for vpos, vtype in variants:
                    # Validate variant position against the node length
                    is_valid_pos = False
                    # Insertion *after* vpos, so vpos can be -1 up to node_len-1
                    if vtype == 'I' and -1 <= vpos < node_len:
                         is_valid_pos = True
                    # Deletion/Mismatch *at* vpos, so vpos must be 0 to node_len-1
                    elif vtype in ('D', 'X') and 0 <= vpos < node_len:
                         is_valid_pos = True

                    if is_valid_pos:
                        # Store the original read offset and sequence associated with this variant
                        reads_by_variant[(vpos, vtype)].append((read_offset, read_seq))
                    #else: # Debugging for invalid variant positions
                    #    print(f"Debug [Worker {pid}][Node {node_id}]: Invalid variant pos/type ({vpos}, {vtype}) vs node_len {node_len}. Read offset {read_offset}, CIGAR {read_cigar}", file=sys.stderr)


        except Exception as var_err:
             print(f"❌ Error [Worker {pid}][Node {node_id}]: Error during variant detection phase: {var_err}", file=sys.stderr)
             traceback.print_exc(file=sys.stderr)
             return node_id, {} # Fail the node if variant detection has issues

        # --- Generate Pileups ---
        try:
            window = 60 # Pileup window size (total width)
            half = window // 2 # Bases before/after variant position
            for (vpos, vtype), reads in reads_by_variant.items():
                if not reads: continue # Skip if no reads somehow associated

                # Create pileup matrix: rows=reads, cols=window, values=base index
                mat = np.full((len(reads), window), BASE_TO_INDEX['N'], dtype=np.uint8)

                for i, (read_offset, read_seq) in enumerate(reads):
                    read_len = len(read_seq)
                    # Calculate start index *in the read* corresponding to start of pileup window
                    # Pileup window on node: [vpos - half, vpos - half + window)
                    # Start of window on node = vpos - half
                    # Start of window in read = (start of window on node) - read_offset
                    start_in_read = vpos - read_offset - half

                    # Fill the matrix row for this read
                    for j in range(window): # j is the column index in the pileup matrix (0 to window-1)
                        # Corresponding base index in the read string
                        read_idx = start_in_read + j
                        if 0 <= read_idx < read_len:
                            base = read_seq[read_idx].upper()
                            # Get index from BASE_TO_INDEX, default to 'N' if base is unexpected
                            mat[i, j] = BASE_TO_INDEX.get(base, BASE_TO_INDEX['N'])
                        # else: Base is outside the read boundaries, matrix already initialized to 'N'

                # Store pileup matrix, converting to list of lists for JSON serialization
                # Use a unique key combining position and type
                pileups[f"{vpos}_{vtype}"] = mat.tolist()

        except Exception as pileup_err:
            print(f"❌ Error [Worker {pid}][Node {node_id}]: Error during pileup generation phase: {pileup_err}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            return node_id, {} # Fail the node if pileup generation has issues

    # --- Catch-all for unexpected errors in the main worker logic ---
    except Exception as e_main:
        print(f"❌❌❌ [Worker {pid}][Node {node_id}] UNCAUGHT MAJOR error during processing: {e_main}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return node_id, {} # Indicate failure

    # Verbose debug:
    # if pileups: print(f"[Worker {pid}][Node {node_id}] Processed successfully. Found {len(pileups)} variants with pileups.")
    # else: print(f"[Worker {pid}][Node {node_id}] Processed successfully. No variants found/stored.")
    return node_id, pileups # Return result on success


# ─────────────────────────────────────────────────────────────────────────────
# Task Generator Function

def generate_tasks(node_index, node_sequences):
    """
    Generator function to yield task tuples (node_id, offset, n_records, sequence)
    one by one, avoiding building a massive list in memory.
    """
    print("    -> Starting task generation...")
    generated_count = 0
    nodes_missing_sequence = 0
    skipped_zero_records = 0
    start_time = time.time()

    # Iterate through the index directly
    for node_id, (offset, n_records) in node_index.items():
        if n_records <= 0:
            skipped_zero_records += 1
            continue # Skip nodes with 0 records specified in index

        sequence = node_sequences.get(node_id) # Get sequence from pre-loaded dict
        if sequence is not None:
            # Yield the tuple needed by process_node_parallel
            yield (node_id, offset, n_records, sequence)
            generated_count += 1
            if generated_count % 500000 == 0: # Print progress every 500k tasks yielded
                 elapsed = time.time() - start_time
                 rate = generated_count / elapsed if elapsed > 0 else 0
                 print(f"      Generated {generated_count:,} tasks... ({rate:.0f} tasks/s)")
        else:
            nodes_missing_sequence += 1 # Count nodes missing from sequence cache/GFA

    # Print summary after generator is exhausted
    print(f"    -> Finished task generation.")
    print(f"    -> Total tasks yielded: {generated_count:,}")
    if nodes_missing_sequence > 0:
        print(f"    -> Nodes skipped (missing sequence): {nodes_missing_sequence:,}")
    if skipped_zero_records > 0:
        print(f"    -> Nodes skipped (zero records in index): {skipped_zero_records:,}")
    # Note: This function doesn't easily return the counts, they are just printed.

# ─────────────────────────────────────────────────────────────────────────────
# Main Execution Logic

def main():
    parser = argparse.ArgumentParser(
        description="Parallel variant-centered pileup generation (Robust Generator Version).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Show default values in help
    )
    # Input/Output Arguments
    parser.add_argument("dat", help=".dat file path (binary read alignment data)")
    parser.add_argument("idx", help=".idx file path (binary index for .dat file)")
    parser.add_argument("output", help="JSON output file path for pileups")
    # Sequence Source Arguments (choose one method)
    parser.add_argument("--gfa", help="GFA graph file path (used to build node sequence cache if cache not loaded)")
    parser.add_argument("--load-cache", help="Load node sequences from this JSON cache file (skip GFA parsing)")
    parser.add_argument("--save-cache", help="Save node sequences to this JSON cache file (required if using --gfa without --load-cache)")
    # Performance Arguments
    parser.add_argument("-w", "--workers", type=int, default=8,
                        help="Number of worker processes. Tune based on I/O tests and CPU usage!")
    # chunksize argument removed as it's not used with submit/as_completed

    # Print help if no arguments are given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    # --- Start Timer ---
    script_start_time = time.time()
    print(f"--- Starting execution at {time.strftime('%Y-%m-%d %H:%M:%S')} ---")

    # --- Input Validation ---
    if not os.path.isfile(args.dat): print(f"❌ Error: DAT file not found: {args.dat}", file=sys.stderr); sys.exit(1)
    if not os.path.isfile(args.idx): print(f"❌ Error: Index file not found: {args.idx}", file=sys.stderr); sys.exit(1)
    if not args.load_cache and not args.gfa: print("❌ Error: Must provide sequence source: --gfa (to build cache) or --load-cache.", file=sys.stderr); sys.exit(1)
    if args.gfa and not args.save_cache and not args.load_cache: print("❌ Error: If using --gfa without --load-cache, must specify --save-cache.", file=sys.stderr); sys.exit(1)
    if args.load_cache and not os.path.isfile(args.load_cache): print(f"❌ Error: Cache file to load not found: {args.load_cache}", file=sys.stderr); sys.exit(1)
    if args.gfa and not os.path.isfile(args.gfa): print(f"❌ Error: GFA file not found: {args.gfa}", file=sys.stderr); sys.exit(1)

    # --- Load Index ---
    start_time = time.time()
    node_index = parse_idx_file(args.idx)
    if not node_index: print("❌ Error: Failed to parse node index.", file=sys.stderr); sys.exit(1)
    print(f"✔ Index parsing took {time.time() - start_time:.2f} seconds.")

    # --- Load or Build Node Sequences ---
    node_sequences = {}
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"🔹 Loading node sequences from cache: {args.load_cache}...")
        start_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                loaded_data = json.load(cf)
                # Convert JSON string keys back to integer node IDs
                node_sequences = {int(k): v for k, v in loaded_data.items()}
            print(f"✔ Loaded {len(node_sequences):,} sequences from cache in {time.time() - start_time:.2f} seconds.")
        except json.JSONDecodeError as jde:
             print(f"❌ Error decoding JSON from cache file {args.load_cache}: {jde}", file=sys.stderr); sys.exit(1)
        except Exception as e:
             print(f"❌ Error loading cache file {args.load_cache}: {e}", file=sys.stderr); sys.exit(1)
    elif args.gfa:
        print(f"🔹 Building node sequence cache from GFA: {args.gfa}...")
        start_time = time.time()
        # Load only sequences for nodes present in the index
        node_sequences = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        print(f"✔ Sequence loading from GFA took {time.time() - start_time:.2f} seconds.")
        if not node_sequences:
             print("⚠️ Warning: No sequences loaded from GFA for nodes in index. Output will be empty.", file=sys.stderr)
             # Optionally exit if this is considered a fatal error
             # sys.exit(1)

        if args.save_cache:
            print(f"🔹 Saving {len(node_sequences):,} node sequences to cache: {args.save_cache}...")
            start_time = time.time()
            try:
                # Ensure directory exists for cache file
                cache_dir = os.path.dirname(args.save_cache)
                if cache_dir and not os.path.exists(cache_dir):
                     os.makedirs(cache_dir, exist_ok=True)
                     print(f"    Created directory: {cache_dir}")

                with open(args.save_cache, 'w') as cf:
                    # Save keys as strings in JSON, as required by the format
                    json.dump({str(k): v for k, v in node_sequences.items()}, cf)
                print(f"✔ Saved cache in {time.time() - start_time:.2f} seconds.")
            except IOError as ioe:
                 print(f"❌ Error saving cache file {args.save_cache}: {ioe}", file=sys.stderr) # Continue even if save fails
            except Exception as e:
                 print(f"❌ Error during cache saving to {args.save_cache}: {e}", file=sys.stderr)

    # --- Prepare and Execute Tasks using Generator ---
    print("🔹 Creating task generator...")
    # Create the generator object (doesn't build the list in memory)
    task_generator = generate_tasks(node_index, node_sequences)

    # Estimate total tasks for initial message (actual count after submission)
    estimated_total_tasks = len(node_index) # Approximation based on index size
    num_workers = min(args.workers, os.cpu_count() or 1) # Use available CPUs if less than requested
    if args.workers > os.cpu_count():
         print(f"⚠️ Warning: Requested workers ({args.workers}) exceeds available CPUs ({os.cpu_count()}). Using {os.cpu_count()} workers.", file=sys.stderr)
         num_workers = os.cpu_count()

    print(f"🔹 Processing estimated ~{estimated_total_tasks:,} nodes using {num_workers} workers...")

    results = {}
    start_proc_time = time.time()
    processed_count = 0
    errors_reported = 0
    submitted_tasks = 0 # Keep track of actual submitted tasks

    try:
        with ProcessPoolExecutor(max_workers=num_workers,
                                 initializer=init_worker,
                                 initargs=(args.dat,)) as executor:

            print("--- Submitting tasks from generator...")
            # Submit tasks *from the generator* and store future mapped to node_id
            # This loop pulls tasks one by one and submits them
            futures_map = {
                executor.submit(process_node_parallel, task_data): task_data[0]
                for task_data in task_generator
            }
            submitted_tasks = len(futures_map) # Get the actual count of tasks submitted
            print(f"--- Submitted {submitted_tasks:,} actual tasks. Waiting for results...")

            # Check if any tasks were actually submitted
            if submitted_tasks == 0:
                 print("✅ Warning: No tasks were generated or submitted. Check index and sequence data.", file=sys.stderr)
                 # No need to proceed further if there's no work
                 sys.exit(0)

            # Process results as they complete
            for future in as_completed(futures_map):
                node_id = futures_map[future] # Get node_id associated with this future
                try:
                    # Get result - Handles exceptions caught inside worker via return value
                    # Re-raises exceptions from worker startup or pickling issues
                    _returned_node_id, pileup_dict = future.result()

                    # Check if worker returned a valid (non-empty) dict
                    if pileup_dict:
                         results[node_id] = pileup_dict
                    else:
                         # Worker returned {} indicating a caught error, already logged by worker
                         errors_reported += 1

                except Exception as exc:
                    # Handle errors during future retrieval itself (e.g., worker crashed hard)
                    print(f'❌❌❌ Main thread error getting result for Node ID {node_id}: {exc}', file=sys.stderr)
                    # Include traceback for these unexpected errors
                    traceback.print_exc(file=sys.stderr)
                    errors_reported += 1

                # --- Progress Update ---
                processed_count += 1
                # Print progress periodically based on actual submitted tasks
                if processed_count % 10000 == 0 or processed_count == submitted_tasks:
                     elapsed = time.time() - start_proc_time
                     rate = processed_count / elapsed if elapsed > 0 else 0
                     error_msg = f"({errors_reported} errors reported)" if errors_reported > 0 else ""
                     # Use submitted_tasks as the denominator for progress
                     print(f"✔ {processed_count}/{submitted_tasks} tasks processed {error_msg}. Rate: {rate:.1f} tasks/s. Elapsed: {elapsed:.2f}s")

    except Exception as pool_exc:
         print(f"\n❌ An error occurred within the ProcessPoolExecutor context: {pool_exc}", file=sys.stderr)
         traceback.print_exc(file=sys.stderr)
         print("⚠️ Attempting to write any partial results obtained...", file=sys.stderr)

    # --- Final Summary ---
    total_proc_time = time.time() - start_proc_time
    print(f"\n✔ Parallel processing finished in {total_proc_time:.2f} seconds.")
    print(f"Processed {processed_count:,} tasks. Stored pileup results for {len(results):,} nodes.")
    if errors_reported > 0:
         print(f"⚠️ Encountered {errors_reported:,} errors during processing (check logs above for details).", file=sys.stderr)
    if processed_count < submitted_tasks:
         print(f"⚠️ Warning: Only {processed_count:,} tasks completed out of {submitted_tasks:,} submitted. Some tasks may have failed to run or complete.", file=sys.stderr)


    # --- Write Output ---
    if results: # Only write if there are results
        print(f"🔹 Writing {len(results):,} results to JSON output: {args.output}")
        start_write_time = time.time()
        try:
            # Ensure output directory exists
            output_dir = os.path.dirname(args.output)
            if output_dir and not os.path.exists(output_dir):
                 os.makedirs(output_dir, exist_ok=True)
                 print(f"    Created output directory: {output_dir}")

            with open(args.output, 'w') as out_f:
                # Save node IDs as strings for JSON compatibility
                json.dump({str(k): v for k, v in results.items()}, out_f, indent=2) # Use indent for readability
            write_elapsed_time = time.time() - start_write_time
            print(f"✔ Output written in {write_elapsed_time:.2f} seconds.")
        except IOError as ioe:
             print(f"❌ Error writing output JSON to {args.output}: {ioe}", file=sys.stderr)
        except Exception as e:
             print(f"❌ Unexpected error writing output JSON: {e}", file=sys.stderr)
             traceback.print_exc(file=sys.stderr)
    else:
         print("ℹ️ No results generated to write to output file.")


    # --- End Timer ---
    total_script_time = time.time() - script_start_time
    print(f"--- Total execution time: {total_script_time:.2f} seconds ---")
    print(f"--- Finished at {time.strftime('%Y-%m-%d %H:%M:%S')} ---")
    print(f"✅ Done. Output saved to {args.output}")


if __name__ == '__main__':
    # Good practice: protect the main execution block
    # Required for multiprocessing on some platforms (like Windows)
    main()