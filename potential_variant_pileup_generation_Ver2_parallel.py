#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np # Using numpy for efficient pileup matrix creation
from collections import defaultdict
import re # Import re at the top level
from concurrent.futures import ProcessPoolExecutor
# from os import cpu_count # Can uncomment to default workers

# ─────────────────────────────────────────────────────────────────────────────
# Constants
RECORD_STRUCT = struct.Struct("<h150s150s20shc") # From original code
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4,
                 'a': 0, 'c': 1, 'g': 2, 't': 3, 'n': 4} # Expanded slightly
N_INDEX = BASE_TO_INDEX['N'] # Convenience

# Default window size
DEFAULT_WINDOW_SIZE = 100

# Global for worker process state (file handle)
worker_dat_file = None

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions

def reverse_complement(sequence):
    """Computes the reverse complement of a DNA sequence."""
    if not isinstance(sequence, str): # Basic check
        return ""
    try:
        complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return sequence.translate(complement_map)[::-1]
    except Exception:
        # Log error or return empty? Returning empty for now.
        return ""

def parse_cigar(cigar_string):
    """Parses a CIGAR string into a list of (length, operation) tuples."""
    if not cigar_string or cigar_string == '*':
        return []
    try:
        # Include common CIGAR ops: M=match/mismatch, X=mismatch, = =match, I=insertion, D=deletion, S=soft clip, N=skip
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr)
        return []

# --- Re-implemented Variant Detection ---
def detect_variants_core(ref_seq, read_seq, start_offset, cigar_ops):
    """
    Core CIGAR-aware variant detection logic. Yields specific variants.

    Yields:
        tuple: (ref_pos, variant_type, ref_base, alt_base)
    """
    ref_len = len(ref_seq)
    read_len = len(read_seq)
    ref_ptr = start_offset
    read_ptr = 0
    last_valid_ref_ptr = start_offset - 1 # Position *before* which an insertion occurs

    try:
        for length, op in cigar_ops:
            if op in ('M', '=', 'X'): # Match, Alignment match, Alignment mismatch
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    current_read_pos = read_ptr + i

                    if current_ref_pos >= ref_len or current_read_pos >= read_len:
                        # print(f"Debug: Bounds exceeded M/= /X Ref:{current_ref_pos}/{ref_len}, Read:{current_read_pos}/{read_len}")
                        return # Alignment goes off sequence ends

                    ref_base = ref_seq[current_ref_pos].upper()
                    read_base = read_seq[current_read_pos].upper()
                    last_valid_ref_ptr = current_ref_pos

                    # Report mismatch ('X') if bases differ and are not 'N'
                    if op != '=' and ref_base != 'N' and read_base != 'N' and ref_base != read_base:
                        yield (current_ref_pos, 'X', ref_base, read_base)

                ref_ptr += length
                read_ptr += length

            elif op == 'I': # Insertion to the reference
                if read_ptr + length > read_len:
                    # print(f"Debug: Bounds exceeded I Read:{read_ptr + length}/{read_len}")
                    return # Alignment goes off sequence ends
                inserted_bases = read_seq[read_ptr : read_ptr + length].upper()
                # Insertion occurs *after* the last valid reference base position
                yield (last_valid_ref_ptr, 'I', '-', inserted_bases)
                read_ptr += length
                # Does not advance ref_ptr or change last_valid_ref_ptr meaningfully here

            elif op == 'D': # Deletion from the reference
                start_del_pos = ref_ptr
                if ref_ptr + length > ref_len:
                    # print(f"Debug: Bounds exceeded D Ref:{ref_ptr + length}/{ref_len}")
                    return # Alignment goes off sequence ends
                deleted_bases = ref_seq[ref_ptr : ref_ptr + length].upper()
                yield (start_del_pos, 'D', deleted_bases, '-')
                ref_ptr += length
                last_valid_ref_ptr = ref_ptr - 1 # Update last valid position to end of deletion
                # Does not advance read_ptr

            elif op == 'S': # Soft clipping (consumes read)
                read_ptr += length

            elif op == 'N': # Skipped region from reference (consumes reference)
                 ref_ptr += length
                 last_valid_ref_ptr = ref_ptr - 1

            elif op in ('H', 'P'): # Hard clipping, Padding (consumes neither)
                 pass # Ignored

    except IndexError:
        print(f"⚠️ Warning: IndexError during variant detection. RefLen={ref_len}, ReadLen={read_len}, Offset={start_offset}, CIGAR='{''.join(map(str,cigar_ops))}', RefPtr={ref_ptr}, ReadPtr={read_ptr}", file=sys.stderr)
    except Exception as e:
        print(f"❌ Error: Unexpected error in detect_variants_core: {e}", file=sys.stderr)

# --- Re-implemented Pileup Row Creation ---
def create_pileup_row(read_seq, read_offset, cigar_ops, window_start, window_end):
    """
    Creates a single pileup row (list of base indices) based on CIGAR alignment
    within a specified reference window.
    """
    window_len = window_end - window_start
    if window_len <= 0: return None # Should not happen if window calculated correctly

    # Initialize pileup row with 'N' index
    pileup_row = [N_INDEX] * window_len
    read_len = len(read_seq)

    ref_ptr = read_offset
    read_ptr = 0
    row_filled = False # Track if any base was actually placed in the window

    try:
        for length, op in cigar_ops:
            # Optimization: If current ref pointer is already past the window end,
            # and we haven't filled anything, the read doesn't overlap the start.
            if ref_ptr >= window_end and not row_filled:
                 return None

            if op in ('M', '=', 'X'): # Consumes both Ref and Read
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    current_read_pos = read_ptr + i

                    if current_read_pos >= read_len: # Check read bounds first
                         # print(f"Debug: Read index out of bounds pileup M/=/X: {current_read_pos}/{read_len}")
                         return pileup_row if row_filled else None # Return what we have

                    # Check if current reference position falls within the window
                    if window_start <= current_ref_pos < window_end:
                        pileup_idx = current_ref_pos - window_start
                        base = read_seq[current_read_pos]
                        pileup_row[pileup_idx] = BASE_TO_INDEX.get(base.upper(), N_INDEX)
                        row_filled = True

                ref_ptr += length
                read_ptr += length

            elif op == 'I': # Insertion (consumes Read only)
                if read_ptr + length > read_len: # Check read bounds
                    # print(f"Debug: Read index out of bounds pileup I: {read_ptr + length}/{read_len}")
                    return pileup_row if row_filled else None
                read_ptr += length
                # Insertions don't occupy reference space in the pileup row

            elif op == 'D': # Deletion (consumes Ref only)
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    if window_start <= current_ref_pos < window_end:
                        pileup_idx = current_ref_pos - window_start
                        # Deletion means no base from *this read* aligns here, keep 'N'.
                        pileup_row[pileup_idx] = N_INDEX
                        row_filled = True
                ref_ptr += length

            elif op == 'S': # Soft clipping (consumes Read only)
                read_ptr += length

            elif op == 'N': # Skipped region (consumes Ref only)
                 for i in range(length):
                    current_ref_pos = ref_ptr + i
                    if window_start <= current_ref_pos < window_end:
                         pileup_idx = current_ref_pos - window_start
                         # Skipped region also means no read base here.
                         pileup_row[pileup_idx] = N_INDEX
                         row_filled = True
                 ref_ptr += length

            elif op in ('H', 'P'): # Consumes nothing
                 pass

        # Return the row only if it overlapped the window meaningfully
        return pileup_row if row_filled else None

    except IndexError:
         print(f"⚠️ Warning: IndexError during pileup creation. ReadLen={read_len}, Offset={read_offset}, CIGAR='{''.join(map(str,cigar_ops))}', RefPtr={ref_ptr}, ReadPtr={read_ptr}, Win=[{window_start}:{window_end}]", file=sys.stderr)
         return pileup_row if row_filled else None # Return partial row if possible
    except Exception as e:
         print(f"❌ Error: Unexpected error in create_pileup_row: {e}", file=sys.stderr)
         return None


# -- File parsing functions from original code (kept as is) --
def parse_idx_file(idx_path):
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

             expected_min_size = 4 + num_nodes * 22 # NodeID(I)+Offset(Q)+BlockSize(I)+NRecs(I)+Padding(H) = 4+8+4+4+2 = 22
             if file_size < expected_min_size:
                 print(f"⚠️ Warning: Index file {idx_path} may be truncated. Expected at least {expected_min_size} bytes for {num_nodes} nodes, found {file_size}.", file=sys.stderr)

             print(f"🔹 Reading index for {num_nodes} nodes...")
             nodes_parsed_count = 0
             for i in range(num_nodes):
                 record_bytes = f.read(22)
                 if len(record_bytes) < 22:
                      print(f"❌ Error: Index file ended prematurely while reading record {i+1}/{num_nodes}.", file=sys.stderr)
                      break
                 # < NodeID(uint32), Offset(uint64), BlockSize(uint32), NumRecords(uint32), Padding(uint16) >
                 node_id, offset, block_size, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                 if n_records > 0: # Only store nodes that have records
                    node_index[node_id] = (offset, n_records)
                    nodes_parsed_count += 1
                 # else: print(f"Debug: Node {node_id} has 0 records in index, skipping.")


         print(f"✔ Parsed {nodes_parsed_count} nodes with records from index file.")
         if nodes_parsed_count == 0 and num_nodes > 0:
              print(f"⚠️ Warning: Header indicated {num_nodes} nodes, but 0 nodes with records were found in the index data.", file=sys.stderr)
         elif len(node_index) != num_nodes and nodes_parsed_count > 0 : # Check if parsed count != header if we found some records
              print(f"⚠️ Warning: Header indicated {num_nodes} nodes, but parsed {nodes_parsed_count} nodes with records.", file=sys.stderr)
         return node_index
     except FileNotFoundError:
         print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr)
         sys.exit(1)
     except Exception as e:
         print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
         sys.exit(1)


def load_node_sequences_from_gfa(gfa_path, target_node_ids):
     node_sequences = {}
     target_node_set = set(target_node_ids) # Use a set for faster lookups
     nodes_found_count = 0
     try:
         with open(gfa_path, 'r') as f:
             print(f"🔹 Reading GFA file: {gfa_path}")
             line_counter = 0
             for line in f:
                 line_counter += 1
                 if line_counter % 5_000_000 == 0: # Progress update
                     print(f"  Checked {line_counter:,} lines... found {nodes_found_count}/{len(target_node_ids)} target sequences.")

                 if not line.startswith('S\t'):
                     continue

                 parts = line.strip().split('\t')
                 if len(parts) < 3:
                     continue # Skip malformed lines

                 try:
                     nid = int(parts[1])
                 except ValueError:
                     # print(f"⚠️ Warning: Skipping S line {line_counter} with non-integer ID: {parts[1]}", file=sys.stderr)
                     continue # Skip nodes with non-integer IDs

                 if nid in target_node_set:
                     seq = parts[2].upper() # Store uppercase
                     if seq != '*' and re.match(r'^[ACGTN]+$', seq): # Basic validation
                         node_sequences[nid] = seq
                         nodes_found_count += 1
                         target_node_set.remove(nid) # Remove from set once found
                         if not target_node_set: # Stop early if all found
                            print(f"✔ Found all {len(target_node_ids)} target sequences after {line_counter:,} lines.")
                            break
                     # else: print(f"Debug: Skipping node {nid} due to invalid sequence '{seq[:10]}...'")


             print(f"✔ Finished checking {line_counter:,} lines in GFA.")
             print(f"✔ Loaded sequences for {len(node_sequences)} target nodes.")
             if target_node_set: # Check if any were not found
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
    # print(f"[Worker {os.getpid()}] Initializing and opening {dat_file_path}")
    try:
        worker_dat_file = open(dat_file_path, 'rb')
    except FileNotFoundError:
         print(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path}", file=sys.stderr)
         sys.exit(1) # Exit worker if file missing
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path}: {e}", file=sys.stderr)
        sys.exit(1) # Exit worker on other errors


def process_node_parallel(task_args):
    """
    Function executed by each worker process.
    Processes a single node to find variants and generate pileups.
    """
    node_id, file_read_offset, n_records, node_sequence, window_size = task_args # Added window_size
    global worker_dat_file

    if worker_dat_file is None:
         print(f"❌ Error [Worker {os.getpid()}]: Worker DAT file handle not initialized for node {node_id}.", file=sys.stderr)
         return node_id, {}

    if not node_sequence:
         # This check should ideally happen before adding task, but good safeguard
         return node_id, {}

    node_len = len(node_sequence)
    half_window = window_size // 2

    # Store reads relevant to potential variants: {variant_key: [(offset, oriented_seq, cigar_ops), ...]}
    reads_by_variant = defaultdict(list)
    processed_read_count = 0
    mapq_filter = 10 # Example filter value

    try:
        # Seek to the correct position for this node's records
        # Assuming the offset from index points directly to the first record (not block metadata)
        worker_dat_file.seek(file_read_offset)

        for record_idx in range(n_records):
            try:
                data = worker_dat_file.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE:
                    print(f"⚠️ Warning [Worker {os.getpid()}]: Short read ({len(data)} bytes) for node {node_id}, record {record_idx+1}/{n_records}. Stopping reads.", file=sys.stderr)
                    break

                # Unpack: offset_in_dat, seq, bq, cigar, mapq, strand
                read_offset_on_node, raw_seq, _, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)

                if mapq < mapq_filter:
                    continue

                try:
                    seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
                    cigar_str = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
                    strand_char = strand_byte.decode('ascii')
                except UnicodeDecodeError:
                    continue # Skip record if decoding fails

                if not seq or not cigar_str or cigar_str == '*':
                    continue

                processed_read_count += 1
                oriented_seq = seq
                current_offset = read_offset_on_node # Use original offset by default

                # Handle reverse strand
                if strand_char == '-':
                    oriented_seq = reverse_complement(seq)
                    if not oriented_seq: continue # Skip if RC failed

                    # Reverse CIGAR operations for reverse strand alignment relative to forward reference
                    # This part is tricky and assumes the stored offset/CIGAR relate to the forward strand.
                    # If the stored CIGAR already accounts for reverse strand, we might not need to reverse it.
                    # Assuming CIGAR needs reversal for consistency with coordinate adjustment:
                    # parsed_cigar_ops = parse_cigar(cigar_str)[::-1] # Reverse order of ops

                    # Let's *assume* the stored CIGAR is relative to the forward reference strand
                    # and the offset is also relative to the forward strand start.
                    # We already reverse complemented the sequence. The CIGAR ops themselves
                    # describe alignment to the forward strand. We don't need to reverse the CIGAR ops order.
                    # The offset adjustment logic seems specific to the data generation process.
                    # Let's use the simpler approach: keep CIGAR as is, use original offset,
                    # but use the reverse-complemented sequence. Variant/pileup functions
                    # will work relative to the forward node sequence.
                    # Let's comment out the complex offset adjustment from original code for now.
                    # read_len = len(seq)
                    # adj_off = node_len - read_offset_on_node - read_len
                    # if adj_off < 0: continue # Skip inconsistent read
                    # current_offset = adj_off

                # Parse CIGAR *once* per read
                cigar_ops = parse_cigar(cigar_str)
                if not cigar_ops: continue # Skip reads with unparseable CIGAR

                # --- Detect Variants for THIS read ---
                variants_found_in_read = list(detect_variants_core(node_sequence, oriented_seq, current_offset, cigar_ops))

                # --- Group read by variants it supports ---
                for vpos, vtype, ref_base, alt_base in variants_found_in_read:
                    # Basic check: ensure variant position is plausible relative to node length
                    # More precise check inside detect_variants_core should handle most cases
                    if (vtype == 'I' and -1 <= vpos < node_len) or \
                       (vtype in ('D', 'X') and 0 <= vpos < node_len):
                        variant_key = f"{vpos}_{vtype}_{ref_base}_{alt_base}"
                        # Store original offset, oriented sequence, and *parsed* CIGAR ops
                        reads_by_variant[variant_key].append((current_offset, oriented_seq, cigar_ops))
                    #else: print(f"Debug: Var {vpos}_{vtype} skipped, pos issue rel node_len {node_len}")


            except struct.error as se:
                 print(f"❌ Error [Worker {os.getpid()}]: Failed to unpack record {record_idx+1} for node {node_id}: {se}. Stopping reads.", file=sys.stderr)
                 break
            except Exception as e_inner:
                 print(f"❌ Error [Worker {os.getpid()}]: Unexpected error processing record {record_idx+1} for node {node_id}: {e_inner}", file=sys.stderr)
                 continue # Continue to next record

    except IOError as ioe:
         print(f"❌ Error [Worker {os.getpid()}]: I/O error reading data for node {node_id} at offset {file_read_offset}: {ioe}", file=sys.stderr)
         return node_id, {} # Return empty result
    except Exception as e_outer:
        print(f"❌ Error [Worker {os.getpid()}]: Unexpected error seeking/reading node {node_id}: {e_outer}", file=sys.stderr)
        return node_id, {}

    # --- Generate Pileups for detected variants ---
    final_pileups = {}
    if not reads_by_variant:
        # print(f"Debug: Node {node_id} - No variants found or no reads passed filters.")
        return node_id, {}

    # print(f"Debug: Node {node_id} - Processed {processed_read_count} reads. Found {len(reads_by_variant)} unique variant events.")

    for variant_key, supporting_reads in reads_by_variant.items():
        if not supporting_reads: continue

        # Extract variant position from key (needed for window calculation)
        try:
            vpos = int(variant_key.split('_')[0])
        except (ValueError, IndexError):
            print(f"⚠️ Warning [Worker {os.getpid()}]: Could not parse position from variant key '{variant_key}' for node {node_id}. Skipping.", file=sys.stderr)
            continue

        # --- Calculate Window Boundaries (Handles short nodes & edges) ---
        if node_len <= window_size:
            # Node is shorter than or equal to window, use the whole node
            window_start = 0
            window_end = node_len
            actual_window_len = node_len
        else:
            # Node is longer, calculate window around vpos, clamp to node boundaries
            # For insertions (key looks like 'pos_I_-_ALT'), vpos is base *before* insertion. Center there.
            # For mismatches/deletions, vpos is the first affected base. Center there.
            center_pos = vpos
            window_start = max(0, center_pos - half_window)
            # Calculate end, ensuring it doesn't exceed node length
            window_end = min(node_len, window_start + window_size)
            # Recalculate start in case window end clamping shortened the effective window near the end
            window_start = max(0, window_end - window_size)
            actual_window_len = window_end - window_start

        if actual_window_len <= 0:
            # print(f"Debug: Node {node_id} VarKey {variant_key}: Zero width window [{window_start}:{window_end}] calculated.")
            continue # Skip if window is invalid

        # --- Create Pileup Matrix for this variant ---
        pileup_rows = []
        for read_offset, read_seq, cigar_ops in supporting_reads:
            row = create_pileup_row(read_seq, read_offset, cigar_ops, window_start, window_end)
            if row: # Only append if the read actually overlapped the window
                pileup_rows.append(row)

        # Only store if we got any valid pileup rows
        if pileup_rows:
            # Convert list of lists to numpy array if desired for analysis later,
            # but keep as list of lists for JSON output.
            # mat = np.array(pileup_rows, dtype=np.uint8)
            # final_pileups[variant_key] = mat.tolist()
            final_pileups[variant_key] = pileup_rows # Store as list of lists
        #else: print(f"Debug: Node {node_id} VarKey {variant_key}: No supporting pileup rows generated for window [{window_start}:{window_end}].")


    # print(f"[Worker {os.getpid()}] Finished node {node_id}. Generated pileups for {len(final_pileups)} variants.")
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
    parser.add_argument("--gfa", required=True, help="GFA graph file path (required for node sequences)")
    # Removed cache logic for simplicity, always load from GFA for now
    # parser.add_argument("--load-cache", help="Load node sequences from this JSON cache file")
    # parser.add_argument("--save-cache", help="Save node sequences to this JSON cache file (used if --gfa is provided)")
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW_SIZE, help="Pileup window size around variant")
    parser.add_argument("-w", "--workers", type=int, default=os.cpu_count(), help="Number of worker processes to use")
    parser.add_argument("-c", "--chunksize", type=int, default=100, help="Approximate number of nodes per worker task")
    args = parser.parse_args()

    # --- Input Validation ---
    if not os.path.isfile(args.dat):
         print(f"❌ Error: DAT file not found: {args.dat}", file=sys.stderr)
         sys.exit(1)
    if not os.path.isfile(args.idx):
         print(f"❌ Error: Index file not found: {args.idx}", file=sys.stderr)
         sys.exit(1)
    if not os.path.isfile(args.gfa):
         print(f"❌ Error: GFA file not found: {args.gfa}", file=sys.stderr)
         sys.exit(1)
    if args.window <= 0 or args.window % 2 != 0:
         print(f"❌ Error: Window size (--window={args.window}) must be a positive even number for symmetric half-window calculation.", file=sys.stderr)
         # Or adjust logic for odd windows if preferred, but even is simpler here.
         # Let's allow odd windows and adjust half_window calc slightly if needed later.
         # Reverting to allow odd:
         if args.window <= 0:
            print(f"❌ Error: Window size (--window={args.window}) must be positive.", file=sys.stderr)
            sys.exit(1)


    # --- Load Index ---
    print("🔹 Parsing index file...")
    start_time = time.time()
    node_index = parse_idx_file(args.idx)
    if not node_index:
         print("❌ Error: Failed to parse node index. Exiting.", file=sys.stderr)
         sys.exit(1)
    print(f"✔ Index parsing took {time.time() - start_time:.2f} seconds.")

    # --- Load Node Sequences (Simplified: Always from GFA) ---
    print(f"🔹 Loading node sequences from GFA: {args.gfa}...")
    start_time = time.time()
    # We only need sequences for nodes present in the index file
    target_ids = node_index.keys()
    node_sequences = load_node_sequences_from_gfa(args.gfa, target_ids)
    print(f"✔ Sequence loading from GFA took {time.time() - start_time:.2f} seconds.")

    if not node_sequences:
         print("⚠️ Warning: No node sequences were loaded. Output will be empty.", file=sys.stderr)
         # Exit if no sequences are loaded for the indexed nodes? Or proceed?
         # Let's exit, as no work can be done.
         print("❌ Exiting as no sequences were available for indexed nodes.")
         sys.exit(1)


    # --- Prepare Tasks for Parallel Processing ---
    print("🔹 Preparing tasks for parallel processing...")
    tasks = []
    nodes_missing_sequence = 0
    skipped_nodes_no_records = 0 # Track nodes skipped because n_records=0 in index parsing

    for node_id, (file_read_offset, n_records) in node_index.items():
        if n_records <= 0: # Should be filtered by parse_idx_file now
            skipped_nodes_no_records += 1
            continue

        sequence = node_sequences.get(node_id)
        if sequence is not None:
            # Add task only if sequence exists and records exist
            tasks.append((node_id, file_read_offset, n_records, sequence, args.window)) # Pass window size
        else:
            # Node in index (with records) but sequence missing from GFA
            nodes_missing_sequence += 1

    if skipped_nodes_no_records > 0:
        print(f"ℹ️ Note: {skipped_nodes_no_records} nodes were listed in index header but had 0 records or were filtered.")
    if nodes_missing_sequence > 0:
        print(f"⚠️ Warning: Skipped {nodes_missing_sequence} nodes present in index (with records) but missing from GFA sequence data.")

    if not tasks:
        print("❌ Error: No tasks to process (no nodes with both index records and sequence data found). Exiting.", file=sys.stderr)
        sys.exit(1)

    total_tasks = len(tasks)
    num_workers = min(args.workers, total_tasks, os.cpu_count() or 1) # Cap workers
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
            # map processes tasks and returns results preserving order (though order doesn't matter for dict assembly)
            # Consider submit for more fine-grained control if needed.
            future_results = executor.map(process_node_parallel, tasks, chunksize=args.chunksize)

            # Iterate through results as they become available
            for node_id, pileup_dict in future_results:
                 nodes_processed_count += 1
                 if pileup_dict: # Store result only if pileups were generated (not empty dict)
                     results[node_id] = pileup_dict

                 # Print progress update periodically
                 if nodes_processed_count % max(1, total_tasks // 20) == 0 or nodes_processed_count == total_tasks: # Aim for ~20 updates
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
            # Save node IDs as strings for compatibility
            json.dump({str(k): v for k, v in results.items()}, out_f, indent=2) # Use indent for readability
        write_elapsed_time = time.time() - start_write_time
        print(f"✔ Output written in {write_elapsed_time:.2f} seconds.")
        print(f"✅ Done. Output saved to {args.output}")
    except IOError as ioe:
         print(f"❌ Error writing output JSON to {args.output}: {ioe}", file=sys.stderr)
         sys.exit(1)
    except Exception as e:
         print(f"❌ Unexpected error writing output JSON: {e}", file=sys.stderr)
         sys.exit(1)


if __name__ == '__main__':
    main()