#!/usr/bin/env python3

import argparse
import logging
import struct
import json
import os
import multiprocessing
import time
import sys
import io
import re
from collections import defaultdict

# --- Constants ---
# Map DNA bases and 'N' to integer indices for the pileup matrix
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4,
                 'a': 0, 'c': 1, 'g': 2, 't': 3, 'n': 4}
INDEX_TO_BASE = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N'}

# Define the structure of binary records in the .dat file
# Adjust field sizes (500s, 100s) if your sequences/CIGARs are longer
# Fields: offset (int), sequence (bytes), qual_placeholder(int), cigar (bytes), mapq (int), strand (char)
RECORD_STRUCT_FORMAT = '<i 500s i 100s i c'
RECORD_STRUCT = struct.Struct(RECORD_STRUCT_FORMAT)
RECORD_SIZE = RECORD_STRUCT.size

# --- Helper Functions ---

def setup_logging(level=logging.INFO):
    """Configures logging."""
    # Include process ID in logging format for parallel debugging
    log_format = '%(asctime)s [%(levelname)s] [%(process)d] %(message)s'
    logging.basicConfig(level=level,
                        format=log_format,
                        datefmt='%Y-%m-%d %H:%M:%S')

def reverse_complement(sequence):
    """Computes the reverse complement of a DNA sequence."""
    if not isinstance(sequence, str):
        logging.warning(f"Attempted to reverse complement non-string: {type(sequence)}. Returning empty string.")
        return ""
    try:
        complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return sequence.translate(complement_map)[::-1]
    except Exception as e:
        logging.error(f"Error in reverse_complement for sequence '{sequence[:50]}...': {e}")
        return ""


def parse_cigar(cigar_string):
    """Parses a CIGAR string into a list of (length, operation) tuples."""
    if not cigar_string or cigar_string == '*':
        return []
    try:
        # Use regex to find all occurrences of number + letter
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        logging.error(f"Error parsing CIGAR string '{cigar_string}': {e}")
        return []

def detect_variants_from_cigar(ref_seq, read_seq, read_offset, cigar_string):
    """
    Detects variants (mismatches, insertions, deletions) based on CIGAR alignment.

    Args:
        ref_seq (str): The reference node sequence.
        read_seq (str): The read sequence (already oriented to reference).
        read_offset (int): The 0-based start position of the alignment on the reference.
        cigar_string (str): The CIGAR string for the alignment.

    Yields:
        tuple: (ref_pos, variant_type, ref_base, alt_base)
               variant_type is 'X' (mismatch), 'I' (insertion), 'D' (deletion).
               For 'I', ref_base is '-', alt_base is inserted sequence. ref_pos is *before* insertion.
               For 'D', ref_base is deleted sequence, alt_base is '-'. ref_pos is start of deletion.
               For 'X', ref_base is from ref_seq, alt_base is from read_seq. ref_pos is mismatch location.
               Bases are returned in uppercase.
    """
    if not ref_seq or not read_seq or read_offset < 0:
        logging.debug(f"detect_variants returning early: ref_seq empty={not ref_seq}, read_seq empty={not read_seq}, read_offset={read_offset}")
        return

    ref_len = len(ref_seq)
    read_len = len(read_seq)
    parsed_cigar = parse_cigar(cigar_string)

    ref_ptr = read_offset
    read_ptr = 0
    # Position *before* which an insertion would occur (0-based)
    # Needs careful handling at the start of alignment/after deletions
    last_valid_ref_ptr = read_offset - 1

    try:
        for length, op in parsed_cigar:
            if op in ('M', '=', 'X'): # Match, Alignment match, Alignment mismatch
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    current_read_pos = read_ptr + i

                    # Check bounds carefully before accessing sequences
                    if current_ref_pos >= ref_len or current_read_pos >= read_len:
                        logging.debug(f"Bounds check failed in M/= /X: Ref {current_ref_pos}/{ref_len}, Read {current_read_pos}/{read_len}, CIGAR {cigar_string}, Offset {read_offset}")
                        # Stop processing this segment, likely an issue with alignment/data
                        return # Use return instead of break to stop the generator entirely for this read

                    ref_base = ref_seq[current_ref_pos].upper()
                    read_base = read_seq[current_read_pos].upper()
                    last_valid_ref_ptr = current_ref_pos # Update last aligned ref base position

                    if op != '=' and ref_base != 'N' and read_base != 'N' and ref_base != read_base:
                         # Report mismatch ('X')
                         yield (current_ref_pos, 'X', ref_base, read_base)

                ref_ptr += length
                read_ptr += length

            elif op == 'I': # Insertion to the reference
                 # Check bounds before accessing read sequence
                 if read_ptr + length > read_len:
                     logging.debug(f"Bounds check failed in I: Read {read_ptr + length}/{read_len}, CIGAR {cigar_string}, Offset {read_offset}")
                     return # Stop generator

                 inserted_bases = read_seq[read_ptr : read_ptr + length].upper()
                 # Insertion occurs *after* the last aligned reference base position
                 yield (last_valid_ref_ptr, 'I', '-', inserted_bases)
                 read_ptr += length
                 # Does not advance ref_ptr or change last_valid_ref_ptr

            elif op == 'D': # Deletion from the reference
                 start_del_pos = ref_ptr # Deletion starts at current ref pointer
                 # Check bounds before accessing reference sequence
                 if ref_ptr + length > ref_len:
                     logging.debug(f"Bounds check failed in D: Ref {ref_ptr + length}/{ref_len}, CIGAR {cigar_string}, Offset {read_offset}")
                     return # Stop generator

                 deleted_bases = ref_seq[ref_ptr : ref_ptr + length].upper()
                 # Deletion occurs *at* the current reference position
                 yield (start_del_pos, 'D', deleted_bases, '-')
                 ref_ptr += length
                 last_valid_ref_ptr = ref_ptr - 1 # Update last aligned position to end of deletion
                 # Does not advance read_ptr

            elif op == 'S': # Soft clipping (consumes read)
                 read_ptr += length
                 # Soft clips are often at ends, doesn't usually affect last_valid_ref_ptr for internal indels

            elif op == 'N': # Skipped region from reference (consumes reference)
                 ref_ptr += length
                 last_valid_ref_ptr = ref_ptr - 1 # Update last aligned position
                 # Does not advance read_ptr

            elif op in ('H', 'P'): # Hard clipping, Padding (consumes neither)
                 pass

    except IndexError:
        # This can happen with malformed CIGAR or inconsistent read/ref lengths
        logging.warning(f"IndexError during CIGAR processing. RefLen={ref_len}, ReadLen={read_len}, Offset={read_offset}, CIGAR='{cigar_string}', RefPtr={ref_ptr}, ReadPtr={read_ptr}")
    except Exception as e:
        logging.error(f"Unexpected error in detect_variants_from_cigar (Offset={read_offset}, CIGAR='{cigar_string}'): {e}", exc_info=False) # Set exc_info=True for trace


def create_pileup_row(read_seq, read_offset, cigar_str, window_start, window_end, ref_len):
    """
    Creates a single row for the pileup matrix corresponding to one read segment
    within the specified reference window.

    Args:
        read_seq (str): The read sequence (oriented to the reference forward strand).
        read_offset (int): 0-based starting position of the read on the reference.
        cigar_str (str): CIGAR string of the alignment.
        window_start (int): 0-based start position of the window on the reference.
        window_end (int): 0-based end position (exclusive) of the window on the reference.
        ref_len (int): Total length of the reference sequence (node_len).

    Returns:
        list[int] or None: A list of base indices (0-4) of length (window_end - window_start),
                           or None if the read does not overlap the window according to CIGAR.
    """
    if not read_seq or window_start >= window_end:
        return None

    window_len = window_end - window_start
    # Initialize pileup row with 'N' index
    pileup_row = [BASE_TO_INDEX['N']] * window_len
    read_len = len(read_seq)
    parsed_cigar = parse_cigar(cigar_str)

    ref_ptr = read_offset
    read_ptr = 0
    row_filled = False # Track if any base was actually placed in the window

    try:
        for length, op in parsed_cigar:
            # Optimization: If current ref pointer is already past the window end,
            # and we haven't filled anything, the read won't overlap.
            if ref_ptr >= window_end and not row_filled:
                 return None

            if op in ('M', '=', 'X'): # Consumes both Ref and Read
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    current_read_pos = read_ptr + i

                    # Check read bounds FIRST before accessing read_seq
                    if current_read_pos >= read_len:
                         logging.debug(f"Read index out of bounds: {current_read_pos}/{read_len} during CIGAR op {op}{length}, CIGAR '{cigar_str}', Offset {read_offset}")
                         # Stop processing this read's pileup if CIGAR/read mismatch
                         return pileup_row if row_filled else None

                    # Check if current reference position falls within the window
                    if window_start <= current_ref_pos < window_end:
                        pileup_idx = current_ref_pos - window_start
                        base = read_seq[current_read_pos]
                        pileup_row[pileup_idx] = BASE_TO_INDEX.get(base.upper(), BASE_TO_INDEX['N']) # Use uppercase base
                        row_filled = True

                ref_ptr += length
                read_ptr += length

            elif op == 'I': # Insertion (consumes Read only)
                # Check read bounds before incrementing pointer
                if read_ptr + length > read_len:
                     logging.debug(f"Read index out of bounds for insertion: {read_ptr + length}/{read_len}, CIGAR '{cigar_str}', Offset {read_offset}")
                     return pileup_row if row_filled else None
                read_ptr += length
                # Insertions don't occupy reference space in the pileup row

            elif op == 'D': # Deletion (consumes Ref only)
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    if window_start <= current_ref_pos < window_end:
                        pileup_idx = current_ref_pos - window_start
                        # Pileup represents alignment to reference, so deleted ref bases
                        # mean the read provides no info here. Keep 'N' (or use specific gap marker).
                        pileup_row[pileup_idx] = BASE_TO_INDEX['N'] # Explicitly mark N for deletion
                        row_filled = True
                ref_ptr += length

            elif op == 'S': # Soft clipping (consumes Read only)
                read_ptr += length

            elif op == 'N': # Skipped region (consumes Ref only)
                 for i in range(length):
                    current_ref_pos = ref_ptr + i
                    if window_start <= current_ref_pos < window_end:
                         pileup_idx = current_ref_pos - window_start
                         # Similar to deletion, skipped reference region means no read info at these ref coords
                         pileup_row[pileup_idx] = BASE_TO_INDEX['N'] # Explicitly mark N for skipped region
                         row_filled = True
                 ref_ptr += length

            elif op in ('H', 'P'): # Consumes nothing
                 pass

        # Return the row only if it overlapped the window meaningfully
        return pileup_row if row_filled else None

    except IndexError:
         # This might indicate CIGAR goes past reported read length
         logging.warning(f"IndexError during pileup creation. ReadLen={read_len}, Offset={read_offset}, CIGAR='{cigar_str}', RefPtr={ref_ptr}, ReadPtr={read_ptr}, Win=[{window_start}:{window_end}]")
         return pileup_row if row_filled else None # Return partial row if possible
    except Exception as e:
         logging.error(f"Unexpected error in create_pileup_row (Offset={read_offset}, CIGAR='{cigar_str}', Win=[{window_start}:{window_end}]): {e}", exc_info=False)
         return None


def process_node_parallel(args_tuple):
    """
    Processes reads associated with a single node to find variants and generate pileups.
    Designed to be run in parallel.
    """
    node_id, node_sequence, dat_file_path, idx_data, window, min_mapq = args_tuple
    if not node_sequence:
         # logging.warning(f"Node {node_id}: Skipping due to empty sequence.") # Too verbose
         return node_id, None

    node_len = len(node_sequence)
    half_window = (window - 1) // 2 # Integer division for centering
    segments = [] # Store (offset, seq, cigar) tuples for reads passing filters

    # --- Determine windowing strategy ---
    use_full_node_as_window = node_len < window
    # logging.debug(f"Node {node_id}: Length={node_len}, Window={window}. Use full node as window: {use_full_node_as_window}")
    # --- End windowing strategy ---

    # Check if index data exists for this node_id (using string key)
    node_id_str = str(node_id) # Ensure key is string if index uses string keys
    if node_id_str not in idx_data:
        # logging.warning(f"No index data found for node {node_id_str}. Skipping.") # Too verbose if many nodes skipped
        return node_id_str, None # Return None if node not in index

    file_offset, num_records = idx_data[node_id_str]

    # --- Read data from .dat file ---
    try:
        with open(dat_file_path, 'rb') as f:
            f.seek(file_offset)
            for i in range(num_records):
                record_start_pos = f.tell()
                data = f.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE:
                    logging.warning(f"Node {node_id_str}: Incomplete record read at offset {record_start_pos}. Expected {RECORD_SIZE}, got {len(data)}. Read {i+1}/{num_records}.")
                    break # Stop reading if record is incomplete

                try:
                    # Unpack data according to the defined structure
                    off, raw_seq, _, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)

                    # Basic filtering
                    if mapq < min_mapq:
                        continue

                    seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
                    cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
                    strand_char = strand_byte.decode('ascii', errors='ignore')

                    if not seq or not cigar or cigar == '*':
                        # logging.debug(f"Node {node_id_str}: Skipping record with empty sequence or CIGAR.")
                        continue

                    # Handle reverse strand
                    current_seq = seq
                    if strand_char == '-':
                        current_seq = reverse_complement(seq)
                        if not current_seq: # Handle potential error in RC
                             # logging.warning(f"Node {node_id_str}: Skipping read at offset {off} due to reverse complement error.") # Can be verbose
                             continue

                    # Validate offset
                    if off < 0: # Offset should not be negative
                       logging.warning(f"Node {node_id_str}: Invalid negative offset {off} encountered. Skipping read.")
                       continue
                    # Rough check if alignment start makes sense (can be refined)
                    # CIGAR processing later will handle exact bounds
                    # if off >= node_len:
                    #    logging.warning(f"Node {node_id_str}: Read offset {off} >= node length {node_len}. Skipping.")
                    #    continue

                    segments.append((off, current_seq, cigar))

                except struct.error as e:
                    logging.error(f"Node {node_id_str}: Error unpacking record at file offset {record_start_pos}: {e}")
                    continue # Skip corrupted record
                except UnicodeDecodeError as e:
                    logging.error(f"Node {node_id_str}: Error decoding sequence/CIGAR/strand at file offset {record_start_pos}: {e}")
                    continue # Skip record with decoding errors
                except Exception as e:
                     logging.error(f"Node {node_id_str}: Unexpected error processing record from offset {record_start_pos}: {e}")
                     continue

    except FileNotFoundError:
        logging.error(f"Data file not found: {dat_file_path}")
        return node_id_str, None
    except Exception as e:
        logging.error(f"Error reading file {dat_file_path} for node {node_id_str}: {e}")
        return node_id_str, None

    if not segments:
        # logging.debug(f"Node {node_id_str}: No valid read segments found after filtering.") # Too verbose
        return node_id_str, None

    # --- Process collected segments ---
    variants_pileups = defaultdict(list)

    for offset, read_seq, cigar in segments:
        try:
            # Detect variants based on CIGAR string and node sequence
            # Pass oriented read_seq here
            detected_vars = list(detect_variants_from_cigar(node_sequence, read_seq, offset, cigar))

            if not detected_vars:
                continue

            processed_windows_for_read = set() # Avoid redundant pileup if read covers multiple variant windows

            for vpos, vtype, ref_base, alt_base in detected_vars:

                # --- Generate the specific variant key ---
                # Includes position, type, ref sequence, and alt sequence
                # Handles multi-base indels uniquely.
                variant_key = f"{vpos}_{vtype}_{ref_base}_{alt_base}"
                # --- End key generation ---

                # Determine window bounds based on node length vs window size
                if use_full_node_as_window:
                    # Use the entire node sequence as the window
                    start_pos = 0
                    end_pos = node_len
                else:
                    # Calculate window around the variant position
                    # Center window on the variant position (vpos)
                    center_pos = vpos
                    start_pos = max(0, center_pos - half_window)
                    end_pos = min(node_len, center_pos + half_window + 1) # Use +1 for exclusive end

                # Ensure start is not negative and end does not exceed node length
                start_pos = max(0, start_pos)
                end_pos = min(node_len, end_pos)

                # Sanity check window validity
                if start_pos >= end_pos:
                    # logging.debug(f"Node {node_id_str}, VarKey {variant_key}: Invalid window [{start_pos}:{end_pos}]")
                    continue # Skip if window is zero or negative width

                window_tuple = (start_pos, end_pos)
                if window_tuple in processed_windows_for_read:
                    # logging.debug(f"Node {node_id_str}, Read offset {offset}: Skipping pileup generation for already processed window {window_tuple}")
                    continue # Already created pileup for this read covering this exact window region

                # Create the pileup row for THIS read segment covering the calculated window
                pileup_row = create_pileup_row(
                    read_seq, offset, cigar, start_pos, end_pos, node_len
                )

                # Only add the row if it was successfully created (e.g., the read covers part of the window)
                if pileup_row:
                    variants_pileups[variant_key].append(pileup_row)
                    processed_windows_for_read.add(window_tuple) # Mark window as processed for this read

        except Exception as e:
            logging.error(f"Node {node_id_str}: Error processing read segment (offset={offset}, cigar='{cigar}') for variants/pileup: {e}", exc_info=False) # Set exc_info=True for full traceback
            continue # Skip this problematic read segment

    # Filter out variants that ended up with no supporting pileup rows
    final_variants_pileups = {k: v for k, v in variants_pileups.items() if v}

    if not final_variants_pileups:
        # logging.debug(f"Node {node_id_str}: No variants with valid pileups found.") # Too verbose
        return node_id_str, None

    logging.info(f"Node {node_id_str}: Processed. Found {len(final_variants_pileups)} variants with pileups.")
    return node_id_str, final_variants_pileups


def load_index(index_file):
    """Loads the index file mapping node IDs to file offsets and record counts."""
    idx_data = {}
    logging.info(f"Loading index from {index_file}...")
    try:
        with open(index_file, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                line = line.strip()
                if not line or line.startswith('#'): continue # Skip empty/comment lines
                try:
                    parts = line.split('\t')
                    if len(parts) != 3:
                         logging.warning(f"Skipping malformed line {line_num} in index file (expected 3 tab-separated fields): {line}")
                         continue
                    node_id_str, offset_str, count_str = parts
                    idx_data[node_id_str] = (int(offset_str), int(count_str))
                except ValueError:
                    logging.warning(f"Skipping malformed line {line_num} in index file (numeric conversion error): {line}")
                    continue
    except FileNotFoundError:
        logging.error(f"Index file not found: {index_file}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error loading index file {index_file}: {e}")
        sys.exit(1)
    logging.info(f"Loaded index for {len(idx_data)} nodes.")
    return idx_data


def load_node_sequences(sequence_file):
    """
    Loads node sequences from a GFA or a simple TSV cache file.
    Expects TSV format: node_id\tsequence
    Expects GFA format: Extracts S lines (Segment lines)
    """
    node_sequences = {}
    logging.info(f"Loading node sequences from {sequence_file}...")
    file_ext = os.path.splitext(sequence_file)[1].lower()
    is_gfa = file_ext == ".gfa"

    try:
        with open(sequence_file, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                line = line.strip()
                if not line: continue

                try:
                    if is_gfa:
                        if line.startswith('S'):
                            parts = line.split('\t')
                            if len(parts) >= 3:
                                node_id_str = parts[1]
                                seq = parts[2]
                                if seq != '*' and re.match(r'^[ACGTNacgtn]+$', seq): # Basic sequence validation
                                     node_sequences[node_id_str] = seq.upper() # Store uppercase
                                else:
                                     logging.warning(f"Skipping GFA Segment line {line_num} with invalid sequence: Node {node_id_str}, Seq '{seq[:20]}...'")
                            else:
                                logging.warning(f"Skipping malformed GFA Segment line {line_num}: {line}")
                    else: # Assume TSV format
                        parts = line.split('\t')
                        if len(parts) == 2:
                            node_id_str = parts[0]
                            seq = parts[1]
                            if re.match(r'^[ACGTNacgtn]+$', seq): # Basic sequence validation
                                node_sequences[node_id_str] = seq.upper() # Store uppercase
                            else:
                                logging.warning(f"Skipping TSV line {line_num} with invalid sequence: Node {node_id_str}, Seq '{seq[:20]}...'")
                        else:
                             logging.warning(f"Skipping malformed TSV line {line_num} (expected 2 tab-separated fields): {line}")
                except Exception as e:
                     logging.warning(f"Error processing line {line_num}: '{line[:50]}...' - {e}")
                     continue # Skip problematic line

    except FileNotFoundError:
        logging.error(f"Sequence file not found: {sequence_file}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error loading sequence file {sequence_file}: {e}")
        sys.exit(1)

    if not node_sequences:
         logging.error(f"No node sequences loaded from {sequence_file}. Check file format and content.")
         sys.exit(1)

    logging.info(f"Loaded sequences for {len(node_sequences)} nodes.")
    return node_sequences


# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate pileup matrices for variants from aligned reads (.dat format).")
    parser.add_argument("sequence_file", help="Path to the node sequence file (GFA or node_id\\tsequence TSV).")
    parser.add_argument("dat_file", help="Path to the binary alignment data file (.dat).")
    parser.add_argument("index_file", help="Path to the index file for the .dat file (node_id\\toffset\\tcount TSV).")
    parser.add_argument("output_json", help="Path to the output JSON file.")
    parser.add_argument("-w", "--window", type=int, default=75, help="Window size around the variant for pileup (default: 75). Must be odd.")
    parser.add_argument("-q", "--min_mapq", type=int, default=10, help="Minimum mapping quality (MAPQ) to consider a read (default: 10).")
    parser.add_argument("-p", "--processes", type=int, default=multiprocessing.cpu_count(), help="Number of parallel processes to use (default: number of CPU cores).")
    parser.add_argument("-v", "--verbose", action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.INFO, help="Enable verbose logging.")

    args = parser.parse_args()

    setup_logging(args.loglevel)
    logging.info(f"Script started at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Current time: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}") # Add timezone if possible
    logging.info(f"Command: {' '.join(sys.argv)}")
    logging.info(f"Args: {args}")


    # Validate window size
    if args.window <= 0 or args.window % 2 == 0:
        logging.error("Window size must be a positive odd integer.")
        sys.exit(1)

    start_time = time.time()

    # Load necessary data
    node_sequences = load_node_sequences(args.sequence_file)
    idx_data = load_index(args.index_file)

    # --- Prepare arguments for parallel processing ---
    tasks = []
    nodes_in_index = set(idx_data.keys())
    nodes_in_seqs = set(node_sequences.keys())

    # Find nodes present in both index and sequence file
    common_nodes = nodes_in_index.intersection(nodes_in_seqs)
    # Find nodes missing in one or the other (for reporting)
    missing_seqs = nodes_in_index - nodes_in_seqs
    missing_index = nodes_in_seqs - nodes_in_index

    if missing_seqs:
         logging.warning(f"{len(missing_seqs)} nodes found in index but not in sequence file. Examples: {list(missing_seqs)[:5]}")
    if missing_index:
         logging.warning(f"{len(missing_index)} nodes found in sequence file but not in index. Examples: {list(missing_index)[:5]}")

    logging.info(f"Preparing tasks for {len(common_nodes)} nodes common to index and sequence file.")

    for node_id_str in common_nodes:
        node_seq = node_sequences.get(node_id_str) # Get sequence using string ID
        # Small sanity check, though common_nodes should guarantee presence
        if node_seq:
             task_args = (
                 node_id_str,         # Use string ID consistent with keys
                 node_seq,
                 args.dat_file,
                 idx_data,            # Pass the whole index dict (might be large, consider alternatives if memory is extreme)
                 args.window,
                 args.min_mapq
             )
             tasks.append(task_args)
        else: # Should not happen if common_nodes logic is correct
             logging.error(f"Internal consistency error: Node {node_id_str} sequence unexpectedly missing during task preparation.")

    if not tasks:
         logging.error("No common nodes found between index and sequence files. No tasks to run. Exiting.")
         sys.exit(1)


    # --- Run parallel processing ---
    all_results = {}
    processed_count = 0
    nodes_with_pileups = 0
    logging.info(f"Starting parallel processing with {args.processes} processes for {len(tasks)} tasks...")

    # Using imap_unordered for potentially better performance with uneven task times
    try:
        with multiprocessing.Pool(processes=args.processes) as pool:
            results_iterator = pool.imap_unordered(process_node_parallel, tasks)

            for i, result in enumerate(results_iterator):
                processed_count = i + 1
                if result: # Check if result is not None
                    node_id_res, node_pileups = result
                    if node_pileups: # Only add if pileups were actually generated
                        all_results[node_id_res] = node_pileups
                        nodes_with_pileups += 1

                # Log progress periodically without flooding logs
                if processed_count % max(1, len(tasks) // 20000) == 0 or processed_count == len(tasks): # Log ~20 times + end
                    logging.info(f"Processed {processed_count}/{len(tasks)} nodes... ({nodes_with_pileups} nodes have pileups so far)")

    except Exception as e:
        logging.error(f"An error occurred during parallel processing: {e}", exc_info=True)
        # Consider cleanup or partial results saving here if needed
        sys.exit(1)
    finally:
        # Ensure pool is closed properly even if errors occur above
        # The 'with' statement handles this automatically.
        pass


    logging.info(f"Finished processing {processed_count} nodes.")
    logging.info(f"Found variants with pileups in {nodes_with_pileups} nodes.")

    # --- Write results to JSON ---
    logging.info(f"Writing results to {args.output_json}...")
    try:
        # Use a compact separator for potentially smaller file size
        with open(args.output_json, 'w') as outfile:
            json.dump(all_results, outfile, separators=(',', ':'))
            # Use indent=2 for human readability instead of separators:
            # json.dump(all_results, outfile, indent=2)
    except IOError as e:
        logging.error(f"Error writing output JSON file {args.output_json}: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred while writing JSON output: {e}")
        sys.exit(1)


    end_time = time.time()
    logging.info(f"Script finished successfully in {end_time - start_time:.2f} seconds.")
    logging.info(f"Output written to {args.output_json}")