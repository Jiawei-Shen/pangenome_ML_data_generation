#!/usr/bin/env python3
import json
import struct
import argparse
import sys
import time
import copy
import subprocess
import shutil  # For checking tool availability
import re  # For extracting chromosome name


# ─────────────────────────────────────────────────────────────────────────────
# Functions for reading IDX file (UNCHANGED)
# ─────────────────────────────────────────────────────────────────────────────
def load_index(idx_path):
    """
    Loads the node index from a .idx file into a dictionary.
    """
    node_index = {}
    try:
        with open(idx_path, "rb") as f:
            blocks_num_bytes = f.read(4)
            if not blocks_num_bytes or len(blocks_num_bytes) < 4:
                print(f"Error: Could not read blocks_num from {idx_path}. File might be empty or corrupted.",
                      file=sys.stderr)
                return {}
            blocks_num, = struct.unpack("<I", blocks_num_bytes)
            for i in range(blocks_num):
                header_data = f.read(4 + 8 + 4 + 4 + 2)
                if len(header_data) < (4 + 8 + 4 + 4 + 2):
                    # End of file reached prematurely
                    break
                block_id, block_start, block_size, n_records, metadata_len = struct.unpack("<I Q I I H", header_data)
                if metadata_len > 0:
                    skipped_bytes = f.read(metadata_len)
                    if len(skipped_bytes) < metadata_len:
                        # End of file reached prematurely
                        break
                node_index[block_id] = {"start": block_start, "size": block_size, "n_records": n_records}
    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except struct.error as e:
        print(f"Error: Could not unpack data from IDX file {idx_path}. Details: {e}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return {}
    return node_index


# ─────────────────────────────────────────────────────────────────────────────
# Function for reading JSON data, extracting Node IDs, and creating a node map (UNCHANGED)
# ─────────────────────────────────────────────────────────────────────────────
def load_json_data_ids_and_map(json_filepath):
    """
    Loads the main JSON file, extracts all node IDs, and creates a map for quick lookups.
    """
    main_json_data = {}
    node_ids_set = set()
    json_nodes_map = {}
    try:
        with open(json_filepath, 'r') as f_json:
            main_json_data = json.load(f_json)
            if not isinstance(main_json_data, dict):
                print(f"Error: Root JSON content in {json_filepath} is not an object/dictionary.", file=sys.stderr)
                return None, set(), {}
            node_list_from_json = main_json_data.get("nodes")
            if not isinstance(node_list_from_json, list):
                print(f"Warning: 'nodes' key not found or not a list in {json_filepath}.", file=sys.stderr)
                return main_json_data, set(), {}
            for item_index, item in enumerate(node_list_from_json):
                if not isinstance(item, dict): continue
                try:
                    node_id_str = item.get("node_id")
                    if node_id_str is None:
                        continue
                    node_id_int = int(node_id_str)
                    node_ids_set.add(node_id_int)
                    json_nodes_map[node_id_int] = item
                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Node ID format error for item {item_index} in 'nodes' (JSON): '{node_id_str}'. Error: {e}. Skipping.",
                        file=sys.stderr)
            return main_json_data, node_ids_set, json_nodes_map
    except FileNotFoundError:
        print(f"Error: JSON file not found: {json_filepath}", file=sys.stderr)
        return None, set(), {}
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {json_filepath}.", file=sys.stderr)
        return None, set(), {}
    except Exception as e:
        print(f"Unexpected error reading JSON {json_filepath}: {e}", file=sys.stderr)
        return None, set(), {}


# ─────────────────────────────────────────────────────────────────────────────
# Function to extract chromosome from path pattern (UNCHANGED)
# ─────────────────────────────────────────────────────────────────────────────
def extract_chromosome_from_path_pattern(path_pattern):
    """
    Extracts a chromosome name (e.g., 'chr1', 'chrX') from a string pattern.
    """
    if not path_pattern or not isinstance(path_pattern, str):
        return None
    match = re.search(r'(chr([0-9A-Za-z_]+|X|Y|M|MT))', path_pattern)
    if match:
        return match.group(1)
    return None


# ─────────────────────────────────────────────────────────────────────────────
# REVISED Function to check for overlap using bedtools
# ─────────────────────────────────────────────────────────────────────────────
def check_bed_overlap(bedtools_path, bed_filepath, chromosome, start_0based, end_0based):
    """
    Checks if a given genomic region overlaps with regions in a BED file using bedtools.
    Returns True if there is an overlap, False otherwise.
    """
    # bedtools uses 0-based, half-open intervals, which matches the required input format.
    node_region_bed_format = f"{chromosome}\t{start_0based}\t{end_0based}\n"
    command = [bedtools_path, 'intersect', '-a', 'stdin', '-b', bed_filepath]
    try:
        # We pass the node's region as a string to the stdin of the bedtools process.
        process = subprocess.run(
            command,
            input=node_region_bed_format,
            capture_output=True,
            text=True,
            check=False  # Do not raise an exception on a non-zero exit code.
        )

        # bedtools intersect outputs the intersecting line(s) to stdout.
        # If stdout is not empty, it means there was an overlap.
        if process.stdout:
            return True

        # Check for actual tool errors, but ignore messages about no overlaps.
        if process.returncode != 0:
            stderr_output = process.stderr.strip()
            if stderr_output:
                error_keywords = ["error", "fail", "could not open", "malformed"]
                if any(keyword in stderr_output.lower() for keyword in error_keywords):
                    print(f"Error from 'bedtools intersect' for region {chromosome}:{start_0based}-{end_0based}: {stderr_output}",
                          file=sys.stderr)
        return False

    except FileNotFoundError:
        # This error should ideally be caught once before entering the main loop.
        print(f"Error: bedtools command not found at '{bedtools_path}'. BED file operations will fail.", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Unexpected error during 'bedtools intersect' for region {chromosome}:{start_0based}-{end_0based}: {e}", file=sys.stderr)
        return False


# ─────────────────────────────────────────────────────────────────────────────
# REVISED Main processing function
# ─────────────────────────────────────────────────────────────────────────────
def filter_json_nodes_and_write(json_filepath, idx_filepath, output_json_filepath, bed_file=None, txt_output_path=None):
    """
    Main function to coordinate the loading, filtering, and writing of data.
    """
    bedtools_path = None
    if bed_file:
        bedtools_path = shutil.which("bedtools")
        if not bedtools_path:
            print("Warning: --bed_file provided, but bedtools not found in PATH. BED file filtering will be disabled.",
                  file=sys.stderr)
        else:
            print(f"Using bedtools found at: {bedtools_path}")

    print("Step 1: Loading JSON data and creating node map...")
    main_json_structure, target_node_ids_from_json, json_nodes_map = load_json_data_ids_and_map(json_filepath)
    if main_json_structure is None:
        print("Critical error: Could not load main JSON structure. Cannot proceed.", file=sys.stderr)
        return
    print(f"Loaded {len(json_nodes_map)} nodes from JSON for potential processing.")

    chromosome_for_filter = None
    if bed_file and bedtools_path:
        path_pattern = main_json_structure.get("path_name_input_pattern")
        chromosome_for_filter = extract_chromosome_from_path_pattern(path_pattern)
        if not chromosome_for_filter:
            print(
                f"Warning: Could not extract chromosome from 'path_name_input_pattern': '{path_pattern}'. "
                "Filtering will fail for nodes unless they have required coordinate information.",
                file=sys.stderr)
        else:
            print(f"Extracted chromosome '{chromosome_for_filter}' for use in filtering.")

    print("\nStep 2: Loading node index from IDX file...")
    idx_data = load_index(idx_filepath)
    available_idx_node_ids = set(idx_data.keys())
    print(f"Found {len(available_idx_node_ids)} unique node IDs in {idx_filepath}")

    print("\nStep 3: Identifying common node IDs between JSON and IDX...")
    common_node_ids_int = target_node_ids_from_json.intersection(available_idx_node_ids)
    num_initial_common_nodes = len(common_node_ids_int)
    if num_initial_common_nodes == 0:
        print("No common node IDs found. Output 'nodes' list will be empty.")
    else:
        print(f"Found {num_initial_common_nodes} common node IDs to process.")

    # Initialize lists and counters for summary statistics
    ultimate_filtered_nodes_list = []
    total_nodes_queried_with_bedtools = 0
    total_nodes_filtered_out_by_prereq = 0
    total_nodes_filtered_out_by_no_overlap = 0
    total_nodes_kept_after_overlap_check = 0

    batch_node_counter_step4 = 0
    batch_queried_step4 = 0
    batch_kept_step4 = 0
    batch_start_time_step4 = time.time()

    if num_initial_common_nodes > 0:
        print(f"\nStep 4: Processing {num_initial_common_nodes} common nodes for BED file filtering (if applicable)...")
        sorted_common_ids = sorted(list(common_node_ids_int))

        for node_idx_step4, node_id_int in enumerate(sorted_common_ids):
            current_node_item = json_nodes_map.get(node_id_int)
            if not current_node_item:
                sys.stderr.write(f"Internal Warning: Common node ID {node_id_int} not found in JSON map. Skipping.\n")
                continue

            include_node_in_final_output = True

            # If a BED file is provided, perform filtering.
            if bed_file and bedtools_path:
                try:
                    start_0based_str = current_node_item.get("grch38_position_start")
                    length_str = current_node_item.get("length")

                    # Check for prerequisites: chromosome, start, and length must be valid.
                    if start_0based_str is None or length_str is None or chromosome_for_filter is None:
                        include_node_in_final_output = False
                        total_nodes_filtered_out_by_prereq += 1
                    else:
                        start_0based = int(start_0based_str)
                        length = int(length_str)
                        if length <= 0:
                            include_node_in_final_output = False
                            total_nodes_filtered_out_by_prereq += 1
                        else:
                            # BED format is [start, end), so end is start + length.
                            end_0based = start_0based + length
                            total_nodes_queried_with_bedtools += 1
                            batch_queried_step4 += 1

                            # Perform the overlap check using bedtools.
                            has_overlap = check_bed_overlap(bedtools_path, bed_file, chromosome_for_filter,
                                                            start_0based, end_0based)

                            if not has_overlap:
                                include_node_in_final_output = False
                                total_nodes_filtered_out_by_no_overlap += 1
                            else:
                                total_nodes_kept_after_overlap_check += 1
                                batch_kept_step4 += 1

                except (ValueError, TypeError):
                    # This catches errors from int() conversion if values are not valid numbers.
                    include_node_in_final_output = False
                    total_nodes_filtered_out_by_prereq += 1

            if include_node_in_final_output:
                # The node passed all filters (or no filters were applied), so add it to the final list.
                ultimate_filtered_nodes_list.append(current_node_item)

            # --- Batch progress reporting ---
            batch_node_counter_step4 += 1
            if (batch_node_counter_step4 % 1000 == 0) or (node_idx_step4 + 1 == num_initial_common_nodes):
                if batch_node_counter_step4 > 0:
                    batch_time_step4 = time.time() - batch_start_time_step4
                    print(
                        f"  Processed {node_idx_step4 + 1}/{num_initial_common_nodes} common nodes. "
                        f"In last batch of {batch_node_counter_step4}: Queried={batch_queried_step4}, Kept={batch_kept_step4}. "
                        f"Time: {batch_time_step4:.2f}s.")
                    batch_node_counter_step4, batch_queried_step4, batch_kept_step4 = 0, 0, 0
                    batch_start_time_step4 = time.time()

    # --- Final Output Generation ---
    output_json_structure = copy.deepcopy(main_json_structure)
    output_json_structure["nodes"] = ultimate_filtered_nodes_list
    num_ultimate_nodes = len(ultimate_filtered_nodes_list)
    print(f"\nStep 5: Final processing complete. A total of {num_ultimate_nodes} nodes will be written to the output.")

    if bed_file and bedtools_path:
        print("\nBED Filter Summary:")
        print(f"  - Nodes eligible and queried with bedtools: {total_nodes_queried_with_bedtools}")
        print(f"  - Nodes filtered out due to missing prerequisites (coords/length): {total_nodes_filtered_out_by_prereq}")
        print(f"  - Nodes filtered out due to no overlap in BED file: {total_nodes_filtered_out_by_no_overlap}")
        print(f"  - Nodes kept that overlapped with a BED region: {total_nodes_kept_after_overlap_check}")

    if txt_output_path:
        print(f"\nStep 6: Writing {num_ultimate_nodes} filtered node ID(s) to {txt_output_path}...")
        try:
            with open(txt_output_path, 'w') as f_txt:
                if ultimate_filtered_nodes_list:
                    ids_to_write = [str(node.get("node_id")) for node in ultimate_filtered_nodes_list if node.get("node_id") is not None]
                    try: sorted_ids = sorted(ids_to_write, key=int)
                    except ValueError: sorted_ids = sorted(ids_to_write)
                    for node_id_str in sorted_ids: f_txt.write(f"{node_id_str}\n")
            print(f"Successfully wrote node IDs to {txt_output_path}")
        except Exception as e:
            print(f"Error writing filtered node IDs to TXT file {txt_output_path}: {e}", file=sys.stderr)

    print(f"\nStep 7: Writing final JSON (with {num_ultimate_nodes} nodes) to {output_json_filepath}...")
    try:
        with open(output_json_filepath, 'w') as f_out_json:
            json.dump(output_json_structure, f_out_json, indent=4)
        print(f"Successfully wrote resulting JSON to {output_json_filepath}")
    except Exception as e:
        print(f"Error writing output JSON {output_json_filepath}: {e}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# REVISED Command-line interface
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Filter 'nodes' in a JSON file based on an IDX file and, optionally, a BED file of high-confidence regions. "
                    "Outputs a modified JSON and optionally a TXT list of the kept node IDs.",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("json_path", help="Path to the input JSON file.")
    parser.add_argument("idx_path", help="Path to the corresponding .idx file.")
    parser.add_argument("output_json_path", help="Path for the output filtered JSON file.")
    parser.add_argument("--bed_file",
                        help="Optional: Path to a BED file containing regions for filtering.\n"
                             "Nodes will be kept only if their genomic coordinates overlap with a region in this file.\n"
                             "Requires 'bedtools' to be installed and in the system's PATH.",
                        default=None)
    parser.add_argument("--txt", dest="txt_output_path",
                        help="Optional: Path to write the final list of filtered node IDs (one per line).",
                        default=None)
    args = parser.parse_args()

    filter_json_nodes_and_write(
        args.json_path,
        args.idx_path,
        args.output_json_path,
        bed_file=args.bed_file,
        txt_output_path=args.txt_output_path
    )


if __name__ == "__main__":
    main()
