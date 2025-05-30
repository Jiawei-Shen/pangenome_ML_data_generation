import json
import struct
import argparse
import sys
import copy
import subprocess
import shutil  # For checking bcftools availability
import re  # For extracting chromosome name


# ─────────────────────────────────────────────────────────────────────────────
# Functions for reading IDX file (remains the same)
# ─────────────────────────────────────────────────────────────────────────────
def load_index(idx_path):
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
                    print(
                        f"Warning: Truncated block header in {idx_path} at block index {i}. Reached end of file unexpectedly.",
                        file=sys.stderr)
                    break
                block_id, block_start, block_size, n_records, metadata_len = struct.unpack("<I Q I I H", header_data)
                if metadata_len > 0:
                    skipped_bytes = f.read(metadata_len)
                    if len(skipped_bytes) < metadata_len:
                        print(
                            f"Warning: Truncated metadata in {idx_path} for block_id {block_id}. Reached end of file unexpectedly.",
                            file=sys.stderr)
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
# Function for reading JSON data, extracting Node IDs, and creating a node map
# ─────────────────────────────────────────────────────────────────────────────
def load_json_data_ids_and_map(json_filepath):
    """
    Reads the main JSON object, extracts node IDs from the "nodes" list,
    and creates a map of integer node_id to the original node object.
    Returns:
        tuple: (main_json_object, set_of_integer_node_ids, map_of_int_node_id_to_node_object)
               Returns (None, set(), {}) if a critical error occurs.
    """
    main_json_data = {}
    node_ids_set = set()
    json_nodes_map = {}  # Map int(node_id) to node_item object
    try:
        with open(json_filepath, 'r') as f_json:
            main_json_data = json.load(f_json)
            if not isinstance(main_json_data, dict):
                print(f"Error: Root JSON content in {json_filepath} is not an object/dictionary.", file=sys.stderr)
                return None, set(), {}

            node_list_from_json = main_json_data.get("nodes")
            if not isinstance(node_list_from_json, list):
                print(
                    f"Warning: 'nodes' key not found or not a list in {json_filepath}. No nodes to process from JSON.",
                    file=sys.stderr)
                return main_json_data, set(), {}

            for item_index, item in enumerate(node_list_from_json):
                if not isinstance(item, dict):
                    continue
                try:
                    node_id_str = item.get("node_id")
                    if node_id_str is None:
                        print(
                            f"Warning: Item at index {item_index} in 'nodes' list (JSON) is missing 'node_id'. Skipping.",
                            file=sys.stderr)
                        continue
                    node_id_int = int(node_id_str)
                    node_ids_set.add(node_id_int)
                    json_nodes_map[node_id_int] = item
                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Could not convert node_id '{node_id_str}' to int or item format error for item {item_index} in 'nodes' list (JSON). Error: {e}. Skipping.",
                        file=sys.stderr)

            num_items_in_list = len(node_list_from_json)
            num_unique_ids = len(node_ids_set)
            num_mapped_nodes = len(json_nodes_map)
            print(f"Successfully read JSON from {json_filepath}. Found {num_items_in_list} items in 'nodes' list.")
            print(f"Extracted {num_unique_ids} unique node IDs and mapped {num_mapped_nodes} nodes for lookup.")
            return main_json_data, node_ids_set, json_nodes_map

    except FileNotFoundError:
        print(f"Error: JSON file not found at {json_filepath}", file=sys.stderr)
        return None, set(), {}
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {json_filepath}.", file=sys.stderr)
        return None, set(), {}
    except Exception as e:
        print(f"An unexpected error occurred while reading JSON file {json_filepath}: {e}", file=sys.stderr)
        return None, set(), {}


# ─────────────────────────────────────────────────────────────────────────────
# Function to extract chromosome from path pattern (remains the same)
# ─────────────────────────────────────────────────────────────────────────────
def extract_chromosome_from_path_pattern(path_pattern):
    if not path_pattern or not isinstance(path_pattern, str):
        return None
    match = re.search(r'(chr([0-9A-Za-z_]+|X|Y|M|MT))', path_pattern)
    if match:
        return match.group(1)
    return None


# ─────────────────────────────────────────────────────────────────────────────
# Function to run bcftools view to get VCF records for a region (remains the same)
# ─────────────────────────────────────────────────────────────────────────────
def run_bcftools_region_extract(bcftools_path, vcf_filepath, chromosome, start_1based, end_1based):
    results = []
    region = f"{chromosome}:{start_1based}-{end_1based}"
    command = [bcftools_path, 'view', '-r', region, vcf_filepath]
    try:
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        if process.stdout:
            for line in process.stdout.splitlines():
                if not line.startswith('#'):
                    results.append(line)
        if process.returncode != 0:
            stderr_output = process.stderr.strip()
            if stderr_output:
                if "failed" in stderr_output.lower() or "error" in stderr_output.lower() or "could not" in stderr_output.lower() or "non-existent" in stderr_output.lower():
                    print(f"Error running 'bcftools view' for region {region}: {stderr_output}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: bcftools command not found at '{bcftools_path}'. Cannot perform VCF operation.", file=sys.stderr)
        return []
    except Exception as e:
        print(f"An unexpected error occurred during 'bcftools view' for region {region}: {e}", file=sys.stderr)
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Main processing function
# ─────────────────────────────────────────────────────────────────────────────
def filter_json_nodes_and_write(json_filepath, idx_filepath, output_json_filepath, vcf_file=None, txt_output_path=None):
    bcftools_path = None
    if vcf_file:
        bcftools_path = shutil.which("bcftools")
        if not bcftools_path:
            print("Warning: --vcf_file provided, but bcftools command not found in PATH. VCF queries will be skipped.",
                  file=sys.stderr)
        else:
            print(f"Using bcftools found at: {bcftools_path}")

    print("Step 1: Loading main JSON structure, node IDs, and creating node map from input JSON file...")
    main_json_structure, target_node_ids_from_json, json_nodes_map = load_json_data_ids_and_map(json_filepath)

    if main_json_structure is None:
        print("Could not load main JSON structure. Cannot proceed.")
        return
    if not json_nodes_map and target_node_ids_from_json:
        print("Warning: Node IDs were extracted from JSON, but node map creation failed. Check JSON integrity.",
              file=sys.stderr)

    chromosome_for_vcf = None
    if vcf_file and bcftools_path:
        path_pattern = main_json_structure.get("path_name_input_pattern")
        chromosome_for_vcf = extract_chromosome_from_path_pattern(path_pattern)
        if not chromosome_for_vcf:
            print(f"Warning: Could not extract chromosome from 'path_name_input_pattern': '{path_pattern}'.",
                  file=sys.stderr)
            print("Attempting to use 'chr1' as a fallback chromosome for VCF queries. This may not be correct for your data.",
                  file=sys.stderr)
            chromosome_for_vcf = "chr1"
        else:
            print(f"Extracted chromosome '{chromosome_for_vcf}' for VCF queries.")

    print("\nStep 2: Loading node index from IDX file...")
    idx_data = load_index(idx_filepath)
    available_idx_node_ids = set(idx_data.keys())
    print(f"Found {len(available_idx_node_ids)} unique node IDs in {idx_filepath}")

    print("\nStep 3: Identifying common node IDs (filtered nodes)...")
    common_node_ids_int = set()
    if target_node_ids_from_json and available_idx_node_ids:
        common_node_ids_int = target_node_ids_from_json.intersection(available_idx_node_ids)

    num_common_nodes = len(common_node_ids_int)
    if num_common_nodes == 0:
        print("No common node IDs found between JSON and IDX. Output 'nodes' list will be empty, and TXT output (if requested) will be empty.")
    else:
        print(f"Found {num_common_nodes} common (filtered) node ID(s). These will be processed for the output JSON.")

    # --- Moved and modified section for writing FILTERED node IDs to TXT file ---
    if txt_output_path:
        if common_node_ids_int:
            print(f"\nStep 3.5: Writing {len(common_node_ids_int)} filtered node IDs to {txt_output_path}...")
            try:
                with open(txt_output_path, 'w') as f_txt:
                    for node_id in sorted(list(common_node_ids_int)): # Sort for consistent output
                        f_txt.write(f"{node_id}\n")
                print(f"Successfully wrote filtered node IDs to {txt_output_path}")
            except IOError as e:
                print(f"Error: Could not write filtered node IDs to TXT file {txt_output_path}. Details: {e}", file=sys.stderr)
            except Exception as e:
                print(f"An unexpected error occurred while writing filtered node IDs to TXT file {txt_output_path}: {e}", file=sys.stderr)
        else:
            print(f"\nStep 3.5: No filtered node IDs to write to {txt_output_path} (file will be empty or not created if it doesn't exist).")
            # Optionally, create an empty file if that's desired behavior for no filtered nodes
            try:
                with open(txt_output_path, 'w') as f_txt: # Creates or truncates the file
                    pass # Ensure an empty file is created if path is given but no common nodes
                print(f"Created empty TXT file at {txt_output_path} as no filtered nodes were found.")
            except IOError as e:
                 print(f"Error: Could not create/truncate TXT file {txt_output_path}. Details: {e}", file=sys.stderr)

    # --- End of moved and modified section ---

    output_json_structure = copy.deepcopy(main_json_structure)
    filtered_nodes_list_in_json = []
    nodes_queried_with_bcftools_count = 0

    if num_common_nodes > 0: # This check is now equivalent to 'if common_node_ids_int:'
        print(f"\nStep 4: Processing {num_common_nodes} common nodes for inclusion in output and potential VCF query...")
        sorted_common_ids = sorted(list(common_node_ids_int))

        for node_id_int in sorted_common_ids:
            original_node_item = json_nodes_map.get(node_id_int)
            if not original_node_item:
                print(f"Warning: Common node ID {node_id_int} not found in JSON node map. Skipping.", file=sys.stderr)
                continue

            current_node_copy = copy.deepcopy(original_node_item)

            if vcf_file and bcftools_path and chromosome_for_vcf:
                try:
                    start_0based_str = current_node_copy.get("grch38_position_start")
                    length_str = current_node_copy.get("length")

                    if start_0based_str is None or length_str is None:
                        raise KeyError("grch38_position_start or length missing")

                    start_0based = int(start_0based_str)
                    length = int(length_str)

                    if length <= 0:
                        current_node_copy["vcf_query_results"] = []
                    else:
                        query_start_1based = start_0based + 1
                        query_end_1based = start_0based + length
                        vcf_results = run_bcftools_region_extract(bcftools_path, vcf_file, chromosome_for_vcf,
                                                                  query_start_1based, query_end_1based)
                        current_node_copy["vcf_query_results"] = vcf_results
                        nodes_queried_with_bcftools_count += 1

                        if nodes_queried_with_bcftools_count > 0 and nodes_queried_with_bcftools_count % 100 == 0:
                            print(
                                f"    ... {nodes_queried_with_bcftools_count} nodes queried with bcftools (out of {num_common_nodes} common nodes).",
                                file=sys.stderr)
                except KeyError as e:
                    print(
                        f"Warning: Node ID {node_id_int} missing required field ('{e}') for VCF query. Skipping VCF query.",
                        file=sys.stderr)
                    current_node_copy["vcf_query_results"] = []
                except ValueError as e:
                    print(
                        f"Warning: Node ID {node_id_int} has invalid numeric value for position/length ({e}). Skipping VCF query.",
                        file=sys.stderr)
                    current_node_copy["vcf_query_results"] = []

            filtered_nodes_list_in_json.append(current_node_copy)

    output_json_structure["nodes"] = filtered_nodes_list_in_json
    print(
        f"\nStep 5: Writing final JSON structure (with {len(filtered_nodes_list_in_json)} nodes in 'nodes' list) to {output_json_filepath}...")

    if vcf_file:
        if bcftools_path:
            if nodes_queried_with_bcftools_count > 0 and nodes_queried_with_bcftools_count % 100 != 0:
                print(f"    ... completed all {nodes_queried_with_bcftools_count} VCF queries for common nodes.",
                      file=sys.stderr)
            elif nodes_queried_with_bcftools_count == 0 and num_common_nodes > 0:
                print(
                    f"    ... no VCF queries were performed (possibly due to missing fields or zero length for all {num_common_nodes} common nodes).",
                    file=sys.stderr)
            print(
                f"Completed: Attempted VCF queries with bcftools for {nodes_queried_with_bcftools_count} filtered and eligible node(s).")

    try:
        with open(output_json_filepath, 'w') as f_out_json:
            json.dump(output_json_structure, f_out_json, indent=4)
        print(f"Successfully wrote the resulting JSON structure to {output_json_filepath}")
    except IOError as e:
        print(f"Error: Could not write to output JSON file {output_json_filepath}. Details: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred while writing the output JSON file: {e}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# Command-line interface
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Filter 'nodes' list in a JSON object based on IDX file and optionally query VCF for each filtered node. Can also output filtered node IDs to a text file.") # Updated description
    parser.add_argument("json_path", help="Path to the input JSON file (object with a 'nodes' list).")
    parser.add_argument("idx_path", help="Path to the .idx file.")
    parser.add_argument("output_json_path", help="Path to the output JSON file for the modified JSON object.")
    parser.add_argument("--vcf_file",
                        help="Optional: Path to a bgzipped and tabix-indexed VCF file to query with bcftools for each filtered node.",
                        default=None)
    parser.add_argument("--txt",
                        dest="txt_output_path",
                        help="Optional: Path to a .txt file to output the filtered node IDs (IDs common to JSON and IDX) (one ID per line).", # Updated help text
                        default=None)

    args = parser.parse_args()

    filter_json_nodes_and_write(args.json_path, args.idx_path, args.output_json_path, args.vcf_file, args.txt_output_path)


if __name__ == "__main__":
    main()