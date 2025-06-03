#!/usr/bin/env python3
import json
import struct
import argparse
import sys
import copy
import subprocess
import shutil  # For checking bcftools availability
import re  # For extracting chromosome name


# ─────────────────────────────────────────────────────────────────────────────
# Functions for reading IDX file
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
                        print(f"Warning: Item at index {item_index} in 'nodes' (JSON) missing 'node_id'. Skipping.",
                              file=sys.stderr)
                        continue
                    node_id_int = int(node_id_str)
                    node_ids_set.add(node_id_int)
                    json_nodes_map[node_id_int] = item
                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Node ID format error for item {item_index} in 'nodes' (JSON): '{node_id_str}'. Error: {e}. Skipping.",
                        file=sys.stderr)
            print(
                f"Read JSON {json_filepath}: {len(node_list_from_json)} items in 'nodes' list, {len(node_ids_set)} unique IDs, {len(json_nodes_map)} mapped nodes.")
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
# Function to extract chromosome from path pattern
# ─────────────────────────────────────────────────────────────────────────────
def extract_chromosome_from_path_pattern(path_pattern):
    if not path_pattern or not isinstance(path_pattern, str):
        return None
    match = re.search(r'(chr([0-9A-Za-z_]+|X|Y|M|MT))', path_pattern)
    if match:
        return match.group(1)
    return None


# ─────────────────────────────────────────────────────────────────────────────
# Function to run bcftools view to get VCF records for a region
# ─────────────────────────────────────────────────────────────────────────────
def run_bcftools_region_extract(bcftools_path, vcf_filepath, chromosome, start_1based, end_1based):
    results = []
    region = f"{chromosome}:{start_1based}-{end_1based}"
    # MODIFIED: Added -H to suppress VCF header lines
    command = [bcftools_path, 'view', '-H', '-r', region, vcf_filepath]
    try:
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        if process.stdout:
            # MODIFIED: No need to check for '#' since -H is used
            for line in process.stdout.splitlines():
                results.append(line)

        if process.returncode != 0:
            stderr_output = process.stderr.strip()
            if stderr_output:
                error_keywords = ["failed", "error", "could not", "non-existent", "unable to open", "no lines"]
                if any(keyword in stderr_output.lower() for keyword in error_keywords):
                    # "no lines matching" is not necessarily an error if the region is empty,
                    # but other errors should be reported.
                    if not (
                            "no lines matching" in stderr_output.lower() and not results):  # Don't print error for empty region if no actual error message
                        print(f"Error/Warning running 'bcftools view -H' for region {region}: {stderr_output}",
                              file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: bcftools command not found at '{bcftools_path}'. VCF operations will fail.", file=sys.stderr)
        return []  # Return empty list as bcftools could not be run
    except Exception as e:
        print(f"Unexpected error during 'bcftools view -H' for region {region}: {e}", file=sys.stderr)
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Main processing function
# ─────────────────────────────────────────────────────────────────────────────
def filter_json_nodes_and_write(json_filepath, idx_filepath, output_json_filepath, vcf_file=None, txt_output_path=None):
    bcftools_path = None
    if vcf_file:
        bcftools_path = shutil.which("bcftools")
        if not bcftools_path:
            print("Warning: --vcf_file provided, but bcftools not found. VCF queries/filtering disabled.",
                  file=sys.stderr)
        else:
            print(f"Using bcftools at: {bcftools_path}")

    print("Step 1: Loading JSON data and node map...")
    main_json_structure, target_node_ids_from_json, json_nodes_map = load_json_data_ids_and_map(json_filepath)
    if main_json_structure is None:
        print("Critical error: Could not load main JSON. Cannot proceed.");
        return

    chromosome_for_vcf = None
    if vcf_file and bcftools_path:
        path_pattern = main_json_structure.get("path_name_input_pattern")
        chromosome_for_vcf = extract_chromosome_from_path_pattern(path_pattern)
        if not chromosome_for_vcf:
            print(
                f"Warning: Could not extract chromosome from 'path_name_input_pattern': '{path_pattern}'. Using 'chr1' as fallback.",
                file=sys.stderr)
            chromosome_for_vcf = "chr1"  # Fallback, user should be aware
        else:
            print(f"Extracted chromosome '{chromosome_for_vcf}' for VCF queries.")

    print("\nStep 2: Loading node index (IDX)...")
    idx_data = load_index(idx_filepath)
    available_idx_node_ids = set(idx_data.keys())
    print(f"Found {len(available_idx_node_ids)} unique node IDs in {idx_filepath}")

    print("\nStep 3: Identifying common node IDs...")
    common_node_ids_int = target_node_ids_from_json.intersection(available_idx_node_ids)
    num_initial_common_nodes = len(common_node_ids_int)
    if num_initial_common_nodes == 0:
        print("No common node IDs found between JSON and IDX. Output 'nodes' list will be empty.")
    else:
        print(f"Found {num_initial_common_nodes} common node IDs between JSON and IDX.")

    ultimate_filtered_nodes_list = []
    nodes_actually_queried_with_bcftools = 0
    nodes_filtered_out_by_vcf_content = 0  # New counter

    if num_initial_common_nodes > 0:
        print(
            f"\nStep 4: Processing {num_initial_common_nodes} common nodes for filtering and VCF query (if applicable)...")
        sorted_common_ids = sorted(list(common_node_ids_int))

        for node_id_int in sorted_common_ids:
            original_node_item = json_nodes_map.get(node_id_int)
            if not original_node_item:
                print(f"Warning: Common node ID {node_id_int} not in JSON map. Skipping.", file=sys.stderr)
                continue

            current_node_copy = copy.deepcopy(original_node_item)
            include_node_in_final_output = True

            if vcf_file and bcftools_path and chromosome_for_vcf:
                try:
                    start_0based_str = current_node_copy.get("grch38_position_start")
                    length_str = current_node_copy.get("length")

                    if start_0based_str is None or length_str is None:
                        print(
                            f"Info: Node ID {current_node_copy.get('node_id')} missing coordinate/length. Excluding from VCF query and final output.",
                            file=sys.stderr)
                        include_node_in_final_output = False
                    else:
                        start_0based = int(start_0based_str)
                        length = int(length_str)
                        if length <= 0:
                            print(
                                f"Info: Node ID {current_node_copy.get('node_id')} has non-positive length ({length}). Excluding from VCF query and final output.",
                                file=sys.stderr)
                            include_node_in_final_output = False
                        else:
                            query_start_1based = start_0based + 1
                            query_end_1based = start_0based + length
                            vcf_results = run_bcftools_region_extract(bcftools_path, vcf_file, chromosome_for_vcf,
                                                                      query_start_1based, query_end_1based)
                            current_node_copy["vcf_query_results"] = vcf_results
                            nodes_actually_queried_with_bcftools += 1

                            # MODIFIED: Filter if vcf_results is empty after successful query
                            if not vcf_results:
                                print(
                                    f"Info: Node ID {current_node_copy.get('node_id')} had no VCF results for region {chromosome_for_vcf}:{query_start_1based}-{query_end_1based}. Excluding from final output.",
                                    file=sys.stderr)
                                include_node_in_final_output = False
                                nodes_filtered_out_by_vcf_content += 1

                            if nodes_actually_queried_with_bcftools > 0 and nodes_actually_queried_with_bcftools % 100 == 0:
                                print(f"    ... {nodes_actually_queried_with_bcftools} nodes queried with bcftools.",
                                      file=sys.stderr)
                except ValueError as e:
                    print(
                        f"Info: Node ID {current_node_copy.get('node_id')} invalid numeric for position/length ('{e}'). Excluding from VCF query and final output.",
                        file=sys.stderr)
                    include_node_in_final_output = False

            if include_node_in_final_output:
                ultimate_filtered_nodes_list.append(current_node_copy)

    output_json_structure = copy.deepcopy(main_json_structure)
    output_json_structure["nodes"] = ultimate_filtered_nodes_list
    num_ultimate_nodes = len(ultimate_filtered_nodes_list)
    print(f"\nStep 5: Final processing complete. {num_ultimate_nodes} nodes will be in output JSON.")

    if vcf_file and bcftools_path:  # Provide summary of VCF operations
        if nodes_actually_queried_with_bcftools > 0 and nodes_actually_queried_with_bcftools % 100 != 0:  # Final count if not multiple of 100
            print(f"    ... completed all {nodes_actually_queried_with_bcftools} VCF queries.", file=sys.stderr)
        elif nodes_actually_queried_with_bcftools == 0 and num_initial_common_nodes > 0:
            print(f"    ... no VCF queries were performed (e.g., all common nodes failed VCF prerequisite checks).",
                  file=sys.stderr)
        print(f"Completed: VCF queries with bcftools for {nodes_actually_queried_with_bcftools} eligible node(s).")
        if nodes_filtered_out_by_vcf_content > 0:
            print(f"Filtered out {nodes_filtered_out_by_vcf_content} node(s) due to empty VCF query results.")

    if txt_output_path:
        print(f"\nStep 6: Writing {num_ultimate_nodes} ultimate filtered node ID(s) to {txt_output_path}...")
        if ultimate_filtered_nodes_list:
            try:
                with open(txt_output_path, 'w') as f_txt:
                    ids_to_write = [str(node.get("node_id")) for node in ultimate_filtered_nodes_list if
                                    node.get("node_id") is not None]
                    try:
                        sorted_ids = sorted(ids_to_write, key=int)
                    except ValueError:
                        print("Warning: Node IDs in TXT not all numeric, using lexicographical sort.", file=sys.stderr)
                        sorted_ids = sorted(ids_to_write)
                    for node_id_str in sorted_ids: f_txt.write(f"{node_id_str}\n")
                print(f"Successfully wrote ultimate filtered node IDs to {txt_output_path}")
            except Exception as e:
                print(f"Error writing ultimate filtered node IDs to TXT {txt_output_path}: {e}", file=sys.stderr)
        else:
            print(f"No ultimate filtered node IDs to write to {txt_output_path}.")
            try:
                with open(txt_output_path, 'w') as f_txt:
                    pass  # Create empty file
                print(f"Created empty TXT file at {txt_output_path}.")
            except Exception as e:
                print(f"Error creating empty TXT file {txt_output_path}: {e}", file=sys.stderr)

    print(f"\nStep 7: Writing final JSON (with {num_ultimate_nodes} nodes) to {output_json_filepath}...")
    try:
        with open(output_json_filepath, 'w') as f_out_json:
            json.dump(output_json_structure, f_out_json, indent=4)
        print(f"Successfully wrote resulting JSON to {output_json_filepath}")
    except Exception as e:
        print(f"Error writing output JSON {output_json_filepath}: {e}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# Command-line interface
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Filter 'nodes' list in a JSON based on IDX file and VCF content. Outputs modified JSON and optionally a TXT list of filtered node IDs.")
    parser.add_argument("json_path", help="Path to the input JSON file.")
    parser.add_argument("idx_path", help="Path to the .idx file.")
    parser.add_argument("output_json_path", help="Path for the output filtered JSON file.")
    parser.add_argument("--vcf_file",
                        help="Optional: Path to a bgzipped and tabix-indexed VCF file. Nodes will be filtered if they lack VCF prerequisites OR if their region yields no variants from this VCF.",
                        default=None)
    parser.add_argument("--txt", dest="txt_output_path",
                        help="Optional: Path to output filtered node IDs (one per line).", default=None)
    args = parser.parse_args()
    filter_json_nodes_and_write(args.json_path, args.idx_path, args.output_json_path, args.vcf_file,
                                args.txt_output_path)


if __name__ == "__main__":
    main()