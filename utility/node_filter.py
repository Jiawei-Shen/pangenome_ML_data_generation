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
                    # print( # Commented out for less verbosity
                    #     f"Warning: Truncated block header in {idx_path} at block index {i}. Reached end of file unexpectedly.",
                    #     file=sys.stderr)
                    break
                block_id, block_start, block_size, n_records, metadata_len = struct.unpack("<I Q I I H", header_data)
                if metadata_len > 0:
                    skipped_bytes = f.read(metadata_len)
                    if len(skipped_bytes) < metadata_len:
                        # print( # Commented out for less verbosity
                        #     f"Warning: Truncated metadata in {idx_path} for block_id {block_id}. Reached end of file unexpectedly.",
                        #     file=sys.stderr)
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
                # This warning is important if "nodes" is expected but missing/wrong type
                print(f"Warning: 'nodes' key not found or not a list in {json_filepath}.", file=sys.stderr)
                return main_json_data, set(), {}
            for item_index, item in enumerate(node_list_from_json):
                if not isinstance(item, dict): continue
                try:
                    node_id_str = item.get("node_id")
                    if node_id_str is None:
                        continue  # Silently skip items missing node_id for less verbosity
                    node_id_int = int(node_id_str)
                    node_ids_set.add(node_id_int)
                    json_nodes_map[node_id_int] = item
                except (ValueError, TypeError) as e:
                    # This warning can be important for malformed node_ids
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
    command = [bcftools_path, 'view', '-H', '-r', region, vcf_filepath]
    try:
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        if process.stdout:
            for line in process.stdout.splitlines():
                results.append(line)

        if process.returncode != 0:
            stderr_output = process.stderr.strip()
            if stderr_output:
                error_keywords = ["failed", "error", "could not", "non-existent", "unable to open"]
                is_actual_error = any(keyword in stderr_output.lower() for keyword in error_keywords)
                is_no_match_message = "no lines matching" in stderr_output.lower()
                if is_actual_error or (is_no_match_message and results):
                    # Still print actual errors from bcftools, but not the "no lines matching" if results is empty
                    if not (is_no_match_message and not results):
                        print(f"Error/Warning from 'bcftools view -H' for region {region}: {stderr_output}",
                              file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: bcftools command not found at '{bcftools_path}'. VCF operations will fail.", file=sys.stderr)
        return []
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
    print(f"Loaded {len(json_nodes_map)} nodes from JSON for potential processing.")

    chromosome_for_vcf = None
    if vcf_file and bcftools_path:
        path_pattern = main_json_structure.get("path_name_input_pattern")
        chromosome_for_vcf = extract_chromosome_from_path_pattern(path_pattern)
        if not chromosome_for_vcf:
            print(
                f"Warning: Could not extract chromosome from 'path_name_input_pattern': '{path_pattern}'. Using 'chr1' as fallback.",
                file=sys.stderr)
            chromosome_for_vcf = "chr1"
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
        print(f"Found {num_initial_common_nodes} common node IDs between JSON and IDX for potential processing.")

    ultimate_filtered_nodes_list = []
    total_nodes_queried_with_bcftools = 0
    total_nodes_filtered_out_by_vcf_prereq = 0
    total_nodes_filtered_out_by_empty_vcf = 0
    total_nodes_kept_with_vcf_results = 0

    batch_node_counter_step4 = 0
    batch_queried_for_vcf_step4 = 0
    batch_kept_with_vcf_results_step4 = 0
    batch_start_time_step4 = time.time()

    if num_initial_common_nodes > 0:
        print(
            f"\nStep 4: Processing {num_initial_common_nodes} common nodes for filtering and VCF query (if applicable)...")
        sorted_common_ids = sorted(list(common_node_ids_int))

        for node_idx_step4, node_id_int in enumerate(sorted_common_ids):
            original_node_item = json_nodes_map.get(node_id_int)
            if not original_node_item:
                sys.stderr.write(f"Internal Warning: Common node ID {node_id_int} not in JSON map. Skipping.\n")
                continue

            current_node_copy = copy.deepcopy(original_node_item)
            include_node_in_final_output = True
            node_was_eligible_for_vcf_query = False
            node_had_vcf_results = False

            if vcf_file and bcftools_path and chromosome_for_vcf:
                node_was_eligible_for_vcf_query = True
                try:
                    start_0based_str = current_node_copy.get("grch38_position_start")
                    length_str = current_node_copy.get("length")

                    if start_0based_str is None or length_str is None:
                        # sys.stderr.write(f"Info: Node ID {current_node_copy.get('node_id')} missing coordinate/length. Excluding from VCF query and final output.\n") # REMOVED
                        include_node_in_final_output = False
                        node_was_eligible_for_vcf_query = False
                        total_nodes_filtered_out_by_vcf_prereq += 1
                    else:
                        start_0based = int(start_0based_str)
                        length = int(length_str)
                        if length <= 0:
                            # sys.stderr.write(f"Info: Node ID {current_node_copy.get('node_id')} has non-positive length ({length}). Excluding from VCF query and final output.\n") # REMOVED
                            include_node_in_final_output = False
                            node_was_eligible_for_vcf_query = False
                            total_nodes_filtered_out_by_vcf_prereq += 1
                        else:
                            query_start_1based = start_0based + 1
                            query_end_1based = start_0based + length

                            vcf_results = run_bcftools_region_extract(bcftools_path, vcf_file, chromosome_for_vcf,
                                                                      query_start_1based, query_end_1based)
                            current_node_copy["vcf_query_results"] = vcf_results
                            total_nodes_queried_with_bcftools += 1
                            batch_queried_for_vcf_step4 += 1

                            if not vcf_results:
                                # sys.stderr.write(f"Info: Node ID {current_node_copy.get('node_id')} had no VCF results for region {chromosome_for_vcf}:{query_start_1based}-{query_end_1based}. Excluding from final output.\n") # REMOVED
                                include_node_in_final_output = False
                                total_nodes_filtered_out_by_empty_vcf += 1
                            else:
                                node_had_vcf_results = True
                except ValueError as e:  # Catches int conversion errors
                    # sys.stderr.write(f"Info: Node ID {current_node_copy.get('node_id')} invalid numeric for position/length ('{e}'). Excluding from VCF query and final output.\n") # REMOVED
                    include_node_in_final_output = False
                    node_was_eligible_for_vcf_query = False
                    total_nodes_filtered_out_by_vcf_prereq += 1

            if include_node_in_final_output:
                ultimate_filtered_nodes_list.append(current_node_copy)
                if node_was_eligible_for_vcf_query and node_had_vcf_results:
                    batch_kept_with_vcf_results_step4 += 1
                    total_nodes_kept_with_vcf_results += 1

            batch_node_counter_step4 += 1
            if (batch_node_counter_step4 % 1000 == 0 and batch_node_counter_step4 > 0) or \
                    (node_idx_step4 + 1 == num_initial_common_nodes and batch_node_counter_step4 > 0):
                batch_time_step4 = time.time() - batch_start_time_step4
                print(
                    f"  Step 4 Batch: Processed {batch_node_counter_step4} nodes ({node_idx_step4 + 1}/{num_initial_common_nodes} total common). "
                    f"VCF queried in batch: {batch_queried_for_vcf_step4}. Kept with VCF results in batch: {batch_kept_with_vcf_results_step4}. "
                    f"Time: {batch_time_step4:.2f}s.")
                batch_node_counter_step4 = 0
                batch_queried_for_vcf_step4 = 0
                batch_kept_with_vcf_results_step4 = 0
                batch_start_time_step4 = time.time()

    output_json_structure = copy.deepcopy(main_json_structure)
    output_json_structure["nodes"] = ultimate_filtered_nodes_list
    num_ultimate_nodes = len(ultimate_filtered_nodes_list)
    print(f"\nStep 5: Final processing complete. {num_ultimate_nodes} nodes will be in output JSON.")

    if vcf_file and bcftools_path:
        print(f"VCF Query Summary:")
        print(f"  Total nodes eligible and queried with bcftools: {total_nodes_queried_with_bcftools}")
        print(
            f"  Nodes filtered out due to VCF prerequisites (missing coords/length, etc.): {total_nodes_filtered_out_by_vcf_prereq}")
        print(f"  Nodes filtered out due to empty VCF query results: {total_nodes_filtered_out_by_empty_vcf}")
        print(f"  Nodes kept that had non-empty VCF results: {total_nodes_kept_with_vcf_results}")

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
                        sorted_ids = sorted(ids_to_write)
                    for node_id_str in sorted_ids: f_txt.write(f"{node_id_str}\n")
                print(f"Successfully wrote ultimate filtered node IDs to {txt_output_path}")
            except Exception as e:
                print(f"Error writing ultimate filtered node IDs to TXT {txt_output_path}: {e}", file=sys.stderr)
        else:
            print(f"No ultimate filtered node IDs to write to {txt_output_path}.")
            try:
                with open(txt_output_path, 'w') as f_txt:
                    pass
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
