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
# Function for reading JSON data and extracting Node IDs
# ─────────────────────────────────────────────────────────────────────────────
def load_json_structure_and_ids(json_filepath):
    main_json_data = {}
    node_ids_set = set()
    try:
        with open(json_filepath, 'r') as f_json:
            main_json_data = json.load(f_json)
            if not isinstance(main_json_data, dict):
                print(f"Error: Root JSON content in {json_filepath} is not an object/dictionary.", file=sys.stderr)
                return None, set()
            node_list_from_json = main_json_data.get("nodes")
            if not isinstance(node_list_from_json, list):
                print(f"Warning: 'nodes' key not found or not a list in {json_filepath}.", file=sys.stderr)
                return main_json_data, set()
            for item_index, item in enumerate(node_list_from_json):
                if not isinstance(item, dict):
                    # This warning can be verbose if many items are not dicts, consider reducing verbosity if needed
                    # print(f"Warning: Item at index {item_index} in 'nodes' list (JSON) is not a dictionary.", file=sys.stderr)
                    continue
                try:
                    node_id_str = item.get("node_id")
                    if node_id_str is None: continue
                    node_ids_set.add(int(node_id_str))
                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Could not convert node_id '{node_id_str}' to int for item {item_index} in 'nodes' list (JSON). Error: {e}.",
                        file=sys.stderr)
            print(
                f"Successfully read JSON from {json_filepath}. Found {len(node_list_from_json if node_list_from_json else [])} items in 'nodes' list, extracted {len(node_ids_set)} unique node IDs.")
            return main_json_data, node_ids_set
    except FileNotFoundError:
        print(f"Error: JSON file not found at {json_filepath}", file=sys.stderr)
        return None, set()
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {json_filepath}.", file=sys.stderr)
        return None, set()
    except Exception as e:
        print(f"An unexpected error occurred while reading JSON file {json_filepath}: {e}", file=sys.stderr)
        return None, set()


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
def filter_json_nodes_and_write(json_filepath, idx_filepath, output_json_filepath, vcf_file=None):
    bcftools_path = None
    if vcf_file:
        bcftools_path = shutil.which("bcftools")
        if not bcftools_path:
            print("Warning: --vcf_file provided, but bcftools command not found in PATH. VCF queries will be skipped.",
                  file=sys.stderr)
        else:
            print(f"Using bcftools found at: {bcftools_path}")

    print("Step 1: Loading main JSON structure and node IDs from input JSON file...")
    main_json_structure, target_node_ids_from_json = load_json_structure_and_ids(json_filepath)

    if main_json_structure is None:
        print("Could not load main JSON structure. Cannot proceed.")
        return

    chromosome_for_vcf = None
    if vcf_file and bcftools_path:
        path_pattern = main_json_structure.get("path_name_input_pattern")
        chromosome_for_vcf = extract_chromosome_from_path_pattern(path_pattern)
        if not chromosome_for_vcf:
            print(f"Warning: Could not extract chromosome from 'path_name_input_pattern': '{path_pattern}'.",
                  file=sys.stderr)
            print(
                "Attempting to use 'chr1' as a fallback chromosome for VCF queries. This may not be correct for your data.",
                file=sys.stderr)
            chromosome_for_vcf = "chr1"
        else:
            print(f"Extracted chromosome '{chromosome_for_vcf}' for VCF queries.")

    print("\nStep 2: Loading node index from IDX file...")
    idx_data = load_index(idx_filepath)
    available_idx_node_ids = set(idx_data.keys())
    print(f"Found {len(available_idx_node_ids)} unique node IDs in {idx_filepath}")

    print("\nStep 3: Identifying common node IDs for filtering the 'nodes' list...")
    common_node_ids_int = set()
    if target_node_ids_from_json and available_idx_node_ids:
        common_node_ids_int = target_node_ids_from_json.intersection(available_idx_node_ids)

    num_common_nodes = len(common_node_ids_int)
    if num_common_nodes == 0 and target_node_ids_from_json:
        print("No common node IDs found to filter the 'nodes' list.")
    elif target_node_ids_from_json:
        print(f"Found {num_common_nodes} common node ID(s) for filtering the 'nodes' list.")

    output_json_structure = copy.deepcopy(main_json_structure)
    filtered_nodes_list_in_json = []
    nodes_queried_with_bcftools_count = 0  # Counter for VCF queries

    original_nodes_list = main_json_structure.get("nodes")
    if isinstance(original_nodes_list, list):
        print(f"Processing {len(original_nodes_list)} nodes from input JSON for filtering and potential VCF query...")
        for node_item_idx, node_item in enumerate(original_nodes_list):
            current_node_copy = copy.deepcopy(node_item)  # Work with a copy
            if isinstance(current_node_copy, dict):
                node_id_str = current_node_copy.get("node_id")
                if node_id_str is not None:
                    try:
                        node_id_int = int(node_id_str)
                        if node_id_int in common_node_ids_int:  # Node passes primary filter
                            if vcf_file and bcftools_path and chromosome_for_vcf:
                                try:
                                    start_0based_str = current_node_copy.get("grch38_position_start")
                                    length_str = current_node_copy.get("length")

                                    if start_0based_str is None or length_str is None:
                                        raise KeyError("grch38_position_start or length missing")

                                    start_0based = int(start_0based_str)
                                    length = int(length_str)

                                    if length <= 0:
                                        # This node won't be queried, no vcf_query_results field added unless explicitly desired
                                        # print(f"Note: Node ID {node_id_str} has non-positive length ({length}). Skipping VCF query.", file=sys.stderr)
                                        current_node_copy[
                                            "vcf_query_results"] = []  # Add empty list for consistency if desired
                                    else:
                                        query_start_1based = start_0based + 1
                                        query_end_1based = start_0based + length
                                        vcf_results = run_bcftools_region_extract(bcftools_path, vcf_file,
                                                                                  chromosome_for_vcf,
                                                                                  query_start_1based, query_end_1based)
                                        current_node_copy["vcf_query_results"] = vcf_results
                                        nodes_queried_with_bcftools_count += 1

                                        # --- PROGRESS UPDATE ---
                                        if nodes_queried_with_bcftools_count > 0 and nodes_queried_with_bcftools_count % 100 == 0:
                                            print(
                                                f"    ... {nodes_queried_with_bcftools_count} nodes queried with bcftools.",
                                                file=sys.stderr)
                                        # --- END PROGRESS UPDATE ---

                                except KeyError as e:
                                    print(
                                        f"Warning: Node ID {node_id_str} missing required field ('{e}') for VCF query. Skipping VCF query for this node.",
                                        file=sys.stderr)
                                    current_node_copy["vcf_query_results"] = []
                                except ValueError as e:
                                    print(
                                        f"Warning: Node ID {node_id_str} has invalid numeric value for position/length ({e}). Skipping VCF query for this node.",
                                        file=sys.stderr)
                                    current_node_copy["vcf_query_results"] = []
                            # else: VCF querying not enabled or prerequisites not met
                            #  if "vcf_query_results" in current_node_copy: # Clean up if it existed
                            #       del current_node_copy["vcf_query_results"]
                            filtered_nodes_list_in_json.append(current_node_copy)
                    except (ValueError, TypeError):
                        pass
        output_json_structure["nodes"] = filtered_nodes_list_in_json
        print(
            f"\nStep 4: Writing modified JSON structure (with {len(filtered_nodes_list_in_json)} nodes in 'nodes' list) to {output_json_filepath}...")
    else:
        print(
            f"\nStep 4: Input JSON does not have a 'nodes' list or it's not a list. Writing structure to {output_json_filepath}...")

    if vcf_file:
        if bcftools_path:
            print(
                f"Completed: Attempted VCF queries with bcftools for {nodes_queried_with_bcftools_count} filtered and eligible node(s).")
        # else: warning already printed if bcftools not found

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
        description="Filter 'nodes' list in a JSON object based on IDX file and optionally query VCF for each filtered node.")
    parser.add_argument("json_path", help="Path to the input JSON file (object with a 'nodes' list).")
    parser.add_argument("idx_path", help="Path to the .idx file.")
    parser.add_argument("output_json_path", help="Path to the output JSON file for the modified JSON object.")
    parser.add_argument("--vcf_file",
                        help="Optional: Path to a bgzipped and tabix-indexed VCF file to query with bcftools for each filtered node.",
                        default=None)

    args = parser.parse_args()

    filter_json_nodes_and_write(args.json_path, args.idx_path, args.output_json_path, args.vcf_file)


if __name__ == "__main__":
    main()