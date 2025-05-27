import json
import struct
import argparse
import sys
import copy  # For deep copying the main JSON structure


# ─────────────────────────────────────────────────────────────────────────────
# Functions for reading IDX file (remains the same)
# ─────────────────────────────────────────────────────────────────────────────
def load_index(idx_path):
    """
    Loads the node index from the .idx file.

    Args:
        idx_path (str): Path to the .idx file.

    Returns:
        dict: A dictionary where keys are integer node IDs (block_id)
              and values are dictionaries with "start", "size", "n_records".
              Returns an empty dictionary if the file cannot be processed.
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

                node_index[block_id] = {
                    "start": block_start,
                    "size": block_size,
                    "n_records": n_records
                }
    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except struct.error as e:
        print(
            f"Error: Could not unpack data from IDX file {idx_path}. It might be corrupted or not match the expected format. Details: {e}",
            file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading the IDX file {idx_path}: {e}", file=sys.stderr)
        return {}
    return node_index


# ─────────────────────────────────────────────────────────────────────────────
# Function for reading JSON data and extracting Node IDs from the "nodes" list
# ─────────────────────────────────────────────────────────────────────────────
def load_json_structure_and_ids(json_filepath):
    """
    Reads the main JSON object from a file. Expects a "nodes" key
    containing a list of node objects, from which IDs are extracted.

    Args:
        json_filepath (str): Path to the JSON file.

    Returns:
        tuple: (main_json_object, set_of_integer_node_ids_from_nodes_list)
               Returns (None, set()) if a critical error occurs or if the
               "nodes" list is not found or is not a list.
    """
    main_json_data = {}
    node_ids_set = set()
    try:
        with open(json_filepath, 'r') as f_json:
            main_json_data = json.load(f_json)

            if not isinstance(main_json_data, dict):
                print(f"Error: Root JSON content in {json_filepath} is not an object/dictionary as expected.",
                      file=sys.stderr)
                return None, set()

            node_list_from_json = main_json_data.get("nodes")
            if not isinstance(node_list_from_json, list):
                print(f"Error: The key 'nodes' was not found in {json_filepath}, or its value is not a list.",
                      file=sys.stderr)
                # Return the main data structure but an empty set for IDs, as filtering target is missing
                return main_json_data, set()

            for item_index, item in enumerate(node_list_from_json):
                if not isinstance(item, dict):
                    print(
                        f"Warning: Item at index {item_index} within the 'nodes' list in {json_filepath} is not a dictionary, skipping ID extraction.",
                        file=sys.stderr)
                    continue
                try:
                    node_id_str = item.get("node_id")
                    if node_id_str is None:
                        # print(f"Debug: Item at index {item_index} in 'nodes' list (JSON) does not have 'node_id' key.")
                        continue
                    node_ids_set.add(int(node_id_str))
                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Could not convert node_id '{node_id_str}' to an integer for item at index {item_index} in 'nodes' list (JSON). Error: {e}. Skipping.",
                        file=sys.stderr)
                    continue

        print(
            f"Successfully read JSON from {json_filepath}. Found {len(node_list_from_json if node_list_from_json else [])} items in 'nodes' list, extracted {len(node_ids_set)} unique node IDs.")
        return main_json_data, node_ids_set

    except FileNotFoundError:
        print(f"Error: JSON file not found at {json_filepath}", file=sys.stderr)
        return None, set()
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {json_filepath}. Check for syntax errors.", file=sys.stderr)
        return None, set()
    except Exception as e:
        print(f"An unexpected error occurred while reading the JSON file {json_filepath}: {e}", file=sys.stderr)
        return None, set()


# ─────────────────────────────────────────────────────────────────────────────
# Main processing function
# ─────────────────────────────────────────────────────────────────────────────
def filter_json_nodes_and_write(json_filepath, idx_filepath, output_json_filepath):
    """
    Reads a main JSON object, filters its "nodes" list based on IDs from IDX,
    and writes the modified main JSON object (with filtered "nodes") to an output file.
    """
    print("Step 1: Loading main JSON structure and node IDs from input JSON file...")
    main_json_structure, target_node_ids_from_json = load_json_structure_and_ids(json_filepath)

    if main_json_structure is None:  # Critical error during JSON loading
        print("Could not load main JSON structure. Cannot proceed.")
        return
    if not target_node_ids_from_json and isinstance(main_json_structure.get("nodes"), list) and main_json_structure.get(
            "nodes"):
        print(
            "Warning: No parsable 'node_id' fields were found or converted in the 'nodes' list of the JSON, though a 'nodes' list exists.")
    elif not main_json_structure.get("nodes") and not isinstance(main_json_structure.get("nodes"), list):
        print("Warning: 'nodes' key is missing or not a list in the input JSON. Output will reflect this.")

    print("\nStep 2: Loading node index from IDX file...")
    idx_data = load_index(idx_filepath)
    if not idx_data:
        print("No node data loaded from IDX file. Cannot proceed with filtering 'nodes' list.")
        # If idx fails, we could still write the original json, but the goal is filtering.
        # For now, we'll proceed and the intersection will be empty if idx_data is empty.

    available_idx_node_ids = set(idx_data.keys())
    print(f"Found {len(available_idx_node_ids)} unique node IDs in {idx_filepath}")

    print("\nStep 3: Identifying common node IDs for filtering the 'nodes' list...")
    common_node_ids_int = set()
    if target_node_ids_from_json and available_idx_node_ids:
        common_node_ids_int = target_node_ids_from_json.intersection(available_idx_node_ids)

    num_common_nodes = len(common_node_ids_int)
    if num_common_nodes == 0 and target_node_ids_from_json:  # Only print if there were targets
        print("No common node IDs found to filter the 'nodes' list.")
    elif target_node_ids_from_json:
        print(f"Found {num_common_nodes} common node ID(s) for filtering the 'nodes' list.")

    # Prepare the output structure
    # We use deepcopy to ensure that if 'nodes' is not present or not a list in input,
    # we don't accidentally add it if no filtering happens.
    output_json_structure = copy.deepcopy(main_json_structure)

    # Filter the "nodes" list if it exists and is a list
    original_nodes_list = main_json_structure.get("nodes")
    if isinstance(original_nodes_list, list):
        filtered_nodes_list_in_json = []
        if num_common_nodes > 0:  # Only filter if there are common IDs
            for node_item in original_nodes_list:
                if isinstance(node_item, dict):
                    node_id_str = node_item.get("node_id")
                    if node_id_str is not None:
                        try:
                            node_id_int = int(node_id_str)
                            if node_id_int in common_node_ids_int:
                                filtered_nodes_list_in_json.append(node_item)
                        except (ValueError, TypeError):
                            pass
                            # Replace the "nodes" list in the output structure
        output_json_structure["nodes"] = filtered_nodes_list_in_json
        print(
            f"\nStep 4: Writing modified JSON structure (with {len(filtered_nodes_list_in_json)} nodes in 'nodes' list) to {output_json_filepath}...")
    else:
        print(
            f"\nStep 4: Input JSON does not have a 'nodes' list or it's not a list. Writing the original structure (or structure without 'nodes') to {output_json_filepath}...")
        if "nodes" in output_json_structure and not isinstance(output_json_structure["nodes"], list):
            print(
                f"Warning: 'nodes' key exists in input but is not a list. Its original value will be preserved: {type(output_json_structure['nodes'])}")
        elif "nodes" not in output_json_structure:
            print("Warning: 'nodes' key was not found in the input JSON.")

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
        description="Filter the 'nodes' list within an input JSON object based on node IDs present in an IDX file, and output the modified JSON object.")
    parser.add_argument("json_path", help="Path to the input JSON file (expects an object with a 'nodes' list).")
    parser.add_argument("idx_path", help="Path to the .idx file.")
    parser.add_argument("output_json_path", help="Path to the output JSON file for the modified JSON object.")

    args = parser.parse_args()

    filter_json_nodes_and_write(args.json_path, args.idx_path, args.output_json_path)


if __name__ == "__main__":
    main()