import json
import struct
import argparse
import sys


# ─────────────────────────────────────────────────────────────────────────────
# Functions for reading IDX file (adapted from your provided script)
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
                # Read header for each block
                header_data = f.read(4 + 8 + 4 + 4 + 2)  # block_id, block_start, block_size, n_records, metadata_len
                if len(header_data) < (4 + 8 + 4 + 4 + 2):
                    print(
                        f"Warning: Truncated block header in {idx_path} at block index {i}. Reached end of file unexpectedly.",
                        file=sys.stderr)
                    break
                block_id, block_start, block_size, n_records, metadata_len = struct.unpack("<I Q I I H", header_data)

                if metadata_len > 0:
                    skipped_bytes = f.read(metadata_len)  # skip metadata
                    if len(skipped_bytes) < metadata_len:
                        print(
                            f"Warning: Truncated metadata in {idx_path} for block_id {block_id}. Reached end of file unexpectedly.",
                            file=sys.stderr)
                        break  # Stop processing further blocks as file position might be unreliable

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
# Function for reading JSON data and extracting Node IDs
# ─────────────────────────────────────────────────────────────────────────────
def load_json_data_and_ids(json_filepath):
    """
    Reads the full list of JSON objects from a file and also extracts
    a set of integer node IDs from this data.

    Args:
        json_filepath (str): Path to the JSON file.

    Returns:
        tuple: (list_of_original_json_objects, set_of_integer_node_ids)
               Returns (None, set()) if an error occurs.
    """
    original_json_list = []
    node_ids_set = set()
    try:
        with open(json_filepath, 'r') as f_json:
            original_json_list = json.load(f_json)

            if not isinstance(original_json_list, list):
                print(f"Error: JSON content in {json_filepath} is not a list as expected.", file=sys.stderr)
                return None, set()

            for item_index, item in enumerate(original_json_list):
                if not isinstance(item, dict):
                    print(
                        f"Warning: Item at index {item_index} in {json_filepath} is not a dictionary, skipping ID extraction for this item.",
                        file=sys.stderr)
                    continue
                try:
                    node_id_str = item.get("node_id")
                    if node_id_str is None:
                        # This item might not have a node_id, which is fine if the goal is to preserve it if other criteria match
                        # However, for ID-based filtering, we need it. We'll log if it's missing for ID extraction.
                        # print(f"Debug: Item at index {item_index} in {json_filepath} does not have a 'node_id' key for ID set creation.")
                        continue  # Skip adding to node_ids_set if no node_id
                    node_ids_set.add(int(node_id_str))
                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Could not convert node_id '{node_id_str}' to an integer for item at index {item_index} in {json_filepath}. Error: {e}. Skipping ID extraction for this item.",
                        file=sys.stderr)
                    continue

        print(
            f"Successfully read {len(original_json_list)} items and extracted {len(node_ids_set)} unique node IDs from {json_filepath}")
        return original_json_list, node_ids_set

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
def filter_and_write_nodes_to_json(json_filepath, idx_filepath, output_json_filepath):
    """
    Reads JSON data, gets all node IDs from IDX, filters for common nodes,
    and writes the filtered JSON items to an output JSON file.
    """
    print("Step 1: Loading data and node IDs from input JSON file...")
    original_json_data, target_node_ids_from_json = load_json_data_and_ids(json_filepath)

    if original_json_data is None or not target_node_ids_from_json:  # Check if target_node_ids is empty as well
        print("No data or node IDs loaded from input JSON. Cannot proceed.")
        if not target_node_ids_from_json and original_json_data is not None:
            print("Specifically, no parsable 'node_id' fields were found or converted in the JSON.")
        return

    print("\nStep 2: Loading node index from IDX file...")
    idx_data = load_index(idx_filepath)
    if not idx_data:
        print("No node data loaded from IDX file. Cannot proceed.")
        return

    available_idx_node_ids = set(idx_data.keys())
    print(f"Found {len(available_idx_node_ids)} unique node IDs in {idx_filepath}")

    print("\nStep 3: Identifying common node IDs...")
    common_node_ids_int = target_node_ids_from_json.intersection(available_idx_node_ids)

    num_common_nodes = len(common_node_ids_int)
    if num_common_nodes == 0:
        print("No common node IDs found between the JSON file and the IDX file.")
    else:
        print(f"Found {num_common_nodes} common node ID(s).")

    print(f"\nStep 4: Filtering original JSON data and writing to {output_json_filepath}...")
    filtered_json_output = []
    if num_common_nodes > 0:
        for item in original_json_data:
            if isinstance(item, dict):
                node_id_str = item.get("node_id")
                if node_id_str is not None:
                    try:
                        node_id_int = int(node_id_str)
                        if node_id_int in common_node_ids_int:
                            filtered_json_output.append(item)
                    except (ValueError, TypeError):
                        # Already warned during ID extraction, pass here
                        pass

    try:
        with open(output_json_filepath, 'w') as f_out_json:
            json.dump(filtered_json_output, f_out_json, indent=4)  # indent for readability
        print(f"Successfully wrote {len(filtered_json_output)} filtered items to {output_json_filepath}")
        if len(filtered_json_output) != num_common_nodes:
            print(
                f"Note: {num_common_nodes} common IDs were found, but {len(filtered_json_output)} items were written. This could happen if some IDs were duplicated in the input JSON or if items lacked 'node_id'.")
    except IOError as e:
        print(f"Error: Could not write to output JSON file {output_json_filepath}. Details: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred while writing the output JSON file: {e}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# Command-line interface
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Filter items from an input JSON file based on node IDs present in an IDX file, and output the filtered items to a new JSON file.")
    parser.add_argument("json_path", help="Path to the input JSON file containing target items.")
    parser.add_argument("idx_path", help="Path to the .idx file.")
    parser.add_argument("output_json_path", help="Path to the output JSON file for filtered items.")

    args = parser.parse_args()

    filter_and_write_nodes_to_json(args.json_path, args.idx_path, args.output_json_path)


if __name__ == "__main__":
    main()