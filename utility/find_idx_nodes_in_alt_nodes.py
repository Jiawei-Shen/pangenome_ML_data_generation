import argparse
import struct
import sys


def load_index(idx_path):
    """
    Loads node IDs from a binary .idx file.

    The structure of the .idx file is assumed to be:
    - 4 bytes: Number of blocks (unsigned int)
    - For each block:
        - 4 bytes: block_id (unsigned int)
        - 8 bytes: block_start (unsigned long long)
        - 4 bytes: block_size (unsigned int)
        - 4 bytes: n_records (unsigned int)
        - 2 bytes: metadata_len (unsigned short)
        - metadata_len bytes: metadata (to be skipped)

    Args:
        idx_path (str): The path to the .idx file.

    Returns:
        dict: A dictionary mapping node IDs (block_id) to their info.
              Returns an empty dictionary if an error occurs.
    """
    node_index = {}
    try:
        with open(idx_path, "rb") as f:
            # Read the total number of blocks
            blocks_num_bytes = f.read(4)
            if not blocks_num_bytes or len(blocks_num_bytes) < 4:
                print(f"Error: Could not read blocks_num from {idx_path}. File might be empty or corrupted.",
                      file=sys.stderr)
                return {}

            blocks_num, = struct.unpack("<I", blocks_num_bytes)

            # Read each block's header
            for i in range(blocks_num):
                header_data = f.read(4 + 8 + 4 + 4 + 2)
                if len(header_data) < (4 + 8 + 4 + 4 + 2):
                    print(f"Warning: Reached end of file unexpectedly while reading block {i + 1}/{blocks_num}.",
                          file=sys.stderr)
                    break

                block_id, block_start, block_size, n_records, metadata_len = struct.unpack("<I Q I I H", header_data)

                # Skip metadata if it exists
                if metadata_len > 0:
                    skipped_bytes = f.read(metadata_len)
                    if len(skipped_bytes) < metadata_len:
                        print(f"Warning: Reached end of file unexpectedly while skipping metadata for block {i + 1}.",
                              file=sys.stderr)
                        break

                node_index[block_id] = {"start": block_start, "size": block_size, "n_records": n_records}
    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except struct.error as e:
        print(f"Error: Could not unpack data from IDX file {idx_path}. The file may be corrupted. Details: {e}",
              file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return {}

    return node_index


def count_matches(tsv_path, idx_nodes):
    """
    Counts how many alt_nodes from a TSV file are present in a given set of nodes.

    Args:
        tsv_path (str): The path to the input TSV file.
        idx_nodes (set): A set of node IDs loaded from the index file.

    Returns:
        int: The total number of matching nodes found. Returns -1 on file error.
    """
    match_count = 0
    line_num = 0
    try:
        with open(tsv_path, 'r') as f:
            for line in f:
                line_num += 1
                # Skip header or empty lines
                if not line.strip() or line.startswith("chr"):
                    continue

                # Split by any whitespace to handle both tabs and spaces
                parts = line.split()

                # The alt_node is the last element. Check for sufficient columns.
                # Based on the example, we need at least 7 columns.
                if len(parts) < 7:
                    print(f"Warning: Skipping malformed line {line_num}: not enough columns.", file=sys.stderr)
                    continue

                try:
                    # The alt_node is the last element in the row
                    print(parts, parts[-1])
                    alt_node = int(parts[-1])
                    if alt_node in idx_nodes:
                        match_count += 1
                except ValueError:
                    print(f"Warning: Skipping line {line_num} due to non-integer alt_node: '{parts[-1]}'",
                          file=sys.stderr)
                    continue
    except FileNotFoundError:
        print(f"Error: TSV file not found at {tsv_path}", file=sys.stderr)
        return -1  # Return an error code
    except Exception as e:
        print(f"An unexpected error occurred while reading TSV file {tsv_path}: {e}", file=sys.stderr)
        return -1

    return match_count


def main():
    """
    Main function to orchestrate the node counting process.
    """
    parser = argparse.ArgumentParser(
        description="Count how many 'alt_node' entries (last column) from a TSV file exist in a binary .idx file."
    )
    parser.add_argument("tsv_file", help="Path to the input TSV file.")
    parser.add_argument("idx_file", help="Path to the .idx index file.")
    args = parser.parse_args()

    # 1. Load node IDs from the index file
    print(f"Loading nodes from index file: {args.idx_file}...", file=sys.stderr)
    index_data = load_index(args.idx_file)
    if not index_data:
        print("Could not load index file. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Get the keys (which are the node IDs) and convert to a set for O(1) lookups
    idx_nodes = set(index_data.keys())
    print(f"Loaded {len(idx_nodes)} unique nodes from the index.", file=sys.stderr)

    # 2. Count the matches in the TSV file
    print(f"Comparing with alt_nodes from TSV file: {args.tsv_file}...", file=sys.stderr)
    total_matches = count_matches(args.tsv_file, idx_nodes)

    # 3. Output the final result
    if total_matches != -1:
        print("\n--- Results ---")
        print(f"Total matching alternate nodes found: {total_matches}")


if __name__ == "__main__":
    main()
