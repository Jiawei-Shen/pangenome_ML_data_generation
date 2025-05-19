import json
import argparse
import sys

# Define the reverse mapping from index to base
# This should match the BASE_TO_INDEX in your pileup generation script
INDEX_TO_BASE = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N'}


# If your pileup script uses '*' for deletions and maps it to index 4 as well,
# this mapping is fine. If '*' has a distinct index, update accordingly.
# For example, if BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 5}
# then INDEX_TO_BASE would be {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N', 5: '*'}


def load_json_data(json_filepath):
    """Loads data from a JSON file."""
    try:
        with open(json_filepath, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        print(f"❌ Error: JSON file not found at {json_filepath}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"❌ Error: Could not decode JSON from {json_filepath}. Ensure it's a valid JSON file.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error loading JSON file {json_filepath}: {e}", file=sys.stderr)
        sys.exit(1)


def print_pileup_matrix(node_id, variant_key, pileup_matrix_indices):
    """
    Prints a single pileup matrix in a human-readable format.
    """
    print(f"\n--- Node ID: {node_id} | Variant: {variant_key} ---")
    if not pileup_matrix_indices:
        print("  (No reads in pileup for this variant)")
        return

    max_reads_to_display = 20  # Limit display for very deep pileups
    displayed_reads = 0

    for i, row_indices in enumerate(pileup_matrix_indices):
        if displayed_reads >= max_reads_to_display:
            print(f"  ... (and {len(pileup_matrix_indices) - max_reads_to_display} more reads)")
            break

        pileup_row_str = "".join([INDEX_TO_BASE.get(idx, '?') for idx in row_indices])  # '?' for unknown index
        print(f"  Read {i + 1:2d}: {pileup_row_str}")
        displayed_reads += 1

    if not pileup_matrix_indices:  # Should be caught by the first check, but as a safeguard
        print("  Pileup matrix is empty.")


def main():
    parser = argparse.ArgumentParser(
        description="View variant-centered pileups from a JSON output file in a human-readable format.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "pileup_json_file",
        help="Path to the JSON file containing pileup data (output from single_node_pileup_v2.py)."
    )
    parser.add_argument(
        "--max_reads",
        type=int,
        default=20,
        help="Maximum number of reads to display per pileup matrix to avoid overly long output."
    )

    args = parser.parse_args()

    # Update global max_reads_to_display if needed (though it's passed to function now)
    # global max_reads_to_display # Not needed if passed as arg
    # max_reads_to_display = args.max_reads

    pileup_data = load_json_data(args.pileup_json_file)

    if not pileup_data:
        print("No data loaded from JSON file.", file=sys.stderr)
        return

    if not isinstance(pileup_data, dict):
        print(f"❌ Error: Expected JSON data to be a dictionary (with node ID as key), but got {type(pileup_data)}.",
              file=sys.stderr)
        return

    # Since the input JSON is for a single node, it should have one top-level key
    if len(pileup_data) == 0:
        print("ℹ️ The JSON file is empty or does not contain any node data.", file=sys.stderr)
        return

    if len(pileup_data) > 1:
        print(
            "⚠️ Warning: The JSON file contains data for more than one node. This script will display pileups for all of them.",
            file=sys.stderr)

    for node_id_str, variants_dict in pileup_data.items():
        if not isinstance(variants_dict, dict):
            print(
                f"⚠️ Warning: Data for node '{node_id_str}' is not in the expected format (should be a dictionary of variants). Skipping.",
                file=sys.stderr)
            continue

        if not variants_dict:
            print(f"ℹ️ No variants found or pileups generated for node ID: {node_id_str}")
            continue

        print(f"\n=== Displaying Pileups for Node ID: {node_id_str} ===")
        for variant_key, pileup_matrix_indices in variants_dict.items():
            # --- Modification for max_reads display ---
            # print_pileup_matrix(node_id_str, variant_key, pileup_matrix_indices)
            print(f"\n--- Variant: {variant_key} ---")
            if not pileup_matrix_indices:
                print("  (No reads in pileup for this variant)")
                continue

            displayed_reads_count = 0
            for i, row_indices in enumerate(pileup_matrix_indices):
                if displayed_reads_count >= args.max_reads:
                    print(f"  ... (and {len(pileup_matrix_indices) - args.max_reads} more reads)")
                    break

                pileup_row_str = "".join([INDEX_TO_BASE.get(idx, '?') for idx in row_indices])
                print(f"  Read {i + 1:3d}: {pileup_row_str}")  # Adjusted padding for read number
                displayed_reads_count += 1

    print("\n✅ Done viewing pileups.")


if __name__ == "__main__":
    main()
