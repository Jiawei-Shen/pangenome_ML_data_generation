#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import shutil
import sys


# --- Helper Functions ---

def query_vcf_chr1(vcf_file_path):
    """
    Queries a VCF file for variants on chr1 using bcftools.
    Assumes 'bcftools' is in the system PATH.
    Returns a set of (position, ref, alt) tuples for chr1 variants.
    """
    print(f"Querying VCF file '{vcf_file_path}' for chr1 variants using 'bcftools' from PATH...")
    chr1_variants = set()
    cmd = ['bcftools', 'view', '-r', 'chr1', vcf_file_path]
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            print(f"Error: bcftools command failed with exit code {process.returncode}", file=sys.stderr)
            if stderr:
                print(f"bcftools stderr:\n{stderr.strip()}", file=sys.stderr)
            return None

        for line in stdout.splitlines():
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    try:
                        pos = int(fields[1])
                        ref = fields[3].upper()
                        alt_alleles = fields[4].upper().split(',')
                        for alt in alt_alleles:
                            chr1_variants.add((pos, ref, alt))
                    except ValueError:
                        print(f"Warning: Could not parse VCF line: {line.strip()}", file=sys.stderr)

        print(f"Found {len(chr1_variants)} unique ALT variants on chr1 in VCF.")
        return chr1_variants
    except FileNotFoundError:
        print(f"Error: 'bcftools' not found in PATH. Please ensure it's installed and PATH is configured correctly.",
              file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error during VCF query: {e}", file=sys.stderr)
        return None


def load_node_positions(json_file_path):
    """
    Loads node position information from the provided JSON file (expected JSON Lines format).
    Returns a dictionary mapping top_level_node_id (str) to grch38_start_position (int).
    Only records entries where grch38_start_position is a valid integer.
    """
    print(f"Loading node positions from '{json_file_path}'...")
    node_positions = {}
    try:
        with open(json_file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                try:
                    entry = json.loads(line)  # Parse each line as a JSON object
                except json.JSONDecodeError as jde:
                    print(f"Warning: Could not decode JSON line {line_num}: {line}. Error: {jde}", file=sys.stderr)
                    continue

                # Process the entry (which should now be a dictionary)
                if entry.get("grch38_start_position") is not None:  # Ensure key exists and is not null
                    try:
                        node_id_val = entry.get("top_level_node_id")
                        if node_id_val is None:
                            print(
                                f"Warning: Skipping entry on line {line_num} with missing 'top_level_node_id': {entry}",
                                file=sys.stderr)
                            continue
                        node_id = str(node_id_val)

                        start_pos_val = entry.get("grch38_start_position")
                        # Handles if start_pos_val is already int or a string representation of int
                        if isinstance(start_pos_val, (str, int, float)):  # Allow int/float and convert
                            start_pos = int(start_pos_val)
                            node_positions[node_id] = start_pos
                        else:
                            print(
                                f"Warning: Skipping entry on line {line_num} with non-numeric start position type ('{type(start_pos_val)}'): {entry}",
                                file=sys.stderr)
                            continue

                    except (ValueError, TypeError) as e:
                        print(
                            f"Warning: Skipping entry on line {line_num} due to invalid/missing node ID or start position: {entry} (Error: {e})",
                            file=sys.stderr)
                # else: # No start position, or it's null
                # print(f"Debug: Entry on line {line_num} skipped, no valid grch38_start_position: {entry}", file=sys.stderr)

        print(f"Loaded start positions for {len(node_positions)} nodes.")
        return node_positions
    except FileNotFoundError:
        print(f"Error: Node position JSON file not found at '{json_file_path}'", file=sys.stderr)
        return None
    except Exception as e:  # Catch other potential errors like read permissions
        print(f"Error loading node positions: {e}", file=sys.stderr)
        return None


def calculate_genomic_position(node_start_pos, variant_pos_on_node):
    """
    Calculates the genomic position.
    Assumes VCF is 1-based and .pth variant position is 0-based relative to node start.
    Node start from JSON is assumed to be 1-based GRCh38 coordinate.
    """
    return node_start_pos + variant_pos_on_node


# --- Main Script Logic ---

def main():
    parser = argparse.ArgumentParser(description="Classify .pth variant files based on VCF concordance on chr1.")
    parser.add_argument("pth_folder_path",
                        help="Path to the folder containing node subdirectories with .pth files and variant_summary.json.")
    parser.add_argument("vcf_file", help="Path to the VCF file to query (e.g., dbSNP for chr1).")
    parser.add_argument("node_pos_json",
                        help="Path to the JSON file with node GRCh38 position information (JSON Lines format expected).")
    parser.add_argument("--output_folder", default="./classification_results",
                        help="Path to the folder where 'true' and 'false' subfolders and summary will be created.")
    args = parser.parse_args()

    print("Starting .pth file classification process...")

    # 1. Create output directories
    true_folder = os.path.join(args.output_folder, "true")
    false_folder = os.path.join(args.output_folder, "false")
    os.makedirs(true_folder, exist_ok=True)
    os.makedirs(false_folder, exist_ok=True)
    print(f"Output folders created/ensured: '{true_folder}', '{false_folder}'")

    # 2. Query VCF for chr1 variants
    vcf_chr1_variants = query_vcf_chr1(args.vcf_file)
    if vcf_chr1_variants is None:
        print("Exiting due to VCF query failure.", file=sys.stderr)
        sys.exit(1)
    if not vcf_chr1_variants:
        print(
            "Warning: No chr1 variants found in the VCF file. All .pth files will likely be classified as 'false' unless they are not on chr1 (and thus also 'false' by this VCF).",
            file=sys.stderr)

    # 3. Load node position information
    # This will load positions for all nodes, not just chr1
    all_node_start_positions = load_node_positions(args.node_pos_json)
    if all_node_start_positions is None:
        print("Exiting due to failure in loading node positions.", file=sys.stderr)
        sys.exit(1)
    if not all_node_start_positions:
        print("Warning: No node start positions loaded. Cannot calculate genomic positions for .pth files.",
              file=sys.stderr)

    # 4. Iterate through .pth files, calculate positions, classify, and copy
    classification_summary = []
    total_pth_files_processed_from_summaries = 0
    pth_files_copied_true = 0
    pth_files_copied_false = 0
    actual_pth_files_found_and_considered = 0

    print(f"Processing .pth files from base folder: {args.pth_folder_path}")
    for node_id_str in os.listdir(args.pth_folder_path):
        node_dir_path = os.path.join(args.pth_folder_path, node_id_str)
        if not os.path.isdir(node_dir_path):
            continue

        # Check if we have this node's start position (it's no longer filtered by chr1 here)
        if node_id_str not in all_node_start_positions:
            # This means the node_id from the directory name was not in the positions JSON,
            # or its entry was invalid.
            # print(f"Info: Skipping node {node_id_str}: Start position unknown based on node position JSON.")
            continue

        node_grch38_start_pos = all_node_start_positions[node_id_str]  # 1-based

        summary_json_path = os.path.join(node_dir_path, "variant_summary.json")
        if not os.path.isfile(summary_json_path):
            # print(f"Warning: variant_summary.json not found in {node_dir_path}. Skipping this node's .pth files.", file=sys.stderr)
            continue

        try:
            with open(summary_json_path, 'r') as sjf:
                node_summary_data = json.load(sjf)
        except Exception as e:
            print(f"Warning: Could not read or parse {summary_json_path}: {e}. Skipping.", file=sys.stderr)
            continue

        for variant_info in node_summary_data.get("variants", []):
            total_pth_files_processed_from_summaries += 1
            pth_filename = variant_info.get("tensor_file")
            variant_key = variant_info.get("variant_key")
            if not pth_filename or not variant_key:
                print(f"Warning: Missing 'tensor_file' or 'variant_key' in {summary_json_path} for an entry. Skipping.",
                      file=sys.stderr)
                continue

            pth_file_full_path = os.path.join(node_dir_path, pth_filename)
            if not os.path.isfile(pth_file_full_path):
                print(
                    f"Warning: .pth file not found: {pth_file_full_path} (referenced in summary). Skipping this entry.",
                    file=sys.stderr)
                continue

            actual_pth_files_found_and_considered += 1

            try:
                parts = variant_key.split('_')
                variant_pos_on_node = int(parts[0])
                variant_ref_on_node = parts[2].upper()
                variant_alt_on_node = parts[3].upper()
            except (IndexError, ValueError) as e:
                print(
                    f"Warning: Could not parse variant_key '{variant_key}' from {summary_json_path}: {e}. Skipping {pth_filename}.",
                    file=sys.stderr)
                continue

            genomic_variant_pos = calculate_genomic_position(node_grch38_start_pos, variant_pos_on_node)

            # Classification is still based on vcf_chr1_variants.
            # Variants not on chr1 (based on their calculated genomic_variant_pos falling outside typical chr1 range,
            # or more directly, their original node not being chr1 according to node_pos_json if that info was kept)
            # will not be in vcf_chr1_variants and thus be 'false'.
            vcf_tuple_to_check = (genomic_variant_pos, variant_ref_on_node, variant_alt_on_node)
            is_true_variant = vcf_tuple_to_check in vcf_chr1_variants  # This implicitly filters for chr1 positions too

            variant_summary_entry = {
                "pth_file": pth_filename,
                "original_path": pth_file_full_path,
                "node_id": node_id_str,
                "variant_key": variant_key,
                "grch38_calculated_pos_1_based": genomic_variant_pos,
                "ref_allele_from_pth_summary": variant_ref_on_node,
                "alt_allele_from_pth_summary": variant_alt_on_node,
                "classification": "true" if is_true_variant else "false",
                "vcf_match_key": str(vcf_tuple_to_check),
                "found_in_vcf_chr1": is_true_variant  # Clarified that vcf set is for chr1
            }
            classification_summary.append(variant_summary_entry)

            destination_filename = f"{node_id_str}_{pth_filename}"

            if is_true_variant:
                shutil.copy2(pth_file_full_path, os.path.join(true_folder, destination_filename))
                pth_files_copied_true += 1
            else:
                shutil.copy2(pth_file_full_path, os.path.join(false_folder, destination_filename))
                pth_files_copied_false += 1

    # 5. Write the overall classification summary JSON
    summary_output_path = os.path.join(args.output_folder, "classification_summary.json")
    try:
        with open(summary_output_path, 'w') as f_summary:
            json.dump(classification_summary, f_summary, indent=4)
        print(f"Classification summary written to: {summary_output_path}")
    except Exception as e:
        print(f"Error writing classification summary: {e}", file=sys.stderr)

    print("\n--- Classification Stats ---")
    print(f"Total variant entries processed from summaries: {total_pth_files_processed_from_summaries}")
    print(f"Total .pth files found and considered for copying: {actual_pth_files_found_and_considered}")
    print(f"Copied to 'true' folder: {pth_files_copied_true}")
    print(f"Copied to 'false' folder: {pth_files_copied_false}")
    print("Classification finished.")


if __name__ == "__main__":
    if sys.version_info < (3, 6):
        sys.stderr.write("Error: This script requires Python 3.6 or higher.\n")
        sys.exit(1)
    main()