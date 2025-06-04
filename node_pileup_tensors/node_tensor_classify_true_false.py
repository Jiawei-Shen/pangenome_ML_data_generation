#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import shutil
import sys


# --- Helper Functions ---

def query_vcf_for_chromosome(vcf_file_path, chromosome):  # Renamed and added chromosome parameter
    """
    Queries a VCF file for variants on the specified chromosome using bcftools.
    Assumes 'bcftools' is in the system PATH.
    Returns a set of (position, ref, alt) tuples for variants on the specified chromosome.
    """
    print(f"Querying VCF file '{vcf_file_path}' for '{chromosome}' variants using 'bcftools' from PATH...")
    chromosome_variants = set()
    cmd = ['bcftools', 'view', '-r', chromosome, vcf_file_path]  # Use the chromosome parameter
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            print(f"Error: bcftools command failed with exit code {process.returncode} for chromosome '{chromosome}'",
                  file=sys.stderr)
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
                            chromosome_variants.add((pos, ref, alt))
                    except ValueError:
                        print(f"Warning: Could not parse VCF line: {line.strip()}", file=sys.stderr)

        print(f"Found {len(chromosome_variants)} unique ALT variants on '{chromosome}' in VCF.")
        return chromosome_variants
    except FileNotFoundError:
        print(f"Error: 'bcftools' not found in PATH. Please ensure it's installed and PATH is configured correctly.",
              file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error during VCF query for chromosome '{chromosome}': {e}", file=sys.stderr)
        return None


def load_node_positions(json_file_path):
    """
    Loads node position information from the provided JSON file.
    Expects a single JSON object at the top level, with a "nodes" key containing a list of node objects.
    Returns a dictionary mapping node_id (str) to grch38_position_start (int).
    """
    print(f"Loading node positions from '{json_file_path}'...")
    node_positions = {}
    try:
        with open(json_file_path, 'r') as f:
            loaded_json_object = json.load(f)

        if not isinstance(loaded_json_object, dict):
            print(
                f"Error: Node position JSON file '{json_file_path}' is not a JSON object at the top level as expected.",
                file=sys.stderr)
            return None

        node_entries_list = loaded_json_object.get("nodes")

        if not isinstance(node_entries_list, list):
            print(f"Error: The 'nodes' key in '{json_file_path}' does not contain a list of node entries.",
                  file=sys.stderr)
            return None

        for entry_num, entry in enumerate(node_entries_list, 1):
            if not isinstance(entry, dict):
                print(
                    f"Warning: Item {entry_num} in the 'nodes' list of '{json_file_path}' is not a dictionary. Skipping.",
                    file=sys.stderr)
                continue

            start_pos_val = entry.get("grch38_position_start")

            if start_pos_val is not None:
                try:
                    node_id_val = entry.get("node_id")
                    if node_id_val is None:
                        print(f"Warning: Skipping entry {entry_num} in 'nodes' list due to missing 'node_id': {entry}",
                              file=sys.stderr)
                        continue
                    node_id_str = str(node_id_val)

                    start_pos_int = int(start_pos_val)
                    node_positions[node_id_str] = start_pos_int

                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Skipping entry {entry_num} in 'nodes' list due to invalid data type for 'node_id' or 'grch38_position_start': {entry} (Error: {e})",
                        file=sys.stderr)

        print(f"Loaded start positions for {len(node_positions)} nodes from the 'nodes' list.")
        return node_positions
    except FileNotFoundError:
        print(f"Error: Node position JSON file not found at '{json_file_path}'", file=sys.stderr)
        return None
    except json.JSONDecodeError as jde:
        print(f"Error: Could not decode JSON from '{json_file_path}'. Malformed JSON. Error: {jde}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error loading node positions: {e}", file=sys.stderr)
        return None


def calculate_genomic_position(node_start_pos, variant_pos_on_node):
    """
    Calculates the 1-based genomic position.
    Assumes node_start_pos is 1-based GRCh38 coordinate.
    Assumes variant_pos_on_node is 0-based offset from the start of the node.
    """
    # Corrected calculation: 1-based start + 0-based offset = 1-based result
    return node_start_pos + variant_pos_on_node + 1


# --- Main Script Logic ---

def main():
    parser = argparse.ArgumentParser(description="Classify .npy tensor files based on VCF concordance on chr1.")
    parser.add_argument("tensor_folder_path",
                        help="Path to the folder containing node subdirectories with .npy files and variant_summary.json.")
    parser.add_argument("vcf_file", help="Path to the VCF file to query.")
    parser.add_argument("node_pos_json",
                        help="Path to the JSON file with node GRCh38 position information (expects a top-level object with a 'nodes' list).")
    parser.add_argument("--output_folder", default="./classification_results",
                        help="Path to the folder where 'true' and 'false' subfolders and summary will be created.")
    # Added --chr argument from previous request, ensuring it's here.
    parser.add_argument("--chr", "--chromosome", default="chr1",
                        help="Chromosome to query in the VCF and use for classification (e.g., chr1, chrX, 1). Default: chr1.")
    args = parser.parse_args()

    print(f"Starting .npy file classification process for chromosome '{args.chr}'...")

    # 1. Create output directories
    true_folder = os.path.join(args.output_folder, "true")
    false_folder = os.path.join(args.output_folder, "false")
    os.makedirs(true_folder, exist_ok=True)
    os.makedirs(false_folder, exist_ok=True)
    print(f"Output folders created/ensured: '{true_folder}', '{false_folder}'")

    # 2. Query VCF for specified chromosome variants
    # Corrected call using the function name from the previous version that accepted chromosome argument
    vcf_chromosome_variants = query_vcf_for_chromosome(args.vcf_file, args.chr)

    if vcf_chromosome_variants is None:
        print(f"Exiting due to VCF query failure for chromosome '{args.chr}'.", file=sys.stderr)
        sys.exit(1)
    if not vcf_chromosome_variants:
        print(
            f"Warning: No variants found for chromosome '{args.chr}' in the VCF file. All .npy files will likely be classified as 'false'.",
            file=sys.stderr)

    # 3. Load node position information
    all_node_start_positions = load_node_positions(args.node_pos_json)
    if all_node_start_positions is None:
        print("Exiting due to failure in loading node positions.", file=sys.stderr)
        sys.exit(1)
    if not all_node_start_positions:
        print("Warning: No node start positions loaded. Cannot calculate genomic positions for .npy files.",
              file=sys.stderr)

    # 4. Iterate through tensor files, calculate positions, classify, and copy
    classification_summary = []
    total_tensor_files_processed_from_summaries = 0
    tensor_files_copied_true = 0
    tensor_files_copied_false = 0
    actual_tensor_files_found_and_considered = 0
    nodes_processed_for_report = 0

    print(f"Processing .npy files from base folder: {args.tensor_folder_path}")
    for node_id_str in os.listdir(args.tensor_folder_path):
        node_dir_path = os.path.join(args.tensor_folder_path, node_id_str)
        if not os.path.isdir(node_dir_path):
            continue

        if node_id_str not in all_node_start_positions:
            continue

        nodes_processed_for_report += 1

        node_grch38_start_pos = all_node_start_positions[node_id_str]

        summary_json_path = os.path.join(node_dir_path, "variant_summary.json")
        if not os.path.isfile(summary_json_path):
            if nodes_processed_for_report > 0 and nodes_processed_for_report % 100 == 0:
                print(f"\nProgress Report: After processing {nodes_processed_for_report} node directories.")
                print(f"  Total .npy files considered (from summaries): {total_tensor_files_processed_from_summaries}")
                print(f"  Total .npy files physically found & processed: {actual_tensor_files_found_and_considered}")
                print(
                    f"  Copied to 'true': {tensor_files_copied_true}, Copied to 'false': {tensor_files_copied_false}\n")
            continue

        try:
            with open(summary_json_path, 'r') as sjf:
                node_summary_data = json.load(sjf)
        except Exception as e:
            print(f"Warning: Could not read or parse {summary_json_path}: {e}. Skipping.", file=sys.stderr)
            if nodes_processed_for_report > 0 and nodes_processed_for_report % 100 == 0:
                print(f"\nProgress Report: After processing {nodes_processed_for_report} node directories.")
                print(f"  Total .npy files considered (from summaries): {total_tensor_files_processed_from_summaries}")
                print(f"  Total .npy files physically found & processed: {actual_tensor_files_found_and_considered}")
                print(
                    f"  Copied to 'true': {tensor_files_copied_true}, Copied to 'false': {tensor_files_copied_false}\n")
            continue

        for variant_info in node_summary_data.get("variants", []):
            total_tensor_files_processed_from_summaries += 1
            tensor_filename = variant_info.get("tensor_file")
            variant_key = variant_info.get("variant_key")
            if not tensor_filename or not variant_key:
                print(f"Warning: Missing 'tensor_file' or 'variant_key' in {summary_json_path} for an entry. Skipping.",
                      file=sys.stderr)
                continue

            if not tensor_filename.endswith(".npy"):
                print(
                    f"Warning: Expected .npy file in summary, but found '{tensor_filename}'. Skipping entry in {summary_json_path}.",
                    file=sys.stderr)
                continue

            tensor_file_full_path = os.path.join(node_dir_path, tensor_filename)
            if not os.path.isfile(tensor_file_full_path):
                print(
                    f"Warning: Tensor file not found: {tensor_file_full_path} (referenced in summary). Skipping this entry.",
                    file=sys.stderr)
                continue

            actual_tensor_files_found_and_considered += 1

            try:
                parts = variant_key.split('_')
                variant_pos_on_node = int(parts[0])  # This is 0-based
                variant_ref_on_node = parts[2].upper()
                variant_alt_on_node = parts[3].upper()
            except (IndexError, ValueError) as e:
                print(
                    f"Warning: Could not parse variant_key '{variant_key}' from {summary_json_path}: {e}. Skipping {tensor_filename}.",
                    file=sys.stderr)
                continue

            allele_frequency = variant_info.get("alt_allele_frequency", "N/A")

            # Use the corrected calculate_genomic_position
            genomic_variant_pos = calculate_genomic_position(node_grch38_start_pos, variant_pos_on_node)

            vcf_tuple_to_check = (genomic_variant_pos, variant_ref_on_node, variant_alt_on_node)
            print(vcf_tuple_to_check)
            is_true_variant = vcf_tuple_to_check in vcf_chromosome_variants

            variant_summary_entry = {
                "tensor_file": tensor_filename,
                "original_path": tensor_file_full_path,
                "node_id": node_id_str,
                "variant_key": variant_key,
                "grch38_calculated_pos_1_based": genomic_variant_pos,
                "ref_allele_from_tensor_summary": variant_ref_on_node,
                "alt_allele_from_tensor_summary": variant_alt_on_node,
                "allele_frequency": allele_frequency,
                "classification": "true" if is_true_variant else "false",
                "vcf_match_key": str(vcf_tuple_to_check),
                "found_in_vcf_for_queried_chromosome": is_true_variant
            }
            classification_summary.append(variant_summary_entry)

            destination_filename = f"{node_id_str}_{tensor_filename}"

            if is_true_variant:
                shutil.copy2(tensor_file_full_path, os.path.join(true_folder, destination_filename))
                tensor_files_copied_true += 1
            else:
                shutil.copy2(tensor_file_full_path, os.path.join(false_folder, destination_filename))
                tensor_files_copied_false += 1

        if nodes_processed_for_report > 0 and nodes_processed_for_report % 100 == 0:
            print(f"\nProgress Report: After processing {nodes_processed_for_report} node directories.")
            print(f"  Total .npy files considered (from summaries): {total_tensor_files_processed_from_summaries}")
            print(f"  Total .npy files physically found & processed: {actual_tensor_files_found_and_considered}")
            print(f"  Copied to 'true': {tensor_files_copied_true}, Copied to 'false': {tensor_files_copied_false}\n")

    # 5. Write the overall classification summary JSON
    final_summary_output = {
        "queried_chromosome": args.chr,
        "classification_details": classification_summary
    }
    summary_output_path = os.path.join(args.output_folder, "classification_summary.json")
    try:
        with open(summary_output_path, 'w') as f_summary:
            json.dump(final_summary_output, f_summary, indent=4)
        print(f"Classification summary written to: {summary_output_path}")
    except Exception as e:
        print(f"Error writing classification summary: {e}", file=sys.stderr)

    print("\n--- Classification Stats ---")
    print(f"Queried Chromosome for VCF: {args.chr}")
    print(f"Total variant entries processed from summaries: {total_tensor_files_processed_from_summaries}")
    print(f"Total .npy files found and considered for copying: {actual_tensor_files_found_and_considered}")
    print(f"Copied to 'true' folder: {tensor_files_copied_true}")
    print(f"Copied to 'false' folder: {tensor_files_copied_false}")
    print(f"Total node directories processed for reporting: {nodes_processed_for_report}")
    print("Classification finished.")


if __name__ == "__main__":
    if sys.version_info < (3, 6):
        sys.stderr.write("Error: This script requires Python 3.6 or higher.\n")
        sys.exit(1)
    main()