#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import shutil
import sys  # Added for sys.exit and sys.stderr


# --- Helper Functions ---

def query_vcf_chr1(vcf_file_path):
    """
    Queries a VCF file for variants on chr1 using bcftools.
    Assumes 'bcftools' is in the system PATH.
    Returns a set of (position, ref, alt) tuples for chr1 variants.
    """
    print(f"🧬 Querying VCF file '{vcf_file_path}' for chr1 variants using 'bcftools' from PATH...")
    chr1_variants = set()
    # bcftools view -r chr1 <vcf_file_path>
    # We need to handle potential chromosome name variations (e.g., "1" vs "chr1")
    # For simplicity, this example assumes "chr1". You might need to adjust.
    cmd = ['bcftools', 'view', '-r', 'chr1', vcf_file_path]
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE)  # Added stderr
        stdout, stderr = process.communicate()  # Capture stdout and stderr

        if process.returncode != 0:
            print(f"❌ Error: bcftools command failed with exit code {process.returncode}", file=sys.stderr)
            if stderr:
                print(f"bcftools stderr:\n{stderr.strip()}", file=sys.stderr)
            return None

        for line in stdout.splitlines():  # Process captured stdout
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    # VCF: CHROM POS ID REF ALT ...
                    # VCF is 1-based, convert POS to int
                    try:
                        pos = int(fields[1])
                        ref = fields[3].upper()  # Normalize to uppercase
                        alt_alleles = fields[4].upper().split(',')  # Normalize to uppercase
                        for alt in alt_alleles:
                            chr1_variants.add((pos, ref, alt))
                    except ValueError:
                        print(f"⚠️ Warning: Could not parse VCF line: {line.strip()}", file=sys.stderr)

        print(f"✔ Found {len(chr1_variants)} unique ALT variants on chr1 in VCF.")
        return chr1_variants
    except FileNotFoundError:
        print(f"❌ Error: 'bcftools' not found in PATH. Please ensure it's installed and PATH is configured correctly.",
              file=sys.stderr)
        return None
    except Exception as e:
        print(f"❌ Error during VCF query: {e}", file=sys.stderr)
        return None


def load_node_positions(json_file_path):
    """
    Loads node position information from the provided JSON file.
    Returns a dictionary mapping target_node_id (str) to grch38_start_position.
    Only records entries for 'chr1' and where grch38_start_position is not null.
    """
    print(f"🗺️ Loading node positions from '{json_file_path}' for chr1...")
    node_positions = {}
    try:
        with open(json_file_path, 'r') as f:
            data = json.load(f)
        for entry in data:
            if (entry.get("grch38_chromosome_region") == "chr1" and
                    entry.get("grch38_start_position") is not None):
                try:
                    node_id = str(entry.get("target_node_id_queried"))
                    start_pos = int(entry.get("grch38_start_position"))
                    node_positions[node_id] = start_pos
                except (ValueError, TypeError):
                    print(f"⚠️ Warning: Skipping entry with invalid start position or node ID: {entry}",
                          file=sys.stderr)

        print(f"✔ Loaded start positions for {len(node_positions)} nodes on chr1.")
        return node_positions
    except FileNotFoundError:
        print(f"❌ Error: Node position JSON file not found at '{json_file_path}'", file=sys.stderr)
        return None
    except json.JSONDecodeError:
        print(f"❌ Error: Could not decode JSON from '{json_file_path}'", file=sys.stderr)
        return None
    except Exception as e:
        print(f"❌ Error loading node positions: {e}", file=sys.stderr)
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
    parser.add_argument("node_pos_json", help="Path to the JSON file with node GRCh38 position information.")
    parser.add_argument("--output_folder", default="./classification_results",
                        help="Path to the folder where 'true' and 'false' subfolders and summary will be created.")
    args = parser.parse_args()

    print("🚀 Starting .pth file classification process...")

    # 1. Create output directories
    true_folder = os.path.join(args.output_folder, "true")
    false_folder = os.path.join(args.output_folder, "false")
    os.makedirs(true_folder, exist_ok=True)
    os.makedirs(false_folder, exist_ok=True)
    print(f"✔️ Output folders created/ensured: '{true_folder}', '{false_folder}'")

    # 2. Query VCF for chr1 variants
    vcf_chr1_variants = query_vcf_chr1(args.vcf_file)  # Removed bcftools_path argument
    if vcf_chr1_variants is None:
        print("❌ Exiting due to VCF query failure.", file=sys.stderr)
        sys.exit(1)
    if not vcf_chr1_variants:
        print("⚠️ Warning: No chr1 variants found in the VCF file. All .pth files might be classified as 'false'.",
              file=sys.stderr)

    # 3. Load node position information for chr1
    node_start_positions_chr1 = load_node_positions(args.node_pos_json)
    if node_start_positions_chr1 is None:
        print("❌ Exiting due to failure in loading node positions.", file=sys.stderr)
        sys.exit(1)
    if not node_start_positions_chr1:
        print("⚠️ Warning: No chr1 node start positions loaded. Cannot calculate genomic positions for .pth files.",
              file=sys.stderr)

    # 4. Iterate through .pth files, calculate positions, classify, and copy
    classification_summary = []
    total_pth_files_processed_from_summaries = 0
    pth_files_copied_true = 0
    pth_files_copied_false = 0
    actual_pth_files_found_and_considered = 0

    print(f"📁 Processing .pth files from base folder: {args.pth_folder_path}")
    for node_id_str in os.listdir(args.pth_folder_path):
        node_dir_path = os.path.join(args.pth_folder_path, node_id_str)
        if not os.path.isdir(node_dir_path):
            continue

        if node_id_str not in node_start_positions_chr1:
            continue

        node_grch38_start_pos = node_start_positions_chr1[node_id_str]

        summary_json_path = os.path.join(node_dir_path, "variant_summary.json")
        if not os.path.isfile(summary_json_path):
            continue

        try:
            with open(summary_json_path, 'r') as sjf:
                node_summary_data = json.load(sjf)
        except Exception as e:
            print(f"⚠️ Warning: Could not read or parse {summary_json_path}: {e}. Skipping.", file=sys.stderr)
            continue

        for variant_info in node_summary_data.get("variants", []):
            total_pth_files_processed_from_summaries += 1
            pth_filename = variant_info.get("tensor_file")
            variant_key = variant_info.get("variant_key")
            if not pth_filename or not variant_key:
                print(
                    f"⚠️ Warning: Missing 'tensor_file' or 'variant_key' in {summary_json_path} for an entry. Skipping.",
                    file=sys.stderr)
                continue

            pth_file_full_path = os.path.join(node_dir_path, pth_filename)
            if not os.path.isfile(pth_file_full_path):
                print(
                    f"⚠️ Warning: .pth file not found: {pth_file_full_path} (referenced in summary). Skipping this entry.",
                    file=sys.stderr)
                continue

            actual_pth_files_found_and_considered += 1

            try:
                parts = variant_key.split('_')
                variant_pos_on_node = int(parts[0])
                # variant_type = parts[1] # Not directly used for classification logic here
                variant_ref_on_node = parts[2].upper()  # Normalize to uppercase
                variant_alt_on_node = parts[3].upper()  # Normalize to uppercase
            except (IndexError, ValueError) as e:
                print(
                    f"⚠️ Warning: Could not parse variant_key '{variant_key}' from {summary_json_path}: {e}. Skipping {pth_filename}.",
                    file=sys.stderr)
                continue

            genomic_variant_pos = calculate_genomic_position(node_grch38_start_pos, variant_pos_on_node)
            vcf_tuple_to_check = (genomic_variant_pos, variant_ref_on_node, variant_alt_on_node)
            is_true_variant = vcf_tuple_to_check in vcf_chr1_variants

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
                "found_in_vcf": is_true_variant
            }
            classification_summary.append(variant_summary_entry)

            destination_filename = f"{node_id_str}_{pth_filename}"  # Prepend node_id to avoid name clashes

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
        print(f"✔️ Classification summary written to: {summary_output_path}")
    except Exception as e:
        print(f"❌ Error writing classification summary: {e}", file=sys.stderr)

    print("\n--- Classification Stats ---")
    print(f"Total variant entries processed from summaries: {total_pth_files_processed_from_summaries}")
    print(f"Total .pth files found and considered for copying: {actual_pth_files_found_and_considered}")
    print(f"Copied to 'true' folder: {pth_files_copied_true}")
    print(f"Copied to 'false' folder: {pth_files_copied_false}")
    print("🎉 Classification finished.")


if __name__ == "__main__":
    if sys.version_info < (3, 6):  # Ensure Python 3.6+ for f-strings, subprocess.communicate etc.
        sys.stderr.write("❌ This script requires Python 3.6 or higher.\n")
        sys.exit(1)
    main()