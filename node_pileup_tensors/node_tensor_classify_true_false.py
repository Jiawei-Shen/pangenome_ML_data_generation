#!/usr/bin/env python3
import argparse
import json
import os
import sys
import subprocess
import shutil


# --- Helper Functions ---

def query_vcf_chr1(bcftools_path, vcf_file_path):
    """
    Queries a VCF file for variants on chr1 using bcftools.
    Returns a set of (position, ref, alt) tuples for chr1 variants.
    """
    print(f"🧬 Querying VCF file '{vcf_file_path}' for chr1 variants...")
    chr1_variants = set()
    # bcftools view -r chr1 <vcf_file_path>
    # We need to handle potential chromosome name variations (e.g., "1" vs "chr1")
    # For simplicity, this example assumes "chr1". You might need to adjust.
    cmd = [bcftools_path, 'view', '-r', 'chr1', vcf_file_path]
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        for line in process.stdout:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    # VCF: CHROM POS ID REF ALT ...
                    # VCF is 1-based, convert POS to int
                    try:
                        pos = int(fields[1])
                        ref = fields[3]
                        alt_alleles = fields[4].split(',')  # ALT can be comma-separated
                        for alt in alt_alleles:
                            chr1_variants.add((pos, ref, alt))
                    except ValueError:
                        print(f"⚠️ Warning: Could not parse VCF line: {line.strip()}", file=sys.stderr)
        process.wait()
        if process.returncode != 0:
            print(f"❌ Error: bcftools command failed with exit code {process.returncode}", file=sys.stderr)
            return None
        print(f"✔ Found {len(chr1_variants)} unique ALT variants on chr1 in VCF.")
        return chr1_variants
    except FileNotFoundError:
        print(f"❌ Error: bcftools not found at '{bcftools_path}'. Please ensure it's installed and path is correct.",
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
            if (entry.get("status") != "task_execution_error" and
                    entry.get(
                        "grch38_chromosome_region") == "chr1" and  # Assuming exact match, adjust if needed (e.g. "1")
                    entry.get("grch38_start_position") is not None):
                try:
                    node_id = str(entry.get("target_node_id_queried"))  # Ensure string key
                    start_pos = int(entry.get("grch38_start_position"))  # Store as int
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
    # If node_start_pos is 1-based and variant_pos_on_node is 0-based:
    # Genomic position (1-based) = node_start_pos (1-based) + variant_pos_on_node (0-based)
    return node_start_pos + variant_pos_on_node


# --- Main Script Logic ---

def main():
    parser = argparse.ArgumentParser(description="Classify .pth variant files based on VCF concordance on chr1.")
    parser.add_argument("pth_folder_path",
                        help="Path to the folder containing node subdirectories with .pth files and variant_summary.json.")
    parser.add_argument("vcf_file", help="Path to the VCF file to query (e.g., dbSNP for chr1).")
    parser.add_argument("node_pos_json", help="Path to the JSON file with node GRCh38 position information.")
    parser.add_argument("--bcftools_path", default="bcftools",
                        help="Path to the bcftools executable (default: 'bcftools').")
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
    # This set will store (position, REF, ALT) from VCF for quick lookup
    # VCF positions are 1-based.
    vcf_chr1_variants = query_vcf_chr1(args.bcftools_path, args.vcf_file)
    if vcf_chr1_variants is None:
        print("❌ Exiting due to VCF query failure.", file=sys.stderr)
        sys.exit(1)
    if not vcf_chr1_variants:
        print("⚠️ Warning: No chr1 variants found in the VCF file. All .pth files might be classified as 'false'.",
              file=sys.stderr)

    # 3. Load node position information for chr1
    # This dict will be node_id (str) -> grch38_start_pos (int, 1-based)
    node_start_positions_chr1 = load_node_positions(args.node_pos_json)
    if node_start_positions_chr1 is None:
        print("❌ Exiting due to failure in loading node positions.", file=sys.stderr)
        sys.exit(1)
    if not node_start_positions_chr1:
        print("⚠️ Warning: No chr1 node start positions loaded. Cannot calculate genomic positions for .pth files.",
              file=sys.stderr)

    # 4. Iterate through .pth files, calculate positions, classify, and copy
    classification_summary = []
    total_pth_files_processed = 0
    pth_files_true = 0
    pth_files_false = 0

    print(f"📁 Processing .pth files from base folder: {args.pth_folder_path}")
    # The .pth files are in node_id subdirectories, along with variant_summary.json
    for node_id_str in os.listdir(args.pth_folder_path):
        node_dir_path = os.path.join(args.pth_folder_path, node_id_str)
        if not os.path.isdir(node_dir_path):
            continue

        # Check if this node is on chr1 and we have its start position
        if node_id_str not in node_start_positions_chr1:
            # print(f"ℹ️ Skipping node {node_id_str}: Not on chr1 or start position unknown based on node position JSON.")
            continue

        node_grch38_start_pos = node_start_positions_chr1[node_id_str]  # 1-based

        summary_json_path = os.path.join(node_dir_path, "variant_summary.json")
        if not os.path.isfile(summary_json_path):
            # print(f"⚠️ Warning: variant_summary.json not found in {node_dir_path}. Skipping this node's .pth files.", file=sys.stderr)
            continue

        try:
            with open(summary_json_path, 'r') as sjf:
                node_summary_data = json.load(sjf)
        except Exception as e:
            print(f"⚠️ Warning: Could not read or parse {summary_json_path}: {e}. Skipping.", file=sys.stderr)
            continue

        # Ensure the summary is for the correct chromosome via the node position JSON (already filtered)
        # We assume all variants in a chr1 node's summary are also on chr1.

        for variant_info in node_summary_data.get("variants", []):
            total_pth_files_processed += 1
            pth_filename = variant_info.get("tensor_file")
            variant_key = variant_info.get("variant_key")  # e.g., "pos_type_ref_alt"
            if not pth_filename or not variant_key:
                print(
                    f"⚠️ Warning: Missing 'tensor_file' or 'variant_key' in {summary_json_path} for an entry. Skipping.",
                    file=sys.stderr)
                continue

            pth_file_full_path = os.path.join(node_dir_path, pth_filename)
            if not os.path.isfile(pth_file_full_path):
                print(f"⚠️ Warning: .pth file not found: {pth_file_full_path}. Skipping.", file=sys.stderr)
                continue

            # Parse variant_key: "v_pos_v_type_v_ref_defined_v_alt_defined"
            # v_pos is 0-based position on the node
            try:
                parts = variant_key.split('_')
                variant_pos_on_node = int(parts[0])  # 0-based
                # variant_type = parts[1]
                variant_ref_on_node = parts[2]  # From the node sequence at that position (or surrounding for indel)
                variant_alt_on_node = parts[3]  # The ALT allele observed in reads
            except (IndexError, ValueError) as e:
                print(
                    f"⚠️ Warning: Could not parse variant_key '{variant_key}' from {summary_json_path}: {e}. Skipping {pth_filename}.",
                    file=sys.stderr)
                continue

            # Calculate genomic position (1-based for VCF comparison)
            # node_grch38_start_pos is 1-based, variant_pos_on_node is 0-based
            # For SNPs and Deletions, variant_pos_on_node is the first affected base.
            # For Insertions, variant_pos_on_node is the base *before* the insertion.
            # The VCF position for an insertion is typically the base before the inserted sequence.
            # The VCF position for a deletion is the base before the deleted sequence (or first base of del).
            # The VCF position for SNP is the SNP's position.
            # This simple addition should align well for SNPs. For Indels, VCF representation can vary slightly.
            # For this script, we'll use this direct calculation. More complex indel normalization might be needed for perfect matching.

            # Simple SNP/Point variant case:
            # Genomic position of variant (1-based) = node_grch38_start_pos (1-based) + variant_pos_on_node (0-based)
            genomic_variant_pos = calculate_genomic_position(node_grch38_start_pos, variant_pos_on_node)

            # VCF REF/ALT alleles might need to be handled carefully for indels.
            # For now, we assume a simple direct comparison for point mutations.
            # The variant_ref_on_node and variant_alt_on_node are from the *graph node's perspective*.
            # If the node itself has a variation relative to GRCh38 standard reference, this needs careful thought.
            # This script currently assumes variant_ref_on_node corresponds to GRCh38 REF at genomic_variant_pos.

            # For this example, let's assume:
            #   - variant_ref_on_node is the reference allele at the calculated genomic_variant_pos on GRCh38.
            #   - variant_alt_on_node is the alternate allele.
            #   This assumption is strong and might not hold if the graph node itself is a variant.

            vcf_tuple_to_check = (genomic_variant_pos, variant_ref_on_node.upper(), variant_alt_on_node.upper())

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
                "found_in_vcf_details": str(vcf_tuple_to_check) if is_true_variant else "Not found"
            }
            classification_summary.append(variant_summary_entry)

            if is_true_variant:
                shutil.copy2(pth_file_full_path, os.path.join(true_folder, pth_filename))
                pth_files_true += 1
            else:
                shutil.copy2(pth_file_full_path, os.path.join(false_folder, pth_filename))
                pth_files_false += 1

        # Optional: print progress per node directory
        # print(f"  Processed node {node_id_str}")

    # 5. Write the overall classification summary JSON
    summary_output_path = os.path.join(args.output_folder, "classification_summary.json")
    try:
        with open(summary_output_path, 'w') as f_summary:
            json.dump(classification_summary, f_summary, indent=4)
        print(f"✔️ Classification summary written to: {summary_output_path}")
    except Exception as e:
        print(f"❌ Error writing classification summary: {e}", file=sys.stderr)

    print("\n--- Classification Stats ---")
    print(f"Total .pth files processed entries from summaries: {total_pth_files_processed}")
    print(f"Copied to 'true' folder: {pth_files_true}")
    print(f"Copied to 'false' folder: {pth_files_false}")
    print("🎉 Classification finished.")


if __name__ == "__main__":
    # Add a check for Python 3.6+ for f-strings and subprocess improvements if necessary,
    # though current usage is broadly compatible.
    if sys.version_info < (3, 6):
        sys.stderr.write("❌ This script requires Python 3.6 or higher.\n")
        sys.exit(1)
    main()