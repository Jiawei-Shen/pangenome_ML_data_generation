#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import shutil
import sys
import time
from multiprocessing import Pool, cpu_count
from functools import partial


def format_time(seconds):
    """Formats a duration in seconds into HH:MM:SS format."""
    seconds = int(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02}:{minutes:02}:{seconds:02}"


def query_vcf_for_chromosome(vcf_file_path, chromosome):
    """
    Queries a VCF file for a specific chromosome using a streaming approach
    to handle potentially large outputs efficiently.
    """
    print(f"Querying VCF for {chromosome} using bcftools...")
    variants = set()
    # Use Popen to stream output, avoiding high memory use for large VCFs
    cmd = ['bcftools', 'view', '-r', chromosome, vcf_file_path]
    try:
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1) as proc:
            for line in proc.stdout:
                if line.startswith("#"):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    pos = int(fields[1])
                    ref = fields[3].upper()
                    # Handle multiple alternate alleles
                    for alt in fields[4].upper().split(','):
                        variants.add((pos, ref, alt))

        print(f"Loaded {len(variants)} variants from VCF for {chromosome}")
        return variants
    except FileNotFoundError:
        print(f"Error: 'bcftools' not found. Please ensure it's installed and in your PATH.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error running bcftools: {e.stderr}", file=sys.stderr)
        return set()


def load_node_positions(json_path):
    """Loads node position data from a JSON file."""
    print(f"Loading node position data from {json_path}...")
    try:
        with open(json_path) as f:
            data = json.load(f)
        node_dict = {
            str(node.get("node_id")): node.get("grch38_position_start")
            for node in data.get("nodes", [])
            if node.get("node_id") and isinstance(node.get("grch38_position_start"), int)
        }
        print(f"Loaded {len(node_dict)} node positions")
        return node_dict
    except FileNotFoundError:
        print(f"Error: Node position file not found at {json_path}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {json_path}", file=sys.stderr)
        sys.exit(1)


def process_node_directory(node_dir, vcf_variants, node_positions, true_dir, false_dir, use_symlinks):
    """
    Processes a single node directory, classifies its variants, and
    creates a symlink or copy of the tensor files.
    """
    node_id = os.path.basename(node_dir)
    summary_path = os.path.join(node_dir, "variant_summary.json")

    if not os.path.isfile(summary_path):
        return [], 0, 0

    try:
        with open(summary_path) as f:
            summary = json.load(f)
    except (IOError, json.JSONDecodeError) as e:
        print(f"\nWarning: Could not read or parse {summary_path}: {e}", file=sys.stderr)
        return [], 0, 0

    start_pos = node_positions.get(node_id)
    if start_pos is None:
        return [], 0, 0

    local_true_count = 0
    local_false_count = 0
    summary_records = []

    for variant in summary.get("variants_passing_af_filter", []):
        tensor_file = variant.get("tensor_file")
        variant_key = variant.get("variant_key")

        if not all([tensor_file, variant_key, tensor_file.endswith(".npy")]):
            continue

        tensor_path = os.path.join(node_dir, tensor_file)
        if not os.path.isfile(tensor_path):
            continue

        try:
            parts = variant_key.split("_")
            offset, ref, alt = int(parts[0]), parts[2].upper(), parts[3].upper()
        except (ValueError, IndexError):
            continue

        grch38_pos = start_pos + offset + 1
        is_match = (grch38_pos, ref, alt) in vcf_variants

        if is_match:
            dest_dir = true_dir
            local_true_count += 1
        else:
            dest_dir = false_dir
            local_false_count += 1

        destination_path = os.path.join(dest_dir, f"{node_id}_{tensor_file}")
        try:
            if use_symlinks:
                os.symlink(os.path.abspath(tensor_path), destination_path)
            else:
                shutil.copyfile(tensor_path, destination_path)
        except FileExistsError:
            pass  # Skip if file already exists from a previous run
        except OSError as e:
            print(f"\nWarning: Could not create link/copy for {tensor_path}: {e}", file=sys.stderr)

        summary_records.append({
            "node_id": node_id,
            "tensor_file": tensor_file,
            "variant_key": variant_key,
            "genomic_position": grch38_pos,
            "ref": ref,
            "alt": alt,
            "classification": "true" if is_match else "false"
        })

    return summary_records, local_true_count, local_false_count


def main():
    parser = argparse.ArgumentParser(description="Classify variant tensors against a VCF file.")
    parser.add_argument("tensor_folder_path", help="Path to the parent folder containing node tensor directories.")
    parser.add_argument("vcf_file", help="Path to the VCF file for ground truth variants.")
    parser.add_argument("node_pos_json", help="Path to the JSON file with node positions.")
    parser.add_argument("--output_folder", default="./classification_results", help="Folder to save results.")
    parser.add_argument("--chr", default="chr1", help="Chromosome to process.")
    parser.add_argument("--use-symlinks", action="store_true",
                        help="Use symlinks instead of copying files for major speedup.")
    parser.add_argument("-j", "--workers", type=int, default=cpu_count(), help="Number of worker processes to use.")
    args = parser.parse_args()

    # --- Setup ---
    os.makedirs(args.output_folder, exist_ok=True)
    true_dir = os.path.join(args.output_folder, "true")
    false_dir = os.path.join(args.output_folder, "false")
    os.makedirs(true_dir, exist_ok=True)
    os.makedirs(false_dir, exist_ok=True)

    # --- Data Loading ---
    vcf_variants = query_vcf_for_chromosome(args.vcf_file, args.chr)
    node_positions = load_node_positions(args.node_pos_json)

    node_dirs = [d.path for d in os.scandir(args.tensor_folder_path) if d.is_dir()]
    total_nodes = len(node_dirs)
    if total_nodes == 0:
        print("No node directories found to process. Exiting.")
        return

    # --- Parallel Processing ---
    print(f"\nStarting classification of {total_nodes} nodes using {args.workers} workers...")
    worker_func = partial(
        process_node_directory,
        vcf_variants=vcf_variants,
        node_positions=node_positions,
        true_dir=true_dir,
        false_dir=false_dir,
        use_symlinks=args.use_symlinks
    )

    all_summary_records = []
    total_true = 0
    total_false = 0
    start_time = time.monotonic()

    with Pool(processes=args.workers) as pool:
        results = pool.imap_unordered(worker_func, node_dirs)

        for i, (records, true_c, false_c) in enumerate(results):
            processed_count = i + 1
            all_summary_records.extend(records)
            total_true += true_c
            total_false += false_c

            elapsed_time = time.monotonic() - start_time

            # --- Progress and ETA Calculation ---
            progress = processed_count / total_nodes
            rate = processed_count / elapsed_time if elapsed_time > 0 else 0

            eta_str = "..."
            # Display ETA only after a few results for a more stable estimate
            if processed_count > 5 and rate > 0:
                remaining_items = total_nodes - processed_count
                eta_seconds = remaining_items / rate
                eta_str = format_time(eta_seconds)

            elapsed_str = format_time(elapsed_time)
            print(
                f"\rProgress: {processed_count}/{total_nodes} ({progress:.1%}) | "
                f"Elapsed: {elapsed_str} | ETA: {eta_str} | Rate: {rate:.2f} nodes/s  ",
                end=""
            )

    print("\nProcessing complete.")

    # --- Finalization ---
    summary_json_path = os.path.join(args.output_folder, "classification_summary.json")
    print(f"Writing classification summary to {summary_json_path}...")
    with open(summary_json_path, 'w') as f:
        json.dump({"chromosome": args.chr, "results": all_summary_records}, f, indent=2)

    print("\n--- Classification Summary ---")
    print(f"Total nodes processed: {total_nodes}")
    print(f"Total tensors classified: {total_true + total_false}")
    print(f"True (matches in VCF): {total_true}")
    print(f"False (not in VCF): {total_false}")


if __name__ == "__main__":
    if sys.version_info < (3, 6):
        sys.exit("This script requires Python 3.6+.")
    main()