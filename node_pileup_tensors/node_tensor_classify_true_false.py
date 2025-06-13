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

# --- Global variables for worker processes ---
# These will be populated by the init_worker function in each process
vcf_variants = set()
node_positions = {}


def format_time(seconds):
    """Formats a duration in seconds into HH:MM:SS format."""
    seconds = int(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02}:{minutes:02}:{seconds:02}"


def query_vcf_for_chromosome(vcf_file_path, chromosome):
    """Queries a VCF file and returns a set of variants."""
    print(f"[{os.getpid()}] Querying VCF for {chromosome}...")
    local_variants = set()
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
                    for alt in fields[4].upper().split(','):
                        local_variants.add((pos, ref, alt))
        print(f"[{os.getpid()}] Loaded {len(local_variants)} variants from VCF.")
        return local_variants
    except FileNotFoundError:
        print(f"Error: 'bcftools' not found. Please ensure it's installed and in your PATH.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error running bcftools in process {os.getpid()}: {e}", file=sys.stderr)
        return set()


def load_node_positions(json_path):
    """Loads node position data from a JSON file."""
    print(f"[{os.getpid()}] Loading node position data...")
    try:
        with open(json_path) as f:
            data = json.load(f)
        local_node_dict = {
            str(node.get("node_id")): node.get("grch38_position_start")
            for node in data.get("nodes", [])
            if node.get("node_id") and isinstance(node.get("grch38_position_start"), int)
        }
        print(f"[{os.getpid()}] Loaded {len(local_node_dict)} node positions.")
        return local_node_dict
    except Exception as e:
        print(f"Error loading node positions in process {os.getpid()}: {e}", file=sys.stderr)
        sys.exit(1)


def init_worker(vcf_path, vcf_chrom, node_pos_path):
    """
    Initializer for each worker process. Loads read-only data into globals.
    This function runs ONCE per worker process.
    """
    global vcf_variants, node_positions
    print(f"Initializing worker {os.getpid()}...")
    vcf_variants = query_vcf_for_chromosome(vcf_path, vcf_chrom)
    node_positions = load_node_positions(node_pos_path)


def process_node_directory(node_dir, true_dir, false_dir, use_symlinks):
    """
    Processes a single node directory.
    NOTE: It now uses the global 'vcf_variants' and 'node_positions'
    instead of receiving them as arguments.
    """
    node_id = os.path.basename(node_dir)
    summary_path = os.path.join(node_dir, "variant_summary.json")

    if not os.path.isfile(summary_path):
        return [], 0, 0

    # Each worker has its own copy of this data now.
    start_pos = node_positions.get(node_id)
    if start_pos is None:
        return [], 0, 0

    local_true_count = 0
    local_false_count = 0
    summary_records = []

    try:
        with open(summary_path) as f:
            summary = json.load(f)
    except (IOError, json.JSONDecodeError):
        return [], 0, 0  # Silently fail for a single bad file

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

        dest_dir = true_dir if is_match else false_dir
        if is_match:
            local_true_count += 1
        else:
            local_false_count += 1

        destination_path = os.path.join(dest_dir, f"{node_id}_{tensor_file}")
        try:
            if use_symlinks:
                os.symlink(os.path.abspath(tensor_path), destination_path)
            else:
                shutil.copyfile(tensor_path, destination_path)
        except FileExistsError:
            pass
        except OSError:
            pass  # Silently ignore file creation errors

        summary_records.append({
            "node_id": node_id, "tensor_file": tensor_file,
            "variant_key": variant_key, "genomic_position": grch38_pos,
            "ref": ref, "alt": alt,
            "classification": "true" if is_match else "false"
        })

    return summary_records, local_true_count, local_false_count


def main():
    parser = argparse.ArgumentParser(
        description="Classify variant tensors against a VCF file (Robust Parallel Version).")
    # ... (rest of the argument parser setup is identical)
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

    node_dirs = [d.path for d in os.scandir(args.tensor_folder_path) if d.is_dir()]
    total_nodes = len(node_dirs)
    if total_nodes == 0:
        print("No node directories found to process. Exiting.")
        return

    # --- Parallel Processing ---
    print(f"\nStarting classification of {total_nodes} nodes using {args.workers} workers...")

    worker_func = partial(
        process_node_directory,
        true_dir=true_dir,
        false_dir=false_dir,
        use_symlinks=args.use_symlinks
    )

    all_summary_records = []
    total_true = 0
    total_false = 0
    start_time = time.monotonic()

    # This is the key change!
    # We pass the initializer function and its arguments.
    # Each worker will run init_worker before starting its tasks.
    init_args = (args.vcf_file, args.chr, args.node_pos_json)
    with Pool(processes=args.workers, initializer=init_worker, initargs=init_args) as pool:
        # imap_unordered is memory-efficient and provides results as they complete
        results = pool.imap_unordered(worker_func, node_dirs)

        for i, (records, true_c, false_c) in enumerate(results):
            processed_count = i + 1
            all_summary_records.extend(records)
            total_true += true_c
            total_false += false_c

            elapsed_time = time.monotonic() - start_time
            progress = processed_count / total_nodes
            rate = processed_count / elapsed_time if elapsed_time > 0 else 0

            eta_str = "..."
            if processed_count > 5 and rate > 0:
                eta_seconds = (total_nodes - processed_count) / rate
                eta_str = format_time(eta_seconds)

            elapsed_str = format_time(elapsed_time)
            print(
                f"\rProgress: {processed_count}/{total_nodes} ({progress:.1%}) | "
                f"Elapsed: {elapsed_str} | ETA: {eta_str} | Rate: {rate:.2f} nodes/s  ",
                end=""
            )

    print("\n\nProcessing complete.")

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
    # The if __name__ == "__main__" guard is CRITICAL for multiprocessing on Windows/macOS
    main()