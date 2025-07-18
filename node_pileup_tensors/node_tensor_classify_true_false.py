#!/usr/bin/env python3
import argparse
import json
import os
import shutil
import sys
import time
import random
from multiprocessing import Pool, cpu_count
from functools import partial
import pysam  # Use pysam instead of subprocess

# --- Global variables for worker processes ---
vcf_variants = set()
node_positions = {}

def format_time(seconds):
    seconds = int(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02}:{minutes:02}:{seconds:02}"

def get_variant_count(vcf_file_path, chromosome):
    print(f"Counting variants in {os.path.basename(vcf_file_path)} for {chromosome}...")
    count = 0
    try:
        with pysam.VariantFile(vcf_file_path) as vcf_in:
            for record in vcf_in.fetch(chromosome):
                if record.alts:
                    count += len(record.alts)
        return count
    except ValueError as e:
        print(f"Error: Could not read VCF file with pysam. Is it a valid bgzipped VCF with a .tbi index?", file=sys.stderr)
        print(f"Pysam error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while counting variants: {e}", file=sys.stderr)
        return -1

def query_vcf_for_chromosome(vcf_file_path, chromosome):
    local_variants = set()
    try:
        with pysam.VariantFile(vcf_file_path) as vcf_in:
            for record in vcf_in.fetch(chromosome):
                pos = record.pos
                ref = record.ref.upper()
                if record.alts:
                    for alt in record.alts:
                        local_variants.add((pos, ref, alt.upper()))
        return local_variants
    except ValueError as e:
        print(f"FATAL ERROR in worker: Pysam could not read VCF file.", file=sys.stderr)
        print(f"Please ensure '{vcf_file_path}' is a bgzipped VCF with a .tbi index.", file=sys.stderr)
        print(f"Pysam error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception:
        return set()

def load_node_positions(json_path):
    try:
        with open(json_path) as f:
            data = json.load(f)
        return {
            str(node.get("node_id")): node.get("grch38_position_start")
            for node in data.get("nodes", [])
            if node.get("node_id") and isinstance(node.get("grch38_position_start"), int)
        }
    except Exception:
        return {}

def init_worker(vcf_path, vcf_chrom, node_pos_path):
    global vcf_variants, node_positions
    vcf_variants = query_vcf_for_chromosome(vcf_path, vcf_chrom)
    node_positions = load_node_positions(node_pos_path)

def process_node_directory(node_dir, true_dir, false_dir, use_symlinks):
    node_id = os.path.basename(node_dir)
    summary_path = os.path.join(node_dir, "variant_summary.json")
    if not os.path.isfile(summary_path):
        return [], 0, 0
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
        return [], 0, 0

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
        grch38_pos = start_pos + offset
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
        except (FileExistsError, OSError):
            pass
        summary_records.append({
            "node_id": node_id, "tensor_file": tensor_file,
            "variant_key": variant_key, "genomic_position": grch38_pos,
            "ref": ref, "alt": alt,
            "classification": "true" if is_match else "false"
        })
    return summary_records, local_true_count, local_false_count

def organize_classified_data(base_dir, ratios, seed=None):
    if seed is not None:
        random.seed(seed)
    assert len(ratios) == 3 and abs(sum(ratios) - 1.0) < 1e-6, "Ratios must be three numbers summing to 1.0"
    for label in ["true", "false"]:
        src_dir = os.path.join(base_dir, label)
        if not os.path.isdir(src_dir):
            continue
        files = sorted(os.listdir(src_dir))
        random.shuffle(files)
        n = len(files)
        n_train = int(n * ratios[0])
        n_val = int(n * ratios[1])
        n_test = n - n_train - n_val
        split_map = {
            "train": files[:n_train],
            "val": files[n_train:n_train+n_val],
            "test": files[n_train+n_val:]
        }
        for split, split_files in split_map.items():
            split_dir = os.path.join(base_dir, split, label)
            os.makedirs(split_dir, exist_ok=True)
            for f in split_files:
                shutil.move(os.path.join(src_dir, f), os.path.join(split_dir, f))
        shutil.rmtree(src_dir)

def main():
    parser = argparse.ArgumentParser(description="Classify variant tensors against a VCF file using pysam.")
    parser.add_argument("tensor_folder_path")
    parser.add_argument("vcf_file")
    parser.add_argument("node_pos_json")
    parser.add_argument("--output_folder", default="./classification_results")
    parser.add_argument("--chr", default="chr1")
    parser.add_argument("--use-symlinks", action="store_true")
    parser.add_argument("-j", "--workers", type=int, default=cpu_count())
    parser.add_argument("--organize", nargs=3, type=float, metavar=('TRAIN', 'VAL', 'TEST'), default=None,
                        help="Split ratio for organizing true/false data (e.g., 0.6 0.2 0.2)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    args = parser.parse_args()

    os.makedirs(args.output_folder, exist_ok=True)
    true_dir = os.path.join(args.output_folder, "true")
    false_dir = os.path.join(args.output_folder, "false")
    os.makedirs(true_dir, exist_ok=True)
    os.makedirs(false_dir, exist_ok=True)

    total_variants = get_variant_count(args.vcf_file, args.chr)
    if total_variants >= 0:
        print(f"Found a total of {total_variants:,} unique variants to be used as ground truth.\n")
    else:
        print("Could not count variants due to an error.\n")

    node_dirs = [d.path for d in os.scandir(args.tensor_folder_path) if d.is_dir()]
    total_nodes = len(node_dirs)
    if total_nodes == 0:
        print("No node directories found to process. Exiting.")
        return

    print(f"Found {total_nodes} nodes. Initializing {args.workers} worker processes...")
    print("(This may take a moment before the progress bar appears...)\n")

    worker_func = partial(process_node_directory, true_dir=true_dir, false_dir=false_dir, use_symlinks=args.use_symlinks)
    all_summary_records = []
    total_true = 0
    total_false = 0
    start_time = time.monotonic()

    init_args = (args.vcf_file, args.chr, args.node_pos_json)
    with Pool(processes=args.workers, initializer=init_worker, initargs=init_args) as pool:
        results = pool.imap_unordered(worker_func, node_dirs)
        for i, (records, true_c, false_c) in enumerate(results):
            processed_count = i + 1
            all_summary_records.extend(records)
            total_true += true_c
            total_false += false_c
            elapsed_time = time.monotonic() - start_time
            progress = processed_count / total_nodes
            rate = processed_count / elapsed_time if elapsed_time > 0 else 0
            eta_str = format_time((total_nodes - processed_count) / rate) if rate > 0 else "..."
            elapsed_str = format_time(elapsed_time)
            print(
                f"\rProgress: {processed_count}/{total_nodes} ({progress:.1%}) | "
                f"Elapsed: {elapsed_str} | ETA: {eta_str} | Rate: {rate:.2f} nodes/s  ", end=""
            )

    print("\n\nProcessing complete.")
    summary_json_path = os.path.join(args.output_folder, "classification_summary.json")
    print(f"Writing classification summary to {summary_json_path}...")
    with open(summary_json_path, 'w') as f:
        json.dump({"chromosome": args.chr, "results": all_summary_records}, f, indent=2)

    print("\n--- Classification Summary ---")
    print(f"Total nodes processed: {total_nodes}")
    print(f"Total tensors classified: {total_true + total_false}")
    print(f"True (matches in VCF): {total_true}")
    print(f"False (not in VCF): {total_false}")

    if args.organize:
        print("\nOrganizing classified data into train/val/test splits...")
        organize_classified_data(args.output_folder, args.organize, seed=args.seed)
        print("Organization complete.")

if __name__ == "__main__":
    if sys.version_info < (3, 6):
        sys.exit("This script requires Python 3.6+.")
    main()
