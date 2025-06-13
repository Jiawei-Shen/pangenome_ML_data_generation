#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import shutil
import sys

def query_vcf_for_chromosome(vcf_file_path, chromosome):
    print(f"Querying VCF for {chromosome} using bcftools...")
    variants = set()
    cmd = ['bcftools', 'view', '-r', chromosome, vcf_file_path]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in proc.stdout.splitlines():
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                pos = int(fields[1])
                ref = fields[3].upper()
                for alt in fields[4].upper().split(','):
                    variants.add((pos, ref, alt))
        print(f"Loaded {len(variants)} variants from VCF")
        return variants
    except subprocess.CalledProcessError as e:
        print(f"Error running bcftools: {e.stderr}", file=sys.stderr)
        return set()

def load_node_positions(json_path):
    print(f"Loading node position data from {json_path}...")
    with open(json_path) as f:
        data = json.load(f)

    node_dict = {}
    for node in data.get("nodes", []):
        node_id = str(node.get("node_id"))
        pos = node.get("grch38_position_start")
        if node_id and isinstance(pos, int):
            node_dict[node_id] = pos
    print(f"Loaded {len(node_dict)} node positions")
    return node_dict

def calculate_genomic_position(start, offset):
    return start + offset + 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tensor_folder_path")
    parser.add_argument("vcf_file")
    parser.add_argument("node_pos_json")
    parser.add_argument("--output_folder", default="./classification_results")
    parser.add_argument("--chr", default="chr1")
    args = parser.parse_args()

    os.makedirs(args.output_folder, exist_ok=True)
    true_dir = os.path.join(args.output_folder, "true")
    false_dir = os.path.join(args.output_folder, "false")
    os.makedirs(true_dir, exist_ok=True)
    os.makedirs(false_dir, exist_ok=True)

    vcf_variants = query_vcf_for_chromosome(args.vcf_file, args.chr)
    node_positions = load_node_positions(args.node_pos_json)

    total = true_count = false_count = 0
    node_count = 0
    summary_records = []

    for node_id in os.listdir(args.tensor_folder_path):
        node_dir = os.path.join(args.tensor_folder_path, node_id)
        if not os.path.isdir(node_dir):
            continue

        summary_path = os.path.join(node_dir, "variant_summary.json")
        if not os.path.isfile(summary_path):
            continue

        node_count += 1
        print(f"Processing node {node_id} ({node_count})...")

        try:
            with open(summary_path) as f:
                summary = json.load(f)
        except Exception as e:
            print(f"Could not load {summary_path}: {e}", file=sys.stderr)
            continue

        variants = summary.get("variants_passing_af_filter", [])
        for variant in variants:
            total += 1
            tensor_file = variant.get("tensor_file")
            variant_key = variant.get("variant_key")

            if not tensor_file or not variant_key or not tensor_file.endswith(".npy"):
                continue

            tensor_path = os.path.join(node_dir, tensor_file)
            if not os.path.isfile(tensor_path):
                continue

            try:
                parts = variant_key.split("_")
                offset = int(parts[0])
                ref = parts[2].upper()
                alt = parts[3].upper()
            except Exception:
                continue

            start_pos = node_positions.get(node_id)
            if start_pos is None:
                continue

            grch38_pos = calculate_genomic_position(start_pos, offset)
            match = (grch38_pos, ref, alt) in vcf_variants

            dest_dir = true_dir if match else false_dir
            shutil.copy2(tensor_path, os.path.join(dest_dir, f"{node_id}_{tensor_file}"))
            if match:
                true_count += 1
            else:
                false_count += 1

            summary_records.append({
                "node_id": node_id,
                "tensor_file": tensor_file,
                "variant_key": variant_key,
                "genomic_position": grch38_pos,
                "ref": ref,
                "alt": alt,
                "classification": "true" if match else "false"
            })

    summary_json_path = os.path.join(args.output_folder, "classification_summary.json")
    with open(summary_json_path, 'w') as f:
        json.dump({"chromosome": args.chr, "results": summary_records}, f, indent=2)
    print(f"Classification summary written to {summary_json_path}")

    print("\n--- Classification Summary ---")
    print(f"Total nodes processed: {node_count}")
    print(f"Total tensors processed: {total}")
    print(f"True: {true_count}")
    print(f"False: {false_count}")

if __name__ == "__main__":
    if sys.version_info < (3, 6):
        sys.exit("Requires Python 3.6+")
    main()