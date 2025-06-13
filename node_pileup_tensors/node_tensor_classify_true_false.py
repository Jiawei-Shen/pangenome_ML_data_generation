#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import shutil
import sys

def query_vcf_for_chromosome(vcf_file_path, chromosome):
    print(f"Querying VCF file '{vcf_file_path}' for '{chromosome}' variants using 'bcftools'...")
    chromosome_variants = set()
    cmd = ['bcftools', 'view', '-r', chromosome, vcf_file_path]
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print(f"Error: bcftools failed for chromosome '{chromosome}'", file=sys.stderr)
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
        print(f"Found {len(chromosome_variants)} variants on '{chromosome}' in VCF.")
        return chromosome_variants
    except FileNotFoundError:
        print(f"Error: 'bcftools' not found.", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return None

def load_node_positions(json_file_path):
    print(f"Loading node positions from '{json_file_path}'...")
    try:
        with open(json_file_path, 'r') as f:
            data = json.load(f)
        nodes = data.get("nodes", [])
        return {
            str(entry["node_id"]): int(entry["grch38_position_start"])
            for entry in nodes
            if "node_id" in entry and "grch38_position_start" in entry
        }
    except Exception as e:
        print(f"Error loading node positions: {e}", file=sys.stderr)
        return {}

def calculate_genomic_position(start_pos, offset):
    return start_pos + offset + 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tensor_folder_path")
    parser.add_argument("vcf_file")
    parser.add_argument("node_pos_json")
    parser.add_argument("--output_folder", default="./classification_results")
    parser.add_argument("--chr", default="chr1")
    args = parser.parse_args()

    os.makedirs(os.path.join(args.output_folder, "true"), exist_ok=True)
    os.makedirs(os.path.join(args.output_folder, "false"), exist_ok=True)

    vcf_variants = query_vcf_for_chromosome(args.vcf_file, args.chr)
    if vcf_variants is None:
        sys.exit(1)

    node_positions = load_node_positions(args.node_pos_json)

    summary = []
    for node_id in os.listdir(args.tensor_folder_path):
        node_dir = os.path.join(args.tensor_folder_path, node_id)
        if not os.path.isdir(node_dir):
            continue

        summary_path = os.path.join(node_dir, "variant_summary.json")
        if not os.path.isfile(summary_path):
            continue

        try:
            with open(summary_path, 'r') as f:
                summary_json = json.load(f)
        except Exception as e:
            print(f"Warning: Failed to load {summary_path}: {e}", file=sys.stderr)
            continue

        start_pos = node_positions.get(node_id)
        for variant in summary_json.get("variants_passing_af_filter", []):
            tensor_file = variant.get("tensor_file")
            variant_key = variant.get("variant_key")
            if not tensor_file or not variant_key or not tensor_file.endswith(".npy"):
                continue
            tensor_path = os.path.join(node_dir, tensor_file)
            if not os.path.isfile(tensor_path):
                continue

            try:
                parts = variant_key.split('_')
                pos = int(parts[0])
                ref = parts[2].upper()
                alt = parts[3].upper()
            except Exception:
                continue

            if start_pos is None:
                continue
            gpos = calculate_genomic_position(start_pos, pos)
            match = (gpos, ref, alt) in vcf_variants
            dest = os.path.join(args.output_folder, "true" if match else "false", f"{node_id}_{tensor_file}")
            shutil.copy2(tensor_path, dest)

            summary.append({
                "node_id": node_id,
                "tensor_file": tensor_file,
                "variant_key": variant_key,
                "genomic_position": gpos,
                "ref": ref,
                "alt": alt,
                "classification": "true" if match else "false"
            })

    with open(os.path.join(args.output_folder, "classification_summary.json"), 'w') as f:
        json.dump({"chromosome": args.chr, "results": summary}, f, indent=2)

    print(f"Classification complete. Summary written to {args.output_folder}/classification_summary.json")

if __name__ == "__main__":
    main()