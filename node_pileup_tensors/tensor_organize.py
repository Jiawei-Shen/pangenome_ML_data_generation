#!/usr/bin/env python3
import argparse
import json
import os
import shutil
import sys
import time
import random
from multiprocessing import Pool, cpu_count
import pysam

node_positions = {}

copy_counter = 0


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
    except Exception as e:
        print(f"Error loading VCF: {e}", file=sys.stderr)
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
    except Exception as e:
        print(f"Failed to read VCF: {e}", file=sys.stderr)
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


def process_node_directory(args):
    global copy_counter
    node_dir, true_dir, false_dir, use_symlinks, vcf_variants = args
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
                copy_counter += 1
                if copy_counter % 10000 == 0:
                    print(f"Copied {copy_counter} items so far...")
        except (FileExistsError, OSError):
            pass
        summary_records.append({
            "node_id": node_id, "tensor_file": tensor_file,
            "variant_key": variant_key, "genomic_position": grch38_pos,
            "ref": ref, "alt": alt,
            "classification": "true" if is_match else "false"
        })
    return summary_records, local_true_count, local_false_count

# ... (rest of the script remains unchanged)
