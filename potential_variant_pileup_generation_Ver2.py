#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def parse_idx_file(idx_path):
    node_index = {}
    with open(idx_path, "rb") as file:
        num_nodes = struct.unpack("<I", file.read(4))[0]
        for _ in range(num_nodes):
            node_id, offset, block_size, n_records, _ = struct.unpack("<I Q I I H", file.read(22))
            node_index[node_id] = (offset, n_records)
    print(f"✔ Parsed {len(node_index)} nodes from index file.")
    return node_index

def load_node_sequences_from_gfa(gfa_path, target_node_ids):
    node_sequences = {}
    line_count = 0
    with open(gfa_path, "r") as file:
        for line in file:
            line_count += 1
            if line_count % 10_000_000 == 0:
                print(f"✔ Checked {line_count} lines in GFA")
            if not line.startswith("S\t"):
                continue
            parts = line.split("\t", 2)
            try:
                node_id = int(parts[1])
            except:
                continue
            if node_id in target_node_ids:
                node_sequences[node_id] = parts[2].rstrip("\n")
    print(f"✔ Loaded {len(node_sequences)} sequences from GFA.")
    return node_sequences

def decode_cigar(cigar_string):
    import re
    return re.findall(r'(\d+)([MXDI])', cigar_string)

def detect_variants_from_cigar(offset, cigar_string):
    variants = []
    position = offset
    for length_str, op in decode_cigar(cigar_string):
        length = int(length_str)
        if op in ("X", "I", "D"):
            for _ in range(length):
                if op == "I":
                    variants.append((position - 1, "I"))
                elif op == "D":
                    variants.append((position, "D"))
                    position += 1
                else:
                    variants.append((position, "X"))
                    position += 1
        else:
            position += length
    return variants

def process_node(task_info):
    node_id, node_sequence, dat_path, dat_offset, n_records = task_info
    node_length = len(node_sequence)
    reads_by_variant = defaultdict(list)

    with open(dat_path, "rb") as file:
        file.seek(dat_offset)
        file.read(10)  # skip block header

        for _ in range(n_records):
            data = file.read(RECORD_STRUCT.size)
            offset, seq_raw, bq_raw, cigar_raw, read_quality, strand = RECORD_STRUCT.unpack(data)
            sequence = seq_raw.rstrip(b'\x00').decode('ascii', errors='ignore')
            cigar_string = cigar_raw.rstrip(b'\x00').decode('ascii', errors='ignore')
            strand = strand.decode()
            read_length = len(sequence)

            if strand == '-':
                offset = node_length - offset - read_length
                sequence = reverse_complement(sequence)

            for variant_position, variant_type in detect_variants_from_cigar(offset, cigar_string):
                if 0 <= variant_position < node_length:
                    reads_by_variant[(variant_position, variant_type)].append((offset, sequence))

    pileups = {}
    window_size = 60
    half_window = window_size // 2
    for (variant_position, variant_type), reads in reads_by_variant.items():
        pileup_array = np.full((len(reads), window_size), 4, dtype=np.uint8)
        for row_index, (read_offset, sequence) in enumerate(reads):
            relative_start = variant_position - read_offset - half_window
            for column_index in range(window_size):
                seq_index = relative_start + column_index
                if 0 <= seq_index < len(sequence):
                    pileup_array[row_index, column_index] = BASE_TO_INDEX.get(sequence[seq_index].upper(), 4)
        pileups[f"{variant_position}_{variant_type}"] = pileup_array.tolist()

    print(f"Node: {node_id}, pileups")
    return node_id, pileups

def main():
    parser = argparse.ArgumentParser(description="Variant-centered pileup from .dat/.idx + GFA")
    parser.add_argument("dat")
    parser.add_argument("idx")
    parser.add_argument("output", help="JSON output path")
    parser.add_argument("--gfa", help="GFA path (required if no --load-cache)")
    parser.add_argument("--load-cache", help="JSON cache of node sequences")
    parser.add_argument("--save-cache", help="Where to write sequence cache")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    if not args.load_cache and not args.gfa:
        sys.exit("❌ Error: --gfa required if --load-cache is not provided.")
    if args.save_cache and not args.gfa:
        sys.exit("❌ Error: --save-cache requires --gfa.")

    print("🔹 Step 1: Parse .idx")
    node_index = parse_idx_file(args.idx)

    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"✔ Loading node sequences from cache: {args.load_cache}")
        with open(args.load_cache) as file:
            cached = json.load(file)
            node_sequences = {int(k): v for k, v in cached.items()}
    else:
        print("🔹 Step 2: Parse GFA file")
        node_sequences = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        if args.save_cache:
            print(f"✔ Saving cache to: {args.save_cache}")
            with open(args.save_cache, "w") as cache_file:
                json.dump({str(k): v for k, v in node_sequences.items()}, cache_file)

    print("🔹 Step 3: Pileup with multiprocessing")
    task_list = [
        (node_id, node_sequences[node_id], args.dat, offset, record_count)
        for node_id, (offset, record_count) in node_index.items()
        if node_id in node_sequences
    ]

    final_output = {}
    start_time = time.time()
    milestone = 1_000
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for i, (node_id, pileup) in enumerate(executor.map(process_node, task_list), 1):
            final_output[node_id] = pileup
            if i > milestone:
                elapsed = time.time() - start_time
                print(f"✔ Processed {i} nodes — total time: {elapsed:.2f} sec")
                milestone += milestone

    print("🔹 Step 4: Write JSON output")
    with open(args.output, "w") as out_file:
        json.dump(final_output, out_file, indent=2)

    print("✅ Done. Output written to", args.output)

if __name__ == "__main__":
    main()
