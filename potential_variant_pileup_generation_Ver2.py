#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size

BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
INDEX_TO_BASE = ['A', 'C', 'G', 'T', 'N']

# ─────────────────────────────────────────────────────────────────────────────
def reverse_complement(seq):
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement)[::-1]

# ─────────────────────────────────────────────────────────────────────────────
def parse_idx_file(idx_path):
    node_ids = set()
    with open(idx_path, "rb") as f:
        num_nodes = struct.unpack("<I", f.read(4))[0]
        for _ in range(num_nodes):
            node_id, _, _, _, _ = struct.unpack("<I Q I I H", f.read(22))
            node_ids.add(node_id)
    print(f"✔ Loaded {len(node_ids)} node IDs from index file.")
    return node_ids

# ─────────────────────────────────────────────────────────────────────────────
def load_node_sequences_from_gfa(gfa_path, target_node_ids):
    node_sequences = {}
    line_count = 0
    with open(gfa_path, "r") as f:
        for line in f:
            line_count += 1
            if line_count % 10_000_000 == 0:
                print(f"✔ Checked {line_count} lines in GFA file")
            if not line.startswith("S\t"):
                continue
            fields = line.split('\t', 2)
            try:
                node_id = int(fields[1])
            except (IndexError, ValueError):
                continue
            if node_id in target_node_ids:
                node_sequences[node_id] = fields[2].rstrip("\n")
    print(f"✔ Loaded sequences for {len(node_sequences)} nodes from GFA.")
    return node_sequences

# ─────────────────────────────────────────────────────────────────────────────
def decode_cigar(cigar_str):
    import re
    return re.findall(r'(\d+)([MXDI])', cigar_str)

def detect_variants_from_cigar(offset, cigar_str):
    variants = []
    ref_pos = offset
    for length, op in decode_cigar(cigar_str):
        length = int(length)
        if op in ('X', 'I', 'D'):
            for _ in range(length):
                if op == 'I':
                    variants.append((ref_pos - 1, 'I'))
                elif op == 'D':
                    variants.append((ref_pos, 'D'))
                    ref_pos += 1
                else:
                    variants.append((ref_pos, 'X'))
                    ref_pos += 1
        else:
            ref_pos += length
    return variants

# ─────────────────────────────────────────────────────────────────────────────
def process_node(args):
    node_id, node_sequence, dat_path = args
    print(f"▶ Processing node {node_id} (length: {len(node_sequence)} bp)...")
    node_len = len(node_sequence)
    reads_by_variant = defaultdict(list)

    with open(dat_path, "rb") as f:
        # skip .dat header
        f.seek(0, 2)
        total_size = f.tell()
        f.seek(0)
        f.read(10)

        while f.tell() < total_size:
            header = f.read(10)
            if len(header) < 10:
                break
            this_node_id, _, n_records = struct.unpack("<I I H", header)
            for _ in range(n_records):
                data = f.read(RECORD_SIZE)
                offset_val, seq_raw, bq_raw, cigar_raw, rq, strand = RECORD_STRUCT.unpack(data)
                if this_node_id != node_id:
                    continue
                seq = seq_raw.rstrip(b'\x00').decode('ascii', errors='ignore')
                cigar = cigar_raw.rstrip(b'\x00').decode('ascii', errors='ignore')
                strand = strand.decode()
                read_len = len(seq)

                if strand == '-':
                    offset_val = node_len - offset_val - read_len
                    seq = reverse_complement(seq)

                for vpos, vtype in detect_variants_from_cigar(offset_val, cigar):
                    if 0 <= vpos < node_len:
                        reads_by_variant[(vpos, vtype)].append((offset_val, seq))

    # build 2D pileup per variant
    pileups = {}
    window_size = 60
    half = window_size // 2
    for (vpos, vtype), reads in reads_by_variant.items():
        mat = np.full((len(reads), window_size), 4, dtype=np.uint8)
        for i, (offset_val, seq) in enumerate(reads):
            start = vpos - offset_val - half
            for j in range(window_size):
                idx = start + j
                if 0 <= idx < len(seq):
                    mat[i, j] = BASE_TO_INDEX.get(seq[idx].upper(), 4)
        pileups[f"{vpos}_{vtype}"] = mat.tolist()

    print(f"✔ Node {node_id}: {len(pileups)} variant sites processed.")
    return node_id, pileups

# ─────────────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description="Pileup .dat/.idx reads centered on variants")
    p.add_argument("dat", help="Path to .dat file")
    p.add_argument("idx", help="Path to .idx file")
    p.add_argument("output", help="Path to output JSON file")
    p.add_argument("--gfa",    help="Path to GFA file (required if not loading cache)")
    p.add_argument("--load-cache", help="JSON cache of node sequences to load")
    p.add_argument("--save-cache", help="Write JSON cache of node sequences after parsing GFA")
    p.add_argument("--workers", type=int, default=2, help="Parallel processes")
    args = p.parse_args()

    print("🔹 Step 1: Reading index file")
    node_ids = parse_idx_file(args.idx)

    # require GFA if no cache
    if not args.load_cache and not args.gfa:
        sys.exit("Error: --gfa is required if --load-cache is not provided.")
    if args.save_cache and not args.gfa:
        sys.exit("Error: --save-cache requires --gfa to be specified.")

    # load or build node_sequences
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"✔ Loading node sequences from cache: {args.load_cache}")
        cached = json.load(open(args.load_cache))
        node_sequences = {int(k): v for k, v in cached.items()}
        print(f"✔ Loaded {len(node_sequences)} sequences from cache.")
    else:
        print("🔹 Step 2: Loading node sequences from GFA")
        node_sequences = load_node_sequences_from_gfa(args.gfa, node_ids)
        if args.save_cache:
            print(f"✔ Caching sequences to: {args.save_cache}")
            with open(args.save_cache, "w") as cf:
                json.dump({str(k): v for k, v in node_sequences.items()}, cf)
            print("✔ Cache written.")

    print("🔹 Step 3: Detecting variants and piling up")
    tasks = [(nid, node_sequences[nid], args.dat)
             for nid in node_ids if nid in node_sequences]

    final_pileup = {}
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for node_id, pileups in ex.map(process_node, tasks):
            final_pileup[node_id] = pileups

    print(f"🔹 Step 4: Writing output to {args.output}")
    with open(args.output, "w") as out:
        json.dump(final_pileup, out, indent=2)
    print("✅ Finished writing JSON output.")

if __name__ == "__main__":
    main()
