#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import mmap
import numpy as np
from collections import defaultdict

# ─────────────────────────────────────────────────────────────────────────────
# Globals for file access
GLOBAL_DAT = None            # memory-mapped .dat file
GLOBAL_NODE_SEQS = None      # dict of node_id → sequence
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# ─────────────────────────────────────────────────────────────────────────────
def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

# ─────────────────────────────────────────────────────────────────────────────
def parse_idx_file(idx_path):
    node_index = {}
    with open(idx_path, 'rb') as f:
        num_nodes = struct.unpack('<I', f.read(4))[0]
        for _ in range(num_nodes):
            node_id, offset, block_size, n_records, _ = struct.unpack('<I Q I I H', f.read(22))
            node_index[node_id] = (offset, n_records)
    print(f"✔ Parsed {len(node_index)} nodes from index file.")
    return node_index

# ─────────────────────────────────────────────────────────────────────────────
def load_node_sequences_from_gfa(gfa_path, target_node_ids):
    node_sequences = {}
    line_counter = 0
    with open(gfa_path, 'r') as f:
        for line in f:
            line_counter += 1
            if line_counter % 10_000_000 == 0:
                print(f"✔ Checked {line_counter} lines in GFA file")
            if not line.startswith('S\t'):
                continue
            parts = line.split('\t', 2)
            try:
                nid = int(parts[1])
            except ValueError:
                continue
            if nid in target_node_ids:
                node_sequences[nid] = parts[2].rstrip('\n')
    print(f"✔ Loaded {len(node_sequences)} sequences from GFA.")
    return node_sequences

# ─────────────────────────────────────────────────────────────────────────────
def decode_cigar(cigar_string):
    import re
    return re.findall(r'(\d+)([MXDI])', cigar_string)

# ─────────────────────────────────────────────────────────────────────────────
def detect_variants_from_cigar(offset, cigar_string):
    variants = []
    pos = offset
    for length_str, op in decode_cigar(cigar_string):
        length = int(length_str)
        if op in ('X','I','D'):
            for _ in range(length):
                if op == 'I':
                    variants.append((pos-1, 'I'))
                elif op == 'D':
                    variants.append((pos, 'D'))
                    pos += 1
                else:
                    variants.append((pos, 'X'))
                    pos += 1
        else:
            pos += length
    return variants

# ─────────────────────────────────────────────────────────────────────────────
def process_node(node_id, offset, n_records):
    """Extract reads for node_id and build variant-centered pileups."""
    sequence = GLOBAL_NODE_SEQS[node_id]
    node_len = len(sequence)
    reads_by_variant = defaultdict(list)
    # skip block header of 10 bytes
    ptr = offset + 10
    for _ in range(n_records):
        data = GLOBAL_DAT[ptr: ptr + RECORD_SIZE]
        ptr += RECORD_SIZE
        off, raw_seq, raw_bq, raw_cigar, rq, strand = RECORD_STRUCT.unpack(data)
        seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
        cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
        strand_char = strand.decode()
        length = len(seq)
        if strand_char == '-':
            off = node_len - off - length
            seq = seq.translate(str.maketrans('ACGTacgtNn','TGCAtgcaNn'))[::-1]
        for vpos, vtype in detect_variants_from_cigar(off, cigar):
            if 0 <= vpos < node_len:
                reads_by_variant[(vpos, vtype)].append((off, seq))
    pileups = {}
    window = 60
    half = window // 2
    for (vpos, vtype), reads in reads_by_variant.items():
        mat = np.full((len(reads), window), 4, dtype=np.uint8)
        for i, (roff, rseq) in enumerate(reads):
            start = vpos - roff - half
            for j in range(window):
                idx = start + j
                if 0 <= idx < len(rseq):
                    mat[i, j] = BASE_TO_INDEX.get(rseq[idx].upper(), 4)
        pileups[f"{vpos}_{vtype}"] = mat.tolist()
    return node_id, pileups

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Single-threaded variant-centered pileup")
    parser.add_argument("dat", help=".dat file path")
    parser.add_argument("idx", help=".idx file path")
    parser.add_argument("output", help="JSON output path")
    parser.add_argument("--gfa", help="GFA path if no cache")
    parser.add_argument("--load-cache", help="JSON file of node sequences")
    parser.add_argument("--save-cache", help="Where to write node-sequence cache")
    args = parser.parse_args()

    # Validate
    if not args.load_cache and not args.gfa:
        sys.exit("❌ Provide --gfa or --load-cache.")
    if not args.load_cache and not args.save_cache:
        sys.exit("❌ Provide --save-cache when building from GFA.")

    print("🔹 Parsing index file")
    node_index = parse_idx_file(args.idx)

    # Load or build cache
    cache_path = args.load_cache or args.save_cache
    if args.load_cache and os.path.isfile(cache_path):
        print(f"✔ Loading node sequences from {cache_path}")
        with open(cache_path) as cf:
            data = json.load(cf)
        node_sequences = {int(k): v for k, v in data.items()}
    else:
        print("🔹 Building node sequences from GFA")
        node_sequences = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        with open(args.save_cache, 'w') as cf:
            json.dump({str(k): v for k, v in node_sequences.items()}, cf)
        print(f"✔ Saved node-sequence cache to {args.save-cache}")

    # Initialize memory-mapped .dat and sequences
    global GLOBAL_DAT, GLOBAL_NODE_SEQS
    GLOBAL_NODE_SEQS = node_sequences
    dat_file = open(args.dat, 'rb')
    GLOBAL_DAT = mmap.mmap(dat_file.fileno(), 0, access=mmap.ACCESS_READ)

    # Single-threaded processing
    print("🔹 Processing nodes in single-thread mode")
    results = {}
    start = time.time()
    total = len(node_index)
    for count, (node_id, (offset, nrec)) in enumerate(node_index.items(), start=1):
        if node_id not in node_sequences:
            continue
        print(count, node_id, offset, nrec, "\n")
        _, pileup = process_node(node_id, offset, nrec)
        results[node_id] = pileup
        if count % 1000 == 0 or count == total:
            elapsed = time.time() - start
            print(f"✔ {count}/{total} nodes processed — elapsed {elapsed:.2f}s")

    # Write output
    print("🔹 Writing JSON output")
    with open(args.output, 'w') as out_f:
        json.dump(results, out_f, indent=2)
    print(f"✅ Done. Wrote output to {args.output}")

if __name__ == '__main__':
    main()
