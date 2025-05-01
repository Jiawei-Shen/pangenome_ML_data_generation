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
from concurrent.futures import ProcessPoolExecutor, as_completed

# ─────────────────────────────────────────────────────────────────────────────
# Globals for file access (initialized in worker processes)
GLOBAL_DAT = None            # memory-mapped .dat file
DAT_FILE_HANDLE = None       # file handle for .dat file to keep it open
GLOBAL_NODE_SEQS = None      # dict of node_id → sequence

RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# ─────────────────────────────────────────────────────────────────────────────
# Helper functions (inherited by worker processes)

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

# parse index file (only in main process)
def parse_idx_file(idx_path):
    node_index = {}
    with open(idx_path, 'rb') as f:
        num_nodes = struct.unpack('<I', f.read(4))[0]
        for _ in range(num_nodes):
            node_id, offset, block_size, n_records, _ = struct.unpack('<I Q I I H', f.read(22))
            node_index[node_id] = (offset, n_records)
    print(f"✔ Parsed {len(node_index)} nodes from index file.")
    return node_index

# load node sequences (only in main process)
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

# CIGAR parsing and variant detection (in worker processes)
import re

def decode_cigar(cigar_string):
    return re.findall(r'(\d+)([MXDI])', cigar_string)

def detect_variants_from_cigar(offset, cigar_string):
    variants = []
    pos = offset
    for length_str, op in decode_cigar(cigar_string):
        length = int(length_str)
        if op in ('X', 'I', 'D'):
            for _ in range(length):
                if op == 'I':
                    variants.append((pos - 1, 'I'))
                elif op == 'D':
                    variants.append((pos, 'D'))
                    pos += 1
                else:
                    variants.append((pos, 'X'))
                    pos += 1
        else:
            pos += length
    return variants

# worker initializer: memory-map .dat file

def init_worker(dat_path):
    global GLOBAL_DAT, DAT_FILE_HANDLE
    DAT_FILE_HANDLE = open(dat_path, 'rb')
    GLOBAL_DAT = mmap.mmap(DAT_FILE_HANDLE.fileno(), 0, access=mmap.ACCESS_READ)

# process a single node (in worker processes)

def process_node(node_id, offset, n_records):
    """Extract reads and build variant-centered pileups (MAPQ>=10 filter)."""
    sequence = GLOBAL_NODE_SEQS[node_id]
    node_len = len(sequence)
    # read and filter segments
    segments = []
    ptr = offset + 10  # skip block header
    for _ in range(n_records):
        data = GLOBAL_DAT[ptr:ptr + RECORD_SIZE]
        ptr += RECORD_SIZE
        off, raw_seq, raw_bq, raw_cigar, mapq, strand = RECORD_STRUCT.unpack(data)
        if mapq < 10:
            continue
        seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
        cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
        strand_char = strand.decode()
        read_len = len(seq)
        if strand_char == '-':
            off = node_len - off - read_len
            seq = seq.translate(str.maketrans('ACGTacgtNn','TGCAtgcaNn'))[::-1]
        segments.append((off, seq, cigar))
    # group reads by variant
    reads_by_variant = defaultdict(list)
    for off, seq, cigar in segments:
        for vpos, vtype in detect_variants_from_cigar(off, cigar):
            if 0 <= vpos < node_len:
                reads_by_variant[(vpos, vtype)].append((off, seq))
    # build pileup matrices
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
# Main entry point

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parallel variant-centered pileup")
    parser.add_argument("dat", help=".dat file path")
    parser.add_argument("idx", help=".idx file path")
    parser.add_argument("output", help="JSON output path")
    parser.add_argument("--gfa", help="GFA path if no cache")
    parser.add_argument("--load-cache", help="JSON file of node sequences")
    parser.add_argument("--save-cache", help="Where to write node-sequence cache")
    parser.add_argument("--workers", type=int, default=os.cpu_count(),
                        help="Number of parallel worker processes (default: CPU count)")
    args = parser.parse_args()

    if not args.load_cache and not args.gfa:
        sys.exit("❌ Provide --gfa or --load-cache.")
    if not args.load_cache and not args.save_cache:
        sys.exit("❌ Provide --save-cache when building from GFA.")

    print("🔹 Parsing index file")
    node_index = parse_idx_file(args.idx)

    cache_path = args.load_cache or args.save_cache
    if args.load_cache and os.path.isfile(cache_path):
        print(f"✔ Loading node sequences from {cache_path}")
        with open(cache_path) as cf:
            data = json.load(cf)
        GLOBAL_NODE_SEQS = {int(k): v for k, v in data.items()}
    else:
        print("🔹 Building node sequences from GFA")
        GLOBAL_NODE_SEQS = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        with open(args.save_cache, 'w') as cf:
            json.dump({str(k): v for k, v in GLOBAL_NODE_SEQS.items()}, cf)
        print(f"✔ Saved node-sequence cache to {args.save_cache}")

    print(f"🔹 Processing {len(GLOBAL_NODE_SEQS)} nodes in parallel with {args.workers} workers")
    results = {}
    start = time.time()
    total = len(GLOBAL_NODE_SEQS)

    with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker, initargs=(args.dat,)) as executor:
        futures = {executor.submit(process_node, nid, off, nrec): nid
                   for nid, (off, nrec) in node_index.items() if nid in GLOBAL_NODE_SEQS}
        for count, future in enumerate(as_completed(futures), start=1):
            nid = futures[future]
            try:
                _, pileup = future.result()
                results[nid] = pileup
            except Exception as e:
                print(f"❌ Error processing node {nid}: {e}", file=sys.stderr)
            if count % 1000 == 0 or count == total:
                elapsed = time.time() - start
                print(f"✔ {count}/{total} nodes processed — elapsed {elapsed:.2f}s")

    print("🔹 Writing JSON output")
    with open(args.output, 'w') as out_f:
        json.dump(results, out_f, indent=2)
    print(f"✅ Done. Wrote output to {args.output}")
