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
# Globals for worker processes
GLOBAL_DAT = None
GLOBAL_NODE_SEQS = None
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# ─────────────────────────────────────────────────────────────────────────────
def _init_worker(cache_path, dat_path):
    """
    Initialize each worker: load JSON cache into GLOBAL_NODE_SEQS and mmap the .dat file.
    """
    global GLOBAL_DAT, GLOBAL_NODE_SEQS
    # Load node sequences from JSON cache
    with open(cache_path, 'r') as cf:
        seq_cache = json.load(cf)
    GLOBAL_NODE_SEQS = {int(k): v for k, v in seq_cache.items()}
    # Memory-map the .dat file once
    dat_fh = open(dat_path, 'rb')
    GLOBAL_DAT = mmap.mmap(dat_fh.fileno(), 0, access=mmap.ACCESS_READ)

# ─────────────────────────────────────────────────────────────────────────────
def parse_idx_file(idx_path):
    """Return dict of node_id -> (offset, n_records)"""
    node_index = {}
    with open(idx_path, 'rb') as f:
        count = struct.unpack('<I', f.read(4))[0]
        for _ in range(count):
            node_id, offset, block_size, n_records, _ = struct.unpack('<I Q I I H', f.read(22))
            node_index[node_id] = (offset, n_records)
    print(f"✔ Parsed {len(node_index)} nodes from index file.")
    return node_index

# ─────────────────────────────────────────────────────────────────────────────
def load_node_sequences_from_gfa(gfa_path, target_node_ids):
    """Parse GFA and extract sequences for target nodes."""
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
    seq = GLOBAL_NODE_SEQS[node_id]
    node_len = len(seq)
    reads_by_variant = defaultdict(list)
    base = offset + 10  # skip block header
    for _ in range(n_records):
        data = GLOBAL_DAT[base:base+RECORD_SIZE]
        base += RECORD_SIZE
        off, raw_seq, raw_bq, raw_cigar, rq, strand = RECORD_STRUCT.unpack(data)
        sequence = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
        cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
        strand_char = strand.decode()
        length = len(sequence)
        if strand_char == '-':
            off = node_len - off - length
            sequence = sequence.translate(str.maketrans('ACGTacgtNn','TGCAtgcaNn'))[::-1]
        for vpos, vtype in detect_variants_from_cigar(off, cigar):
            if 0 <= vpos < node_len:
                reads_by_variant[(vpos, vtype)].append((off, sequence))
    pileups = {}
    window = 60
    half = window//2
    for (vpos, vtype), reads in reads_by_variant.items():
        mat = np.full((len(reads), window), 4, dtype=np.uint8)
        for i, (roff, rseq) in enumerate(reads):
            start = vpos - roff - half
            for j in range(window):
                idx = start + j
                if 0 <= idx < len(rseq):
                    mat[i,j] = BASE_TO_INDEX.get(rseq[idx].upper(), 4)
        pileups[f"{vpos}_{vtype}"] = mat.tolist()
    return node_id, pileups

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Efficient pileup (.dat/.idx + GFA)")
    parser.add_argument("dat", help=".dat file path")
    parser.add_argument("idx", help=".idx file path")
    parser.add_argument("output", help="JSON output path")
    parser.add_argument("--gfa", help="GFA path if no cache")
    parser.add_argument("--load-cache", help="JSON file of node sequences")
    parser.add_argument("--save-cache", help="Write JSON cache here if building")
    parser.add_argument("--workers", type=int, default=8, help="Number of processes")
    args = parser.parse_args()
    # Ensure cache/GFA
    if not args.load_cache:
        if not args.gfa or not args.save_cache:
            sys.exit("❌ Provide --gfa and --save-cache when no --load-cache")
    # Step 1: parse index
    print("🔹 Parsing .idx file")
    node_index = parse_idx_file(args.idx)
    # Step 2: load/build sequence cache
    cache_path = args.load_cache or args.save_cache
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"✔ Loading sequence cache: {args.load_cache}")
        with open(args.load_cache) as cf:
            data = json.load(cf)
        node_sequences = {int(k):v for k,v in data.items()}
    else:
        print("🔹 Building sequence cache from GFA")
        node_sequences = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        with open(args.save_cache, 'w') as cf:
            json.dump({str(k):v for k,v in node_sequences.items()}, cf)
        print(f"✔ Saved cache to {args.save_cache}")
    # Step 3: multiprocessing pileup
    print("🔹 Starting multiprocessing pileup")
    total = len(node_index)
    results = {}
    start_time = time.time()
    # Launch pool
    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=_init_worker,
        initargs=(cache_path, args.dat)
    ) as executor:
        futures = {executor.submit(process_node, nid, off, cnt): nid
                   for nid,(off,cnt) in node_index.items()
                   if nid in node_sequences}
        processed = 0
        for future in as_completed(futures):
            processed += 1
            nid = futures[future]
            try:
                _, pileup = future.result(timeout=600)
                results[nid] = pileup
            except Exception as e:
                print(f"❌ Error on node {nid}: {e}")
            if processed % 100_000 == 0 or processed == total:
                print(f"✔ {processed}/{total} done — elapsed {time.time()-start_time:.2f}s")
    # Step 4: write output
    print("🔹 Writing JSON output")
    with open(args.output,'w') as out_f:
        json.dump(results,out_f,indent=2)
    print("✅ Completed. Output:", args.output)

if __name__=='__main__':
    main()
