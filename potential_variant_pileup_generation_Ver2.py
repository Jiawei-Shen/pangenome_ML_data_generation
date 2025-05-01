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
from concurrent.futures import ProcessPoolExecutor

# ─────────────────────────────────────────────────────────────────────────────
# Global variables initialized once per worker
GLOBAL_DAT = None            # memory-mapped .dat file
GLOBAL_NODE_SEQS = None      # dict of node_id → sequence
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# ─────────────────────────────────────────────────────────────────────────────
def _init_worker(cache_path, dat_path):
    """
    Worker initializer: load node sequences into GLOBAL_NODE_SEQS and memory-map the .dat file.
    """
    global GLOBAL_DAT, GLOBAL_NODE_SEQS
    # Load node sequences from JSON cache
    with open(cache_path, 'r') as cf:
        seq_cache = json.load(cf)
    GLOBAL_NODE_SEQS = {int(k): v for k, v in seq_cache.items()}
    # Memory-map the .dat file for fast random access
    dat_fh = open(dat_path, 'rb')
    GLOBAL_DAT = mmap.mmap(dat_fh.fileno(), 0, access=mmap.ACCESS_READ)

# ─────────────────────────────────────────────────────────────────────────────
def parse_idx_file(idx_path):
    """
    Parse .idx and return dict of node_id → (offset, n_records).
    """
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
    """
    Parse GFA and extract sequences for target_node_ids.
    """
    node_sequences = {}
    line_count = 0
    with open(gfa_path, 'r') as f:
        for line in f:
            line_count += 1
            if line_count % 10_000_000 == 0:
                print(f"✔ Checked {line_count} lines in GFA file")
            if not line.startswith('S\t'):
                continue
            parts = line.split('\t', 2)
            try:
                nid = int(parts[1])
            except:
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
    """
    Worker function: extract records for node_id at given offset, build variant-centered pileups.
    """
    sequence = GLOBAL_NODE_SEQS[node_id]
    node_len = len(sequence)
    reads_by_variant = defaultdict(list)
    # Block header is 10 bytes; records start at offset + 10
    base = offset + 10
    for _ in range(n_records):
        data = GLOBAL_DAT[base: base + RECORD_SIZE]
        base += RECORD_SIZE
        off, seq_raw, bq_raw, cigar_raw, rq, strand = RECORD_STRUCT.unpack(data)
        seq = seq_raw.rstrip(b'\x00').decode('ascii', 'ignore')
        cigar = cigar_raw.rstrip(b'\x00').decode('ascii', 'ignore')
        strand_char = strand.decode()
        read_len = len(seq)
        # Normalize to + strand
        if strand_char == '-':
            off = node_len - off - read_len
            seq = seq.translate(str.maketrans('ACGTacgtNn','TGCAtgcaNn'))[::-1]
        for vpos, vtype in detect_variants_from_cigar(off, cigar):
            if 0 <= vpos < node_len:
                reads_by_variant[(vpos, vtype)].append((off, seq))
    # Build pileup matrices
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
    p = argparse.ArgumentParser(description="Efficient variant-centered pileup (.dat/.idx + GFA)")
    p.add_argument("dat", help=".dat file path")
    p.add_argument("idx", help=".idx file path")
    p.add_argument("output", help="JSON output path")
    p.add_argument("--gfa", help="GFA path (required if no --load-cache)")
    p.add_argument("--load-cache", help="JSON cache of node sequences")
    p.add_argument("--save-cache", help="Path to write node-sequence cache")
    p.add_argument("--workers", type=int, default=8, help="Number of processes")
    args = p.parse_args()
    # Validate cache/GFA
    if not args.load_cache:
        if not args.gfa or not args.save_cache:
            sys.exit("❌ Must provide --gfa and --save-cache when no --load-cache")
    # Parse index
    print("🔹 Parsing .idx file")
    node_index = parse_idx_file(args.idx)
    # Load or build node-sequence cache
    cache_path = args.load_cache or args.save_cache
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"✔ Loading sequence cache from {args.load_cache}")
    else:
        print("🔹 Building sequence cache from GFA")
        seqs = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        with open(args.save_cache, 'w') as cf:
            json.dump({str(k): v for k, v in seqs.items()}, cf)
        print(f"✔ Saved sequence cache to {args.save_cache}")
    # Prepare task lists
    node_ids = []
    offsets = []
    counts = []
    for nid, (off, nrec) in node_index.items():
        node_ids.append(nid)
        offsets.append(off)
        counts.append(nrec)
    # ─────────────────────────────────────────────────────────────────────────────
    print("🔹 Step 3: Pileup with multiprocessing")
    start_time = time.time()
    total_nodes = len(node_index)
    results = {}
    # Submit tasks and track futures
    futures = {}
    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=_init_worker,
        initargs=(cache_path, args.dat)
    ) as executor:
        for node_id, (off, nrec) in node_index.items():
            if node_id in GLOBAL_NODE_SEQS or os.path.isfile(cache_path):
                future = executor.submit(process_node, node_id, off, nrec)
                futures[future] = node_id
        processed = 0
        for future in as_completed(futures):
            node_id = futures[future]
            processed += 1
            try:
                _, pileup = future.result(timeout=600)
                results[node_id] = pileup
            except Exception as e:
                print(f"❌ Error processing node {node_id}: {e}")
            # Always print progress per node
            print(f"Processed {processed}/{total_nodes} nodes (Node {node_id})")
            if processed % 100_000 == 0:
                elapsed = time.time() - start_time
                print(f"✔ Completed {processed} nodes — elapsed: {elapsed:.2f}s")

    # Write final output
    print("🔹 Step 4: Writing JSON output")
    with open(args.output, 'w') as out_f:
        json.dump(results, out_f, indent=2)
    print("✅ Done. Output:", args.output)
    print("🔹 Writing JSON output")
    with open(args.output, 'w') as out_f:
        json.dump(results, out_f, indent=2)
    print("✅ Done. Output:", args.output)

if __name__ == '__main__':
    main()
