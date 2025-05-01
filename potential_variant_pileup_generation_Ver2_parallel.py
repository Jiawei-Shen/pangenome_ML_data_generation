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
# Globals for file access (populated in each worker)
GLOBAL_DAT = None
GLOBAL_NODE_SEQS = None
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

def init_worker(dat_path, node_seqs):
    """Initializer for each worker: mmap the .dat and grab the sequences dict."""
    global GLOBAL_DAT, GLOBAL_NODE_SEQS
    dat_f = open(dat_path, 'rb')
    GLOBAL_DAT = mmap.mmap(dat_f.fileno(), 0, access=mmap.ACCESS_READ)
    GLOBAL_NODE_SEQS = node_seqs

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def decode_cigar(cigar_string):
    import re
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

def process_node(node_id, offset, n_records):
    """Worker task: same as before, but uses GLOBAL_DAT & GLOBAL_NODE_SEQS."""
    seq = GLOBAL_NODE_SEQS[node_id]
    node_len = len(seq)

    # Phase 1: read & MQ-filter
    segments = []
    ptr = offset + 10
    for _ in range(n_records):
        data = GLOBAL_DAT[ptr:ptr + RECORD_SIZE]
        ptr += RECORD_SIZE
        off, raw_seq, raw_bq, raw_cigar, mapq, strand = RECORD_STRUCT.unpack(data)
        if mapq < 10:
            continue
        s = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
        cig = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
        strand_char = strand.decode()
        if strand_char == '-':
            off = node_len - off - len(s)
            s = reverse_complement(s)
        segments.append((off, s, cig))

    # Phase 2: detect variants
    reads_by_variant = defaultdict(list)
    for off, s, cig in segments:
        for vpos, vtype in detect_variants_from_cigar(off, cig):
            if 0 <= vpos < node_len:
                reads_by_variant[(vpos, vtype)].append((off, s))

    # Phase 3: build pileups
    pileups = {}
    window, half = 60, 30
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

def parse_idx_file(idx_path):
    node_index = {}
    with open(idx_path, 'rb') as f:
        num_nodes = struct.unpack('<I', f.read(4))[0]
        for _ in range(num_nodes):
            nid, off, bs, nr, _ = struct.unpack('<I Q I I H', f.read(22))
            node_index[nid] = (off, nr)
    print(f"✔ Parsed {len(node_index)} nodes from index file.")
    return node_index

def load_node_sequences_from_gfa(gfa_path, target_node_ids):
    node_sequences = {}
    with open(gfa_path, 'r') as f:
        for lc, line in enumerate(f, start=1):
            if lc % 10_000_000 == 0:
                print(f"✔ Checked {lc} lines in GFA")
            if not line.startswith('S\t'): continue
            _, sid, seq = line.rstrip('\n').split('\t', 2)
            try:
                nid = int(sid)
            except ValueError:
                continue
            if nid in target_node_ids:
                node_sequences[nid] = seq
    print(f"✔ Loaded {len(node_sequences)} sequences from GFA.")
    return node_sequences

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dat", help=".dat file path")
    parser.add_argument("idx", help=".idx file path")
    parser.add_argument("output", help="JSON output path")
    parser.add_argument("--gfa", help="GFA path if no cache")
    parser.add_argument("--load-cache", help="JSON file of node sequences")
    parser.add_argument("--save-cache", help="Where to write node-sequence cache")
    args = parser.parse_args()

    if not args.load_cache and not args.gfa:
        sys.exit("❌ Provide --gfa or --load-cache.")
    if not args.load_cache and not args.save_cache:
        sys.exit("❌ Provide --save-cache when building from GFA.")

    # 1) parse idx & load sequences
    node_index = parse_idx_file(args.idx)
    cache = args.load_cache or args.save_cache
    if args.load_cache and os.path.isfile(cache):
        print(f"✔ Loading node sequences from {cache}")
        with open(cache) as cf:
            data = json.load(cf)
        seqs = {int(k): v for k, v in data.items()}
    else:
        print("🔹 Building node sequences from GFA")
        seqs = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        with open(args.save_cache, 'w') as cf:
            json.dump({str(k): v for k, v in seqs.items()}, cf)
        print(f"✔ Saved node sequence cache to {args.save_cache}")

    # 2) dispatch parallel jobs
    print("🔹 Spawning workers and memory-mapping .dat in each")
    start = time.time()
    results = {}
    total = len(node_index)

    with ProcessPoolExecutor(
            initializer=init_worker,
            initargs=(args.dat, seqs)
        ) as exe:

        # submit one future per node
        future_to_nid = {
            exe.submit(process_node, nid, off, nrec): nid
            for nid, (off, nrec) in node_index.items()
            if nid in seqs
        }

        completed = 0
        for fut in as_completed(future_to_nid):
            nid = future_to_nid[fut]
            try:
                _, pileup = fut.result()
                results[nid] = pileup
            except Exception as e:
                print(f"⚠ Error processing node {nid}: {e}", file=sys.stderr)
            completed += 1
            print(completed)
            if completed % 1000 == 0 or completed == total:
                elapsed = time.time() - start
                print(f"✔ {completed}/{total} nodes processed — elapsed {elapsed:.2f}s")

    # 3) write out
    print("🔹 Writing JSON output")
    with open(args.output, 'w') as out_f:
        json.dump(results, out_f, indent=2)
    print(f"✅ Done. Wrote output to {args.output}")

if __name__ == "__main__":
    main()
