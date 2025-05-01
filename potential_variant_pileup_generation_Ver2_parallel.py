#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import math

# ─────────────────────────────────────────────────────────────────────────────
# Globals (inherited by worker processes via fork)
GLOBAL_DAT_FILE = None        # file handle for .dat reads
GLOBAL_NODE_SEQS = None       # dict of node_id → sequence
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
CHUNK_SIZE = 1000

# ─────────────────────────────────────────────────────────────────────────────
def init_worker(dat_path, node_seqs):
    """Initializer for each worker: open the .dat file and grab sequences dict."""
    global GLOBAL_DAT_FILE, GLOBAL_NODE_SEQS
    GLOBAL_DAT_FILE = open(dat_path, 'rb')
    GLOBAL_NODE_SEQS = node_seqs

# ─────────────────────────────────────────────────────────────────────────────
def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

# ─────────────────────────────────────────────────────────────────────────────
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

# ─────────────────────────────────────────────────────────────────────────────
def process_node(node_id, offset, n_records):
    """Extract reads and build variant-centered pileups (MAPQ>=10)."""
    f = GLOBAL_DAT_FILE
    seq = GLOBAL_NODE_SEQS[node_id]
    node_len = len(seq)

    # Phase 1: read & filter by mapping quality
    segments = []
    f.seek(offset + 10)
    for _ in range(n_records):
        data = f.read(RECORD_SIZE)
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

# ─────────────────────────────────────────────────────────────────────────────
def process_node_chunk(chunk):
    """Process a chunk of nodes, returns dict of results for the chunk."""
    out = {}
    for nid, off, nrec in chunk:
        if nid not in GLOBAL_NODE_SEQS:
            continue
        _, pu = process_node(nid, off, nrec)
        out[nid] = pu
    return out

# ─────────────────────────────────────────────────────────────────────────────
def chunked_iterable(items, size):
    for i in range(0, len(items), size):
        yield items[i:i+size]

# ─────────────────────────────────────────────────────────────────────────────
def parse_idx_file(idx_path):
    node_index = {}
    with open(idx_path, 'rb') as f:
        num = struct.unpack('<I', f.read(4))[0]
        for _ in range(num):
            nid, off, bs, nr, _ = struct.unpack('<I Q I I H', f.read(22))
            node_index[nid] = (off, nr)
    print(f"✔ Parsed {len(node_index)} nodes", file=sys.stderr, flush=True)
    return node_index

# ─────────────────────────────────────────────────────────────────────────────
def load_node_sequences_from_gfa(gfa_path, target_ids):
    seqs = {}
    with open(gfa_path, 'r') as f:
        for lc, line in enumerate(f, start=1):
            if lc % 10_000_000 == 0:
                print(f"✔ Checked {lc} lines in GFA", file=sys.stderr, flush=True)
            if not line.startswith('S\t'): continue
            _, sid, seq = line.rstrip('\n').split('\t',2)
            try:
                nid = int(sid)
            except ValueError:
                continue
            if nid in target_ids:
                seqs[nid] = seq
    print(f"✔ Loaded {len(seqs)} sequences", file=sys.stderr, flush=True)
    return seqs

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dat")
    parser.add_argument("idx")
    parser.add_argument("output")
    parser.add_argument("--gfa")
    parser.add_argument("--load-cache")
    parser.add_argument("--save-cache")
    args = parser.parse_args()

    if not args.load_cache and not args.gfa:
        sys.exit("❌ Provide --gfa or --load-cache.")
    if not args.load_cache and not args.save_cache:
        sys.exit("❌ Provide --save-cache when building from GFA.")

    node_index = parse_idx_file(args.idx)
    cache = args.load_cache or args.save_cache
    if args.load_cache and os.path.isfile(cache):
        print(f"✔ Loading node sequences from {cache}", file=sys.stderr, flush=True)
        with open(cache) as cf:
            data = json.load(cf)
        seqs = {int(k): v for k,v in data.items()}
    else:
        print("🔹 Building node sequences from GFA", file=sys.stderr, flush=True)
        seqs = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        with open(args.save_cache, 'w') as cf:
            json.dump({str(k):v for k,v in seqs.items()}, cf)
        print(f"✔ Saved cache to {args.save_cache}", file=sys.stderr, flush=True)

    # build flat list of items
    items = [(nid, off, nrec)
             for nid,(off,nrec) in node_index.items()
             if nid in seqs]
    total_nodes = len(items)
    total_chunks = math.ceil(total_nodes/CHUNK_SIZE)
    print(f"🔹 Dispatching {total_chunks} chunks for {total_nodes} nodes…", file=sys.stderr, flush=True)

    start = time.time()
    results = {}
    completed = 0

    with ProcessPoolExecutor(initializer=init_worker, initargs=(args.dat, seqs)) as exe:
        chunks = list(chunked_iterable(items, CHUNK_SIZE))
        futures = [exe.submit(process_node_chunk, c) for c in chunks]
        for fut in as_completed(futures):
            chunk_res = fut.result()
            results.update(chunk_res)
            completed += len(chunk_res)
            if completed % CHUNK_SIZE == 0 or completed == total_nodes:
                elapsed = time.time() - start
                print(f"✔ {completed}/{total_nodes} nodes — {elapsed:.1f}s", file=sys.stderr, flush=True)

    print("🔹 Writing JSON output", file=sys.stderr, flush=True)
    with open(args.output, 'w') as out:
        json.dump(results, out, indent=2)
    print(f"✅ Done (wrote {len(results)} nodes).", file=sys.stderr, flush=True)

if __name__ == "__main__":
    main()
