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
from multiprocessing import get_context

# ─────────────────────────────────────────────────────────────────────────────
# Globals for file access (inherited via fork)
GLOBAL_DAT = None            # memory-mapped .dat file
GLOBAL_NODE_SEQS = {}        # loaded before forking
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# ─────────────────────────────────────────────────────────────────────────────
def init_worker(dat_path):
    """Initializer: mmap the .dat once. GLOBAL_NODE_SEQS inherited via fork."""
    global GLOBAL_DAT
    fd = open(dat_path, 'rb')
    GLOBAL_DAT = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
    # GLOBAL_NODE_SEQS remains available from parent via fork

# ─────────────────────────────────────────────────────────────────────────────
def reverse_complement(sequence):
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement)[::-1]

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
    seq = GLOBAL_NODE_SEQS[node_id]
    node_len = len(seq)
    # Phase 1: read & filter MAPQ>=10
    segs = []
    ptr = offset + 10
    for _ in range(n_records):
        data = GLOBAL_DAT[ptr:ptr+RECORD_SIZE]
        ptr += RECORD_SIZE
        off, raw_seq, raw_bq, raw_cigar, mapq, strand = RECORD_STRUCT.unpack(data)
        if mapq < 10:
            continue
        s = raw_seq.rstrip(b"\x00").decode('ascii', errors='ignore')
        cig = raw_cigar.rstrip(b"\x00").decode('ascii', errors='ignore')
        if strand.decode() == '-':
            off = node_len - off - len(s)
            s = reverse_complement(s)
        segs.append((off, s, cig))
    # Phase 2: group by variant
    reads_by_variant = defaultdict(list)
    for off, s, cig in segs:
        for vpos, vtype in detect_variants_from_cigar(off, cig):
            if 0 <= vpos < node_len:
                reads_by_variant[(vpos, vtype)].append((off, s))
    # Phase 3: build pileup
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
def parse_idx_file(idx_path):
    idx = {}
    with open(idx_path, 'rb') as f:
        num = struct.unpack('<I', f.read(4))[0]
        for _ in range(num):
            nid, off, bs, nr, _ = struct.unpack('<I Q I I H', f.read(22))
            idx[nid] = (off, nr)
    print(f"✔ Parsed {len(idx)} nodes", file=sys.stderr, flush=True)
    return idx

# ─────────────────────────────────────────────────────────────────────────────
def load_node_sequences_from_gfa(gfa_path, target_ids):
    seqs = {}
    with open(gfa_path, 'r') as f:
        for lc, line in enumerate(f, start=1):
            if lc % 10_000_000 == 0:
                print(f"✔ Scanned {lc} GFA lines", file=sys.stderr, flush=True)
            if not line.startswith('S\t'): continue
            _, sid, s = line.rstrip().split('\t', 2)
            try:
                nid = int(sid)
            except:
                continue
            if nid in target_ids:
                seqs[nid] = s
    print(f"✔ Loaded {len(seqs)} sequences", file=sys.stderr, flush=True)
    return seqs

# ─────────────────────────────────────────────────────────────────────────────
def main():
    global GLOBAL_NODE_SEQS
    p = argparse.ArgumentParser()
    p.add_argument('dat')
    p.add_argument('idx')
    p.add_argument('output')
    p.add_argument('--gfa')
    p.add_argument('--load-cache')
    p.add_argument('--save-cache')
    args = p.parse_args()

    if not args.load_cache and not args.gfa:
        sys.exit("❌ Provide --gfa or --load-cache.")
    if not args.load_cache and not args.save_cache:
        sys.exit("❌ Provide --save-cache when building from GFA.")

    # index + sequences
    node_index = parse_idx_file(args.idx)
    cache = args.load_cache or args.save_cache
    if args.load_cache and os.path.isfile(cache):
        print(f"✔ Loading node sequences from {cache}", file=sys.stderr, flush=True)
        with open(cache) as cf:
            data = json.load(cf)
        GLOBAL_NODE_SEQS = {int(k):v for k,v in data.items()}
    else:
        print("🔹 Building sequences from GFA", file=sys.stderr, flush=True)
        GLOBAL_NODE_SEQS = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        with open(args.save_cache,'w') as cf:
            json.dump({str(k):v for k,v in GLOBAL_NODE_SEQS.items()}, cf)
        print(f"✔ Saved cache to {args.save_cache}", file=sys.stderr, flush=True)

    # parallel dispatch
    dat_path, idx_map = args.dat, node_index
    total = len(idx_map)
    print(f"🔹 Dispatching {total} nodes…", file=sys.stderr, flush=True)

    start = time.time()
    results = {}

    ctx = get_context('fork')
    with ProcessPoolExecutor(
            max_workers=os.cpu_count(),
            mp_context=ctx,
            initializer=init_worker,
            initargs=(dat_path,)
        ) as exe:

        # submit tasks
        futures = {exe.submit(process_node, nid, off, nrec): nid
                   for nid,(off,nrec) in idx_map.items()
                   if nid in GLOBAL_NODE_SEQS}

        print(f"✔ Scheduled {len(futures)} tasks", file=sys.stderr, flush=True)

        # collect
        for done_count, fut in enumerate(as_completed(futures), start=1):
            nid = futures[fut]
            try:
                _, pu = fut.result()
                results[nid] = pu
            except Exception as e:
                print(f"⚠ Node {nid} error: {e}", file=sys.stderr, flush=True)

            # progress
            if done_count % 1000 == 0 or done_count == total:
                elapsed = time.time() - start
                print(f"✔ {done_count}/{total} done — {elapsed:.2f}s", file=sys.stderr, flush=True)

    # output
    print("🔹 Writing JSON…", file=sys.stderr, flush=True)
    with open(args.output,'w') as of:
        json.dump(results, of, indent=2)
    print(f"✅ Done writing {args.output}", file=sys.stderr, flush=True)

if __name__ == '__main__':
    main()
