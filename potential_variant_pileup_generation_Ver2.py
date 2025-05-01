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

BASE_TO_INDEX = {'A':0,'C':1,'G':2,'T':3,'N':4}

# ─────────────────────────────────────────────────────────────────────────────
def reverse_complement(seq):
    comp = str.maketrans("ACGTacgtNn","TGCAtgcaNn")
    return seq.translate(comp)[::-1]

# ─────────────────────────────────────────────────────────────────────────────
def parse_idx_file(idx_path):
    """
    Returns: dict node_id -> (offset_in_dat, n_records)
    """
    node_index = {}
    with open(idx_path,"rb") as f:
        num = struct.unpack("<I", f.read(4))[0]
        for _ in range(num):
            node_id, offset, block_size, n_records, _ = struct.unpack("<I Q I I H", f.read(22))
            node_index[node_id] = (offset, n_records)
    print(f"✔ Parsed {len(node_index)} nodes from index file.")
    return node_index

# ─────────────────────────────────────────────────────────────────────────────
def load_node_sequences_from_gfa(gfa_path, target_node_ids):
    node_sequences = {}
    line_count = 0
    with open(gfa_path,"r") as f:
        for line in f:
            line_count += 1
            if line_count % 10_000_000 == 0:
                print(f"✔ Checked {line_count} lines in GFA")
            if not line.startswith("S\t"):
                continue
            parts = line.split("\t",2)
            try:
                nid = int(parts[1])
            except:
                continue
            if nid in target_node_ids:
                node_sequences[nid] = parts[2].rstrip("\n")
    print(f"✔ Loaded {len(node_sequences)} sequences from GFA.")
    return node_sequences

# ─────────────────────────────────────────────────────────────────────────────
def decode_cigar(cigar):
    import re
    return re.findall(r'(\d+)([MXDI])', cigar)

def detect_variants_from_cigar(offset, cigar):
    variants=[]
    pos=offset
    for length,op in decode_cigar(cigar):
        l=int(length)
        if op in ("X","I","D"):
            for i in range(l):
                if op=="I":
                    variants.append((pos-1,"I"))
                elif op=="D":
                    variants.append((pos,"D"))
                    pos+=1
                else:
                    variants.append((pos,"X"))
                    pos+=1
        else:
            pos+=l
    return variants

# ─────────────────────────────────────────────────────────────────────────────
def process_node(args):
    nid, seq, dat_path, dat_offset, n_records = args
    print(f"▶ Node {nid}: reading {n_records} records at offset {dat_offset}")
    node_len = len(seq)

    # read exactly this node’s block
    reads_by_variant = defaultdict(list)
    with open(dat_path,"rb") as f:
        f.seek(dat_offset)
        # block header: 4B node_id,4B n_records,2B zero
        header = f.read(10)
        _, recs, _ = struct.unpack("<I I H",header)
        assert recs==n_records

        for _ in range(n_records):
            data = f.read(RECORD_SIZE)
            off, seq_raw, bq_raw, cigar_raw, rq, strand = RECORD_STRUCT.unpack(data)
            s = seq_raw.rstrip(b'\x00').decode("ascii","ignore")
            cigar = cigar_raw.rstrip(b'\x00').decode("ascii","ignore")
            st = strand.decode()
            read_len = len(s)

            if st=="-":
                off = node_len - off - read_len
                s = reverse_complement(s)

            for vpos, vtype in detect_variants_from_cigar(off,cigar):
                if 0<=vpos<node_len:
                    reads_by_variant[(vpos,vtype)].append((off,s))

    # build pileup matrices
    pileups={}
    W=60
    H=W//2
    for (vpos,vtype), reads in reads_by_variant.items():
        M = np.full((len(reads),W),4,dtype=np.uint8)
        for i,(off,s) in enumerate(reads):
            start = vpos - off - H
            for j in range(W):
                idx = start+j
                if 0<=idx<len(s):
                    M[i,j] = BASE_TO_INDEX.get(s[idx].upper(),4)
        pileups[f"{vpos}_{vtype}"] = M.tolist()

    print(f"✔ Node {nid}: {len(pileups)} variant sites")
    return nid,pileups

# ─────────────────────────────────────────────────────────────────────────────
def main():
    p=argparse.ArgumentParser(
       description="Variant-centered pileup from .dat/.idx + GFA")
    p.add_argument("dat")
    p.add_argument("idx")
    p.add_argument("output",help="JSON output path")
    p.add_argument("--gfa", help="GFA path (required if no --load-cache)")
    p.add_argument("--load-cache", help="JSON cache of node sequences")
    p.add_argument("--save-cache", help="Where to write sequence cache")
    p.add_argument("--workers",type=int,default=8)
    args=p.parse_args()

    if not args.load_cache and not args.gfa:
        sys.exit("Error: --gfa required if no --load-cache")
    if args.save_cache and not args.gfa:
        sys.exit("Error: --save-cache requires --gfa")

    print("Step 1: parse .idx")
    node_index = parse_idx_file(args.idx)

    # load or build sequences
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"✔ Loading sequences from cache {args.load_cache}")
        C=json.load(open(args.load_cache))
        node_seqs={int(k):v for k,v in C.items()}
    else:
        print("Step 2: parse GFA")
        node_seqs=load_node_sequences_from_gfa(args.gfa,node_index.keys())
        if args.save_cache:
            print(f"✔ Saving cache to {args.save_cache}")
            with open(args.save_cache,"w") as cf:
                json.dump({str(k):v for k,v in node_seqs.items()},cf)

    print("Step 3: pileup in parallel")
    tasks=[(nid,node_seqs[nid],args.dat,offset,nr)
           for nid,(offset,nr) in node_index.items()
           if nid in node_seqs]

    final={}
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for nid,p in ex.map(process_node,tasks):
            final[nid]=p

    print("Step 4: write JSON")
    with open(args.output,"w") as out:
        json.dump(final,out,indent=2)

    print("✅ Done")

if __name__=="__main__":
    main()
