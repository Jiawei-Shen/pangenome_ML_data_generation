#!/usr/bin/env python3
"""
Parallel pileup of read segments in .dat/.idx format,
filter by mapping quality, and retain per-base quality matrix.

Outputs, per node:
  - "pileup": list of strings (rows of bases, '.' for gap)
  - "bq_matrix": list of lists of ints (same shape, -1 for gap)
  - "allele_frequencies": per-column base frequencies

Usage:
  python pileup_dat_v2.py --dat reads.dat --idx reads.idx --out res.json \
       --min_mapq 20 --threads 8
"""
import argparse, struct, json
import numpy as np
from pathlib import Path
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
RECORD_STRUCT = struct.Struct("<h150s150shc")  # offset, seq, bq, rq, strand
RECORD_SIZE = RECORD_STRUCT.size

def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]

def load_idx(idx_path):
    d = {}
    with open(idx_path, "rb") as f:
        nb, = struct.unpack("<I", f.read(4))
        for _ in range(nb):
            bid, off, sz, nrec, _ = struct.unpack("<I Q I I H", f.read(22))
            d[bid] = (off, sz, nrec)
    return d

def load_node_segments(dat_path, node_id, idx, min_mapq):
    off, _, nrec = idx[node_id]
    segs = []
    with open(dat_path, "rb") as f:
        f.seek(off + 4+4+2)
        for _ in range(nrec):
            raw = f.read(RECORD_SIZE)
            node_off, sb, bqb, rq, strand_c = RECORD_STRUCT.unpack(raw)
            if rq < min_mapq:
                continue
            seq = sb.rstrip(b'\x00').decode('ascii')
            bq = list(bqb[:len(seq)])
            strand = strand_c.decode()
            segs.append({"offset":node_off, "sequence":seq,
                         "bq":bq, "strand":strand})
    return segs

def pileup_with_qualities(segs):
    if not segs:
        return np.empty((0,0),dtype="<U1"), np.empty((0,0),dtype=int)
    # compute width
    w = max(s["offset"]+len(s["sequence"]) for s in segs)
    pile_chars = []
    pile_bq    = []
    for s in segs:
        off = s["offset"]; seq = s["sequence"]; bq = s["bq"]
        if s["strand"] == "-":
            seq = reverse_complement(seq)
            bq  = bq[::-1]
            off = w - off - len(seq)
        # find or create a row
        placed = False
        for row_c, row_q in zip(pile_chars, pile_bq):
            # check if this row can accommodate
            if all((row_c[off+i] == '.') for i in range(len(seq))):
                for i,ch in enumerate(seq):
                    row_c[off+i] = ch.upper() if ch.upper() in "ACGT" else '.'
                    row_q[off+i] = bq[i]
                placed = True
                break
        if not placed:
            new_c = ['.']*w
            new_q = [-1]*w
            for i,ch in enumerate(seq):
                new_c[off+i] = ch.upper() if ch.upper() in "ACGT" else '.'
                new_q[off+i] = bq[i]
            pile_chars.append(new_c)
            pile_bq.append(new_q)
    return np.array(pile_chars,dtype="<U1"), np.array(pile_bq,dtype=int)

def allele_freq(mat):
    h,w=mat.shape; af={}
    for c in range(w):
        bases=[mat[r,c] for r in range(h) if mat[r,c] in "ACGT"]
        if not bases: continue
        cnt=Counter(bases); tot=sum(cnt.values())
        af[c]={b:round(v/tot,4) for b,v in cnt.items()}
    return af

def process_node(nid, dat, idx, min_mapq):
    segs = load_node_segments(dat, nid, idx, min_mapq)
    mat_c, mat_q = pileup_with_qualities(segs)
    af = allele_freq(mat_c)
    return nid, {
        "pileup": ["".join(row) for row in mat_c.tolist()],
        "bq_matrix": mat_q.tolist(),
        "allele_frequencies": {str(c):v for c,v in af.items()}
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--dat",required=True)
    p.add_argument("--idx",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--min_mapq",type=int,default=10)
    p.add_argument("--threads",type=int,default=4)
    args=p.parse_args()

    idx=load_idx(args.idx)
    results={}
    with ProcessPoolExecutor(max_workers=args.threads) as exe:
        fs={exe.submit(process_node,n, args.dat, idx, args.min_mapq):n
            for n in idx}
        for f in fs:
            nid,info=f.result()
            results[str(nid)] = info

    with open(args.out,"w") as fo:
        json.dump(results, fo, indent=2)

if __name__=="__main__":
    main()
