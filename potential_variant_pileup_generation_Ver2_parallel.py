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
import re

RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# ─────────────────────────────────────────────────────────────────────────────
def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]

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

# ─────────────────────────────────────────────────────────────────────────────
def process_node_batch(node_batch, dat_path, node_seqs):
    results = {}
    with open(dat_path, 'rb') as dat_file:
        for node_id, offset, n_records in node_batch:
            if node_id not in node_seqs:
                continue
            sequence = node_seqs[node_id]
            node_len = len(sequence)
            segments = []
            dat_file.seek(offset + 10)
            for _ in range(n_records):
                data = dat_file.read(RECORD_SIZE)
                off, raw_seq, _, raw_cigar, mapq, strand = RECORD_STRUCT.unpack(data)
                if mapq < 10:
                    continue
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
                cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
                strand_char = strand.decode()
                if strand_char == '-':
                    off = node_len - off - len(seq)
                    seq = reverse_complement(seq)
                segments.append((off, seq, cigar))

            reads_by_variant = defaultdict(list)
            for off, seq, cigar in segments:
                for vpos, vtype in detect_variants_from_cigar(off, cigar):
                    if 0 <= vpos < node_len:
                        reads_by_variant[(vpos, vtype)].append((off, seq))

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

            results[node_id] = pileups
    return results

# ─────────────────────────────────────────────────────────────────────────────
def parse_idx_file(idx_path):
    index = []
    with open(idx_path, 'rb') as f:
        num_nodes = struct.unpack('<I', f.read(4))[0]
        for _ in range(num_nodes):
            nid, offset, _, nrec, _ = struct.unpack('<I Q I I H', f.read(22))
            index.append((nid, offset, nrec))
    print(f"✔ Parsed {len(index)} nodes.")
    return index

def load_node_sequences(path):
    with open(path) as f:
        data = json.load(f)
    return {int(k): v for k, v in data.items()}

def chunk_list(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i + size]

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dat", help=".dat file path")
    parser.add_argument("idx", help=".idx file path")
    parser.add_argument("output", help="JSON output path")
    parser.add_argument("--load-cache", required=True, help="JSON file of node sequences")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=1000)
    args = parser.parse_args()

    print("🔹 Loading index and sequences")
    node_index = parse_idx_file(args.idx)
    node_seqs = load_node_sequences(args.load_cache)

    node_chunks = list(chunk_list(node_index, args.chunk_size))
    total_chunks = len(node_chunks)
    print(f"🔹 Dispatching {total_chunks} chunks using {args.workers} workers…")

    start = time.time()
    with open(args.output, 'w') as outf:
        outf.write("{\n")
        first = True
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(process_node_batch, chunk, args.dat, node_seqs)
                for chunk in node_chunks
            ]
            for i, fut in enumerate(as_completed(futures), start=1):
                batch_res = fut.result()
                for nid, pileup in batch_res.items():
                    if not first:
                        outf.write(",\n")
                    first = False
                    json.dump(str(nid), outf)
                    outf.write(": ")
                    json.dump(pileup, outf)
                print(f"✔ {i}/{total_chunks} chunks done ({i * args.chunk_size} nodes)", file=sys.stderr)
        outf.write("\n}\n")

    print(f"✅ Done. Elapsed time: {time.time() - start:.2f}s")

if __name__ == "__main__":
    main()
