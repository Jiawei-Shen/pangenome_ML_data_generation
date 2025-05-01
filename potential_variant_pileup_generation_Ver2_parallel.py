#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import sqlite3
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import math

# ─────────────────────────────────────────────────────────────────────────────
# Globals (set up in main and inherited via fork)
GLOBAL_DAT_FILE = None       # open file handle for .dat
GLOBAL_DB = None             # SQLite connection for node sequences
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# ─────────────────────────────────────────────────────────────────────────────
def init_worker(dat_path, db_path):
    """Initializer for each worker: open .dat file and sqlite DB."""
    global GLOBAL_DAT_FILE, GLOBAL_DB
    GLOBAL_DAT_FILE = open(dat_path, 'rb')
    GLOBAL_DB = sqlite3.connect(db_path, uri=True, check_same_thread=False)

# ─────────────────────────────────────────────────────────────────────────────
def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]

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
def process_node_chunk(chunk, window=60):
    """
    Process a list of nodes: [(node_id, offset, n_records), ...]
    Opens its own sqlite cursor and dat file handle from globals.
    """
    cursor = GLOBAL_DB.cursor()
    # Fetch all sequences needed for this chunk
    node_ids = [nid for nid, _, _ in chunk]
    placeholders = ",".join(['?'] * len(node_ids))
    cursor.execute(f"SELECT node_id, seq FROM sequences WHERE node_id IN ({placeholders})", node_ids)
    seq_map = {nid: seq for nid, seq in cursor.fetchall()}

    out = {}
    for node_id, offset, n_records in chunk:
        seq = seq_map.get(node_id)
        if seq is None:
            continue
        node_len = len(seq)
        # Phase 1: read & filter MQ>=10
        segments = []
        # skip header bytes
        ptr = offset + 10
        for _ in range(n_records):
            GLOBAL_DAT_FILE.seek(ptr)
            data = GLOBAL_DAT_FILE.read(RECORD_SIZE)
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
        half = window // 2
        for (vpos, vtype), reads in reads_by_variant.items():
            mat = np.full((len(reads), window), 4, dtype=np.uint8)
            for i, (roff, rseq) in enumerate(reads):
                start = vpos - roff - half
                for j in range(window):
                    idx = start + j
                    if 0 <= idx < len(rseq):
                        mat[i, j] = BASE_TO_INDEX.get(rseq[idx].upper(), 4)
            out[f"{node_id}_{vpos}_{vtype}"] = mat.tolist()
    return out

# ─────────────────────────────────────────────────────────────────────────────
def parse_idx_file(idx_path):
    node_index = {}
    with open(idx_path, 'rb') as f:
        num_nodes = struct.unpack('<I', f.read(4))[0]
        for _ in range(num_nodes):
            nid, off, bs, nr, _ = struct.unpack('<I Q I I H', f.read(22))
            node_index[nid] = (off, nr)
    print(f"✔ Parsed {len(node_index)} nodes", file=sys.stderr, flush=True)
    return node_index

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Parallel variant pileup without mmap")
    parser.add_argument("dat", help=".dat file path")
    parser.add_argument("idx", help=".idx file path")
    parser.add_argument("output", help="JSON output path")
    parser.add_argument("--load-cache", help="JSON file of node sequences")
    parser.add_argument("--save-cache", help="Where to write node-sequence cache (JSON)")
    args = parser.parse_args()

    if not args.load_cache:
        sys.exit("❌ Provide --load-cache (prebuilt JSON cache)")
    # 1) parse idx
    node_index = parse_idx_file(args.idx)

    # 2) load JSON cache and build sqlite DB
    cache_path = args.load_cache
    print(f"✔ Loading node sequences from {cache_path}", file=sys.stderr, flush=True)
    with open(cache_path) as cf:
        seq_data = json.load(cf)
    # create SQLite DB for fast, memory-efficient lookups
    db_path = cache_path + ".db"
    if not os.path.isfile(db_path):
        print(f"🔹 Building SQLite DB at {db_path}…", file=sys.stderr, flush=True)
        conn = sqlite3.connect(db_path)
        conn.execute("CREATE TABLE sequences(node_id INTEGER PRIMARY KEY, seq TEXT)")
        conn.executemany(
            "INSERT INTO sequences(node_id, seq) VALUES (?, ?)",
            [(int(k), v) for k, v in seq_data.items()]
        )
        conn.commit()
        conn.close()
        print("✔ SQLite DB ready", file=sys.stderr, flush=True)
    # drop in-memory JSON to free RAM
    del seq_data

    # 3) prepare chunks
    items = [ (nid, off, nrec) for nid, (off, nrec) in node_index.items() ]
    total_nodes = len(items)
    CHUNK_SIZE = 1000
    chunks = [ items[i:i+CHUNK_SIZE] for i in range(0, total_nodes, CHUNK_SIZE) ]
    total_chunks = len(chunks)
    print(f"🔹 Dispatching {total_chunks} chunks for {total_nodes} nodes…", file=sys.stderr, flush=True)

    # 4) parallel execution
    start = time.time()
    results = {}
    with ProcessPoolExecutor(
        initializer=init_worker,
        initargs=(args.dat, db_path)
    ) as exe:
        futures = [ exe.submit(process_node_chunk, chunk) for chunk in chunks ]
        completed_nodes = 0
        for fut in as_completed(futures):
            chunk_res = fut.result()
            results.update(chunk_res)
            completed_nodes += len(chunk_res)
            if completed_nodes % (CHUNK_SIZE * 10) == 0 or completed_nodes == total_nodes:
                elapsed = time.time() - start
                print(f"✔ {completed_nodes}/{total_nodes} nodes done — {elapsed:.1f}s", file=sys.stderr, flush=True)

    # 5) write output
    print("🔹 Writing JSON output…", file=sys.stderr, flush=True)
    with open(args.output, 'w') as out_f:
        json.dump(results, out_f)
    print(f"✅ Done. Wrote output to {args.output}", file=sys.stderr, flush=True)

if __name__ == '__main__':
    main()
