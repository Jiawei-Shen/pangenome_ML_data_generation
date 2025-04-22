#!/usr/bin/env python3
import argparse
import gzip
import pickle
import sqlite3
import zlib
import vg_pb2
import gc
import time
from collections import defaultdict

MAX_CACHE_SIZE = 1_000_000  # Flush when a single node gets this many reads

def read_varint(stream):
    value, shift = 0, 0
    while True:
        b = stream.read(1)
        if not b:
            raise EOFError
        byte = b[0]
        value |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            break
        shift += 7
    return value

def file_is_gzip(path):
    with open(path, "rb") as f:
        return f.read(2) == b'\x1f\x8b'

def gam_record_iter(path, tag="GAM"):
    open_func = gzip.open if file_is_gzip(path) else open
    with open_func(path, "rb") as f:
        while True:
            try:
                group_count = read_varint(f)
                if group_count == 0:
                    continue
                tag_len = read_varint(f)
                group_tag = f.read(tag_len).decode()
                if group_tag != tag:
                    for _ in range(group_count - 1):
                        skip_len = read_varint(f)
                        f.seek(skip_len, 1)
                    continue
                for _ in range(group_count - 1):
                    msg_size = read_varint(f)
                    yield f.read(msg_size)
            except:
                break

def process_alignment(raw_msg, wanted_nodes, chrom_filter):
    segs_by_node = defaultdict(list)
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw_msg)

    if chrom_filter and not any(p.name == chrom_filter for p in alignment.refpos):
        return segs_by_node, 0

    read_seq = alignment.sequence
    read_qual = alignment.quality
    mapq = alignment.mapping_quality
    read_off = 0
    read_count = 1

    for mapping in alignment.path.mapping:
        node_id = mapping.position.node_id
        if node_id not in wanted_nodes:
            for edit in mapping.edit:
                read_off += max(edit.from_length, len(edit.sequence))
            continue

        node_off = mapping.position.offset
        strand = "-" if mapping.position.is_reverse else "+"
        parts = []
        qual_bytes = bytearray()

        for edit in mapping.edit:
            if edit.from_length:
                frag = read_seq[read_off:read_off + edit.from_length]
                parts.append(frag.lower() if edit.sequence else frag)
                qual_bytes.extend(read_qual[read_off:read_off + edit.from_length])
                read_off += edit.from_length
            elif edit.sequence:
                l = len(edit.sequence)
                frag = read_seq[read_off:read_off + l]
                parts.append(frag.lower())
                qual_bytes.extend(read_qual[read_off:read_off + l])
                read_off += l

        segs_by_node[node_id].append({
            "offset": node_off,
            "seq": "".join(parts),
            "bq": bytes(qual_bytes),
            "rq": mapq,
            "strand": strand
        })

    return segs_by_node, read_count

def flush_node_to_sqlite(conn, nid, node_cache):
    segs_to_flush = node_cache[nid]
    if not segs_to_flush:
        return
    # Load existing blob if present
    cur = conn.execute("SELECT blob FROM segments WHERE node_id=?", (nid,))
    row = cur.fetchone()
    if row:
        old_list = pickle.loads(zlib.decompress(row[0]))
        old_list.extend(segs_to_flush)
    else:
        old_list = segs_to_flush
    new_blob = zlib.compress(pickle.dumps(old_list, protocol=5))
    conn.execute("INSERT OR REPLACE INTO segments VALUES (?, ?)", (nid, new_blob))
    node_cache[nid].clear()

def flush_all(conn, node_cache):
    for nid in list(node_cache.keys()):
        flush_node_to_sqlite(conn, nid, node_cache)
    conn.commit()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gam", help="GAM file path")
    parser.add_argument("stats", help="Pickle file with node stats")
    parser.add_argument("sqlite", help="Output SQLite path")
    parser.add_argument("--chr", default="", help="Optional chromosome name filter")
    parser.add_argument("--milestone", type=int, default=1_000_000)
    args = parser.parse_args()

    with open(args.stats, "rb") as f:
        stats = pickle.load(f)
    wanted_nodes = {
        int(nid) for nid, v in stats.items()
        if v["not_perfect"] > 1 and v["not_perfect"] / (v["not_perfect"] + v["perfect"]) > 0.1
    }
    total_nodes = len(stats)
    nodes_with_unperfect = sum(1 for s in stats.values() if s["not_perfect"] > 0)
    print("\nNode‑level overview")
    print(f"  Total nodes               : {total_nodes}")
    print(f"  Nodes with ≥1 un‑perfect  : {nodes_with_unperfect} "
          f"({nodes_with_unperfect/total_nodes*100:.2f}%)")
    print(f"  Nodes passing filter      : {len(wanted_nodes)} "
          f"({len(wanted_nodes)/total_nodes*100:.2f}%)\n")

    del stats
    gc.collect()

    conn = sqlite3.connect(args.sqlite)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("CREATE TABLE IF NOT EXISTS segments (node_id INTEGER PRIMARY KEY, blob BLOB)")

    total_reads = 0
    node_cache = defaultdict(list)
    start = time.perf_counter()
    milestone = args.milestone

    for raw_msg in gam_record_iter(args.gam):
        segs, count = process_alignment(raw_msg, wanted_nodes, args.chr)
        total_reads += count

        for nid, items in segs.items():
            node_cache[nid].extend(items)
            if len(node_cache[nid]) >= MAX_CACHE_SIZE:
                flush_node_to_sqlite(conn, nid, node_cache)

        if total_reads >= milestone:
            conn.commit()
            elapsed = time.perf_counter() - start
            print(f"{total_reads} reads processed | {elapsed:.1f} s")
            milestone += args.milestone

    flush_all(conn, node_cache)
    conn.close()
    print(f"✅ Done. Total reads: {total_reads}")

if __name__ == "__main__":
    main()
