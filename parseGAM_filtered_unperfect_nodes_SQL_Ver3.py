#!/usr/bin/env python3
import argparse
import gzip
import sqlite3
import vg_pb2
import gc
import time
import pickle
from collections import defaultdict

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
    entries = []
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw_msg)

    if chrom_filter and not any(p.name == chrom_filter for p in alignment.refpos):
        return [], 0

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
        qual_bytes = []

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

        entries.append((
            node_id,
            node_off,
            strand,
            mapq,
            "".join(parts),
            bytes(qual_bytes).hex()
        ))

    return entries, read_count

def flush_entries(conn, entries):
    if entries:
        conn.execute("BEGIN TRANSACTION")
        conn.executemany("""
            INSERT INTO segments (node_id, offset, strand, rq, seq, bq)
            VALUES (?, ?, ?, ?, ?, ?)
        """, entries)
        conn.execute("COMMIT")

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

    filtered_nodes = {
        int(nid) for nid, v in stats.items()
        if v["not_perfect"] > 1 and v["not_perfect"] / (v["perfect"] + v["not_perfect"]) > 0.1
    }
    total_nodes = len(stats)
    nodes_with_unperfect = sum(1 for v in stats.values() if v["not_perfect"] > 0)

    print("\nNode‑level overview")
    print(f"  Total nodes               : {total_nodes}")
    print(f"  Nodes with ≥1 un‑perfect  : {nodes_with_unperfect} "
          f"({nodes_with_unperfect / total_nodes * 100:.2f} %)")
    print(f"  Nodes passing filter      : {len(filtered_nodes)} "
          f"({len(filtered_nodes) / total_nodes * 100:.2f} %)\n")

    del stats
    gc.collect()

    conn = sqlite3.connect(args.sqlite)
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = OFF")  # Optional for speed
    conn.execute("DROP TABLE IF EXISTS segments")
    conn.execute("""
        CREATE TABLE segments (
            node_id INTEGER,
            offset INTEGER,
            strand TEXT,
            rq INTEGER,
            seq TEXT,
            bq TEXT
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_node ON segments(node_id)")

    total_reads = 0
    batch_entries = []
    milestone = args.milestone
    start = time.perf_counter()

    for raw_msg in gam_record_iter(args.gam):
        entries, count = process_alignment(raw_msg, filtered_nodes, args.chr)
        batch_entries.extend(entries)
        total_reads += count

        if total_reads >= milestone:
            flush_entries(conn, batch_entries)
            batch_entries.clear()
            elapsed = time.perf_counter() - start
            print(f"{total_reads} reads processed | {elapsed:.1f}s")
            milestone += args.milestone

    flush_entries(conn, batch_entries)
    conn.close()
    print(f"✅ Done. Total reads: {total_reads}")

if __name__ == "__main__":
    main()
