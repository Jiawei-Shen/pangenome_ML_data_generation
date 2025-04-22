#!/usr/bin/env python3
import argparse
import gzip
import sqlite3
import pickle
import time
import gc
import vg_pb2

def read_varint(stream):
    value, shift = 0, 0
    while True:
        b = stream.read(1)
        if not b:
            raise EOFError
        value |= (b[0] & 0x7F) << shift
        if not (b[0] & 0x80):
            break
        shift += 7
    return value

def file_is_gzip(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"

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
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw_msg)

    if chrom_filter and not any(p.name == chrom_filter for p in alignment.refpos):
        return []

    read_seq = alignment.sequence
    read_qual = alignment.quality
    mapq = alignment.mapping_quality
    read_off = 0
    rows = []

    for mapping in alignment.path.mapping:
        node_id = mapping.position.node_id
        if node_id not in wanted_nodes:
            for edit in mapping.edit:
                read_off += max(edit.from_length, len(edit.sequence))
            continue

        node_off = mapping.position.offset
        strand = "-" if mapping.position.is_reverse else "+"
        parts, qual_bytes = [], []

        for edit in mapping.edit:
            if edit.from_length:
                frag = read_seq[read_off : read_off + edit.from_length]
                parts.append(frag.lower() if edit.sequence else frag)
                qual_bytes.extend(read_qual[read_off : read_off + edit.from_length])
                read_off += edit.from_length
            elif edit.sequence:
                l = len(edit.sequence)
                frag = read_seq[read_off : read_off + l]
                parts.append(frag.lower())
                qual_bytes.extend(read_qual[read_off : read_off + l])
                read_off += l

        rows.append((
            node_id,
            node_off,
            strand,
            "".join(parts),
            bytes(qual_bytes).hex(),
            mapq
        ))

    return rows

def flush_to_sqlite(conn, buffer):
    conn.execute("BEGIN")
    conn.executemany("""
        INSERT INTO segments (node_id, offset, strand, seq, bq, rq)
        VALUES (?, ?, ?, ?, ?, ?)
    """, buffer)
    conn.commit()

def run_pipeline(gam_path, stats_path, sqlite_path, milestone_step, chrom_filter):
    print(f"Loading stats: {stats_path}")
    with open(stats_path, "rb") as f:
        stats = pickle.load(f)

    wanted_nodes = {
        int(nid) for nid, v in stats.items()
        if v["not_perfect"] > 1 and v["not_perfect"] / (v["not_perfect"] + v["perfect"]) > 0.1
    }

    total_nodes = len(stats)
    nodes_with_unperfect = sum(1 for v in stats.values() if v["not_perfect"] > 0)

    print("\nNode-level overview")
    print(f"  Total nodes               : {total_nodes}")
    print(f"  Nodes with ≥1 un-perfect  : {nodes_with_unperfect} "
          f"({nodes_with_unperfect / total_nodes * 100:.2f} %)")
    print(f"  Nodes passing filter      : {len(wanted_nodes)} "
          f"({len(wanted_nodes) / total_nodes * 100:.2f} %)\n")

    del stats
    gc.collect()

    conn = sqlite3.connect(sqlite_path)
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("DROP TABLE IF EXISTS segments")
    conn.execute("""
        CREATE TABLE segments (
            node_id INTEGER,
            offset INTEGER,
            strand TEXT,
            seq TEXT,
            bq TEXT,
            rq INTEGER
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_node ON segments(node_id)")

    total_reads = 0
    buffer = []
    milestone = milestone_step
    start = time.perf_counter()

    for raw_msg in gam_record_iter(gam_path):
        rows = process_alignment(raw_msg, wanted_nodes, chrom_filter)
        buffer.extend(rows)
        total_reads += 1

        if total_reads >= milestone:
            flush_to_sqlite(conn, buffer)
            buffer.clear()
            elapsed = time.perf_counter() - start
            print(f"{total_reads} reads processed | {elapsed:.1f}s")
            milestone += milestone_step

    flush_to_sqlite(conn, buffer)
    conn.close()
    print(f"✅ Done. Total reads: {total_reads}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Pickle file with node stats")
    parser.add_argument("sqlite_output", help="Path to output SQLite file")
    parser.add_argument("--milestone", type=int, default=1_000_000, help="Milestone for flushing")
    parser.add_argument("--chr", default="", help="Optional chromosome filter")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        sqlite_path=args.sqlite_output,
        milestone_step=args.milestone,
        chrom_filter=args.chr
    )

if __name__ == "__main__":
    main()
