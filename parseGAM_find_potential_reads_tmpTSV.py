#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import sys
import tempfile
import vg_pb2
import json
import sqlite3
import concurrent.futures
from google.protobuf.json_format import MessageToDict
import base64



def read_varint(stream):
    result = 0
    shift = 0
    while True:
        b = stream.read(1)
        if len(b) == 0:
            raise EOFError("Unexpected EOF while reading varint")
        byte = b[0]
        result |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            break
        shift += 7
    return result


def is_gzipped(filename):
    with open(filename, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def parse_gam_file_groups(filename, expected_tag="GAM"):
    open_func = gzip.open if is_gzipped(filename) else open
    with open_func(filename, 'rb') as f:
        group_number = 0
        while True:
            try:
                group_count = read_varint(f)
            except EOFError:
                break

            group_number += 1
            if group_count == 0:
                continue

            try:
                tag_size = read_varint(f)
                tag = f.read(tag_size).decode("utf-8")
            except (EOFError, UnicodeDecodeError):
                break

            if tag != expected_tag:
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                        f.seek(msg_size, os.SEEK_CUR)
                    except EOFError:
                        break
                continue

            messages = []
            for _ in range(group_count - 1):
                try:
                    msg_size = read_varint(f)
                    msg_bytes = f.read(msg_size)
                    if len(msg_bytes) == msg_size:
                        messages.append(msg_bytes)
                except EOFError:
                    break

            yield (group_number, messages)


def is_on_chromosome(alignment, chrom_name):
    return any(rp.name == chrom_name for rp in alignment.refpos)


import json

def process_group_serialized(args):
    group_number, messages, chrom_name = args
    total = 0
    perfect = 0
    not_perfect = 0
    records = []

    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)

        if chrom_name and not is_on_chromosome(alignment, chrom_name):
            continue

        total += 1
        if alignment.identity == 1.0:
            perfect += 1
            continue

        not_perfect += 1
        read_name = alignment.name
        read_seq = alignment.sequence
        read_qual = alignment.quality.decode('latin1')  # already decoded from Protobuf
        read_offset = 0

        for mapping in alignment.path.mapping:
            if not mapping.position.node_id:
                continue

            node_id = mapping.position.node_id
            offset = mapping.position.offset if mapping.HasField("position") else 0

            node_seq = ""
            node_qual = ""

            for edit in mapping.edit:
                aligned_len = edit.from_length

                if aligned_len > 0:
                    # Get the matching segment from the read
                    seg_seq = read_seq[read_offset:read_offset + aligned_len]
                    seg_qual = read_qual[read_offset:read_offset + aligned_len]

                    # If it's a substitution, mark it with lowercase
                    if edit.sequence:
                        seg_seq = seg_seq.lower()

                    node_seq += seg_seq
                    node_qual += seg_qual
                    read_offset += aligned_len

                elif edit.sequence:
                    # Insertion: take and lowercase inserted read segment
                    insert_len = len(edit.sequence)
                    seg_seq = read_seq[read_offset:read_offset + insert_len].lower()
                    seg_qual = read_qual[read_offset:read_offset + insert_len]

                    node_seq += seg_seq
                    node_qual += seg_qual
                    read_offset += insert_len

            record_data = {
                "name": read_name,
                "node_id": node_id,
                "offset": offset,
                "sequence": node_seq,
                "quality": node_qual
            }

            records.append(f"{node_id}\t{json.dumps(record_data)}")

    return group_number, records, total, perfect, not_perfect



def process_groups_pipeline(filename, threads, max_pending=10, chrom_name=""):
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        pending_futures = []
        for group in parse_gam_file_groups(filename):
            args = (group[0], group[1], chrom_name)
            future = executor.submit(process_group_serialized, args)
            pending_futures.append(future)

            if len(pending_futures) >= max_pending:
                done, not_done = concurrent.futures.wait(
                    pending_futures,
                    return_when=concurrent.futures.FIRST_COMPLETED
                )
                for completed in done:
                    yield completed.result()
                pending_futures = list(not_done)

        for future in concurrent.futures.as_completed(pending_futures):
            yield future.result()


def setup_sqlite_db(db_path):
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute("CREATE TABLE IF NOT EXISTS pileup (node_id INTEGER, read_json TEXT)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_node ON pileup (node_id)")
    return conn


def export_sqlite_to_json(conn, output_path):
    cursor = conn.cursor()
    cursor.execute("SELECT node_id, json_group_array(read_json) FROM pileup GROUP BY node_id")

    with open(output_path, "w") as out_f:
        out_f.write("{\n")
        first = True
        for node_id, json_array in cursor:
            if not first:
                out_f.write(",\n")
            first = False
            out_f.write(f'  "{node_id}": {json_array}')
        out_f.write("\n}\n")


def main():
    parser = argparse.ArgumentParser(description="Parse a GAM file, group by node ID, and filter non-perfect alignments using SQLite.")
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=4, help="Number of worker processes (default: 4)")
    parser.add_argument("--max_pending", type=int, default=16, help="Max number of groups in flight (default: 16)")
    parser.add_argument("--chr", default="", help="Chromosome name to filter on (default: \"\")")
    parser.add_argument("--output_json", default="reads_by_node.json", help="Output JSON with node_id -> reads")
    parser.add_argument("--tmp_dir", default='./tmp', help="Directory to store temporary SQLite DB (default: system temp)")
    parser.add_argument("--milestone", type=int, default=100_000_000,
                        help="Number of reads between progress updates (default: 100000000)")

    args = parser.parse_args()

    if args.tmp_dir:
        os.makedirs(args.tmp_dir, exist_ok=True)

    start_time = time.perf_counter()
    total_count = 0
    perfect_count = 0
    not_perfect_count = 0
    milestone = args.milestone

    db_path = os.path.join(args.tmp_dir or ".", "pileup.sqlite")
    conn = setup_sqlite_db(db_path)
    cursor = conn.cursor()

    for group_number, records, t, p, np in process_groups_pipeline(
            args.filename, args.threads, args.max_pending, chrom_name=args.chr):

        total_count += t
        perfect_count += p
        not_perfect_count += np

        to_insert = []
        for record in records:
            node_id_str, read_json = record.split("\t", 1)
            to_insert.append((int(node_id_str), read_json))

        cursor.executemany("INSERT INTO pileup (node_id, read_json) VALUES (?, ?)", to_insert)

        if total_count >= milestone:
            elapsed = time.perf_counter() - start_time
            print("\nEarly Stop Summary:")
            print(f"  Total reads processed: {total_count}")
            print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
            print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
            print(f"Elapsed time: {elapsed:.2f} seconds.")
            milestone += args.milestone

    conn.commit()
    export_sqlite_to_json(conn, args.output_json)
    conn.close()
    os.remove(db_path)

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_count}")
    print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
    print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
    print(f"Elapsed time: {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()
