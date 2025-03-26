#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import sys
import subprocess
import tempfile
import vg_pb2
import json
import concurrent.futures
from google.protobuf.json_format import MessageToDict


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
        aln_dict = MessageToDict(alignment)

        for mapping in alignment.path.mapping:
            if mapping.position.node_id:
                node_id = mapping.position.node_id
                record = f"{node_id}\t{json.dumps(aln_dict)}"
                records.append(record)

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


def build_json_from_sorted_tsv(sorted_tsv_path, output_json_path):
    with open(sorted_tsv_path, "r") as infile, open(output_json_path, "w") as out_f:
        out_f.write("{\n")
        current_node = None
        current_reads = []

        for line in infile:
            node_id_str, json_read = line.strip().split("\t", 1)
            if node_id_str != current_node:
                if current_node is not None:
                    json_entry = f'  "{current_node}": {json.dumps(current_reads)},\n'
                    out_f.write(json_entry)
                current_node = node_id_str
                current_reads = [json.loads(json_read)]
            else:
                current_reads.append(json.loads(json_read))

        if current_node is not None:
            json_entry = f'  "{current_node}": {json.dumps(current_reads)}\n'
            out_f.write(json_entry)
        out_f.write("}\n")


def main():
    parser = argparse.ArgumentParser(description="Parse a GAM file, group by node ID, and filter non-perfect alignments.")
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=4, help="Number of worker processes (default: 4)")
    parser.add_argument("--max_pending", type=int, default=16, help="Max number of groups in flight (default: 16)")
    parser.add_argument("--chr", default="", help="Chromosome name to filter on (default: \"\")")
    parser.add_argument("--output_json", default="reads_by_node.json", help="Output JSON with node_id -> reads")
    parser.add_argument("--tmp_dir", default="./tmp", help="Directory to store temporary files (default: system temp)")
    parser.add_argument("--milestone", type=int, default=100_000_000,
                        help="Number of reads between progress updates (default: 100000000)")
    parser.add_argument("--sort_buffer", default="4G",
                        help="Buffer size for external sort (e.g. 1G, 512M). Default: 4G")

    args = parser.parse_args()

    start_time = time.perf_counter()
    total_count = 0
    perfect_count = 0
    not_perfect_count = 0
    milestone = args.milestone

    with tempfile.NamedTemporaryFile(mode='w+', delete=False, dir=args.tmp_dir) as temp_out:
        temp_path = temp_out.name

        for group_number, records, t, p, np in process_groups_pipeline(
                args.filename, args.threads, args.max_pending, chrom_name=args.chr):

            total_count += t
            perfect_count += p
            not_perfect_count += np

            for record in records:
                temp_out.write(record + "\n")

            if total_count >= milestone:
                elapsed = time.perf_counter() - start_time
                print("\nEarly Stop Summary:")
                print(f"  Total reads processed: {total_count}")
                print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
                print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
                print(f"Elapsed time: {elapsed:.2f} seconds.")
                milestone += args.milestone

    # Sort using external Unix sort with TMPDIR set
    sorted_path = temp_path + ".sorted"
    env = os.environ.copy()
    if args.tmp_dir:
        env["TMPDIR"] = args.tmp_dir

    sort_cmd = f"sort -k1,1n --buffer-size={args.sort_buffer} {temp_path} > {sorted_path}"
    subprocess.run(sort_cmd, shell=True, check=True, env=env)

    # Convert sorted TSV to final JSON
    build_json_from_sorted_tsv(sorted_path, args.output_json)

    # Clean up
    os.remove(temp_path)
    os.remove(sorted_path)

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_count}")
    print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
    print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
    print(f"Elapsed time: {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()
