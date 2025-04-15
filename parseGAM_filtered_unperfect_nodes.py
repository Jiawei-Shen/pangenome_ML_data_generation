#!/usr/bin/env python3
import argparse
import pickle
import gzip
import os
import concurrent.futures
import vg_pb2
import json
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
        while True:
            try:
                group_count = read_varint(f)
            except EOFError:
                break

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
                        f.seek(msg_size, 1)
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

            yield messages


def process_group_extract(args):
    messages, filtered_nodes = args
    result = {}

    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)

        read_seq = alignment.sequence
        base_quality = base64.b64decode(alignment.quality).decode('latin1') if alignment.quality else ""
        read_quality = alignment.mapping_quality
        read_offset = 0

        for mapping in alignment.path.mapping:
            node_id = mapping.position.node_id
            if node_id not in filtered_nodes:
                for edit in mapping.edit:
                    read_offset += max(edit.from_length, len(edit.sequence))
                continue

            offset = mapping.position.offset
            node_seq = ""
            node_bqual = ""

            for edit in mapping.edit:
                aligned_len = edit.from_length

                if aligned_len > 0:
                    seg_seq = read_seq[read_offset:read_offset + aligned_len]
                    seg_qual = base_quality[read_offset:read_offset + aligned_len]
                    if edit.sequence:
                        seg_seq = seg_seq.lower()
                    node_seq += seg_seq
                    node_bqual += seg_qual
                    read_offset += aligned_len

                elif edit.sequence:
                    insert_len = len(edit.sequence)
                    seg_seq = read_seq[read_offset:read_offset + insert_len].lower()
                    seg_qual = base_quality[read_offset:read_offset + insert_len]
                    node_seq += seg_seq
                    node_bqual += seg_qual
                    read_offset += insert_len

            if node_id not in result:
                result[node_id] = []

            result[node_id].append({
                "offset": offset,
                "sequence": node_seq,
                "base_quality": node_bqual,
                "read_quality": read_quality
            })

    return result


def merge_node_results(results_list):
    merged = {}
    for partial in results_list:
        for node_id, segments in partial.items():
            if node_id not in merged:
                merged[node_id] = []
            merged[node_id].extend(segments)
    return merged


def extract_node_segments_parallel(gam_file, node_stats_pickle, output_prefix, threads=4, max_pending=16):
    with open(node_stats_pickle, "rb") as f:
        node_stats = pickle.load(f)

    filtered_nodes = {
        int(node_id)
        for node_id, stats in node_stats.items()
        if stats["not_perfect"] > 1 and stats["not_perfect"] / (stats["not_perfect"] + stats["perfect"]) > 0.1
    }
    print(f"Filtered {len(filtered_nodes)} node IDs meeting criteria.")

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for messages in parse_gam_file_groups(gam_file):
            futures.append(executor.submit(process_group_extract, (messages, filtered_nodes)))
            while len(futures) >= max_pending:
                done, futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    yield fut.result()

        for fut in concurrent.futures.as_completed(futures):
            yield fut.result()


def main():
    parser = argparse.ArgumentParser(description="Extract read segments mapping to nodes with high not-perfect ratio.")
    parser.add_argument("gam_file", help="Input GAM file")
    parser.add_argument("node_stats_pickle", help="Pickle file of node -> {perfect, not_perfect}")
    parser.add_argument("output_prefix", help="Prefix for output files (.json and .pkl)")
    parser.add_argument("--threads", type=int, default=4, help="Number of worker threads (default: 4)")
    parser.add_argument("--max_pending", type=int, default=16, help="Maximum number of futures pending (default: 16)")
    args = parser.parse_args()

    all_results = []
    for partial_result in extract_node_segments_parallel(
            args.gam_file, args.node_stats_pickle, args.output_prefix, args.threads, args.max_pending):
        all_results.append(partial_result)

    merged = merge_node_results(all_results)

    json_path = args.output_prefix + ".json"
    pkl_path = args.output_prefix + ".pkl"

    with open(json_path, "w") as f:
        json.dump(merged, f, indent=2)
    print(f"JSON saved to {json_path}")

    with open(pkl_path, "wb") as f:
        pickle.dump(merged, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Pickle saved to {pkl_path}")


if __name__ == "__main__":
    main()
