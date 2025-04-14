#!/usr/bin/env python3
import argparse
import pickle
import time
import vg_pb2
import concurrent.futures
import os


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


def parse_gam_file_groups(filename, expected_tag="GAM"):
    with open(filename, 'rb') as f:
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
                # Skip group content
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

            yield group_number, messages


def process_group_serialized(args):
    group_number, messages = args
    node_counts = {}

    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)

        node_edit_info = {}  # node_id -> list of is_perfect_chunk

        for mapping in alignment.path.mapping:
            node_id = mapping.position.node_id
            if not node_id:
                continue

            is_perfect = True
            for edit in mapping.edit:
                if edit.sequence or edit.from_length == 0:
                    is_perfect = False
                    break

            node_id = str(node_id)
            if node_id not in node_edit_info:
                node_edit_info[node_id] = []
            node_edit_info[node_id].append(is_perfect)

        for node_id, chunk_results in node_edit_info.items():
            if node_id not in node_counts:
                node_counts[node_id] = {"perfect": 0, "not_perfect": 0}
            if all(chunk_results):
                node_counts[node_id]["perfect"] += 1
            else:
                node_counts[node_id]["not_perfect"] += 1

    return group_number, node_counts


def merge_node_counts(all_counts):
    merged = {}
    for node_counts in all_counts:
        for node_id, stats in node_counts.items():
            if node_id not in merged:
                merged[node_id] = {"perfect": 0, "not_perfect": 0}
            merged[node_id]["perfect"] += stats["perfect"]
            merged[node_id]["not_perfect"] += stats["not_perfect"]
    return merged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_gam", help="Input GAM file (varint+GAM-tagged format)")
    parser.add_argument("output_pickle", help="Output .pkl file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
    args = parser.parse_args()

    start_time = time.time()
    print("Reading alignments and launching threads...")

    grouped_messages = [(i, group) for i, group in parse_gam_file_groups(args.input_gam)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = list(executor.map(process_group_serialized, grouped_messages))

    all_node_counts = [node_counts for _, node_counts in results]
    merged_node_counts = merge_node_counts(all_node_counts)

    print(f"Writing pickle output to {args.output_pickle}...")
    with open(args.output_pickle, "wb") as f:
        pickle.dump(merged_node_counts, f, protocol=pickle.HIGHEST_PROTOCOL)

    elapsed = time.time() - start_time
    print(f"Done in {elapsed:.2f} seconds")


if __name__ == "__main__":
    main()
