#!/usr/bin/env python3
import argparse
import pickle
import time
import vg_pb2
import concurrent.futures


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


def parse_gam_file_groups(filename, group_size=1000):
    open_func = open
    if is_gzipped(filename):
        raise ValueError("Compressed GAM (.gz) not supported without gzip")
    with open_func(filename, 'rb') as f:
        while True:
            group = []
            try:
                for _ in range(group_size):
                    msg_len = read_varint(f)
                    msg = f.read(msg_len)
                    if not msg:
                        raise EOFError()
                    group.append(msg)
            except EOFError:
                if group:
                    yield group
                break
            yield group


def process_group_serialized(args):
    group_number, messages = args
    node_counts = {}

    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)

        node_edit_info = {}  # node_id -> list of is_perfect_chunk (per mapping)

        for mapping in alignment.path.mapping:
            node_id = mapping.position.node_id
            if not node_id:
                continue

            is_perfect = True
            for edit in mapping.edit:
                # If any edit has a sequence (i.e., substitution or insertion) or from_length == 0 (insertion)
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

            # All mappings to the node must be perfect for the read to count as "perfect" on that node
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
    parser.add_argument("input_gam", help="Input GAM file (must be uncompressed)")
    parser.add_argument("output_pickle", help="Output .pkl file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
    args = parser.parse_args()

    start_time = time.time()
    print("Reading alignments and launching threads...")

    grouped_messages = [(i, group) for i, group in enumerate(parse_gam_file_groups(args.input_gam))]
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
