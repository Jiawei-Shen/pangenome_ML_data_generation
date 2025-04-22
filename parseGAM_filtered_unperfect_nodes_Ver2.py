#!/usr/bin/env python3
import argparse
import gzip
import json
import pickle
import time
import gc
import vg_pb2

# ─────────────────────────────────────────────────────────────────────────────
# Utility: Read varint from stream
def read_varint(stream):
    value = 0
    shift_amount = 0
    while True:
        byte_pair = stream.read(1)
        if not byte_pair:
            raise EOFError
        byte_value = byte_pair[0]
        value |= (byte_value & 0x7F) << shift_amount
        if not (byte_value & 0x80):
            return value
        shift_amount += 7

# Utility: Check if a file is gzipped
def file_is_gzip(path):
    with open(path, "rb") as file_handle:
        return file_handle.read(2) == b"\x1f\x8b"

# Generator: Yield groups of GAM records from file
def gam_group_iter(path, tag="GAM"):
    open_func = gzip.open if file_is_gzip(path) else open
    with open_func(path, "rb") as file_handle:
        while True:
            try:
                group_count = read_varint(file_handle)
            except EOFError:
                break
            if group_count == 0:
                continue
            try:
                tag_length = read_varint(file_handle)
                group_tag = file_handle.read(tag_length).decode()
            except (EOFError, UnicodeDecodeError):
                break
            if group_tag != tag:
                for _ in range(group_count - 1):
                    try:
                        skip_length = read_varint(file_handle)
                        file_handle.seek(skip_length, 1)
                    except EOFError:
                        break
                continue
            message_list = []
            for _ in range(group_count - 1):
                try:
                    message_size = read_varint(file_handle)
                    message_bytes = file_handle.read(message_size)
                except EOFError:
                    break
                if len(message_bytes) == message_size:
                    message_list.append(message_bytes)
            if message_list:
                yield message_list

# ─────────────────────────────────────────────────────────────────────────────
# Process one group of alignments
def process_group(message_list, wanted_nodes, chrom_filter):
    segment_dict = {}
    read_count = 0

    for raw_message in message_list:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(raw_message)

        if chrom_filter and not any(pos.name == chrom_filter for pos in alignment.refpos):
            continue

        read_count += 1
        read_sequence = alignment.sequence
        read_quality_bytes = alignment.quality
        mapping_quality = alignment.mapping_quality
        read_offset = 0

        for mapping in alignment.path.mapping:
            node_id = mapping.position.node_id

            if node_id not in wanted_nodes:
                for edit in mapping.edit:
                    read_offset += max(edit.from_length, len(edit.sequence))
                continue

            node_offset = mapping.position.offset
            strand_char = "-" if mapping.position.is_reverse else "+"
            sequence_parts = []
            quality_bytes = bytearray()

            for edit in mapping.edit:
                if edit.from_length:
                    fragment = read_sequence[read_offset : read_offset + edit.from_length]
                    sequence_parts.append(fragment.lower() if edit.sequence else fragment)
                    quality_bytes.extend(read_quality_bytes[read_offset : read_offset + edit.from_length])
                    read_offset += edit.from_length
                elif edit.sequence:
                    insertion_length = len(edit.sequence)
                    fragment = read_sequence[read_offset : read_offset + insertion_length]
                    sequence_parts.append(fragment.lower())
                    quality_bytes.extend(read_quality_bytes[read_offset : read_offset + insertion_length])
                    read_offset += insertion_length

            # Record this read‑segment for the current node
            segment_dict.setdefault(node_id, []).append({
                "offset": node_offset,                   # Start within node
                "sequence": "".join(sequence_parts),       # Bases (insertions lower‑cased)
                "base_quality_hex": quality_bytes.hex(),           # Per‑base Phred scores, hex-encoded
                "read_quality": mapping_quality,               # MAPQ from the read
                "strand": strand_char                    # "+" or "-" on node
            })

            # segment_dict.setdefault(node_id, []).append({
            #     "offset": node_offset,  # Start within node
            #     "sequence": "".join(sequence_parts),  # Bases (insertions lower‑cased)
            #     "base_quality_hex": quality_bytes.hex(),  # Per‑base Phred scores, hex-encoded
            #     "read_quality": mapping_quality,  # MAPQ from the read
            #     "strand": strand_char  # "+" or "-" on node
            # })

            segment_dict.setdefault(node_id, []).append({
                "offset": node_offset,  # offset: Start within node
                "seq": "".join(sequence_parts),  # sequence: Bases (insertions lower‑cased)
                "bq": quality_bytes.hex(),  # base_quality_hex: Per‑base Phred scores, hex-encoded
                "rq": mapping_quality,  # read_quality: MAPQ from the read
                "strand": strand_char  # strand: "+" or "-" on node
            })

    return segment_dict, read_count

# ─────────────────────────────────────────────────────────────────────────────
# Merge multiple segment dictionaries
def merge_partials(partial_list):
    merged_dict = {}
    total_reads = 0
    for segment_dict, read_count in partial_list:
        total_reads += read_count
        for node_id, segment_list in segment_dict.items():
            merged_dict.setdefault(node_id, []).extend(segment_list)
    return merged_dict, total_reads

# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
def run_pipeline(gam_path, stats_path, output_prefix, output_format, milestone_step, chrom_filter):
    print(f"Loading stats: {stats_path}")
    with open(stats_path, "rb") as stats_file:
        stats_data = pickle.load(stats_file)

    wanted_nodes = {
        int(node_id)
        for node_id, stat in stats_data.items()
        if stat["not_perfect"] > 1 and stat["not_perfect"] / (stat["perfect"] + stat["not_perfect"]) > 0.10
    }

    del stats_data
    gc.collect()  # Release memory

    print(f"Filtered {len(wanted_nodes)} nodes from stats file.")

    partial_results = []
    total_reads = 0
    next_milestone = milestone_step
    start_time = time.perf_counter()

    for message_batch in gam_group_iter(gam_path):
        result = process_group(message_batch, wanted_nodes, chrom_filter)
        partial_results.append(result)
        total_reads += result[1]

        if total_reads >= next_milestone:
            elapsed = time.perf_counter() - start_time
            print(f"{total_reads} reads processed | {elapsed:.1f} seconds")
            next_milestone += milestone_step

    merged_segments, _ = merge_partials(partial_results)

    if output_format in ("pkl", "both"):
        pickle_path = output_prefix + ".pkl"
        with open(pickle_path, "wb") as pickle_file:
            pickle.dump(merged_segments, pickle_file, pickle.HIGHEST_PROTOCOL)
        print(f"Wrote Pickle: {pickle_path}")

    if output_format in ("json", "both"):
        json_path = output_prefix + ".json"
        with open(json_path, "w") as json_file:
            json.dump({str(k): v for k, v in merged_segments.items()}, json_file)
        print(f"Wrote JSON: {json_path}")

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads}")
    print(f"  Nodes included: {len(merged_segments)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")

# ─────────────────────────────────────────────────────────────────────────────
# Command-line interface
def main():
    parser = argparse.ArgumentParser(description="Extract read-segments by node from a GAM file.")
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the node stats pickle file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--fmt", choices=["json", "pkl", "both"], default="json", help="Output format")
    parser.add_argument("--milestone", type=int, default=10_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        output_format=args.fmt,
        milestone_step=args.milestone,
        chrom_filter=args.chr
    )

if __name__ == "__main__":
    main()
