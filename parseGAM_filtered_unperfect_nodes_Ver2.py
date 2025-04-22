#!/usr/bin/env python3
import argparse
import gzip
import pickle
import time
import gc
import json
import shelve
from collections import defaultdict
import vg_pb2

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

def file_is_gzip(path):
    with open(path, "rb") as file_handle:
        return file_handle.read(2) == b"\x1f\x8b"

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
            except (EOFError, UnicodeDecodeError):
                break

def process_alignment(raw_message, wanted_nodes, chrom_filter):
    segment_dict = defaultdict(list)
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw_message)

    if chrom_filter and not any(pos.name == chrom_filter for pos in alignment.refpos):
        return segment_dict, 0

    read_sequence = alignment.sequence
    read_quality_bytes = alignment.quality
    mapping_quality = alignment.mapping_quality
    read_offset = 0
    reads = 1

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
                frag = read_sequence[read_offset: read_offset + edit.from_length]
                sequence_parts.append(frag.lower() if edit.sequence else frag)
                quality_bytes.extend(read_quality_bytes[read_offset: read_offset + edit.from_length])
                read_offset += edit.from_length
            elif edit.sequence:
                ins_len = len(edit.sequence)
                frag = read_sequence[read_offset: read_offset + ins_len]
                sequence_parts.append(frag.lower())
                quality_bytes.extend(read_quality_bytes[read_offset: read_offset + ins_len])
                read_offset += ins_len

        segment_dict[node_id].append({
            "offset": node_offset,
            "seq": "".join(sequence_parts),
            "bq": bytes(quality_bytes),
            "rq": mapping_quality,
            "strand": strand_char
        })

    return segment_dict, reads

def flush_segment_dict(shelf, buffer_dict):
    for node_id, segs in buffer_dict.items():
        key = str(node_id)
        if key in shelf:
            old = pickle.loads(shelf[key])
            old.extend(segs)
        else:
            old = segs
        shelf[key] = pickle.dumps(old, protocol=pickle.HIGHEST_PROTOCOL)
    buffer_dict.clear()

def run_pipeline(gam_path, stats_path, output_db, milestone_step, chrom_filter):
    print(f"Loading stats: {stats_path}")
    with open(stats_path, "rb") as stats_file:
        stats_data = pickle.load(stats_file)

    wanted_nodes = {
        int(node_id)
        for node_id, stat in stats_data.items()
        if stat["not_perfect"] > 1 and stat["not_perfect"] / (stat["perfect"] + stat["not_perfect"]) > 0.10
    }

    del stats_data
    gc.collect()
    print(f"Filtered {len(wanted_nodes)} nodes from stats file.")

    segment_buffer = defaultdict(list)
    total_reads = 0
    next_milestone = milestone_step
    start_time = time.perf_counter()

    with shelve.open(output_db, writeback=False) as db:
        for raw_msg in gam_record_iter(gam_path):
            seg_dict, count = process_alignment(raw_msg, wanted_nodes, chrom_filter)
            total_reads += count
            for node_id, segs in seg_dict.items():
                segment_buffer[node_id].extend(segs)

            if total_reads >= next_milestone:
                flush_segment_dict(db, segment_buffer)
                elapsed = time.perf_counter() - start_time
                print(f"{total_reads} reads processed | {elapsed:.1f} seconds")
                next_milestone += milestone_step

        flush_segment_dict(db, segment_buffer)

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")

def main():
    parser = argparse.ArgumentParser(description="Stream and store GAM read-segments to a shelf database.")
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the node stats pickle file")
    parser.add_argument("output_shelf", help="Output shelf DB file prefix")
    parser.add_argument("--milestone", type=int, default=10_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_db=args.output_shelf,
        milestone_step=args.milestone,
        chrom_filter=args.chr
    )

if __name__ == "__main__":
    main()
