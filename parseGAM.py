#!/usr/bin/env python3
import argparse
import io
import os
import sys
import concurrent.futures

# Try to use pysam for BGZF support.
try:
    import pysam

    USE_PYSAM = True
except ImportError:
    import gzip

    USE_PYSAM = False

import vg_pb2  # Import the generated protobuf module


def read_varint(stream):
    """Read a varint from the stream."""
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


def open_file(filename):
    """
    Open the file for reading in binary mode.
    If pysam is available, use pysam.BGZFFile to properly handle BGZF-compressed files.
    Otherwise, if the file is gzipped, use gzip.open.
    """
    if USE_PYSAM:
        return pysam.BGZFFile(filename, "rb")
    else:
        # Fallback to gzip if pysam is not installed.
        with open(filename, 'rb') as f:
            magic = f.read(2)
        if magic == b'\x1f\x8b':
            import gzip
            return gzip.open(filename, "rb")
        else:
            return open(filename, "rb")


def parse_gam_file_groups(filename, expected_tag="GAM"):
    """
    Parse the GAM file using the framing format and yield groups.
    Each group is a tuple (group_number, messages) where messages is a list
    of raw protobuf bytes for each Alignment.
    """
    f = open_file(filename)
    # Wrap in BufferedReader for faster I/O.
    f = io.BufferedReader(f)

    group_number = 0
    try:
        while True:
            try:
                group_count = read_varint(f)
            except EOFError:
                break  # End of file reached
            group_number += 1
            print(f"Group {group_number}: {group_count} messages")

            if group_count == 0:
                print("Empty group.")
                continue

            try:
                tag_size = read_varint(f)
            except EOFError:
                print("Unexpected EOF when reading type tag size.")
                break
            tag_bytes = f.read(tag_size)
            if len(tag_bytes) != tag_size:
                print("Unexpected EOF when reading type tag bytes.")
                break
            tag = tag_bytes.decode("utf-8")
            print(f"  Type tag: {tag}")

            if tag != expected_tag:
                print("  Skipping group with unexpected tag.")
                # Skip the remaining messages in this group.
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                    except EOFError:
                        break
                    f.seek(msg_size, os.SEEK_CUR)
                continue
            else:
                messages = []
                for i in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                    except EOFError:
                        print("Unexpected EOF when reading message size.")
                        break
                    msg_bytes = f.read(msg_size)
                    if len(msg_bytes) != msg_size:
                        print("Unexpected EOF when reading message bytes.")
                        break
                    messages.append(msg_bytes)
                yield (group_number, messages)
    finally:
        f.close()


def process_group(group_tuple):
    """
    Process a group of raw message bytes by parsing each one into a vg_pb2.Alignment.
    Returns a tuple (group_number, list of parsed alignments).
    """
    group_number, messages = group_tuple
    alignments = []
    for i, msg_bytes in enumerate(messages):
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)
        alignments.append(alignment)
    return (group_number, alignments)


def main():
    parser = argparse.ArgumentParser(
        description="Parse a GAM file using vg_pb2 in parallel and print alignments.")
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of worker threads to use (default: 4)")
    args = parser.parse_args()

    groups = list(parse_gam_file_groups(args.filename))
    print(f"Found {len(groups)} groups with expected tag.")

    # Use a ThreadPoolExecutor to process groups in parallel.
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_group = {executor.submit(process_group, group): group for group in groups}
        for future in concurrent.futures.as_completed(future_to_group):
            group_number, alignments = future.result()
            print(f"Processed Group {group_number}: Parsed {len(alignments)} alignments")
            for alignment in alignments:
                print(f"  Alignment name: {alignment.name}")


if __name__ == "__main__":
    main()
