#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import sys
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


def is_gzipped(filename):
    """Check if the file is gzipped by reading its magic number."""
    with open(filename, 'rb') as f:
        magic = f.read(2)
    return magic == b'\x1f\x8b'


def parse_gam_file(filename, expected_tag="GAM"):
    """
    Parse a GAM file with the framing format.

    Each group starts with a varint count and a type tag. If the tag matches
    the expected_tag (default "GAM"), then the subsequent messages in the group
    are parsed as vg_pb2.Alignment objects.
    """
    # Open the file; if gzipped, use gzip.open
    if is_gzipped(filename):
        f = gzip.open(filename, 'rb')
    else:
        f = open(filename, 'rb')

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
                # Skip remaining messages in this group
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                    except EOFError:
                        break
                    f.seek(msg_size, os.SEEK_CUR)
                continue
            else:
                # For each remaining message in the group, decode it as an Alignment.
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
                    alignment = vg_pb2.Alignment()
                    alignment.ParseFromString(msg_bytes)
                    print(f"   Message {i + 1}: {msg_size} bytes, Alignment name: {alignment.name}")
                    yield alignment
    finally:
        f.close()


def main():
    parser = argparse.ArgumentParser(description="Parse a GAM file using vg_pb2 and print alignments.")
    parser.add_argument("filename", help="Path to the GAM file")
    args = parser.parse_args()

    count = 0
    for alignment in parse_gam_file(args.filename):
        # Process the alignment as needed. Here we simply print the full message.
        print(f"\rProgress: {count}", end="", flush=True)
        count += 1
        # print("Parsed Alignment:")
        # print(alignment)


if __name__ == "__main__":
    # main()
    start_time = time.perf_counter()
    # Your code here
    main()  # Simulating delay
    end_time = time.perf_counter()

    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.6f} seconds")