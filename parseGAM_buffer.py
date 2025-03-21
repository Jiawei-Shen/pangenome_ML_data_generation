#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import io
import concurrent.futures
import vg_pb2  # Import the generated protobuf module

def is_gzipped(filename):
    """Check if the file is gzipped by reading its magic number."""
    with open(filename, 'rb') as f:
        magic = f.read(2)
    return magic == b'\x1f\x8b'

def read_varint(stream):
    """Read a varint from a stream and return its value."""
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

def scan_file_groups(filename, expected_tag="GAM"):
    """
    Scan the GAM file sequentially and record the metadata for each group.
    Each group is recorded as a tuple:
         (group_number, group_offset, group_size)
    The group_offset is the file position at the start of the group,
    and group_size is the number of bytes that belong to that group.
    """
    groups = []
    # Use BGZF access if gzipped and random access is needed.
    if is_gzipped(filename):
        try:
            import pysam
            f = pysam.BGZFFile(filename, "rb")
        except ImportError:
            raise RuntimeError("File is gzipped but pysam is not installed.")
    else:
        f = open(filename, "rb")
    try:
        group_number = 0
        while True:
            group_start = f.tell()
            try:
                group_count = read_varint(f)
            except EOFError:
                break  # Reached end of file
            # If group_count is 0, skip it.
            if group_count == 0:
                continue
            group_number += 1
            # Read type tag
            try:
                tag_size = read_varint(f)
            except EOFError:
                break
            tag_bytes = f.read(tag_size)
            if len(tag_bytes) != tag_size:
                break
            tag = tag_bytes.decode("utf-8")
            # Process (or skip) messages in the group based on tag.
            if tag != expected_tag:
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                    except EOFError:
                        break
                    f.seek(msg_size, os.SEEK_CUR)
            else:
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                    except EOFError:
                        break
                    f.seek(msg_size, os.SEEK_CUR)
            group_end = f.tell()
            group_size = group_end - group_start
            groups.append((group_number, group_start, group_size))
    finally:
        f.close()
    return groups

def process_group_from_file(filename, group_meta, expected_tag="GAM"):
    """
    Open the file, seek to the group's offset, read the group's bytes,
    and parse the group into vg_pb2.Alignment objects.
    Returns (group_number, list of alignments).
    """
    group_number, offset, size = group_meta
    # Open file for random access (use pysam.BGZFFile if gzipped)
    if is_gzipped(filename):
        import pysam
        f = pysam.BGZFFile(filename, "rb")
    else:
        f = open(filename, "rb")
    try:
        f.seek(offset)
        data = f.read(size)
    finally:
        f.close()
    stream = io.BytesIO(data)
    try:
        group_count = read_varint(stream)
    except EOFError:
        return (group_number, [])
    try:
        tag_size = read_varint(stream)
    except EOFError:
        return (group_number, [])
    tag_bytes = stream.read(tag_size)
    if len(tag_bytes) != tag_size:
        return (group_number, [])
    tag = tag_bytes.decode("utf-8")
    alignments = []
    if tag != expected_tag:
        # If tag doesn't match, skip this group.
        return (group_number, [])
    for _ in range(group_count - 1):
        try:
            msg_size = read_varint(stream)
        except EOFError:
            break
        msg_bytes = stream.read(msg_size)
        if len(msg_bytes) != msg_size:
            break
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)
        alignments.append(alignment)
    return (group_number, alignments)

def main():
    parser = argparse.ArgumentParser(
        description="Multi-threaded parsing of a GAM file by processing groups concurrently using random access.")
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of worker threads (default: 4)")
    args = parser.parse_args()

    start_time = time.perf_counter()
    groups = scan_file_groups(args.filename, expected_tag="GAM")
    print(f"Found {len(groups)} groups in the file.")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_group = {executor.submit(process_group_from_file, args.filename, meta, "GAM"): meta for meta in groups}
        for future in concurrent.futures.as_completed(future_to_group):
            group_number, alignments = future.result()
            print(f"Processed Group {group_number}: Parsed {len(alignments)} alignments.")
            for aln in alignments:
                print(f"  Alignment name: {aln.name}")
    end_time = time.perf_counter()
    print(f"Elapsed time: {end_time - start_time:.6f} seconds")

if __name__ == "__main__":
    main()
