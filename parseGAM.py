#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import threading
import queue
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


def parse_gam_file_groups(filename, expected_tag="GAM"):
    """
    Generator that traverses the GAM file group by group.

    Each group starts with a varint count and a type tag. If the tag matches
    expected_tag, yield a tuple:
       (group_number, list of raw message bytes for each alignment in the group)
    """
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
            print(f"[Producer] Reading Group {group_number}: {group_count} messages")

            if group_count == 0:
                continue

            try:
                tag_size = read_varint(f)
            except EOFError:
                print("[Producer] Unexpected EOF when reading type tag size.")
                break
            tag_bytes = f.read(tag_size)
            if len(tag_bytes) != tag_size:
                print("[Producer] Unexpected EOF when reading type tag bytes.")
                break
            tag = tag_bytes.decode("utf-8")
            print(f"[Producer]  Type tag: {tag}")

            if tag != expected_tag:
                print(f"[Producer]  Skipping group with unexpected tag: {tag}")
                # Skip remaining messages in this group.
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                    except EOFError:
                        break
                    f.seek(msg_size, os.SEEK_CUR)
                continue
            else:
                messages = []
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                    except EOFError:
                        print("[Producer] Unexpected EOF when reading message size.")
                        break
                    msg_bytes = f.read(msg_size)
                    if len(msg_bytes) != msg_size:
                        print("[Producer] Unexpected EOF when reading message bytes.")
                        break
                    messages.append(msg_bytes)
                yield (group_number, messages)
    finally:
        f.close()


def process_group(group_tuple):
    """
    Parse each raw message in the group into a vg_pb2.Alignment.
    Returns a tuple (group_number, list of parsed alignments).
    """
    group_number, messages = group_tuple
    alignments = []
    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)
        alignments.append(alignment)
    return (group_number, alignments)


def group_producer(filename, q, expected_tag="GAM"):
    """
    Read groups from the GAM file and put each group into the queue.
    After finishing, put a sentinel (None) into the queue.
    """
    for group in parse_gam_file_groups(filename, expected_tag):
        q.put(group)
    q.put(None)  # Sentinel to signal completion


def group_consumer(q):
    """
    Continuously get groups from the queue, process them, and output results.
    Exits when it sees the sentinel (None).
    """
    while True:
        group = q.get()
        if group is None:
            # Propagate the sentinel to other consumers.
            q.put(None)
            q.task_done()
            break
        group_number, alignments = process_group(group)
        print(f"[Consumer] Processed Group {group_number}: Parsed {len(alignments)} alignments")
        for alignment in alignments:
            print(f"  Alignment name: {alignment.name}")
        q.task_done()


def main():
    parser = argparse.ArgumentParser(
        description="Pipeline multi-threaded parsing of a GAM file by group using vg_pb2.")
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=32,
                        help="Number of consumer threads to use (default: 32)")
    parser.add_argument("--max_pending", type=int, default=100,
                        help="Maximum number of groups to buffer (default: 100)")
    args = parser.parse_args()

    start_time = time.perf_counter()
    q = queue.Queue(maxsize=args.max_pending)

    # Start the producer thread.
    prod_thread = threading.Thread(target=group_producer, args=(args.filename, q, "GAM"))
    prod_thread.start()

    # Start consumer threads.
    consumers = []
    for i in range(args.threads):
        t = threading.Thread(target=group_consumer, args=(q,))
        t.start()
        consumers.append(t)

    # Wait for the producer to finish.
    prod_thread.join()
    # Wait until the queue is fully processed.
    q.join()
    # Wait for all consumers to finish.
    for t in consumers:
        t.join()

    end_time = time.perf_counter()
    print(f"Elapsed time: {end_time - start_time:.6f} seconds")


if __name__ == "__main__":
    main()
