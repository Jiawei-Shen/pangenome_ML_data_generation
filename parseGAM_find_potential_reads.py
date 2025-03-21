#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import sys
import concurrent.futures
import vg_pb2  # Import the generated protobuf module

# For dumping protobufs to JSON
import json
from google.protobuf.json_format import MessageToDict


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
    Generator that traverses the GAM file group by group, yielding tuples:
    (group_number, list_of_raw_message_bytes).
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
            if group_count == 0:
                # No messages in this group; move on
                continue

            # Read the type tag
            try:
                tag_size = read_varint(f)
            except EOFError:
                break
            tag_bytes = f.read(tag_size)
            if len(tag_bytes) != tag_size:
                break
            tag = tag_bytes.decode("utf-8")

            if tag != expected_tag:
                # Skip messages in groups we don't care about
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                    except EOFError:
                        break
                    f.seek(msg_size, os.SEEK_CUR)
                continue

            # If it's the correct tag, read all message bytes
            messages = []
            for _ in range(group_count - 1):
                try:
                    msg_size = read_varint(f)
                except EOFError:
                    break
                msg_bytes = f.read(msg_size)
                if len(msg_bytes) != msg_size:
                    break
                messages.append(msg_bytes)

            yield (group_number, messages)
    finally:
        f.close()


def process_group(group_tuple):
    """
    Parse each raw message in the group into a vg_pb2.Alignment.
    Returns (group_number, [parsed_alignments]).
    """
    group_number, messages = group_tuple
    alignments = []
    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)
        alignments.append(alignment)
    return (group_number, alignments)


def process_groups_pipeline(filename, threads, max_pending=10):
    """
    Streams through the GAM file, reading each group, then processing
    in a thread pool. Yields (group_number, [alignments]) as they complete.
    """
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        pending_futures = []
        for group in parse_gam_file_groups(filename):
            future = executor.submit(process_group, group)
            pending_futures.append(future)

            # If too many futures are queued, wait for at least one to finish
            if len(pending_futures) >= max_pending:
                done, not_done = concurrent.futures.wait(
                    pending_futures,
                    return_when=concurrent.futures.FIRST_COMPLETED
                )
                for completed in done:
                    yield completed.result()
                pending_futures = list(not_done)

        # Wait for any remaining groups
        for future in concurrent.futures.as_completed(pending_futures):
            yield future.result()


def is_on_chromosome(alignment, chrom_name):
    """
    Return True if the alignment has any refpos with name == chrom_name.
    """
    for rp in alignment.refpos:
        if rp.name == chrom_name:
            return True
    return False


def main():
    parser = argparse.ArgumentParser(
        description="Faster parsing of a GAM file by group using vg_pb2."
    )
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of worker threads (default: 4)")
    parser.add_argument("--max_pending", type=int, default=10,
                        help="Max number of groups in flight (default: 10)")
    parser.add_argument("--chromosome", default="chr22",
                        help="Chromosome name to filter on (default: chr22)")
    parser.add_argument("--output_json", default="not_perfect_on_chrom.json",
                        help="Where to save the not-perfect reads on the chromosome")

    args = parser.parse_args()

    start_time = time.perf_counter()

    # Counters
    total_count = 0
    perfect_count = 0
    not_perfect_count = 0

    chunk_size = 100000
    not_perfect_buffer = []

    with open(args.output_json, "w") as out_f:
        # Stream and parse the GAM groups
        for group_number, alignments in process_groups_pipeline(
                args.filename,
                args.threads,
                args.max_pending
        ):
            for aln in alignments:
                if is_on_chromosome(aln, args.chromosome) or args.chromosome == "":
                    total_count += 1
                    if aln.identity < 1.0:
                        not_perfect_count += 1
                        # Store in a buffer
                        alignment_dict = MessageToDict(aln)
                        not_perfect_buffer.append(alignment_dict)
                        # If the buffer is large, flush to disk
                        if len(not_perfect_buffer) >= chunk_size:
                            for item in not_perfect_buffer:
                                json.dump(item, out_f)
                                out_f.write("\n")
                            not_perfect_buffer.clear()

                            print(f"Total reads processed: {total_count}")
                            print(f"  Perfect reads: {perfect_count} "
                                  f"({(perfect_count / total_count * 100) if total_count else 0:.2f}% of total)")
                            print(f"  Not-perfect reads: {not_perfect_count} "
                                  f"({(not_perfect_count / total_count * 100) if total_count else 0:.2f}% of total)\n")
                            # print(f"Elapsed time: {elapsed:.2f} seconds.")
                    else:
                        perfect_count += 1

        # After all alignments, flush any leftover
        for item in not_perfect_buffer:
            json.dump(item, out_f)
            out_f.write("\n")
        not_perfect_buffer.clear()

    end_time = time.perf_counter()
    elapsed = end_time - start_time

    # Final summary
    print(f"Total reads processed: {total_count}")
    print(f"  Perfect reads: {perfect_count} "
          f"({(perfect_count / total_count * 100) if total_count else 0:.2f}% of total)")
    print(f"  Not-perfect reads: {not_perfect_count} "
          f"({(not_perfect_count / total_count * 100) if total_count else 0:.2f}% of total)")
    print(f"Elapsed time: {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()
