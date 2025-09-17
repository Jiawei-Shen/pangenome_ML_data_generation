#!/usr/bin/env python3
import argparse
import gzip
import pickle
import struct
import time
import gc
import os
from collections import defaultdict
from multiprocessing.pool import ThreadPool
import vg_pb2  # Assumes you have the compiled protobuf Python file


# ... (Code from Segment class to worker functions remains the same) ...
# ─────────────────────────────────────────────────────────────────────────────
# Segment container
class Segment:
    __slots__ = ('offset', 'seq', 'bq', 'cigar', 'rq', 'strand')

    def __init__(self, offset, seq, bq, cigar, rq, strand):
        self.offset = offset
        self.seq = seq
        self.bq = bq
        self.cigar = cigar
        self.rq = rq
        self.strand = strand


# ... (rest of unchanged code) ...

def file_is_gzip(path):
    with open(path, "rb") as f: return f.read(2) == b"\x1f\x8b"


# ... (rest of unchanged code) ...

# ─────────────────────────────────────────────────────────────────────────────
# PARALLEL FLUSH WORKER FUNCTIONS

g_fd = -1


def init_flush_worker(fd):
    global g_fd
    g_fd = fd


def flush_worker(job):
    buf, write_pos = job
    os.pwrite(g_fd, buf, write_pos)


# ─────────────────────────────────────────────────────────────────────────────
# Main Pipeline (MODIFIED for Batched Processing)
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter, use_existing, num_flush_workers):
    # ... (initialization code is the same) ...
    if use_existing:
        block_infos, dat_path, wanted_nodes = load_existing_output_files(output_prefix)
    else:
        block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)

    # BUFFER_SEGMENTS = 400_000_000
    BUFFER_SEGMENTS = 1_000_000
    next_milestone, total_reads, total_segments = milestone_step, 0, 0
    start_time = time.perf_counter()

    dat_fh = open(dat_path, "r+b")
    segment_buffer = defaultdict(list)

    def flush_segment_buffer_parallel():
        nonlocal total_segments
        if not segment_buffer:
            return

        if not hasattr(os, 'pwrite'):
            raise NotImplementedError("Parallel flushing requires os.pwrite(), not available on this OS.")

        print(f"Preparing {len(segment_buffer)} node blocks for batched parallel flush...")
        flush_start_time = time.perf_counter()
        dat_fd = dat_fh.fileno()

        items_to_process = list(segment_buffer.items())
        segment_buffer.clear()

        # Helper to yield batches of items
        def chunker(seq, size):
            return (seq[pos:pos + size] for pos in range(0, len(seq), size))

        # Tune this based on your system. A few hundred to a few thousand is a good start.
        BATCH_SIZE = 1024

        # Create the pool once for the entire flush operation
        with ThreadPool(processes=num_flush_workers, initializer=init_flush_worker, initargs=(dat_fd,)) as pool:
            for item_batch in chunker(items_to_process, BATCH_SIZE):
                jobs_batch = []
                # SERIAL: Prepare one small batch of jobs
                for nid, segs in item_batch:
                    info = block_infos.get(nid)
                    if not info: continue

                    base_offset = info["offset"] + BLOCK_HDR_SIZE
                    write_pos = base_offset + info["current_pos"] * info["record_size"]
                    info["current_pos"] += len(segs)

                    R, C = info["max_read_len"], info["max_cigar_len"]
                    rec_pack = make_record_struct(R, C)
                    buf = bytearray(rec_pack.size * len(segs))
                    off = 0
                    for s in segs:
                        rec_pack.pack_into(
                            buf, off, int(s.offset), s.seq.ljust(R, b'\x00')[:R],
                            s.bq.ljust(R, b'\x00')[:R], s.cigar.ljust(C, b'\x00')[:C],
                            int(s.rq), s.strand if s.strand in (b'+', b'-') else b'+'
                        )
                        off += rec_pack.size
                    jobs_batch.append((buf, write_pos))

                # DISPATCH: Send the current batch to the workers
                # The workers will process this batch while the main thread
                # prepares the next one.
                pool.map(flush_worker, jobs_batch)

        total_segments = 0
        print(
            f"Parallel flush of {len(items_to_process)} blocks complete in {time.perf_counter() - flush_start_time:.2f} seconds.")

    # ... (The rest of the main loop and main function are unchanged) ...
    for raw_msg in gam_record_iter(gam_path):
        segs_by_node = process_alignment(raw_msg, wanted_nodes, chrom_filter)
        total_reads += 1

        for nid, segs in segs_by_node.items():
            segment_buffer[nid].extend(segs)
            total_segments += len(segs)

        if total_segments >= BUFFER_SEGMENTS:
            flush_segment_buffer_parallel()

        if total_reads >= next_milestone:
            elapsed = time.perf_counter() - start_time
            print(f"{total_reads:,} reads processed | {elapsed:.1f} seconds")
            next_milestone += milestone_step

    flush_segment_buffer_parallel()
    dat_fh.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads:,}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="GAM segment extractor with parallel buffer flushing.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # ... (arguments are the same) ...
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the stats pickle with node maxima")
    parser.add_argument("output_prefix", help="Prefix for output files (.dat, .idx)")
    parser.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    parser.add_argument("--use-existing", action="store_true", help="Reuse existing initialized output files")
    parser.add_argument("--flush-workers", type=int, default=4,
                        help="Number of workers for parallel file writing (I/O-bound)")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr,
        use_existing=args.use_existing,
        num_flush_workers=args.flush_workers
    )