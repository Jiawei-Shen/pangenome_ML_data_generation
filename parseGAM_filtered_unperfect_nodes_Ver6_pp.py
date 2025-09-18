#!/usr/bin/env python3
# Option B: parent packs buffers; workers only do os.pwrite at fixed offsets.

import argparse
import gzip
import pickle
import struct
import time
import gc
import os
from collections import defaultdict
import multiprocessing as mp
import vg_pb2  # compiled protobuf for VG Alignment

# ─────────────────────────────────────────────────────────────────────────────
# Segment container
class Segment:
    __slots__ = ("offset", "seq", "bq", "cigar", "rq", "strand")

    def __init__(self, offset, seq, bq, cigar, rq, strand):
        self.offset = offset     # int
        self.seq = seq           # bytes (unpadded)
        self.bq = bq             # bytes (unpadded)
        self.cigar = cigar       # bytes (unpadded)
        self.rq = rq             # int (MAPQ)
        self.strand = strand     # b'+' or b'-'


# ─────────────────────────────────────────────────────────────────────────────
# File layout
GLOBAL_MAGIC = b"MYFMT\x01"
GLOBAL_VER_PACK = struct.Struct("<BBI16s")  # major, minor, nblocks, reserved16
GLOBAL_MAJOR, GLOBAL_MINOR = 0, 5
GLOBAL_HEADER_SIZE = len(GLOBAL_MAGIC) + GLOBAL_VER_PACK.size
BLOCK_HDR_PACK = struct.Struct("<I I H I I")  # node_id, n_records, flags, maxR, maxC
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size


def make_record_struct(max_read_len: int, max_cigar_len: int) -> struct.Struct:
    # <h seq[R] bq[R] cigar[C] h c>
    return struct.Struct(f"<h{max_read_len}s{max_read_len}s{max_cigar_len}shc")


def record_size(max_read_len: int, max_cigar_len: int) -> int:
    return make_record_struct(max_read_len, max_cigar_len).size


# ─────────────────────────────────────────────────────────────────────────────
# GAM parsing
def read_varint(stream):
    value, shift_amount = 0, 0
    while True:
        b = stream.read(1)
        if not b:
            raise EOFError("EOF while reading varint")
        v = b[0]
        value |= (v & 0x7F) << shift_amount
        if not (v & 0x80):
            return value
        shift_amount += 7


def file_is_gzip(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def gam_record_iter(path, tag="GAM"):
    open_func = gzip.open if file_is_gzip(path) else open
    with open_func(path, "rb") as f:
        while True:
            try:
                group_count = read_varint(f)
            except EOFError:
                break
            if group_count == 0:
                continue
            try:
                tag_len = read_varint(f)
                group_tag = f.read(tag_len).decode()
            except (EOFError, UnicodeDecodeError):
                break
            if group_tag != tag:
                for _ in range(group_count - 1):
                    try:
                        f.seek(read_varint(f), 1)
                    except EOFError:
                        break
                continue
            for _ in range(group_count - 1):
                try:
                    size = read_varint(f)
                    yield f.read(size)
                except EOFError:
                    break


def process_alignment(raw_message, wanted_nodes, chrom_filter):
    """
    Returns: dict[node_id] -> [Segment, ...] for segments belonging to wanted_nodes
    """
    out = {}
    aln = vg_pb2.Alignment()
    aln.ParseFromString(raw_message)

    if aln.mapping_quality <= 10:
        return out
    if chrom_filter and not any(getattr(rp, "name", "") == chrom_filter for rp in aln.refpos):
        return out

    seq = aln.sequence
    bqs = aln.quality
    mapq = aln.mapping_quality
    read_off = 0

    for mapping in aln.path.mapping:
        nid = mapping.position.node_id
        if nid not in wanted_nodes:
            # Skip over edits to keep read_off correct
            for e in mapping.edit:
                read_off += e.to_length
            continue

        node_off = mapping.position.offset
        strand = b"-" if mapping.position.is_reverse else b"+"
        seq_parts = []
        bq_parts = bytearray()
        cigar_parts = []

        for e in mapping.edit:
            fL, tL = e.from_length, e.to_length
            if fL == tL:
                cigar_parts.append(f"{fL}M" if not e.sequence else f"{fL}X")
            elif fL > 0 and tL == 0:
                cigar_parts.append(f"{fL}D")
            elif fL == 0 and tL > 0:
                cigar_parts.append(f"{tL}I")

            if tL > 0:
                seq_parts.append(seq[read_off : read_off + tL].upper())
                bq_parts.extend(bqs[read_off : read_off + tL])
                read_off += tL

        seg = Segment(
            offset=node_off,
            seq=("".join(seq_parts)).encode(),
            bq=bytes(bq_parts),
            cigar=("".join(cigar_parts)).encode(),
            rq=mapq,
            strand=strand,
        )
        out.setdefault(nid, []).append(seg)

    return out


# ─────────────────────────────────────────────────────────────────────────────
# Output initialization
def initialize_output_files(stats_path, output_prefix):
    with open(stats_path, "rb") as fh:
        stats = pickle.load(fh)

    wanted_nodes = set()
    node_counts = {}
    maxima = {}

    for node_id_key, stat in stats.items():
        nid = int(node_id_key)
        perfect = int(stat.get("perfect", 0))
        not_perfect = int(stat.get("not_perfect", 0))
        total = perfect + not_perfect
        if total > 0 and not_perfect > 1 and not_perfect / total > 0.05:
            wanted_nodes.add(nid)
            node_counts[nid] = total
            R = int(stat.get("max_read_length", 1) or 1)
            C = int(stat.get("max_cigar_length", 1) or 1)
            maxima[nid] = (max(1, R), max(1, C))

    print(f"Filtered {len(wanted_nodes)} nodes from {len(stats)} total.")
    del stats
    gc.collect()

    block_infos = {}
    current_offset = GLOBAL_HEADER_SIZE
    for nid in sorted(wanted_nodes):  # deterministic block order
        nrec = node_counts[nid]
        R, C = maxima[nid]
        rec_sz = record_size(R, C)
        blk_sz = BLOCK_HDR_SIZE + nrec * rec_sz
        block_infos[nid] = {
            "offset": current_offset,
            "n_records": nrec,
            "current_pos": 0,
            "max_read_len": R,
            "max_cigar_len": C,
            "record_size": rec_sz,
            "block_size": blk_sz,
        }
        current_offset += blk_sz

    # Create and preallocate .dat
    dat_path = output_prefix + ".dat"
    with open(dat_path, "wb") as f:
        f.write(GLOBAL_MAGIC)
        f.write(GLOBAL_VER_PACK.pack(GLOBAL_MAJOR, GLOBAL_MINOR, len(block_infos), b"\x00" * 16))
        for nid, info in block_infos.items():
            f.write(BLOCK_HDR_PACK.pack(nid, info["n_records"], 0, info["max_read_len"], info["max_cigar_len"]))
        # Pre-allocate the payload region
        if current_offset > GLOBAL_HEADER_SIZE:
            try:
                # Best: fallocate (may not exist on all FS)
                os.posix_fallocate(f.fileno(), 0, current_offset)
            except Exception:
                # Fallback: seek and write a single zero at end-1
                f.seek(current_offset - 1)
                f.write(b"\x00")

    # Write .idx
    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx:
        idx.write(struct.pack("<I", len(block_infos)))
        for nid, info in block_infos.items():
            idx.write(
                struct.pack(
                    "<I Q I I H I I",
                    nid,
                    info["offset"],
                    info["block_size"],
                    info["n_records"],
                    0,
                    info["max_read_len"],
                    info["max_cigar_len"],
                )
            )

    return block_infos, dat_path, wanted_nodes


def load_existing_output_files(output_prefix):
    raise NotImplementedError("--use-existing loader not implemented in this example.")


# ─────────────────────────────────────────────────────────────────────────────
# Worker: only pwrite
g_fd = -1  # per-process global FD


def init_flush_worker(fd):
    global g_fd
    g_fd = fd


def flush_worker(job):
    """job: (write_pos:int, buf:bytes)"""
    write_pos, buf = job
    os.pwrite(g_fd, buf, write_pos)


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter, use_existing, num_flush_workers):
    if use_existing:
        block_infos, dat_path, wanted_nodes = load_existing_output_files(output_prefix)
    else:
        block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)

    BUFFER_SEGMENTS = 2_000_000  # total segments batched before a flush
    next_milestone = milestone_step
    total_reads = 0
    total_segments = 0
    start_time = time.perf_counter()

    # Open once; we won't use high-level writes further, only the FD
    dat_fh = open(dat_path, "r+b")
    dat_fd = dat_fh.fileno()
    try:
        os.set_inheritable(dat_fd, True)
    except Exception:
        pass

    # Prefer "fork" on Linux; fall back if unavailable
    try:
        mpc = mp.get_context("fork")
    except ValueError:
        mpc = mp.get_context()

    pool = mpc.Pool(
        processes=num_flush_workers,
        initializer=init_flush_worker,
        initargs=(dat_fd,),
    )

    segment_buffer = defaultdict(list)

    def flush_segment_buffer_parallel():
        nonlocal total_segments
        if not segment_buffer:
            return
        if not hasattr(os, "pwrite"):
            raise NotImplementedError("os.pwrite() not available on this OS.")

        items = list(segment_buffer.items())
        segment_buffer.clear()

        jobs = []
        for nid, segs in items:
            if not segs:
                continue
            info = block_infos.get(nid)
            if not info:
                continue

            # Overflow guard
            if info["current_pos"] + len(segs) > info["n_records"]:
                raise RuntimeError(
                    f"Block overflow for node {nid}: "
                    f"{info['current_pos']} + {len(segs)} > {info['n_records']}"
                )

            base_offset = info["offset"] + BLOCK_HDR_SIZE
            R, C = info["max_read_len"], info["max_cigar_len"]
            rec_pack = make_record_struct(R, C)
            rec_sz = rec_pack.size

            # Parent packs buffer
            n = len(segs)
            buf = bytearray(rec_sz * n)
            off = 0
            for s in segs:
                rec_pack.pack_into(
                    buf,
                    off,
                    int(s.offset),
                    s.seq.ljust(R, b"\x00")[:R],
                    s.bq.ljust(R, b"\x00")[:R],
                    s.cigar.ljust(C, b"\x00")[:C],
                    int(s.rq),
                    s.strand if s.strand in (b"+", b"-") else b"+",
                )
                off += rec_sz

            write_pos = base_offset + info["current_pos"] * rec_sz
            info["current_pos"] += n

            jobs.append((write_pos, bytes(buf)))  # send immutable bytes

        if jobs:
            # Larger chunksize reduces scheduling overhead
            pool.map(flush_worker, jobs, chunksize=32)
        total_segments = 0

    try:
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
                print(f"{total_reads:,} reads processed | {elapsed:.1f} s")
                next_milestone += milestone_step

        # final drain
        flush_segment_buffer_parallel()
    finally:
        pool.close()
        pool.join()
        dat_fh.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads:,}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
def main():
    p = argparse.ArgumentParser(
        description="GAM segment extractor with parallel pwrite (parent packs).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("gam_path", help="Path to the GAM file")
    p.add_argument("stats_pickle", help="Path to the stats pickle with node maxima")
    p.add_argument("output_prefix", help="Prefix for output files (.dat, .idx)")
    p.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval")
    p.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    p.add_argument("--use-existing", action="store_true", help="Reuse existing initialized output files")
    p.add_argument("--flush-workers", type=int, default=4, help="Workers for I/O (pwrite)")
    args = p.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr,
        use_existing=args.use_existing,
        num_flush_workers=args.flush_workers,
    )


if __name__ == "__main__":
    main()
