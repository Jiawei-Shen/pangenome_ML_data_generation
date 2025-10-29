#!/usr/bin/env python3
import argparse
import gzip
import pickle
import struct
import time
import gc
import os
from collections import defaultdict

import vg_pb2
import fast_writer  # Segment, BlockTable, flush_entire_buffer_parallel_dict

# ─────────────────────────────────────────────────────────────────────────────
# File layout
GLOBAL_MAGIC = b"MYFMT\x01"
GLOBAL_VER_PACK = struct.Struct("<BBI16s")
GLOBAL_MAJOR, GLOBAL_MINOR = 0, 5
GLOBAL_HEADER_SIZE = len(GLOBAL_MAGIC) + GLOBAL_VER_PACK.size
BLOCK_HDR_PACK = struct.Struct("<I I H I I")
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size

def make_record_struct(max_read_len: int, max_cigar_len: int) -> struct.Struct:
    return struct.Struct(f"<h{max_read_len}s{max_read_len}s{max_cigar_len}shc")

def record_size(max_read_len: int, max_cigar_len: int) -> int:
    return make_record_struct(max_read_len, max_cigar_len).size

# ─────────────────────────────────────────────────────────────────────────────
# GAM parsing helpers
def read_varint(stream):
    value, shift_amount = 0, 0
    while True:
        b = stream.read(1)
        if not b: raise EOFError("EOF while reading varint")
        v = b[0]
        value |= (v & 0x7F) << shift_amount
        if not (v & 0x80): return value
        shift_amount += 7

def file_is_gzip(path):
    with open(path, "rb") as f: return f.read(2) == b"\x1f\x8b"

def gam_record_iter(path, tag="GAM"):
    open_func = gzip.open if file_is_gzip(path) else open
    with open_func(path, "rb") as f:
        while True:
            try:
                group_count = read_varint(f)
            except EOFError:
                break
            if group_count == 0: continue
            try:
                tag_len = read_varint(f)
                group_tag = f.read(tag_len).decode()
            except (EOFError, UnicodeDecodeError):
                break
            if group_tag != tag:
                for _ in range(group_count - 1): f.seek(read_varint(f), 1)
                continue
            for _ in range(group_count - 1):
                try:
                    yield f.read(read_varint(f))
                except EOFError:
                    break

def process_alignment(raw_message, wanted_nodes, chrom_filter):
    """Return dict[nid] -> list[fast_writer.Segment]."""
    segment_dict = {}
    aln = vg_pb2.Alignment()
    aln.ParseFromString(raw_message)

    if aln.mapping_quality <= 10 or (chrom_filter and not any(pos.name == chrom_filter for pos in aln.refpos)):
        return segment_dict

    read_sequence, read_quality, mapq, read_offset = aln.sequence, aln.quality, aln.mapping_quality, 0

    for mapping in aln.path.mapping:
        nid = mapping.position.node_id
        if nid not in wanted_nodes:
            for e in mapping.edit: read_offset += e.to_length
            continue

        node_offset = mapping.position.offset
        strand_byte = (b"-" if mapping.position.is_reverse else b"+")[0]

        seq_parts, bq_parts, cigar_parts = [], bytearray(), []
        for e in mapping.edit:
            fL, tL = e.from_length, e.to_length
            if fL == tL:
                cigar_parts.append(f"{fL}M" if not e.sequence else f"{fL}X")
            elif fL > 0 and tL == 0:
                cigar_parts.append(f"{fL}D")
            elif fL == 0 and tL > 0:
                cigar_parts.append(f"{tL}I")
            if tL > 0:
                seq_parts.append(read_sequence[read_offset: read_offset + tL].upper())
                bq_parts.extend(read_quality[read_offset: read_offset + tL])
                read_offset += tL

        seg = fast_writer.Segment(
            offset=int(node_offset),
            seq="".join(seq_parts).encode(),
            bq=bytes(bq_parts),
            cigar="".join(cigar_parts).encode(),
            rq=int(mapq),
            strand=int(strand_byte),
        )
        segment_dict.setdefault(nid, []).append(seg)

    return segment_dict

# ─────────────────────────────────────────────────────────────────────────────
# Output init
def initialize_output_files(stats_path, output_prefix, alt_threshold):
    with open(stats_path, "rb") as fh:
        stats_data = pickle.load(fh)

    wanted_nodes, node_counts, maxima = set(), {}, {}
    for node_id_key, stat in stats_data.items():
        nid = int(node_id_key)
        perfect = int(stat.get("perfect", 0))
        not_perfect = int(stat.get("not_perfect", 0))
        if (perfect + not_perfect) > 0 and not_perfect >= 1 and not_perfect / (perfect + not_perfect) > alt_threshold:
            wanted_nodes.add(nid)
            node_counts[nid] = perfect + not_perfect
            R = int(stat.get("max_read_length", 1) or 1)
            C = int(stat.get("max_cigar_length", 1) or 1)
            maxima[nid] = (max(1, R), max(1, C))

    print(f"Filtered {len(wanted_nodes)} nodes from {len(stats_data)} total.")
    del stats_data
    gc.collect()

    block_infos = {}
    current_offset = GLOBAL_HEADER_SIZE
    for nid in sorted(list(wanted_nodes)):  # deterministic layout
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

    dat_path = output_prefix + ".dat"
    with open(dat_path, "wb") as f:
        # 1) write global header
        f.write(GLOBAL_MAGIC)
        f.write(GLOBAL_VER_PACK.pack(GLOBAL_MAJOR, GLOBAL_MINOR, len(block_infos), b'\x00' * 16))

        # 2) true pre-allocation of the entire file
        final_size = current_offset
        try:
            os.posix_fallocate(f.fileno(), 0, final_size)  # best: real allocation, not sparse
        except AttributeError:
            os.ftruncate(f.fileno(), final_size)           # fallback if posix_fallocate missing
        except OSError:
            os.ftruncate(f.fileno(), final_size)           # fallback if FS doesn't support fallocate

        # 3) write each block header at its declared offset
        for nid, info in block_infos.items():
            f.seek(info["offset"], os.SEEK_SET)
            f.write(BLOCK_HDR_PACK.pack(
                nid,
                info["n_records"],
                0,  # flags
                info["max_read_len"],
                info["max_cigar_len"],
            ))

    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx:
        idx.write(struct.pack("<I", len(block_infos)))
        for nid, info in block_infos.items():
            idx.write(struct.pack(
                "<I Q I I H I I",
                nid,
                info["offset"],
                info["block_size"],
                info["n_records"],
                0,
                info["max_read_len"],
                info["max_cigar_len"],
            ))

    return block_infos, dat_path, wanted_nodes

# ─────────────────────────────────────────────────────────────────────────────
# Windowed flush helper: bound each C++ call to a max address span
def _flush_in_windows(dat_fd, segbuf, state, max_span_bytes, num_threads):
    """
    segbuf: dict[nid] -> list[Segment]
    state:  fast_writer.BlockTable (current_pos updated by C++ after each flush)
    """
    if not segbuf:
        return

    # pull current block infos (includes current_pos)
    infos = state.to_py_dict()

    # Build (nid, start, end, segs) for this batch
    entries = []
    for nid, segs in segbuf.items():
        if not segs:
            continue
        bi = infos.get(nid)
        if bi is None:
            continue
        start = bi["offset"] + BLOCK_HDR_SIZE + bi["current_pos"] * bi["record_size"]
        end   = start + len(segs) * bi["record_size"]
        entries.append((nid, start, end, segs))

    if not entries:
        return

    entries.sort(key=lambda x: x[1])

    def _flush_window(win_entries):
        if not win_entries:
            return
        chunk = {nid: segs for (nid, _s, _e, segs) in win_entries}
        fast_writer.flush_entire_buffer_parallel_dict(
            dat_fd,
            chunk,
            state,
            BLOCK_HDR_SIZE,
            num_threads=num_threads,
            sort_by_offset=True,
        )

    # Greedy pack into windows by contiguous address span
    window = []
    win_start = entries[0][1]
    win_end   = entries[0][2]
    for item in entries:
        nid, start, end, segs = item
        # If adding this entry would exceed the span budget, flush current window
        if end - win_start > max_span_bytes and window:
            _flush_window(window)
            window.clear()
            win_start = start
            win_end   = end
        else:
            win_end = max(win_end, end)
        window.append(item)

    if window:
        _flush_window(window)

# ─────────────────────────────────────────────────────────────────────────────
# Pipeline
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter, use_existing, num_threads, buffer_segments, alt_threshold, max_span_bytes):
    if use_existing:
        raise NotImplementedError("--use-existing is not implemented yet.")
    block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix, alt_threshold)

    # Persistent C++ state (tracks current_pos)
    state = fast_writer.BlockTable(block_infos)

    BUFFER_SEGMENTS = buffer_segments

    next_milestone, total_reads, total_segments = milestone_step, 0, 0
    start_time = time.perf_counter()

    dat_fh = open(dat_path, "r+b", buffering=0)
    dat_fd = dat_fh.fileno()

    segment_buffer = defaultdict(list)

    for raw_msg in gam_record_iter(gam_path):
        segs_by_node = process_alignment(raw_msg, wanted_nodes, chrom_filter)
        total_reads += 1

        for nid, segs in segs_by_node.items():
            segment_buffer[nid].extend(segs)
            total_segments += len(segs)

        if total_segments >= BUFFER_SEGMENTS:
            print(f"Buffer full. Flushing {len(segment_buffer)} node blocks in windows (threads={num_threads})...")
            t0 = time.perf_counter()
            _flush_in_windows(dat_fd, segment_buffer, state, max_span_bytes, num_threads)
            segment_buffer.clear()
            total_segments = 0
            print(f"Flush complete in {time.perf_counter() - t0:.2f} s.")

        if total_reads >= next_milestone:
            elapsed = time.perf_counter() - start_time
            print(f"{total_reads:,} reads processed | {elapsed:.1f} s")
            next_milestone += milestone_step

    if segment_buffer:
        print(f"Final flush of {len(segment_buffer)} node blocks in windows (threads={num_threads})...")
        t0 = time.perf_counter()
        _flush_in_windows(dat_fd, segment_buffer, state, max_span_bytes, num_threads)
        segment_buffer.clear()
        print(f"Final flush done in {time.perf_counter() - t0:.2f} s.")

    dat_fh.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads:,}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} s")

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="GAM segment extractor with C++ pack+write (parallel mmap with range-msync).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("gam_path")
    parser.add_argument("stats_pickle")
    parser.add_argument("output_prefix")
    parser.add_argument("--buffer", type=int, default=50_000_000,
                        help="Flush when accumulated segments reach this count.")
    parser.add_argument("--alt", type=float, default=0.05, help="ALT fraction threshold (0–1)")
    parser.add_argument("--milestone", type=int, default=1_000_000)
    parser.add_argument("--chr", default="")
    parser.add_argument("--use-existing", action="store_true")
    parser.add_argument("--threads", type=int, default=4, help="C++ writer threads")
    parser.add_argument("--max-span", type=str, default="2G",
                        help="Max file-offset span per flush window, e.g., 1G, 2G, 512M.")
    args = parser.parse_args()

    # parse --max-span like "2G" / "512M"
    unit = args.max_span[-1].upper() if args.max_span[-1].isalpha() else ""
    val = int(args.max_span[:-1]) if unit else int(args.max_span)
    mult = {"K": 1024, "M": 1024**2, "G": 1024**3}.get(unit, 1)
    max_span_bytes = val * mult

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr,
        use_existing=args.use_existing,
        num_threads=args.threads,
        buffer_segments=args.buffer,
        alt_threshold=args.alt,
        max_span_bytes=max_span_bytes,
    )

if __name__ == "__main__":
    main()
