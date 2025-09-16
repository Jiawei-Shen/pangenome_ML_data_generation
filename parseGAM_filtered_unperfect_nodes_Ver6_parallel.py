#!/usr/bin/env python3
import argparse
import gzip
import pickle
import struct
import time
import gc
import mmap
import os
from collections import defaultdict
import vg_pb2
import multiprocessing
from multiprocessing import Process, Queue, Value, Lock
from ctypes import c_ulonglong


# ─────────────────────────────────────────────────────────────────────────────
# Data structures (Segment, Struct definitions) - UNCHANGED
# ─────────────────────────────────────────────────────────────────────────────
class Segment:
    __slots__ = ('offset', 'seq', 'bq', 'cigar', 'rq', 'strand')

    def __init__(self, offset, seq, bq, cigar, rq, strand):
        self.offset, self.seq, self.bq, self.cigar, self.rq, self.strand = offset, seq, bq, cigar, rq, strand


GLOBAL_MAGIC = b"MYFMT\x01"
GLOBAL_VER_PACK = struct.Struct("<BBI16s")
GLOBAL_MAJOR_EXPECTED, GLOBAL_MINOR_EXPECTED = 0, 5
GLOBAL_HEADER_SIZE = len(GLOBAL_MAGIC) + GLOBAL_VER_PACK.size
BLOCK_HDR_PACK = struct.Struct("<I I H I I")
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size


def make_record_struct(max_read_len: int, max_cigar_len: int) -> struct.Struct:
    if max_read_len <= 0 or max_cigar_len <= 0:
        raise ValueError("max_read_len and max_cigar_len must be > 0")
    return struct.Struct(f"<h{max_read_len}s{max_read_len}s{max_cigar_len}shc")


def record_size(max_read_len: int, max_cigar_len: int) -> int:
    return make_record_struct(max_read_len, max_cigar_len).size


# ─────────────────────────────────────────────────────────────────────────────
# GAM Parsing and CIGAR builder - UNCHANGED
# ─────────────────────────────────────────────────────────────────────────────
def read_varint(stream):
    value, shift_amount = 0, 0
    while True:
        byte_pair = stream.read(1)
        if not byte_pair: raise EOFError("EOF while reading varint")
        byte_value = byte_pair[0]
        value |= (byte_value & 0x7F) << shift_amount
        if not (byte_value & 0x80): return value
        shift_amount += 7


def file_is_gzip(path):
    with open(path, "rb") as f: return f.read(2) == b"\x1f\x8b"


def process_alignment(raw_message, wanted_nodes, chrom_filter):
    segment_dict = defaultdict(list)
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw_message)
    if alignment.mapping_quality <= 10: return segment_dict
    if chrom_filter and not any(pos.name == chrom_filter for pos in alignment.refpos): return segment_dict
    read_sequence, read_quality, mapping_quality, read_offset = alignment.sequence, alignment.quality, alignment.mapping_quality, 0
    for mapping in alignment.path.mapping:
        node_id = mapping.position.node_id
        current_read_offset_for_node = 0
        for edit in mapping.edit: current_read_offset_for_node += edit.to_length
        if node_id not in wanted_nodes:
            read_offset += current_read_offset_for_node
            continue
        node_offset, strand_char = mapping.position.offset, b"-" if mapping.position.is_reverse else b"+"
        sequence_parts, quality_parts, cigar_parts = [], bytearray(), []
        for edit in mapping.edit:
            from_len, to_len, edit_len = edit.from_length, edit.to_length, len(edit.sequence)
            if from_len == to_len:
                cigar_parts.append(f"{from_len}{'X' if edit_len else 'M'}")
            elif from_len > 0 and to_len == 0:
                cigar_parts.append(f"{from_len}D")
            elif from_len == 0 and to_len > 0:
                cigar_parts.append(f"{to_len}I")
            else:
                raise ValueError(f"Unexpected edit: from={from_len}, to={to_len}, seq={edit_len}")
            sequence_parts.append(read_sequence[read_offset: read_offset + to_len].upper())
            quality_parts.extend(read_quality[read_offset: read_offset + to_len])
            read_offset += to_len
        seg = Segment(node_offset, "".join(sequence_parts).encode(), bytes(quality_parts),
                      "".join(cigar_parts).encode(), mapping_quality, strand_char)
        segment_dict[node_id].append(seg)
    return segment_dict


# ─────────────────────────────────────────────────────────────────────────────
# File initialization and validation - UNCHANGED
# ─────────────────────────────────────────────────────────────────────────────
def initialize_output_files(stats_path, output_prefix):
    with open(stats_path, "rb") as stats_file:
        stats_data = pickle.load(stats_file)
    block_infos, wanted_nodes, current_offset, total_nodes = {}, set(), GLOBAL_HEADER_SIZE, 0
    for node_id_key, stat in stats_data.items():
        total_nodes += 1;
        node_id = int(node_id_key);
        perfect, not_perfect = int(stat.get("perfect", 0)), int(stat.get("not_perfect", 0))
        if (perfect + not_perfect) > 0 and not_perfect > 1 and not_perfect / (perfect + not_perfect) > 0.05:
            wanted_nodes.add(node_id);
            n_records = perfect + not_perfect
            R, C = int(stat.get("max_read_length", 1) or 1), int(stat.get("max_cigar_length", 1) or 1)
            if R <= 0 or C <= 0: raise RuntimeError(f"Invalid maxima for node {node_id}: R={R}, C={C}")
            rec_sz = record_size(R, C);
            blk_sz = BLOCK_HDR_SIZE + n_records * rec_sz
            block_infos[node_id] = {"offset": current_offset, "n_records": n_records, "max_read_len": R,
                                    "max_cigar_len": C, "record_size": rec_sz, "block_size": blk_sz}
            current_offset += blk_sz
    print(f"Filtered {len(wanted_nodes)}/{total_nodes} nodes ({len(wanted_nodes) / max(total_nodes, 1):.2%}).")
    dat_path = output_prefix + ".dat"
    with open(dat_path, "wb") as f:
        f.write(GLOBAL_MAGIC);
        f.write(GLOBAL_VER_PACK.pack(GLOBAL_MAJOR_EXPECTED, GLOBAL_MINOR_EXPECTED, len(block_infos), b'\x00' * 16))
        for node_id, info in block_infos.items():
            f.write(BLOCK_HDR_PACK.pack(node_id, info["n_records"], 0, info["max_read_len"], info["max_cigar_len"]))
            f.seek(info["n_records"] * info["record_size"], os.SEEK_CUR)  # Pre-allocate space
    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx_file:
        idx_file.write(struct.pack("<I", len(block_infos)))
        for node_id, info in block_infos.items():
            idx_file.write(
                struct.pack("<I Q I I H I I", node_id, info["offset"], info["block_size"], info["n_records"], 0,
                            info["max_read_len"], info["max_cigar_len"]))
    return block_infos, dat_path, wanted_nodes


def _validate_latest_dat_header(dat_path):
    with open(dat_path, "rb") as df:
        if df.read(len(GLOBAL_MAGIC)) != GLOBAL_MAGIC: raise RuntimeError("Invalid .dat magic")
        major, minor, _, _ = GLOBAL_VER_PACK.unpack(df.read(GLOBAL_VER_PACK.size))
        if (major, minor) != (GLOBAL_MAJOR_EXPECTED, GLOBAL_MINOR_EXPECTED): raise RuntimeError(
            f"Unsupported .dat version {major}.{minor}")


def _validate_latest_idx_header_and_size(idx_path):
    with open(idx_path, "rb") as f:
        data = f.read()
    if len(data) < 4: raise RuntimeError("Corrupt .idx: too small")
    (count,) = struct.unpack("<I", data[:4])
    if len(data) != 4 + count * 30: raise RuntimeError(f"Unsupported .idx size")


def load_existing_output_files(output_prefix):
    idx_path, dat_path = f"{output_prefix}.idx", f"{output_prefix}.dat"
    if not (os.path.exists(idx_path) and os.path.exists(dat_path)): raise FileNotFoundError(
        f"{idx_path} and {dat_path}")
    _validate_latest_idx_header_and_size(idx_path);
    _validate_latest_dat_header(dat_path)
    entries = []
    with open(idx_path, "rb") as f:
        (count,) = struct.unpack("<I", f.read(4))
        for _ in range(count): entries.append(struct.unpack("<I Q I I H I I", f.read(30)))
    block_infos, wanted_nodes = {}, set()
    with open(dat_path, "rb") as df:
        _, _, dat_count, _ = GLOBAL_VER_PACK.unpack(df.read(GLOBAL_HEADER_SIZE)[len(GLOBAL_MAGIC):])
        if dat_count != len(entries): raise RuntimeError(f".dat/.idx count mismatch")
        for node_id, offset, block_size, n_records, _, R_idx, C_idx in entries:
            df.seek(offset)
            nid2, nrec2, _, R2, C2 = BLOCK_HDR_PACK.unpack(df.read(BLOCK_HDR_SIZE))
            if nid2 != node_id or nrec2 != n_records or R2 <= 0 or C2 <= 0: raise RuntimeError(
                f".dat/.idx mismatch for node {node_id}")
            block_infos[node_id] = {"offset": offset, "n_records": n_records, "max_read_len": R2, "max_cigar_len": C2,
                                    "record_size": record_size(R2, C2), "block_size": block_size}
            wanted_nodes.add(node_id)
    print(f"Reusing existing output with {len(block_infos)} node blocks.")
    return block_infos, dat_path, wanted_nodes


# ─────────────────────────────────────────────────────────────────────────────
# READER PROCESS: Reads GAM and queues raw messages for workers - UNCHANGED
# ─────────────────────────────────────────────────────────────────────────────
def reader_process(gam_path, queue, total_reads_counter, milestone_step, num_workers):
    open_func = gzip.open if file_is_gzip(gam_path) else open
    batch, batch_size = [], 1000
    start_time, next_milestone = time.perf_counter(), milestone_step

    with open_func(gam_path, "rb") as f:
        while True:
            try:
                # This simplified GAM iterator assumes all chunks are Alignments
                read_varint(f)  # group_count
                read_varint(f)  # tag_len
                f.read(3)  # "GAM" tag
                msg_size = read_varint(f)
                msg = f.read(msg_size)
                if not msg: break
                batch.append(msg)
                if len(batch) >= batch_size:
                    queue.put(batch)
                    with total_reads_counter.get_lock():
                        total_reads_counter.value += len(batch)
                    batch = []
                if total_reads_counter.value >= next_milestone:
                    elapsed = time.perf_counter() - start_time
                    print(f"Reader: {total_reads_counter.value} reads queued | Total Elapsed: {elapsed:.1f}s")
                    next_milestone += milestone_step
            except EOFError:
                break
    if batch:
        queue.put(batch)
        with total_reads_counter.get_lock(): total_reads_counter.value += len(batch)
    for _ in range(num_workers): queue.put(None)  # Sentinel values to signal workers to stop


# ─────────────────────────────────────────────────────────────────────────────
# WORKER PROCESS: Processes messages and writes to .dat file via mmap
# ─────────────────────────────────────────────────────────────────────────────
def worker_process(queue, dat_path, block_infos, wanted_nodes, chrom_filter, shared_positions, locks):
    """
    **MODIFIED**: This worker now reports its flush statistics.
    """
    dat_fh = open(dat_path, "r+b")
    fileno = dat_fh.fileno()
    pid = os.getpid()  # Get Process ID once for reporting

    for batch in iter(queue.get, None):
        local_segment_buffer = defaultdict(list)
        for raw_msg in batch:
            segment_dict = process_alignment(raw_msg, wanted_nodes, chrom_filter)
            for node_id, segs in segment_dict.items():
                local_segment_buffer[node_id].extend(segs)

        if not local_segment_buffer:
            continue

        # --- NEW: Start timing and counting before the flush ---
        flush_start_time = time.perf_counter()
        total_segments_in_batch = sum(len(segs) for segs in local_segment_buffer.values())

        # This is the main flush loop
        for node_id, segs in local_segment_buffer.items():
            if not segs: continue

            info = block_infos[node_id]
            R, C, rec_size = info["max_read_len"], info["max_cigar_len"], info["record_size"]
            rec_pack = make_record_struct(R, C)

            # Atomically reserve a block of records to write to
            num_segs = len(segs)
            with locks[node_id]:
                start_record_idx = shared_positions[node_id].value
                shared_positions[node_id].value += num_segs

            # Calculate absolute byte offset for this reserved block
            block_records_start_offset = info["offset"] + BLOCK_HDR_SIZE
            write_start_offset = block_records_start_offset + (start_record_idx * rec_size)

            # Write the data using mmap
            with mmap.mmap(fileno, length=num_segs * rec_size, access=mmap.ACCESS_WRITE,
                           offset=write_start_offset) as mm:
                current_pos_in_map = 0
                for seg in segs:
                    rec_pack.pack_into(
                        mm, current_pos_in_map,
                        int(seg.offset), seg.seq.ljust(R, b'\x00')[:R], seg.bq.ljust(R, b'\x00')[:R],
                        seg.cigar.ljust(C, b'\x00')[:C], int(seg.rq), seg.strand if seg.strand in (b'+', b'-') else b'+'
                    )
                    current_pos_in_map += rec_size

        # --- NEW: Report flush stats after writing is complete ---
        flush_duration = time.perf_counter() - flush_start_time
        print(f"  [Worker {pid}] Flushed {total_segments_in_batch} segments in {flush_duration:.4f} seconds.")

    dat_fh.close()


# ─────────────────────────────────────────────────────────────────────────────
# Main Pipeline Orchestrator - UNCHANGED
# ─────────────────────────────────────────────────────────────────────────────
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter, use_existing, num_workers):
    if use_existing:
        block_infos, dat_path, wanted_nodes = load_existing_output_files(output_prefix)
    else:
        block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)

    start_time = time.perf_counter()

    manager = multiprocessing.Manager()
    locks = manager.dict({node_id: Lock() for node_id in block_infos})
    shared_positions = manager.dict({node_id: Value(c_ulonglong, 0) for node_id in block_infos})

    work_queue = Queue(maxsize=num_workers * 4)
    total_reads_counter = Value('L', 0)

    reader = Process(target=reader_process,
                     args=(gam_path, work_queue, total_reads_counter, milestone_step, num_workers))
    reader.start()

    workers = []
    for _ in range(num_workers):
        p = Process(target=worker_process,
                    args=(work_queue, dat_path, block_infos, wanted_nodes, chrom_filter, shared_positions, locks))
        workers.append(p)
        p.start()

    reader.join()
    for p in workers:
        p.join()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads_counter.value}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point - UNCHANGED
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Parallel GAM segment extractor using mmap.")
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the node stats pickle file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    parser.add_argument("--use-existing", action="store_true", help="Reuse existing initialized output")
    parser.add_argument("--workers", type=int, default=os.cpu_count(), help="Number of worker processes")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr,
        use_existing=args.use_existing,
        num_workers=args.workers
    )


if __name__ == "__main__":
    main()