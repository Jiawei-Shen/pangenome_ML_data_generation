#!/usr/bin/env python3
import argparse
import gzip
import pickle
import struct
import time
import gc
import os
from collections import defaultdict
import vg_pb2      # Assumes you have the compiled protobuf Python file
import fast_writer # Our new C++ module

# ─────────────────────────────────────────────────────────────────────────────
# Segment container
class Segment:
    __slots__ = ('offset', 'seq', 'bq', 'cigar', 'rq', 'strand')

    def __init__(self, offset, seq, bq, cigar, rq, strand):
        self.offset = offset
        self.seq = seq      # bytes (unpadded)
        self.bq = bq        # bytes (unpadded)
        self.cigar = cigar  # bytes (unpadded)
        self.rq = rq        # int (MAPQ)
        self.strand = strand[0] # Pass the integer value of the byte

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
# GAM parsing and Alignment Processing (Serial Helpers)
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
        strand_char = b"-" if mapping.position.is_reverse else b"+"
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

        seg = Segment(offset=node_offset, seq="".join(seq_parts).encode(), bq=bytes(bq_parts),
                      cigar="".join(cigar_parts).encode(), rq=mapq, strand=strand_char)
        segment_dict.setdefault(nid, []).append(seg)
    return segment_dict

# ─────────────────────────────────────────────────────────────────────────────
# File Initialization
def initialize_output_files(stats_path, output_prefix):
    with open(stats_path, "rb") as fh:
        stats_data = pickle.load(fh)

    wanted_nodes, node_counts, maxima = set(), {}, {}
    for node_id_key, stat in stats_data.items():
        nid = int(node_id_key)
        perfect = int(stat.get("perfect", 0))
        not_perfect = int(stat.get("not_perfect", 0))
        if (perfect + not_perfect) > 0 and not_perfect > 1 and not_perfect / (perfect + not_perfect) > 0.05:
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
    for nid in sorted(list(wanted_nodes)): # Sort for deterministic layout
        nrec = node_counts[nid]
        R, C = maxima[nid]
        rec_sz = record_size(R, C)
        blk_sz = BLOCK_HDR_SIZE + nrec * rec_sz
        block_infos[nid] = {"offset": current_offset, "n_records": nrec, "current_pos": 0,
                            "max_read_len": R, "max_cigar_len": C, "record_size": rec_sz, "block_size": blk_sz}
        current_offset += blk_sz

    dat_path = output_prefix + ".dat"
    with open(dat_path, "wb") as f:
        f.write(GLOBAL_MAGIC)
        f.write(GLOBAL_VER_PACK.pack(GLOBAL_MAJOR, GLOBAL_MINOR, len(block_infos), b'\x00' * 16))
        for nid, info in block_infos.items():
            f.write(BLOCK_HDR_PACK.pack(nid, info["n_records"], 0, info["max_read_len"], info["max_cigar_len"]))
        # Pre-allocate the file
        if current_offset > GLOBAL_HEADER_SIZE:
            f.seek(current_offset - 1)
            f.write(b'\x00')

    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx:
        idx.write(struct.pack("<I", len(block_infos)))
        for nid, info in block_infos.items():
            idx.write(struct.pack("<I Q I I H I I", nid, info["offset"], info["block_size"],
                                  info["n_records"], 0, info["max_read_len"], info["max_cigar_len"]))

    return block_infos, dat_path, wanted_nodes

def load_existing_output_files(output_prefix):
    # This function would need to be implemented to read the .idx file
    # and reconstruct the block_infos dictionary if you use --use-existing.
    # For now, it's a placeholder.
    raise NotImplementedError("--use-existing is not fully implemented yet.")

# ─────────────────────────────────────────────────────────────────────────────
# Main Pipeline
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter, use_existing):
    if use_existing:
        print("Reusing existing .dat/.idx...")
        block_infos, dat_path, wanted_nodes = load_existing_output_files(output_prefix)
    else:
        print("Initializing output files from PKL maxima...")
        block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)
        print(f"Output file created: {dat_path}")

    BUFFER_SEGMENTS = 1_000_000
    next_milestone, total_reads, total_segments = milestone_step, 0, 0
    start_time = time.perf_counter()

    dat_fh = open(dat_path, "r+b")
    dat_fd = dat_fh.fileno() # Get the file descriptor once
    segment_buffer = defaultdict(list)

    # Single-threaded flush function using the C++ module
    def flush_segment_buffer():
        nonlocal total_segments
        if not segment_buffer:
            return

        print(f"Flushing {len(segment_buffer)} node blocks to disk...")
        flush_start_time = time.perf_counter()

        # Process each node's buffered segments in the main thread
        for nid, segs in segment_buffer.items():
            info = block_infos.get(nid)
            if not info: continue

            # 1. Python calculates the metadata (where to write, padding sizes)
            base_offset = info["offset"] + BLOCK_HDR_SIZE
            write_pos = base_offset + info["current_pos"] * info["record_size"]
            info["current_pos"] += len(segs)

            # 2. C++ does the heavy lifting: serialization and writing
            try:
                fast_writer.flush_node_buffer(
                    dat_fd,
                    write_pos,
                    segs,
                    info["max_read_len"],
                    info["max_cigar_len"]
                )
            except Exception as e:
                print(f"Error during C++ flush for node {nid}: {e}")
                # Decide how to handle errors, e.g., exit or continue
                raise

        segment_buffer.clear()
        total_segments = 0
        print(f"Flush complete in {time.perf_counter() - flush_start_time:.2f} seconds.")

    # Main processing loop
    for raw_msg in gam_record_iter(gam_path):
        segs_by_node = process_alignment(raw_msg, wanted_nodes, chrom_filter)
        total_reads += 1

        for nid, segs in segs_by_node.items():
            segment_buffer[nid].extend(segs)
            total_segments += len(segs)

        if total_segments >= BUFFER_SEGMENTS:
            flush_segment_buffer()

        if total_reads >= next_milestone:
            elapsed = time.perf_counter() - start_time
            print(f"{total_reads:,} reads processed | {elapsed:.1f} seconds")
            next_milestone += milestone_step

    # Final flush before exiting
    flush_segment_buffer()
    dat_fh.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads:,}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")

# ─────────────────────────────────────────────────────────────────────────────
# Entry Point
def main():
    parser = argparse.ArgumentParser(
        description="GAM segment extractor with C++-accelerated writing.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the stats pickle with node maxima")
    parser.add_argument("output_prefix", help="Prefix for output files (.dat, .idx)")
    parser.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    parser.add_argument("--use-existing", action="store_true", help="Reuse existing initialized output files")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr,
        use_existing=args.use_existing
    )

if __name__ == "__main__":
    main()