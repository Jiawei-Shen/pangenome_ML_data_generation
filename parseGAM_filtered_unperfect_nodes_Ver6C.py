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


# ─────────────────────────────────────────────────────────────────────────────
# Segment container (stores raw, unpadded bytes for seq/bq/cigar; padded at pack time)
class Segment:
    __slots__ = ('offset', 'seq', 'bq', 'cigar', 'rq', 'strand')
    def __init__(self, offset, seq, bq, cigar, rq, strand):
        self.offset = offset
        self.seq    = seq      # bytes (unpadded)
        self.bq     = bq       # bytes (unpadded)
        self.cigar  = cigar    # bytes (unpadded)
        self.rq     = rq       # int (MAPQ)
        self.strand = strand   # b'+' or b'-'


# ─────────────────────────────────────────────────────────────────────────────
# File layout (STRICT latest format only: no node_length; per-block maxima)

# Global header: magic + version block
GLOBAL_MAGIC = b"MYFMT\x01"                              # 6 bytes: "MYFMT" + 0x01
GLOBAL_VER_PACK = struct.Struct("<BBI16s")               # major, minor, block_count, reserved[16]
GLOBAL_MAJOR_EXPECTED, GLOBAL_MINOR_EXPECTED = 0, 5      # STRICT: require version 0.5
GLOBAL_HEADER_SIZE = len(GLOBAL_MAGIC) + GLOBAL_VER_PACK.size  # 6 + 22 = 28

# Per-node block header (NO node_length)
# <I I H I I> -> node_id (u32), n_records (u32), flags (u16=0),
#                max_read_length (u32), max_cigar_length (u32)
BLOCK_HDR_PACK = struct.Struct("<I I H I I")
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size  # 18

def make_record_struct(max_read_len: int, max_cigar_len: int) -> struct.Struct:
    """
    Build the struct for a single segment record for given maxima.
    Layout: <h {R}s {R}s {C}s h c
      - i16 offset
      - seq[R] bytes
      - bq[R] bytes
      - cigar[C] bytes (ASCII, null-padded to C)
      - i16 rq (MAPQ)
      - char strand ('+' / '-')
    """
    if max_read_len <= 0 or max_cigar_len <= 0:
        raise ValueError("max_read_len and max_cigar_len must be > 0 for latest format")
    return struct.Struct(f"<h{max_read_len}s{max_read_len}s{max_cigar_len}shc")

def record_size(max_read_len: int, max_cigar_len: int) -> int:
    return make_record_struct(max_read_len, max_cigar_len).size


# ─────────────────────────────────────────────────────────────────────────────
# GAM parsing

def read_varint(stream):
    value = 0
    shift_amount = 0
    while True:
        byte_pair = stream.read(1)
        if not byte_pair:
            raise EOFError("EOF while reading varint")
        byte_value = byte_pair[0]
        value |= (byte_value & 0x7F) << shift_amount
        if not (byte_value & 0x80):
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
                    skip_len = read_varint(f)
                    f.seek(skip_len, 1)
                continue
            for _ in range(group_count - 1):
                try:
                    msg_size = read_varint(f)
                    yield f.read(msg_size)
                except EOFError:
                    break


# ─────────────────────────────────────────────────────────────────────────────
# CIGAR builder (unchanged semantics)
def build_cigar(mapping_edits):
    cigar_parts = []
    for edit in mapping_edits:
        from_len = edit.from_length
        to_len   = edit.to_length
        edit_len = len(edit.sequence)

        if from_len == to_len:
            if edit_len == 0:
                cigar_parts.append(f"{from_len}M")  # match
            else:
                cigar_parts.append(f"{from_len}X")  # substitution
        elif from_len > to_len:
            cigar_parts.append(f"{from_len - to_len}D")      # deletion
        elif from_len < to_len:
            cigar_parts.append(f"{to_len - from_len}I")      # insertion
        else:
            raise ValueError(f"Unexpected edit: from_length={from_len}, to_length={to_len}")
    return "".join(cigar_parts)


# ─────────────────────────────────────────────────────────────────────────────
# Alignment → Segment (same logic as your original)
def process_alignment(raw_message, wanted_nodes, chrom_filter):
    segment_dict = {}
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw_message)

    # Filter out low mapping quality
    if alignment.mapping_quality <= 10:
        return segment_dict

    if chrom_filter and not any(pos.name == chrom_filter for pos in alignment.refpos):
        return segment_dict

    read_sequence = alignment.sequence
    read_quality  = alignment.quality
    mapping_quality = alignment.mapping_quality
    read_offset = 0

    for mapping in alignment.path.mapping:
        node_id = mapping.position.node_id

        # Advance read_offset even if node not wanted, to keep position in read correct
        if node_id not in wanted_nodes:
            for edit in mapping.edit:
                read_offset += edit.to_length
            continue

        node_offset = mapping.position.offset
        strand_char = b"-" if mapping.position.is_reverse else b"+"

        sequence_parts = []
        quality_parts = bytearray()
        cigar_parts = []

        for edit in mapping.edit:
            from_len = edit.from_length
            to_len   = edit.to_length
            edit_len = len(edit.sequence)

            # CIGAR part (same semantics)
            if from_len == to_len:
                if edit_len == 0:
                    cigar_parts.append(f"{from_len}M")
                else:
                    cigar_parts.append(f"{from_len}X")
            elif from_len > 0 and to_len == 0:
                cigar_parts.append(f"{from_len}D")
            elif from_len == 0 and to_len > 0:
                cigar_parts.append(f"{to_len}I")
            else:
                raise ValueError(f"Unexpected edit: from_length={from_len}, sequence_length={edit_len}")

            # Append sequence/quality by to_len
            edit_length = to_len
            sequence_fragment = read_sequence[read_offset: read_offset + edit_length]
            quality_fragment  = read_quality[read_offset: read_offset + edit_length]
            sequence_parts.append(sequence_fragment.upper())
            quality_parts.extend(quality_fragment)
            read_offset += edit_length

        cigar_string = "".join(cigar_parts).encode()
        seq_bytes = "".join(sequence_parts).encode()
        bq_bytes  = bytes(quality_parts)

        seg = Segment(
            offset=node_offset,
            seq=seq_bytes,          # unpadded; will be padded to max_read_length on write
            bq=bq_bytes,            # unpadded; will be padded to max_read_length on write
            cigar=cigar_string,     # unpadded; will be padded to max_cigar_length on write
            rq=mapping_quality,
            strand=strand_char
        )
        segment_dict.setdefault(node_id, []).append(seg)

    return segment_dict


# ─────────────────────────────────────────────────────────────────────────────
# STRICT header validation helpers (latest only)

def _validate_latest_dat_header(dat_path):
    with open(dat_path, "rb") as df:
        magic = df.read(len(GLOBAL_MAGIC))
        if magic != GLOBAL_MAGIC:
            raise RuntimeError(f"Invalid .dat magic (got {magic!r}, want {GLOBAL_MAGIC!r})")
        major, minor, _count, _ = GLOBAL_VER_PACK.unpack(df.read(GLOBAL_VER_PACK.size))
        if (major, minor) != (GLOBAL_MAJOR_EXPECTED, GLOBAL_MINOR_EXPECTED):
            raise RuntimeError(f"Unsupported .dat version {major}.{minor}; require {GLOBAL_MAJOR_EXPECTED}.{GLOBAL_MINOR_EXPECTED}")

def _validate_latest_idx_header_and_size(idx_path):
    with open(idx_path, "rb") as f:
        data = f.read()
    if len(data) < 4:
        raise RuntimeError("Corrupt .idx: too small to contain count")
    (count,) = struct.unpack("<I", data[:4])
    expected_size = 4 + count * 30  # STRICT latest entry size
    if len(data) != expected_size:
        raise RuntimeError(f"Unsupported .idx layout/size: got {len(data)} bytes; expected {expected_size} for {count} entries of 30B")


# ─────────────────────────────────────────────────────────────────────────────
# Initialize outputs (read per-node maxima directly from stats PKL; no node_length)
def initialize_output_files(stats_path, output_prefix):
    with open(stats_path, "rb") as stats_file:
        stats_data = pickle.load(stats_file)

    block_infos = {}
    wanted_nodes = set()

    current_offset = GLOBAL_HEADER_SIZE
    total_nodes = 0

    for node_id_key, stat in stats_data.items():
        total_nodes += 1
        node_id = int(node_id_key)

        perfect = int(stat.get("perfect", 0))
        not_perfect = int(stat.get("not_perfect", 0))

        # pull maxima from PKL (fallback to 1 if missing/invalid)
        R = int(stat.get("max_read_length", 1) or 1)
        C = int(stat.get("max_cigar_length", 1) or 1)
        if R <= 0 or C <= 0:
            raise RuntimeError(f"Invalid maxima in stats for node {node_id}: R={R}, C={C}")

        # your selection rule (unchanged)
        if (perfect + not_perfect) > 0 and not_perfect > 1 and not_perfect / (perfect + not_perfect) > 0.05:
            wanted_nodes.add(node_id)
            n_records = perfect + not_perfect

            rec_sz = record_size(R, C)
            blk_sz = BLOCK_HDR_SIZE + n_records * rec_sz

            block_infos[node_id] = {
                "offset": current_offset,   # absolute from file start
                "n_records": n_records,
                "current_pos": 0,
                "max_read_len": R,
                "max_cigar_len": C,
                "record_size": rec_sz,
                "block_size": blk_sz,
            }
            current_offset += blk_sz

    print(f"Filtered {len(wanted_nodes)} nodes from {total_nodes} total nodes "
          f"({(len(wanted_nodes) / max(total_nodes,1)):.2%} selected).")
    del stats_data
    gc.collect()

    dat_path = output_prefix + ".dat"

    # Write .dat file: global header + blocks (STRICT latest format)
    with open(dat_path, "wb") as f:
        # Global header
        f.write(GLOBAL_MAGIC)
        f.write(GLOBAL_VER_PACK.pack(GLOBAL_MAJOR_EXPECTED, GLOBAL_MINOR_EXPECTED, len(block_infos), b'\x00' * 16))

        # Blocks
        for node_id, info in block_infos.items():
            n_records = info["n_records"]
            R = info["max_read_len"]
            C = info["max_cigar_len"]

            # Block header with maxima
            f.write(BLOCK_HDR_PACK.pack(node_id, n_records, 0, R, C))

            # Preallocate blank records for the block
            rec_pack = make_record_struct(R, C)
            blank = rec_pack.pack(
                0,
                b'\x00' * R,
                b'\x00' * R,
                b'\x00' * C,
                0,
                b'+'
            )
            for _ in range(n_records):
                f.write(blank)

    # Write .idx: store maxima (STRICT latest format: 30B/entry)
    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx_file:
        # count
        idx_file.write(struct.pack("<I", len(block_infos)))
        # entries: node_id (u32), offset (u64), block_size (u32), n_records (u32),
        #          flags (u16), max_read_len (u32), max_cigar_len (u32)
        for node_id, info in block_infos.items():
            idx_file.write(struct.pack(
                "<I Q I I H I I",
                node_id,
                info["offset"],
                info["block_size"],
                info["n_records"],
                0,
                info["max_read_len"],
                info["max_cigar_len"]
            ))

    return block_infos, dat_path, wanted_nodes


# ─────────────────────────────────────────────────────────────────────────────
# Reuse existing .dat/.idx (STRICT latest format only)
def load_existing_output_files(output_prefix):
    idx_path = output_prefix + ".idx"
    dat_path = output_prefix + ".dat"

    if not (os.path.exists(idx_path) and os.path.exists(dat_path)):
        raise FileNotFoundError(f"Expected existing files: {idx_path} and {dat_path}")

    # Strictly validate headers/sizes first
    _validate_latest_idx_header_and_size(idx_path)
    _validate_latest_dat_header(dat_path)

    # Parse .idx (fixed 30-byte entries)
    entries = []
    with open(idx_path, "rb") as f:
        (count,) = struct.unpack("<I", f.read(4))
        for _ in range(count):
            data = f.read(30)
            if len(data) != 30:
                raise RuntimeError("Corrupt .idx: truncated entry")
            node_id, offset, block_size, n_records, flags, R_idx, C_idx = struct.unpack("<I Q I I H I I", data)
            entries.append((node_id, offset, block_size, n_records, flags, R_idx, C_idx))

    # Verify .dat header & block_count; trust maxima from .dat and validate strictly
    with open(dat_path, "rb") as df:
        magic = df.read(len(GLOBAL_MAGIC))
        if magic != GLOBAL_MAGIC:
            raise RuntimeError("Invalid .dat magic/version")
        major, minor, dat_count, _ = GLOBAL_VER_PACK.unpack(df.read(GLOBAL_VER_PACK.size))
        if (major, minor) != (GLOBAL_MAJOR_EXPECTED, GLOBAL_MINOR_EXPECTED):
            raise RuntimeError(f"Unsupported .dat version {major}.{minor}; require {GLOBAL_MAJOR_EXPECTED}.{GLOBAL_MINOR_EXPECTED}")
        if dat_count != len(entries):
            raise RuntimeError(f".dat block_count ({dat_count}) != .idx count ({len(entries)}) for latest format")

        block_infos = {}
        wanted_nodes = set()
        for node_id, offset, block_size, n_records, flags, _R_idx, _C_idx in entries:
            df.seek(offset, os.SEEK_SET)
            hdr = df.read(BLOCK_HDR_SIZE)
            if len(hdr) != BLOCK_HDR_SIZE:
                raise RuntimeError(f"Corrupt .dat: cannot read block header at {offset}")
            nid2, nrec2, flg2, R2, C2 = struct.unpack("<I I H I I", hdr)
            if nid2 != node_id or nrec2 != n_records:
                raise RuntimeError(f".dat/.idx mismatch for node {node_id} (idx n={n_records}, dat n={nrec2})")
            if R2 <= 0 or C2 <= 0:
                raise RuntimeError(f"Invalid maxima in .dat for node {node_id}: R={R2}, C={C2}")
            rec_sz = record_size(R2, C2)
            block_infos[node_id] = {
                "offset": offset,
                "n_records": n_records,
                "current_pos": 0,       # assumes empty blocks; not resuming partial writes
                "max_read_len": R2,
                "max_cigar_len": C2,
                "record_size": rec_sz,
                "block_size": block_size,
            }
            wanted_nodes.add(node_id)

    print(f"Reusing existing output with {len(block_infos)} node blocks from {idx_path}")
    return block_infos, dat_path, wanted_nodes


# ─────────────────────────────────────────────────────────────────────────────
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter, use_existing):
    if use_existing:
        print("Reusing existing .dat/.idx (strict latest format)...")
        block_infos, dat_path, wanted_nodes = load_existing_output_files(output_prefix)
    else:
        print("Initializing output files (reading maxima from PKL, strict latest format)...")
        block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)
        print(f"Output file created: {dat_path}")

    BUFFER_SEGMENTS = 400_000_000  # number of segments buffered before flushing


    next_milestone = milestone_step
    total_reads = 0
    total_segments = 0
    start_time = time.perf_counter()

    dat_fh = open(dat_path, "r+b")
    segment_buffer = defaultdict(list)

    def flush_segment_buffer():
        nonlocal total_segments
        if not segment_buffer:
            return

        count_i = 0
        for node_id, segs in segment_buffer.items():
            if not segs:
                continue

            info = block_infos[node_id]
            base_offset = info["offset"] + BLOCK_HDR_SIZE  # start of records for this block

            R = info["max_read_len"]
            C = info["max_cigar_len"]
            rec_pack = make_record_struct(R, C)
            rec_size = info["record_size"]

            # Absolute byte range we’re about to write for this node
            start = base_offset + info["current_pos"] * rec_size
            nbytes = len(segs) * rec_size

            # Map only the needed window, aligned to the system allocation granularity
            page = mmap.ALLOCATIONGRANULARITY
            aligned_off = (start // page) * page
            rel0 = start - aligned_off
            map_len = rel0 + nbytes
            if map_len <= 0:
                continue  # nothing to write

            fileno = dat_fh.fileno()
            with mmap.mmap(fileno, length=map_len, offset=aligned_off, access=mmap.ACCESS_WRITE) as mm:
                rel = rel0
                for seg in segs:
                    rec_pack.pack_into(
                        mm, rel,
                        int(seg.offset),
                        seg.seq.ljust(R, b'\x00')[:R],
                        seg.bq.ljust(R, b'\x00')[:R],
                        seg.cigar.ljust(C, b'\x00')[:C],
                        int(seg.rq),
                        seg.strand if seg.strand in (b'+', b'-') else b'+'
                    )
                    rel += rec_size
                mm.flush()  # ensure OS sees the writes
            count_i += 1
            print(f"\rProgress: {count_i})", end="", flush=True)

            info["current_pos"] += len(segs)

        segment_buffer.clear()
        total_segments = 0

    for raw_msg in gam_record_iter(gam_path):
        segment_dict = process_alignment(raw_msg, wanted_nodes, chrom_filter)
        total_reads += 1

        for node_id, segs in segment_dict.items():
            segment_buffer[node_id].extend(segs)
            total_segments += len(segs)

        if total_segments >= BUFFER_SEGMENTS:
            flush_segment_buffer()
            segment_buffer = defaultdict(list)

        if total_reads >= next_milestone:
            elapsed = time.perf_counter() - start_time
            print(f"{total_reads} reads processed | {elapsed:.1f} seconds")
            next_milestone += milestone_step

    flush_segment_buffer()
    dat_fh.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")


# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="GAM segment extractor with per-block maxima (STRICT latest .dat/.idx format only); maxima pulled from stats PKL."
    )
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the node stats pickle file (must include max_read_length/max_cigar_length)")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    parser.add_argument("--use-existing", action="store_true",
                        help="Reuse existing initialized output (output_prefix.dat/.idx) — latest format only")
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
