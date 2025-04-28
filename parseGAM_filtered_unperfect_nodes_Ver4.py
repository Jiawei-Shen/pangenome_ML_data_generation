#!/usr/bin/env python3
import argparse
import gzip
import pickle
import struct
import time
import gc
import psutil
import tracemalloc
from collections import defaultdict
import vg_pb2

# ─────────────────────────────────────────────────────────────────────────────
class Segment:
    __slots__ = ('offset', 'seq', 'bq', 'rq', 'strand')
    def __init__(self, offset, seq, bq, rq, strand):
        self.offset = offset
        self.seq    = seq
        self.bq     = bq
        self.rq     = rq
        self.strand = strand

RECORD_STRUCT = struct.Struct("<h150s150shc")
RECORD_SIZE = RECORD_STRUCT.size

# ─────────────────────────────────────────────────────────────────────────────
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
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"

def gam_record_iter(path, tag="GAM"):
    open_func = gzip.open if file_is_gzip(path) else open
    # test = b'\n~TCCCTGAGGTGGTGGCGGAGGTGGTGGAGGGGCGGAGGGCGGAGCACCGTAGCCCCCTCTGGCCCGACTCGGGGCGGCCCGATTGCCCCGGTCCCAGCAGCCCTCCAGGGCCTCCAGGCCCCGGCC\x12o\x12\x11\n\x07\x08\xd3\xd7\xd8*\x10\x0f\x12\x04\x08\x11\x10\x11(\x01\x12\x0f\n\x05\x08\xd4\xd7\xd8*\x12\x04\x08 \x10 (\x02\x12\x1e\n\x05\x08\xd5\xd7\xd8*\x12\x04\x08\x19\x10\x19\x12\x07\x08\x01\x10\x01\x1a\x01C\x12\x04\x08\x06\x10\x06(\x03\x12\x0f\n\x05\x08\xd6\xd7\xd8*\x12\x04\x08 \x10 (\x04\x12\x18\n\x05\x08\xd7\xd7\xd8*\x12\x04\x08\x0c\x10\x0c\x12\x07\x08\x01\x10\x01\x1a\x01C(\x05\x1a\x14ERR903030.252281075 "~\x1e \x1f!""\x10\x1f\x0f$\x1f\x0f\x1d\x0f%!\x1d\x0e\x0e\x0e\x1b$\x1f\x0f\x0e\x0f\x18$&\x18"& \x1d\x0e\x1b\x0e\x0e\x0e\x1b\x0e\x19\x0e""!&&\x0e\x19\x19\x19\x1f\x0f\x0f\x19"\x0e$\x0f\x0f$\x0f\x1f#\x0e"\x1e\x18"\x0e\x0e\x19\x0e\x0e\x19\x0e\x18\x0e\x0e""&\x1b"#"\x0f\r"\x1a"\x1a&!\x0e\x0e\x19\x16\r\r\x17\x15\x1f\r"#"\x15#\r\x15\r\x19\x1f\x17\x0e\r\x17""\x02\x02\x02\x02\x02(<0~Z\x16\x1a\x14ERR903030.252281075 \x81\x01\xdf\xf7}\xdf\xf7}\xef?\x8a\x01\x12\n\x05chr22 \xc9\xfb\xff\xff\xff\xff\xff\xff\xff\x01\x9a\x01\x0c\x10\xf5\xf1\x94\x08*\x05chr22\xf1\x01\x00\x00\x00\x00\x00\x80I@\x82\x02\x171927:443.135:148.44:0:1\x99\x02\x00\x00\x00\x00\x00\x14\x93@'
    # while True:
    #     yield test
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


def process_alignment(raw_message, wanted_nodes, chrom_filter, alignment):
    segment_dict = {}
    alignment.Clear()
    alignment.ParseFromString(raw_message)

    if chrom_filter and not any(pos.name == chrom_filter for pos in alignment.refpos):
        return segment_dict, 0

    read_sequence = alignment.sequence
    read_quality_bytes = alignment.quality
    mapping_quality = alignment.mapping_quality
    read_offset = 0
    reads = 1

    for mapping in alignment.path.mapping:
        node_id = mapping.position.node_id

        if node_id not in wanted_nodes:
            for edit in mapping.edit:
                read_offset += max(edit.from_length, len(edit.sequence))
            continue

        node_offset = mapping.position.offset
        strand_char = b"-" if mapping.position.is_reverse else b"+"
        sequence_parts = []
        quality_bytes = bytearray()

        for edit in mapping.edit:
            if edit.from_length:
                frag = read_sequence[read_offset: read_offset + edit.from_length]
                sequence_parts.append(frag.lower() if edit.sequence else frag)
                quality_bytes.extend(read_quality_bytes[read_offset: read_offset + edit.from_length])
                read_offset += edit.from_length
            elif edit.sequence:
                ins_len = len(edit.sequence)
                frag = read_sequence[read_offset: read_offset + ins_len]
                sequence_parts.append(frag.lower())
                quality_bytes.extend(read_quality_bytes[read_offset: read_offset + ins_len])
                read_offset += ins_len

        seg = Segment(
            node_offset,
            "".join(sequence_parts).encode().ljust(150, b'\x00')[:150],
            bytes(quality_bytes).ljust(150, b'\x00')[:150],
            mapping_quality,
            strand_char
        )
        segment_dict.setdefault(node_id, []).append(seg)

    return segment_dict, reads


# ─────────────────────────────────────────────────────────────────────────────
def initialize_output_files(stats_path, output_prefix):
    with open(stats_path, "rb") as stats_file:
        stats_data = pickle.load(stats_file)

    block_infos = {}
    wanted_nodes = set()

    current_offset = 0
    total_nodes = 0

    for node_id_str, stat in stats_data.items():
        total_nodes += 1
        node_id = int(node_id_str)
        perfect = stat["perfect"]
        not_perfect = stat["not_perfect"]
        if not_perfect > 1 and not_perfect / (perfect + not_perfect) > 0.10:
            wanted_nodes.add(node_id)
            n_records = perfect + not_perfect

            block_infos[node_id] = {
                "offset": current_offset,
                "n_records": n_records,
                "current_pos": 0
            }
            current_offset += 4 + 4 + 2 + n_records * RECORD_SIZE

    print(f"Filtered {len(wanted_nodes)} nodes from {total_nodes} total nodes "
          f"({len(wanted_nodes) / total_nodes:.2%} selected).")
    del stats_data
    gc.collect()

    dat_path = output_prefix + ".dat"
    #
    # with open(dat_path, "wb") as f:
    #     f.write(b"MYFMT\1")
    #     f.write(struct.pack("<BBI16s", 0, 1, len(block_infos), b'\x00' * 16))
    #
    #     # 直接顺序写所有block
    #     blank_record = RECORD_STRUCT.pack(0, b'\x00'*150, b'\x00'*150, 0, b'+')
    #
    #     for node_id, info in block_infos.items():
    #         # 构建整个block的内容 (block header + n_records blank)
    #         block_header = struct.pack("<I I H", node_id, info["n_records"], 0)
    #         records_blob = blank_record * info["n_records"]
    #         full_block = block_header + records_blob
    #         f.write(full_block)
    #
    # # idx文件（照旧，单线程）
    # idx_path = output_prefix + ".idx"
    # with open(idx_path, "wb") as idx_file:
    #     idx_file.write(struct.pack("<I", len(block_infos)))
    #     for node_id, info in block_infos.items():
    #         idx_file.write(struct.pack("<I Q I I H", node_id, info["offset"],
    #                                    4 + 4 + 2 + info["n_records"] * RECORD_SIZE, info["n_records"], 0))

    return block_infos, dat_path, wanted_nodes

# ─────────────────────────────────────────────────────────────────────────────
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter):
    print(f"Initializing output files...")
    block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)
    print(f"Output file created: {dat_path}")

    # BUFFER_SEGMENTS = 30_000_000  # 每累积2000万条segments flush一次
    BUFFER_SEGMENTS = 100_000  # 每累积2000万条segments flush一次

    next_milestone = milestone_step
    total_reads = 0
    total_segments = 0
    start_time = time.perf_counter()

    dat_fh = open(dat_path, "r+b")
    segment_buffer = defaultdict(list)
    alignment = vg_pb2.Alignment()

    process = psutil.Process()
    baseline_memory = process.memory_info().rss
    print(f"[Info] baseline_memory: {baseline_memory / 1024 / 1024:.2f} MB")

    # 在run_pipeline开始时加：
    tracemalloc.start()
    flush_counter = 0

    def flush_segment_buffer():
        nonlocal total_segments, baseline_memory, segment_buffer, flush_counter
        for node_id, segs in segment_buffer.items():
            if not segs:
                continue
            info = block_infos[node_id]
            base_offset = info["offset"] + 4 + 4 + 2

            batch_blob = bytearray()
            for seg in segs:
                batch_blob += RECORD_STRUCT.pack(seg.offset, seg.seq, seg.bq, seg.rq, seg.strand)

            pos = base_offset + info["current_pos"] * RECORD_SIZE
            dat_fh.seek(pos)
            dat_fh.write(batch_blob)

            info["current_pos"] += len(segs)

        # flush结束
        segment_buffer = defaultdict(list)  # 重新new一个
        total_segments = 0
        flush_counter += 1

        # 强制GC
        gc.collect()

        # 记录内存
        current_memory = process.memory_info()
        print(
            f"[Info] Memory usage after flush #{flush_counter}: RSS={current_memory.rss / 1024 / 1024:.2f} MB, VMS={current_memory.vms / 1024 / 1024:.2f} MB")

        # 每隔5次flush做一次tracemalloc snapshot
        if flush_counter % 5 == 0:
            snapshot = tracemalloc.take_snapshot()
            top_stats = snapshot.statistics('lineno')

            print("[tracemalloc] Top 5 memory allocations:")
            for stat in top_stats[:5]:
                print(stat)

    for raw_msg in gam_record_iter(gam_path):
        segment_dict, read_count = process_alignment(raw_msg, wanted_nodes, chrom_filter, alignment)
        total_reads += read_count

        for node_id, segs in segment_dict.items():
            segment_buffer[node_id].extend(segs)
            total_segments += len(segs)
            print(raw_msg)

        if total_segments >= BUFFER_SEGMENTS:
            flush_segment_buffer()


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
    parser = argparse.ArgumentParser(description="GAM segment extractor into binary blocks (single-threaded).")
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the node stats pickle file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--milestone", type=int, default=10_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr
    )

if __name__ == "__main__":
    main()
