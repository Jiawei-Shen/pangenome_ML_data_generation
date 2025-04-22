#!/usr/bin/env python3
"""
Aggregate perfect / un‑perfect read‑segments per VG node
(with strand) and save as JSON / Pickle / both – single‑threaded.
"""
import argparse, gzip, json, pickle, time, base64
import vg_pb2                     # protobufs generated from vg.proto

# ───────── helpers ──────────────────────────────────────────────────────────
def read_varint(stream):
    value = shift = 0
    while True:
        b = stream.read(1)
        if not b:
            raise EOFError
        byte = b[0]
        value |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            return value
        shift += 7

def file_is_gzip(path: str) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"

def parse_gam_groups(path: str, tag="GAM"):
    opener = gzip.open if file_is_gzip(path) else open
    with opener(path, "rb") as fh:
        while True:
            try:
                group_count = read_varint(fh)
            except EOFError:
                break
            if group_count == 0:
                continue
            try:
                tag_len   = read_varint(fh)
                group_tag = fh.read(tag_len).decode()
            except (EOFError, UnicodeDecodeError):
                break

            # skip groups that aren’t GAM
            if group_tag != tag:
                for _ in range(group_count - 1):
                    try:
                        skip = read_varint(fh)
                        fh.seek(skip, 1)
                    except EOFError:
                        break
                continue

            msgs = []
            for _ in range(group_count - 1):
                try:
                    size = read_varint(fh)
                    msg  = fh.read(size)
                except EOFError:
                    break
                if len(msg) == size:
                    msgs.append(msg)
            if msgs:
                yield msgs

# ───────── worker (now local call) ─────────────────────────────────────────
def process_group(message_list, wanted_nodes, chrom_filter):
    node_segments, read_count = {}, 0
    for raw in message_list:
        aln = vg_pb2.Alignment(); aln.ParseFromString(raw)

        if chrom_filter and not any(rp.name == chrom_filter for rp in aln.refpos):
            continue

        read_count += 1
        seq, qual_bytes, mapq = aln.sequence, aln.quality, aln.mapping_quality
        read_offset = 0

        for mapping in aln.path.mapping:
            nid = mapping.position.node_id
            if nid not in wanted_nodes:
                for e in mapping.edit:
                    read_offset += max(e.from_length, len(e.sequence))
                continue

            node_offset  = mapping.position.offset
            strand_char  = "-" if mapping.position.is_reverse else "+"

            part_seq, part_q = [], bytearray()
            for edit in mapping.edit:
                if edit.from_length:
                    frag = seq[read_offset : read_offset + edit.from_length]
                    part_seq.append(frag.lower() if edit.sequence else frag)
                    part_q.extend(qual_bytes[read_offset : read_offset + edit.from_length])
                    read_offset += edit.from_length
                elif edit.sequence:          # insertion
                    ins_len = len(edit.sequence)
                    part_seq.append(seq[read_offset : read_offset + ins_len].lower())
                    part_q.extend(qual_bytes[read_offset : read_offset + ins_len])
                    read_offset += ins_len

            node_segments.setdefault(nid, []).append({
                "offset"      : node_offset,
                "sequence"    : "".join(part_seq),
                "base_quality": base64.b64encode(part_q).decode(),
                "read_quality": mapq,
                "strand"      : strand_char
            })
    return node_segments, read_count

def merge_partial(parts):
    merged, total = {}, 0
    for part, cnt in parts:
        total += cnt
        for nid, segs in part.items():
            merged.setdefault(nid, []).extend(segs)
    return merged, total

# ───────── pipeline ────────────────────────────────────────────────────────
def run_pipeline(gam_path, stats_pkl, prefix, fmt, milestone, chrom):
    print(f"Loading stats: {stats_pkl}")
    with open(stats_pkl, "rb") as fh:
        stats = pickle.load(fh)

    total_nodes          = len(stats)
    nodes_with_unperfect = sum(1 for s in stats.values() if s["not_perfect"] > 0)
    filtered_nodes = {
        int(nid)
        for nid, s in stats.items()
        if s["not_perfect"] > 1 and s["not_perfect"] / (s["perfect"] + s["not_perfect"]) > 0.1
    }

    print("\nNode‑level overview")
    print(f"  Total nodes               : {total_nodes}")
    print(f"  Nodes with ≥1 un‑perfect  : {nodes_with_unperfect} "
          f"({nodes_with_unperfect/total_nodes*100:.2f} %)")
    print(f"  Nodes passing filter      : {len(filtered_nodes)} "
          f"({len(filtered_nodes)/total_nodes*100:.2f} %)\n")

    partials, reads_total, next_mstone = [], 0, milestone
    start = time.perf_counter()

    for batch in parse_gam_groups(gam_path):
        part = process_group(batch, filtered_nodes, chrom)
        partials.append(part)
        reads_total += part[1]

        if reads_total >= next_mstone:
            print(f"Milestone {reads_total} reads | {time.perf_counter()-start:.1f}s")
            next_mstone += milestone

    merged, _ = merge_partial(partials)

    # ─── write out ──────────────────────────────────────────────────────────
    print("\n--- pre‑write summary ---")
    print(f"  Reads processed : {reads_total}")
    print(f"  Nodes saved     : {len(merged)}")
    print(f"  Elapsed so far  : {time.perf_counter()-start:.2f}s")

    if fmt in ("json", "both"):
        jfile = prefix + ".json"
        with open(jfile, "w") as jf:
            json.dump({str(k): v for k, v in merged.items()}, jf, indent=2)
        print(f"  JSON written    : {jfile}")
    if fmt in ("pkl", "both"):
        pfile = prefix + ".pkl"
        with open(pfile, "wb") as pf:
            pickle.dump(merged, pf, pickle.HIGHEST_PROTOCOL)
        print(f"  Pickle written  : {pfile}")

    print("\n--- final summary ---")
    print(f"  Reads processed : {reads_total}")
    print(f"  Nodes saved     : {len(merged)}")
    print(f"  Total time      : {time.perf_counter()-start:.2f}s\n")

# ───────── CLI ──────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Extract node‑mapped read segments with strand and high un‑perfect ratio (single‑thread)."
    )
    ap.add_argument("gam_file")
    ap.add_argument("stats_pickle")
    ap.add_argument("output_prefix")
    ap.add_argument("--save-format", choices=["json", "pkl", "both"],
                    default="json", help="File format(s) to write")
    ap.add_argument("--milestone", type=int, default=100_000_000,
                    help="Progress update frequency (reads)")
    ap.add_argument("--chr", default="", help="Restrict to chromosome name in refpos")
    args = ap.parse_args()

    run_pipeline(
        args.gam_file, args.stats_pickle, args.output_prefix,
        args.save_format, args.milestone, args.chr
    )

if __name__ == "__main__":
    main()
