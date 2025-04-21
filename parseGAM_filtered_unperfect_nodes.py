#!/usr/bin/env python3
import argparse, gzip, pickle, json, time, concurrent.futures, vg_pb2, base64

# ---------- helpers ----------------------------------------------------------
def read_varint(stream):
    res = shift = 0
    while True:
        b = stream.read(1)
        if not b:
            raise EOFError
        byte = b[0]
        res |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            return res
        shift += 7


def is_gz(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def parse_gam_groups(path, tag="GAM"):
    opener = gzip.open if is_gz(path) else open
    with opener(path, "rb") as f:
        while True:
            try:
                n = read_varint(f)
            except EOFError:
                break
            if n == 0:
                continue
            try:
                tlen = read_varint(f)
                t = f.read(tlen).decode()
            except (EOFError, UnicodeDecodeError):
                break
            if t != tag:
                for _ in range(n - 1):
                    try:
                        f.seek(read_varint(f), 1)
                    except EOFError:
                        break
                continue
            msgs = []
            for _ in range(n - 1):
                try:
                    sz = read_varint(f)
                    blk = f.read(sz)
                except EOFError:
                    break
                if len(blk) == sz:
                    msgs.append(blk)
            yield msgs


# ---------- worker -----------------------------------------------------------
def proc_group(args):
    msgs, wanted = args
    out, cnt = {}, 0
    for mb in msgs:
        aln = vg_pb2.Alignment()
        aln.ParseFromString(mb)
        cnt += 1
        seq = aln.sequence
        qbytes = aln.quality  # raw per‑base bytes
        mapq = aln.mapping_quality
        roff = 0
        for m in aln.path.mapping:
            nid = m.position.node_id
            if nid not in wanted:
                for e in m.edit:
                    roff += max(e.from_length, len(e.sequence))
                continue
            off = m.position.offset
            seg_seq, seg_q = [], bytearray()
            for e in m.edit:
                if e.from_length:
                    seg_seq.append(
                        seq[roff : roff + e.from_length].lower()
                        if e.sequence
                        else seq[roff : roff + e.from_length]
                    )
                    seg_q.extend(qbytes[roff : roff + e.from_length])
                    roff += e.from_length
                elif e.sequence:  # insertion
                    ins = len(e.sequence)
                    seg_seq.append(seq[roff : roff + ins].lower())
                    seg_q.extend(qbytes[roff : roff + ins])
                    roff += ins
            out.setdefault(nid, []).append(
                {
                    "offset": off,
                    "sequence": "".join(seg_seq),
                    "base_quality": base64.b64encode(seg_q).decode(),
                    "read_quality": mapq,
                }
            )
    return out, cnt


def merge(parts):
    merged, reads = {}, 0
    for part, cnt in parts:
        reads += cnt
        for nid, segs in part.items():
            merged.setdefault(nid, []).extend(segs)
    return merged, reads


# ---------- pipeline ---------------------------------------------------------
def run(gam, stats_pkl, prefix, fmt, threads, max_pending, milestone):
    print(f"Loading stats pickle: {stats_pkl}")
    with open(stats_pkl, "rb") as f:
        stats = pickle.load(f)

    total_nodes = len(stats)
    unperfect_nodes = sum(1 for s in stats.values() if s["not_perfect"] > 0)

    filtered_nodes = {
        int(nid)
        for nid, s in stats.items()
        if s["not_perfect"] > 1
        and s["not_perfect"] / (s["perfect"] + s["not_perfect"]) > 0.1
    }

    print("\nNode‑level overview")
    print(f"  Total nodes in pickle      : {total_nodes}")
    print(
        f"  Nodes with ≥1 un‑perfect   : {unperfect_nodes} "
        f"({unperfect_nodes/total_nodes*100:.2f} %)"
    )
    print(
        f"  Nodes passing filter       : {len(filtered_nodes)} "
        f"({len(filtered_nodes)/total_nodes*100:.2f} %)\n"
    )

    parts, total, nxt, t0 = [], 0, milestone, time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        futs = []
        for msgs in parse_gam_groups(gam):
            futs.append(ex.submit(proc_group, (msgs, filtered_nodes)))
            while len(futs) >= max_pending:
                done, not_done = concurrent.futures.wait(
                    futs, return_when=concurrent.futures.FIRST_COMPLETED
                )
                for ft in done:
                    parts.append(ft.result())
                    total += ft.result()[1]
                futs = list(not_done)
                if total >= nxt:
                    print(
                        f"Milestone {total} reads  |  elapsed {time.perf_counter()-t0:.1f}s"
                    )
                    nxt += milestone
        for ft in concurrent.futures.as_completed(futs):
            parts.append(ft.result())
            total += ft.result()[1]
            if total >= nxt:
                print(
                    f"Milestone {total} reads  |  elapsed {time.perf_counter()-t0:.1f}s"
                )
                nxt += milestone

    merged, _ = merge(parts)

    # -------- pre‑write summary ---------------------------------------------
    print("\n--- pre‑write summary ---")
    print(f"  Reads processed : {total}")
    print(f"  Nodes saved     : {len(merged)}")
    print(f"  Elapsed so far  : {time.perf_counter()-t0:.2f}s")

    if fmt in ("json", "both"):
        jfile = prefix + ".json"
        with open(jfile, "w") as jf:
            json.dump(merged, jf, indent=2)
        print(f"  JSON written    : {jfile}")
    if fmt in ("pkl", "both"):
        pfile = prefix + ".pkl"
        with open(pfile, "wb") as pf:
            pickle.dump(merged, pf, pickle.HIGHEST_PROTOCOL)
        print(f"  Pickle written  : {pfile}")

    # -------- final summary -------------------------------------------------
    print("\n--- final summary ---")
    print(f"  Total reads processed : {total}")
    print(f"  Total nodes saved     : {len(merged)}")
    print(f"  Total time            : {time.perf_counter()-t0:.2f}s\n")


# ---------- CLI --------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Extract node‑mapped read segments with high un‑perfect ratio"
    )
    ap.add_argument("gam_file")
    ap.add_argument("stats_pickle")
    ap.add_argument("output_prefix")
    ap.add_argument(
        "--save-format",
        choices=["json", "pkl", "both"],
        default="both",
        help="Which file type(s) to write (default both)",
    )
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--max_pending", type=int, default=16)
    ap.add_argument("--milestone", type=int, default=100_000_000)
    args = ap.parse_args()

    run(
        args.gam_file,
        args.stats_pickle,
        args.output_prefix,
        args.save_format,
        args.threads,
        args.max_pending,
        args.milestone,
    )


if __name__ == "__main__":
    main()
