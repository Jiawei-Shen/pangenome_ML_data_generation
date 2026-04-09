#!/usr/bin/env python3
"""
Check all records in a changed-representation TSV efficiently.

Features:
- no pandas dependency
- handles hap1 / hap2
- handles multiple chromosomes
- scans each BAM only once
- maps assembly positions to hg38 positions from BAM
- optional graph node lookup from XG using vg + jq
- writes a full TSV and a summary JSON

Input TSV expected columns include at least:
#CHROM, POS, REF, ALT, FULL_ANNOTATION

Example FULL_ANNOTATION:
    chr1_hap1:1195496-C-CA
    chr7_hap2:12345-AT-A

Usage:
python check_changed_repr_variants_nopandas.py \
  --tsv changed_representation_variants.tsv \
  --hap1-bam HG008N_curatedv6_250714_polished6.2.hap1.bam \
  --hap2-bam HG008N_curatedv6_250714_polished6.2.hap2.bam \
  --xg /path/to/graph.xg \
  --out-prefix changed_repr_check \
  --graph-flank 2
"""

import argparse
import csv
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from typing import Any, Dict, List, Optional, Tuple

import pysam


# BAM CIGAR op codes
CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_P = 6
CIGAR_EQ = 7
CIGAR_X = 8

CIGAR_CODE_TO_CHAR = {
    CIGAR_M: "M",
    CIGAR_I: "I",
    CIGAR_D: "D",
    CIGAR_N: "N",
    CIGAR_S: "S",
    CIGAR_H: "H",
    CIGAR_P: "P",
    CIGAR_EQ: "=",
    CIGAR_X: "X",
}


FULL_ANNOT_RE = re.compile(
    r'^(?P<contig>[^:]+):(?P<pos>\d+)-(?P<asm_ref>[^-]+)-(?P<asm_alt>.+)$'
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True)
    p.add_argument("--hap1-bam", required=True)
    p.add_argument("--hap2-bam", required=True)
    p.add_argument("--xg", default=None)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--graph-flank", type=int, default=2)
    p.add_argument("--context-bp", type=int, default=20)
    return p.parse_args()


def leading_clipped_bases(cigartuples: Optional[List[Tuple[int, int]]]) -> int:
    if not cigartuples:
        return 0
    total = 0
    i = 0
    while i < len(cigartuples) and cigartuples[i][0] in (CIGAR_H, CIGAR_S):
        total += cigartuples[i][1]
        i += 1
    return total


def leading_hard_clipped_bases(cigartuples: Optional[List[Tuple[int, int]]]) -> int:
    if not cigartuples:
        return 0
    total = 0
    i = 0
    while i < len(cigartuples) and cigartuples[i][0] == CIGAR_H:
        total += cigartuples[i][1]
        i += 1
    return total


def aligned_query_span_1based(rec: pysam.AlignedSegment) -> Optional[Tuple[int, int]]:
    cigartuples = rec.cigartuples
    if cigartuples is None:
        return None
    qstart = leading_clipped_bases(cigartuples) + 1
    q_aligned_len = sum(length for op, length in cigartuples if op in (CIGAR_M, CIGAR_I, CIGAR_EQ, CIGAR_X))
    qend = qstart + q_aligned_len - 1
    return qstart, qend


def get_query_sequence_offset(rec: pysam.AlignedSegment) -> int:
    return leading_hard_clipped_bases(rec.cigartuples)


def get_query_context(
    rec: pysam.AlignedSegment,
    query_pos_1based: int,
    context_bp: int,
) -> Tuple[str, str, str]:
    seq = rec.query_sequence
    if seq is None:
        return ("", "", "")

    hard_offset = get_query_sequence_offset(rec)

    # keep the user's corrected convention: shift +1
    seq_pos_1based = query_pos_1based - hard_offset + 1

    if seq_pos_1based < 1 or seq_pos_1based > len(seq):
        return ("", "", "")

    idx0 = seq_pos_1based - 1
    prev_seq = seq[max(0, idx0 - context_bp): idx0]
    target_base = seq[idx0]
    next_seq = seq[idx0 + 1: idx0 + 1 + context_bp]
    return (prev_seq, target_base, next_seq)


def build_op_intervals(rec: pysam.AlignedSegment) -> List[Dict[str, Any]]:
    cigartuples = rec.cigartuples
    if cigartuples is None:
        return []

    ref_name = rec.reference_name
    ref_pos = rec.reference_start + 1
    q_pos = leading_clipped_bases(cigartuples) + 1

    intervals = []
    for op, length in cigartuples:
        op_char = CIGAR_CODE_TO_CHAR.get(op, "?")

        query_consuming = op in (CIGAR_M, CIGAR_I, CIGAR_S, CIGAR_EQ, CIGAR_X)
        ref_consuming = op in (CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_EQ, CIGAR_X)

        q_start = q_pos if query_consuming and op != CIGAR_S else None
        q_end = (q_pos + length - 1) if query_consuming and op != CIGAR_S else None

        r_start = ref_pos if ref_consuming else None
        r_end = (ref_pos + length - 1) if ref_consuming else None

        insertion_anchor_ref = None
        insertion_anchor_between = None
        if op == CIGAR_I:
            left_anchor = ref_pos - 1
            right_anchor = ref_pos
            insertion_anchor_ref = left_anchor
            insertion_anchor_between = (left_anchor, right_anchor)

        intervals.append({
            "op": op,
            "op_char": op_char,
            "length": length,
            "query_start": q_start,
            "query_end": q_end,
            "ref_start": r_start,
            "ref_end": r_end,
            "insertion_anchor_ref": insertion_anchor_ref,
            "insertion_anchor_between": insertion_anchor_between,
            "ref_name": ref_name,
        })

        if query_consuming:
            q_pos += length
        if ref_consuming:
            ref_pos += length

    return intervals


def find_target_interval(intervals: List[Dict[str, Any]], query_pos_1based: int) -> Optional[Dict[str, Any]]:
    for it in intervals:
        qs = it["query_start"]
        qe = it["query_end"]
        if qs is not None and qe is not None and qs <= query_pos_1based <= qe:
            return it
    return None


def query_pos_to_ref_pos(rec: pysam.AlignedSegment, query_pos_1based: int) -> Optional[Tuple[str, int]]:
    hard_offset = get_query_sequence_offset(rec)
    local_query_pos_1based = query_pos_1based - hard_offset + 1
    if local_query_pos_1based < 1:
        return None

    target_qpos0 = local_query_pos_1based - 1
    for qpos0, rpos0 in rec.get_aligned_pairs(matches_only=False):
        if qpos0 is not None and qpos0 == target_qpos0:
            if rpos0 is None:
                return None
            return rec.reference_name, rpos0 + 1
    return None


def parse_full_annotation(s: str) -> Dict[str, Any]:
    m = FULL_ANNOT_RE.match(s)
    if not m:
        return {
            "asm_contig": "",
            "asm_pos": "",
            "asm_ref": "",
            "asm_alt": "",
            "hap": "",
        }
    contig = m.group("contig")
    asm_pos = int(m.group("pos"))
    asm_ref = m.group("asm_ref")
    asm_alt = m.group("asm_alt")
    hap = "hap1" if "_hap1" in contig else ("hap2" if "_hap2" in contig else "")
    return {
        "asm_contig": contig,
        "asm_pos": asm_pos,
        "asm_ref": asm_ref,
        "asm_alt": asm_alt,
        "hap": hap,
    }


def graph_lookup_nodes(xg: str, chrom: str, pos: int, flank: int) -> str:
    start = max(1, pos - flank)
    end = pos + flank
    region = f"GRCh38#0#{chrom}:{start}-{end}"

    cmd = (
        f"vg find -x {shell_quote(xg)} -p {shell_quote(region)} "
        f"| vg view -j - "
        f"| jq -c '[.node[] | {{id, sequence}}]'"
    )
    try:
        out = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.DEVNULL).strip()
        return out if out else "[]"
    except subprocess.CalledProcessError:
        return "[]"


def shell_quote(s: str) -> str:
    return "'" + s.replace("'", "'\"'\"'") + "'"


def load_tsv(tsv_path: str) -> List[Dict[str, Any]]:
    rows = []
    with open(tsv_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            row = dict(row)
            row["_row_index"] = i
            parsed = parse_full_annotation(row.get("FULL_ANNOTATION", ""))
            row.update(parsed)
            rows.append(row)
    return rows


def group_targets_by_hap_and_contig(rows: List[Dict[str, Any]]) -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
    grouped: Dict[str, Dict[str, List[Dict[str, Any]]]] = {"hap1": defaultdict(list), "hap2": defaultdict(list)}
    for row in rows:
        hap = row.get("hap", "")
        contig = row.get("asm_contig", "")
        asm_pos = row.get("asm_pos", "")
        if hap in grouped and contig and isinstance(asm_pos, int):
            grouped[hap][contig].append(row)

    for hap in grouped:
        for contig in grouped[hap]:
            grouped[hap][contig].sort(key=lambda r: r["asm_pos"])
    return grouped


def scan_bam_for_targets(
    bam_path: str,
    grouped_targets: Dict[str, List[Dict[str, Any]]],
    context_bp: int,
) -> Dict[int, Dict[str, Any]]:
    results: Dict[int, Dict[str, Any]] = {}

    bam = pysam.AlignmentFile(bam_path, "rb")
    targets_by_contig = grouped_targets

    for rec in bam.fetch(until_eof=True):
        qname = rec.query_name
        if qname not in targets_by_contig:
            continue

        qspan = aligned_query_span_1based(rec)
        if qspan is None:
            continue
        qstart, qend = qspan

        target_list = targets_by_contig[qname]
        if not target_list:
            continue

        intervals = build_op_intervals(rec)

        for row in target_list:
            idx = row["_row_index"]
            if idx in results:
                continue

            qpos = row["asm_pos"]
            if not (qstart <= qpos <= qend):
                continue

            target_item = find_target_interval(intervals, qpos)
            if target_item is None:
                continue

            prev_seq, target_base, next_seq = get_query_context(rec, qpos, context_bp)
            mapped = query_pos_to_ref_pos(rec, qpos)

            result = {
                "bam_query_name": rec.query_name,
                "bam_flag": rec.flag,
                "bam_flag_label": flag_labels(rec.flag),
                "bam_reference": rec.reference_name,
                "bam_reference_start_1based": rec.reference_start + 1,
                "bam_mapping_quality": rec.mapping_quality,
                "query_span_1based": f"{qstart}-{qend}",
                "query_base": target_base,
                "query_context_left_20bp": prev_seq,
                "query_context_right_20bp": next_seq,
                "covering_cigar_fragment": f"{target_item['length']}{target_item['op_char']}",
                "covering_cigar_op": target_item["op_char"],
                "covering_qspan_start": target_item["query_start"],
                "covering_qspan_end": target_item["query_end"],
                "covering_rspan_start": target_item["ref_start"],
                "covering_rspan_end": target_item["ref_end"],
                "map_kind": "",
                "mapped_ref_chrom": "",
                "mapped_ref_pos": "",
                "insertion_left_anchor": "",
                "insertion_right_anchor": "",
            }

            if target_item["op"] == CIGAR_I:
                left_anchor, right_anchor = target_item["insertion_anchor_between"]
                result["map_kind"] = "insertion_anchor"
                result["mapped_ref_chrom"] = rec.reference_name
                result["mapped_ref_pos"] = left_anchor
                result["insertion_left_anchor"] = left_anchor
                result["insertion_right_anchor"] = right_anchor
            elif mapped is None:
                result["map_kind"] = "unaligned_or_insertion"
            else:
                result["map_kind"] = "exact"
                result["mapped_ref_chrom"] = mapped[0]
                result["mapped_ref_pos"] = mapped[1]

            results[idx] = result

    bam.close()
    return results


def flag_labels(flag: int) -> str:
    labels = []
    if flag & 0x100:
        labels.append("SECONDARY")
    if flag & 0x800:
        labels.append("SUPPLEMENTARY")
    if not labels:
        labels.append("PRIMARY")
    return ",".join(labels)


def safe_int(x: Any) -> Optional[int]:
    try:
        return int(x)
    except Exception:
        return None


def build_output_rows(
    rows: List[Dict[str, Any]],
    hap1_results: Dict[int, Dict[str, Any]],
    hap2_results: Dict[int, Dict[str, Any]],
    xg: Optional[str],
    graph_flank: int,
) -> List[Dict[str, Any]]:
    out_rows = []

    for row in rows:
        idx = row["_row_index"]
        hap = row.get("hap", "")
        bam_result = hap1_results.get(idx, {}) if hap == "hap1" else hap2_results.get(idx, {})

        out = dict(row)
        out.update(bam_result)

        tsv_chrom = row.get("#CHROM", "")
        tsv_pos = safe_int(row.get("POS", ""))

        mapped_chrom = out.get("mapped_ref_chrom", "")
        mapped_pos = safe_int(out.get("mapped_ref_pos", ""))

        out["same_chrom_as_tsv"] = ""
        out["pos_delta_vs_tsv"] = ""
        if mapped_chrom and tsv_chrom:
            out["same_chrom_as_tsv"] = str(mapped_chrom == tsv_chrom)
        if mapped_pos is not None and tsv_pos is not None:
            out["pos_delta_vs_tsv"] = mapped_pos - tsv_pos

        if xg and tsv_chrom and tsv_pos is not None:
            out["graph_nodes_json"] = graph_lookup_nodes(xg, tsv_chrom, tsv_pos, graph_flank)
        else:
            out["graph_nodes_json"] = ""

        out_rows.append(out)

    return out_rows


def write_tsv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        with open(path, "w") as f:
            f.write("")
        return

    fieldnames = []
    seen = set()
    for row in rows:
        for k in row.keys():
            if k not in seen and k != "_row_index":
                seen.add(k)
                fieldnames.append(k)

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            clean = {k: v for k, v in row.items() if k != "_row_index"}
            writer.writerow(clean)


def write_summary(path: str, rows: List[Dict[str, Any]]) -> None:
    c_map_kind = Counter()
    c_hap = Counter()
    c_type_change = Counter()
    same_chrom_true = 0
    same_chrom_false = 0
    exact_match_pos = 0

    for r in rows:
        c_map_kind[r.get("map_kind", "missing")] += 1
        c_hap[r.get("hap", "unknown")] += 1
        c_type_change[r.get("TYPE_CHANGE", "")] += 1

        sc = r.get("same_chrom_as_tsv", "")
        if sc == "True":
            same_chrom_true += 1
        elif sc == "False":
            same_chrom_false += 1

        delta = r.get("pos_delta_vs_tsv", "")
        if str(delta) == "0":
            exact_match_pos += 1

    summary = {
        "n_rows": len(rows),
        "hap_counts": dict(c_hap),
        "map_kind_counts": dict(c_map_kind),
        "type_change_counts": dict(c_type_change),
        "same_chrom_true": same_chrom_true,
        "same_chrom_false": same_chrom_false,
        "exact_position_match_count": exact_match_pos,
    }

    with open(path, "w") as f:
        json.dump(summary, f, indent=2)


def main() -> None:
    args = parse_args()

    rows = load_tsv(args.tsv)
    grouped = group_targets_by_hap_and_contig(rows)

    hap1_results = scan_bam_for_targets(
        bam_path=args.hap1_bam,
        grouped_targets=grouped["hap1"],
        context_bp=args.context_bp,
    )
    hap2_results = scan_bam_for_targets(
        bam_path=args.hap2_bam,
        grouped_targets=grouped["hap2"],
        context_bp=args.context_bp,
    )

    out_rows = build_output_rows(
        rows=rows,
        hap1_results=hap1_results,
        hap2_results=hap2_results,
        xg=args.xg,
        graph_flank=args.graph_flank,
    )

    out_tsv = args.out_prefix + ".full.tsv"
    out_json = args.out_prefix + ".summary.json"

    write_tsv(out_tsv, out_rows)
    write_summary(out_json, out_rows)

    print(f"[DONE] Wrote: {out_tsv}")
    print(f"[DONE] Wrote: {out_json}")


if __name__ == "__main__":
    main()