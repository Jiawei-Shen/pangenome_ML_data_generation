#!/usr/bin/env python3
"""
Validate/project changed-representation variants against hap1/hap2 BAMs and optionally
collect graph nodes around the GRCh38 position.

What it does for each TSV row:
1. Parse FULL_ANNOTATION, e.g. chr1_hap1:1195496-C-CA
2. Use the appropriate BAM (hap1 or hap2) to map the assembly query position to GRCh38
   - returns exact mapped base position when available
   - returns insertion anchor when the query base falls in an insertion/query-only segment
3. Compare the mapped/anchored reference position to the TSV #CHROM/POS
4. Optionally query graph nodes around the TSV GRCh38 position with vg find

Designed to be efficient:
- BAM is scanned only once per hap file
- assembly positions are grouped by query contig
- aligned-pair maps are built only for alignments that overlap requested positions
- graph queries are deduplicated by unique (chrom, pos, flank)

Example:
python check_changed_repr_variants.py \
  --tsv /mnt/data/changed_representation_variants.tsv \
  --hap1-bam /scratch/jshen/tmp/rebuild_hprc_d9_plus_HG008N_curatedv6_250714_polished6.2_vg1660_268142/HG008N_curatedv6_250714_polished6.2.hap1.bam \
  --hap2-bam /scratch/jshen/tmp/rebuild_hprc_d9_plus_HG008N_curatedv6_250714_polished6.2_vg1660_268142/HG008N_curatedv6_250714_polished6.2.hap2.bam \
  --xg /scratch/jshen/data/HG008_GIAB/indexes_HG008N_curatedv6_250714_polished6.2_rebuild_hprc_d9_vg1660/hprc-v1.1-mc-grch38.d9.plus_HG008N_curatedv6_250714_polished6_2.xg \
  --out-prefix /scratch/jshen/tmp/changed_repr_check \
  --graph-flank 2 \
  --workers 8
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import re
import subprocess
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd
import pysam

# CIGAR op codes
CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_EQ = 7
CIGAR_X = 8

FULL_ANNOT_RE = re.compile(r"^(?P<contig>[^:]+):(?P<pos>\d+)-(?P<asm_ref>[^-]+)-(?P<asm_alt>.+)$")


@dataclass(frozen=True)
class QuerySite:
    row_idx: int
    qname: str
    qpos: int
    hap: str
    chrom: str
    pos: int
    ref: str
    alt: str
    full_annotation: str


@dataclass
class SiteResult:
    row_idx: int
    qname: str
    qpos: int
    hap: str
    graph_chrom: str
    graph_pos: int
    mapped_ref_chrom: Optional[str]
    mapped_ref_pos: Optional[int]
    map_kind: str  # exact_match / insertion_anchor / not_found / off_target
    insertion_left_anchor: Optional[int]
    insertion_right_anchor: Optional[int]
    query_base: Optional[str]
    query_context_left_20bp: Optional[str]
    query_context_right_20bp: Optional[str]
    covering_cigar_fragment: Optional[str]
    covering_cigar_op: Optional[str]
    covering_qspan_start: Optional[int]
    covering_qspan_end: Optional[int]
    covering_rspan_start: Optional[int]
    covering_rspan_end: Optional[int]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True)
    p.add_argument("--hap1-bam", required=True)
    p.add_argument("--hap2-bam", required=True)
    p.add_argument("--xg", default=None, help="Optional .xg for graph node lookup")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--graph-flank", type=int, default=0)
    p.add_argument("--context-bp", type=int, default=20)
    p.add_argument("--workers", type=int, default=4)
    return p.parse_args()


def parse_full_annotation(value: str) -> Tuple[str, int, str, str, str]:
    m = FULL_ANNOT_RE.match(value)
    if not m:
        raise ValueError(f"Could not parse FULL_ANNOTATION: {value}")
    qname = m.group("contig")
    qpos = int(m.group("pos"))
    asm_ref = m.group("asm_ref")
    asm_alt = m.group("asm_alt")
    if qname.endswith("_hap1"):
        hap = "hap1"
    elif qname.endswith("_hap2"):
        hap = "hap2"
    else:
        raise ValueError(f"FULL_ANNOTATION contig does not end with _hap1/_hap2: {value}")
    return qname, qpos, asm_ref, asm_alt, hap


def load_sites(tsv_path: str) -> Tuple[pd.DataFrame, List[QuerySite]]:
    df = pd.read_csv(tsv_path, sep="\t", dtype={"#CHROM": str, "REF": str, "ALT": str, "FULL_ANNOTATION": str})
    sites: List[QuerySite] = []
    for idx, row in df.iterrows():
        qname, qpos, _asm_ref, _asm_alt, hap = parse_full_annotation(row["FULL_ANNOTATION"])
        sites.append(
            QuerySite(
                row_idx=int(idx),
                qname=qname,
                qpos=qpos,
                hap=hap,
                chrom=str(row["#CHROM"]),
                pos=int(row["POS"]),
                ref=str(row["REF"]),
                alt=str(row["ALT"]),
                full_annotation=str(row["FULL_ANNOTATION"]),
            )
        )
    return df, sites


def leading_hard_clipped_bases(cigartuples: Optional[List[Tuple[int, int]]]) -> int:
    if not cigartuples:
        return 0
    i = 0
    total = 0
    while i < len(cigartuples) and cigartuples[i][0] == CIGAR_H:
        total += cigartuples[i][1]
        i += 1
    return total


def leading_total_clipped_bases(cigartuples: Optional[List[Tuple[int, int]]]) -> int:
    if not cigartuples:
        return 0
    i = 0
    total = 0
    while i < len(cigartuples) and cigartuples[i][0] in (CIGAR_H, CIGAR_S):
        total += cigartuples[i][1]
        i += 1
    return total


def aligned_query_span(rec: pysam.AlignedSegment) -> Optional[Tuple[int, int]]:
    if rec.cigartuples is None:
        return None
    qstart = leading_total_clipped_bases(rec.cigartuples) + 1
    qlen = sum(length for op, length in rec.cigartuples if op in (CIGAR_M, CIGAR_I, CIGAR_EQ, CIGAR_X))
    if qlen <= 0:
        return None
    return qstart, qstart + qlen - 1


def build_op_intervals(rec: pysam.AlignedSegment) -> List[Dict[str, Any]]:
    cigartuples = rec.cigartuples
    if cigartuples is None:
        return []

    ref_name = rec.reference_name
    ref_pos = rec.reference_start + 1  # 1-based
    q_pos = leading_total_clipped_bases(cigartuples) + 1
    items: List[Dict[str, Any]] = []

    for op, length in cigartuples:
        q_consume = op in (CIGAR_M, CIGAR_I, CIGAR_S, CIGAR_EQ, CIGAR_X)
        r_consume = op in (CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_EQ, CIGAR_X)
        op_char = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 7: "=", 8: "X"}.get(op, "?")

        q_start = q_pos if q_consume and op != CIGAR_S else None
        q_end = q_pos + length - 1 if q_consume and op != CIGAR_S else None
        r_start = ref_pos if r_consume else None
        r_end = ref_pos + length - 1 if r_consume else None

        left_anchor = right_anchor = None
        if op == CIGAR_I:
            left_anchor = ref_pos - 1
            right_anchor = ref_pos

        items.append(
            {
                "op": op,
                "op_char": op_char,
                "length": length,
                "query_start": q_start,
                "query_end": q_end,
                "ref_start": r_start,
                "ref_end": r_end,
                "ref_name": ref_name,
                "left_anchor": left_anchor,
                "right_anchor": right_anchor,
            }
        )

        if q_consume:
            q_pos += length
        if r_consume:
            ref_pos += length

    return items


def find_covering_interval(items: List[Dict[str, Any]], qpos: int) -> Optional[Dict[str, Any]]:
    for item in items:
        qs = item["query_start"]
        qe = item["query_end"]
        if qs is not None and qe is not None and qs <= qpos <= qe:
            return item
    return None


def query_context_from_record(rec: pysam.AlignedSegment, qpos: int, context_bp: int) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    seq = rec.query_sequence
    if seq is None:
        return None, None, None
    hard = leading_hard_clipped_bases(rec.cigartuples)
    seq_pos = qpos - hard + 1  # user requested +1 correction from prior debugging
    if seq_pos < 1 or seq_pos > len(seq):
        return None, None, None
    idx0 = seq_pos - 1
    left = seq[max(0, idx0 - context_bp): idx0]
    base = seq[idx0]
    right = seq[idx0 + 1: idx0 + 1 + context_bp]
    return left, base, right


def map_query_pos_in_record(rec: pysam.AlignedSegment, qpos: int) -> Tuple[Optional[str], Optional[int], str]:
    """
    Map global 1-based query coordinate to global 1-based reference coordinate.
    Returns (chrom, pos, kind), where kind is exact_match or query_only.
    Uses the corrected +1 query-sequence convention from debugging.
    """
    hard = leading_hard_clipped_bases(rec.cigartuples)
    local_qpos1 = qpos - hard + 1
    if local_qpos1 < 1:
        return None, None, "not_in_record"
    target_q0 = local_qpos1 - 1

    for q0, r0 in rec.get_aligned_pairs(matches_only=False):
        if q0 is not None and q0 == target_q0:
            if r0 is None:
                return None, None, "query_only"
            return rec.reference_name, r0 + 1, "exact_match"
    return None, None, "not_in_record"


def process_hap_bam(bam_path: str, sites: List[QuerySite], context_bp: int) -> List[SiteResult]:
    needed_by_qname: Dict[str, List[QuerySite]] = defaultdict(list)
    for s in sites:
        needed_by_qname[s.qname].append(s)
    for qname in needed_by_qname:
        needed_by_qname[qname].sort(key=lambda x: x.qpos)

    results: Dict[int, SiteResult] = {}

    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        for rec in bam.fetch(until_eof=True):
            qname = rec.query_name
            if qname not in needed_by_qname:
                continue
            span = aligned_query_span(rec)
            if span is None:
                continue
            qstart, qend = span

            group = needed_by_qname[qname]
            q_positions = [s.qpos for s in group]
            lo = bisect_left(q_positions, qstart)
            i = lo
            if i >= len(group):
                continue
            if group[i].qpos > qend:
                continue

            intervals = build_op_intervals(rec)
            while i < len(group) and group[i].qpos <= qend:
                site = group[i]
                if site.row_idx in results:
                    i += 1
                    continue

                item = find_covering_interval(intervals, site.qpos)
                if item is None:
                    i += 1
                    continue

                left_ctx, base, right_ctx = query_context_from_record(rec, site.qpos, context_bp)
                mapped_chrom, mapped_pos, kind = map_query_pos_in_record(rec, site.qpos)
                left_anchor = item.get("left_anchor")
                right_anchor = item.get("right_anchor")

                if item["op"] == CIGAR_I:
                    map_kind = "insertion_anchor"
                    mapped_chrom = rec.reference_name
                    mapped_pos = left_anchor
                elif kind == "exact_match":
                    map_kind = "exact_match"
                elif kind == "query_only":
                    map_kind = "query_only_nonI"
                else:
                    map_kind = "not_found"

                results[site.row_idx] = SiteResult(
                    row_idx=site.row_idx,
                    qname=site.qname,
                    qpos=site.qpos,
                    hap=site.hap,
                    graph_chrom=site.chrom,
                    graph_pos=site.pos,
                    mapped_ref_chrom=mapped_chrom,
                    mapped_ref_pos=mapped_pos,
                    map_kind=map_kind,
                    insertion_left_anchor=left_anchor,
                    insertion_right_anchor=right_anchor,
                    query_base=base,
                    query_context_left_20bp=left_ctx,
                    query_context_right_20bp=right_ctx,
                    covering_cigar_fragment=f"{item['length']}{item['op_char']}",
                    covering_cigar_op=item["op_char"],
                    covering_qspan_start=item["query_start"],
                    covering_qspan_end=item["query_end"],
                    covering_rspan_start=item["ref_start"],
                    covering_rspan_end=item["ref_end"],
                )
                i += 1
    finally:
        bam.close()

    out: List[SiteResult] = []
    for s in sites:
        if s.row_idx in results:
            out.append(results[s.row_idx])
        else:
            out.append(
                SiteResult(
                    row_idx=s.row_idx,
                    qname=s.qname,
                    qpos=s.qpos,
                    hap=s.hap,
                    graph_chrom=s.chrom,
                    graph_pos=s.pos,
                    mapped_ref_chrom=None,
                    mapped_ref_pos=None,
                    map_kind="not_found",
                    insertion_left_anchor=None,
                    insertion_right_anchor=None,
                    query_base=None,
                    query_context_left_20bp=None,
                    query_context_right_20bp=None,
                    covering_cigar_fragment=None,
                    covering_cigar_op=None,
                    covering_qspan_start=None,
                    covering_qspan_end=None,
                    covering_rspan_start=None,
                    covering_rspan_end=None,
                )
            )
    return out


def vg_find_one(xg: str, chrom: str, pos: int, flank: int) -> Tuple[str, int, int, List[Dict[str, Any]], Optional[str]]:
    start = max(1, pos - flank)
    end = pos + flank
    region = f"GRCh38#0#{chrom}:{start}-{end}"
    cmd = f"vg find -x {subprocess.list2cmdline([xg])} -p '{region}' | vg view -j -"
    try:
        proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        obj = json.loads(proc.stdout)
        nodes = [{"id": n.get("id"), "sequence": n.get("sequence")} for n in obj.get("node", [])]
        return chrom, pos, flank, nodes, None
    except subprocess.CalledProcessError as e:
        return chrom, pos, flank, [], e.stderr.strip() or str(e)
    except json.JSONDecodeError as e:
        return chrom, pos, flank, [], f"JSON decode failed: {e}"


def collect_graph_nodes(xg: str, sites: Iterable[QuerySite], flank: int, workers: int) -> Dict[Tuple[str, int, int], Dict[str, Any]]:
    keys = sorted({(s.chrom, s.pos, flank) for s in sites})
    out: Dict[Tuple[str, int, int], Dict[str, Any]] = {}
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(vg_find_one, xg, chrom, pos, flank) for chrom, pos, flank in keys]
        for fut in cf.as_completed(futs):
            chrom, pos, flank, nodes, err = fut.result()
            out[(chrom, pos, flank)] = {
                "graph_nodes_json": json.dumps(nodes, separators=(",", ":")),
                "graph_node_count": len(nodes),
                "graph_query_error": err,
            }
    return out


def main() -> None:
    args = parse_args()
    out_prefix = Path(args.out_prefix)

    df, sites = load_sites(args.tsv)
    hap1_sites = [s for s in sites if s.hap == "hap1"]
    hap2_sites = [s for s in sites if s.hap == "hap2"]

    hap1_results = process_hap_bam(args.hap1_bam, hap1_sites, args.context_bp)
    hap2_results = process_hap_bam(args.hap2_bam, hap2_sites, args.context_bp)
    all_results = hap1_results + hap2_results

    res_df = pd.DataFrame([r.__dict__ for r in all_results]).sort_values("row_idx")
    merged = df.copy()
    merged["row_idx"] = range(len(merged))
    merged = merged.merge(res_df, on="row_idx", how="left", validate="one_to_one")

    merged["pos_delta_vs_tsv"] = merged.apply(
        lambda r: (int(r["mapped_ref_pos"]) - int(r["POS"])) if pd.notna(r["mapped_ref_pos"]) else pd.NA,
        axis=1,
    )
    merged["same_chrom_as_tsv"] = merged.apply(
        lambda r: (str(r["mapped_ref_chrom"]) == str(r["#CHROM"])) if pd.notna(r["mapped_ref_chrom"]) else pd.NA,
        axis=1,
    )

    if args.xg:
        graph_map = collect_graph_nodes(args.xg, sites, args.graph_flank, args.workers)
        merged["graph_nodes_json"] = [graph_map[(s.chrom, s.pos, args.graph_flank)]["graph_nodes_json"] for s in sites]
        merged["graph_node_count"] = [graph_map[(s.chrom, s.pos, args.graph_flank)]["graph_node_count"] for s in sites]
        merged["graph_query_error"] = [graph_map[(s.chrom, s.pos, args.graph_flank)]["graph_query_error"] for s in sites]

    merged.to_csv(f"{out_prefix}.full.tsv", sep="\t", index=False)

    summary = {
        "n_rows": int(len(merged)),
        "hap1_rows": int(len(hap1_sites)),
        "hap2_rows": int(len(hap2_sites)),
        "map_kind_counts": merged["map_kind"].value_counts(dropna=False).to_dict(),
        "same_chrom_true": int((merged["same_chrom_as_tsv"] == True).sum()),
        "same_chrom_false": int((merged["same_chrom_as_tsv"] == False).sum()),
        "exact_pos_match": int((merged["pos_delta_vs_tsv"] == 0).sum()),
        "pos_delta_min": None if merged["pos_delta_vs_tsv"].dropna().empty else int(merged["pos_delta_vs_tsv"].dropna().min()),
        "pos_delta_max": None if merged["pos_delta_vs_tsv"].dropna().empty else int(merged["pos_delta_vs_tsv"].dropna().max()),
    }
    with open(f"{out_prefix}.summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"[DONE] wrote {out_prefix}.full.tsv")
    print(f"[DONE] wrote {out_prefix}.summary.json")


if __name__ == "__main__":
    main()
