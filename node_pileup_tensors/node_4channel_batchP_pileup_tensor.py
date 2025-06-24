#!/usr/bin/env python3
"""
variant_tensor_generator.py

Generate 4‑channel (bases, base‑quality, mismatch flag, mapping quality) tensors
centred on candidate variants extracted from a custom *.dat* alignment file.

Change log (Jun 2025)
─────────────────────
* **Added allele‑site base‑quality filter**
  ‑ New CLI arg `--min_allele_bq` (float, default 15.0).
  ‑ For every candidate variant we now collect **base‑qualities of the ALT‑supporting
    bases at the variant position only** (`alt_bqs_at_site`).  The variant is kept
    only if the **mean** of those values is ≥ `--min_allele_bq`.
* `get_allele_from_read_at_node_pos()` now also returns the **Phred quality** of
  the allele (or the mean quality for an insertion).  Call sites updated.
"""

import argparse
import struct
import json
import os
import sys
import time
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import torch

# ────────────────────────────────────────────────────────────────────────────
# Constants ------------------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE   = RECORD_STRUCT.size

BASE_TO_INDEX = {
    "A": 2, "C": 3, "G": 5, "T": 7,
    "N": 1, "*": 9, "_PADDING_": 0
}
PADDING_BASE_INDEX = 0

INDEX_TO_BASE_FOR_VIEW = {
    2: "A", 3: "C", 5: "G", 7: "T", 1: "N", 9: "*", 0: " "
}

TENSOR_WINDOW_SIZE      = 100
TENSOR_MAX_READ_ROWS    = 200
DEFAULT_QUALITY_PADDING = 0
DEFAULT_MAPPING_QUALITY_PADDING = -1
MISMATCH_CHANNEL_REF_ROW_VALUE  = 0
MISMATCH_COMPARISON_PADDING_VALUE = -1

# Globals initialised in each worker
worker_dat_file = None
worker_base_output_dir = None

# ────────────────────────────────────────────────────────────────────────────
# Helper functions -----------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(table)[::-1]


def decode_cigar_to_int_ops(cigar: str):
    if not cigar or cigar == "*":
        return []
    try:
        return [(int(l), op) for l, op in re.findall(r"(\d+)([MIDNSHPX=])", cigar)]
    except Exception as e:
        sys.stderr.write(f"Warning: bad CIGAR '{cigar}': {e}\n")
        return []


# --- NEW: allele + base‑quality extractor -----------------------------------

def get_allele_from_read_at_node_pos(read_offset_on_node: int,
                                     read_seq: str,
                                     read_qual: str,
                                     cigar_ops,
                                     target_node_pos: int,
                                     node_seq: str,
                                     expected_var_type=None,
                                     expected_ref_allele_for_indel=None):
    """Return (**allele**, **Phred quality**) observed in *this read* at
    *target_node_pos*.
    For SNV: quality of that single base.
    For insertion: average quality of inserted segment.
    For deletion: returns ("*", None) when the deletion matches; quality N/A.
    Returns (None, None) if the read does not cover the site.
    """
    node_pos = read_offset_on_node
    read_pos = 0

    for length, op in cigar_ops:
        # Match / mismatch ----------------------------------------------------
        if op in ("M", "=", "X"):
            if node_pos <= target_node_pos < node_pos + length:
                if expected_var_type in ("I", "D"):
                    return "REF_STATE_FOR_INDEL", None
                offset = target_node_pos - node_pos
                idx    = read_pos + offset
                if idx < len(read_seq):
                    base = read_seq[idx].upper()
                    bq   = ord(read_qual[idx]) - 33 if idx < len(read_qual) else 0
                    return base, bq
                return None, None
            node_pos += length
            read_pos += length

        # Insertion -----------------------------------------------------------
        elif op == "I":
            if expected_var_type == "I" and (node_pos - 1) == target_node_pos:
                sub_seq  = read_seq[read_pos: read_pos + length]
                sub_qual = read_qual[read_pos: read_pos + length]
                if sub_seq:
                    mean_bq = sum(ord(q) - 33 for q in sub_qual) / len(sub_qual)
                    return sub_seq.upper(), mean_bq
                return None, None
            read_pos += length

        # Deletion ------------------------------------------------------------
        elif op == "D":
            if node_pos <= target_node_pos < node_pos + length:
                if expected_var_type == "I":
                    return "OTHER_FOR_INDEL", None
                if expected_var_type == "D":
                    del_seq = node_seq[node_pos: node_pos + length]
                    if del_seq == expected_ref_allele_for_indel:
                        return "*", None
                    return "OTHER_FOR_INDEL", None
                return "*", None
            node_pos += length

        # Soft‑clip -----------------------------------------------------------
        elif op == "S":
            read_pos += length

        # Skipped region ------------------------------------------------------
        elif op == "N":
            node_pos += length

        # Early exit ----------------------------------------------------------
        if node_pos > target_node_pos + 1 and not (
            expected_var_type == "I" and (node_pos - 1) <= target_node_pos):
            break

    return None, None


# --- Detect candidate variants (unchanged) ----------------------------------

def detect_variants_from_cigar(offset_on_node, cigar_ops, read_seq, node_seq):
    variants = []
    node_p, read_p = offset_on_node, 0
    for L, op in cigar_ops:
        if op in ("M", "=", "X"):
            for i in range(L):
                n, r = node_p + i, read_p + i
                if n < len(node_seq) and r < len(read_seq):
                    if op != "=" and node_seq[n].upper() != read_seq[r].upper():
                        variants.append((n, "X", read_seq[r].upper(), node_seq[n].upper()))
            node_p += L; read_p += L
        elif op == "I":
            ins_seq = read_seq[read_p: read_p + L].upper()
            anchor  = node_p - 1 if node_p > 0 else 0
            ref_base= node_seq[anchor].upper() if anchor < len(node_seq) else "*"
            variants.append((anchor, "I", ins_seq, ref_base))
            read_p += L
        elif op == "D":
            del_seq = node_seq[node_p: node_p + L].upper()
            if del_seq:
                variants.append((node_p, "D", "*", del_seq))
            node_p += L
        elif op == "S":
            read_p += L
        elif op == "N":
            node_p += L
    return variants

# (Other helper functions for window visualisation and tensor‑row assembly are
# unchanged and omitted here for brevity – retain them from the previous
# version.)
# ---------------------------------------------------------------------------

# ────────────────────────────────────────────────────────────────────────────
# Worker initialisation ------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────

def init_worker(dat_path, base_out):
    global worker_dat_file, worker_base_output_dir
    worker_dat_file         = open(dat_path, "rb")
    worker_base_output_dir  = base_out

# ────────────────────────────────────────────────────────────────────────────
# Worker: process a single node ---------------------------------------------
# ────────────────────────────────────────────────────────────────────────────

def process_single_node_for_pileup(task):
    (node_id, offset, n_records, node_seq,
     min_af, min_alt_count, min_bq) = task

    global worker_dat_file, worker_base_output_dir
    node_len = len(node_seq)
    aligned  = []

    # --- read all segments for this node -----------------------------------
    worker_dat_file.seek(offset + 10)  # skip 10‑byte per‑node header
    for _ in range(n_records):
        rec = worker_dat_file.read(RECORD_SIZE)
        if len(rec) < RECORD_SIZE:
            break
        off, raw_seq, raw_qual, raw_cigar, mapq, strand = RECORD_STRUCT.unpack(rec)
        if mapq < 10:
            continue
        seq   = raw_seq.rstrip(b"\0").decode("ascii", "replace")
        qual  = raw_qual.rstrip(b"\0").decode("ascii", "replace")
        cigar = raw_cigar.rstrip(b"\0").decode("ascii", "replace")
        strand= strand.decode()

        if not seq or len(seq) != len(qual):
            continue
        cigar_ops = decode_cigar_to_int_ops(cigar)
        if not cigar_ops and cigar != "*":
            continue
        if strand == "-":
            seq   = reverse_complement(seq)
            qual  = qual[::-1]
            cigar_ops = list(reversed(cigar_ops))
            span  = sum(l for l, op in cigar_ops if op in "MDN=X")
            off   = node_len - off - span
        aligned.append({
            "off": off,
            "seq": seq,
            "qual": qual,
            "cigar": cigar_ops,
            "mapq": mapq,
            "strand": strand
        })

    if not aligned:
        return node_id, {}, 0

    # --- discover candidate variants ---------------------------------------
    cand = defaultdict(lambda: {"alt_count":0, "ref_count":0, "other":0,
                                "cov":0, "alt_bqs":[]})

    for seg in aligned:
        for pos, vtype, alt, ref in detect_variants_from_cigar(seg["off"], seg["cigar"], seg["seq"], node_seq):
            key = (pos, vtype, ref, alt)
            cand.setdefault(key)  # ensure key exists

    # --- gather evidence ----------------------------------------------------
    for key in cand:
        pos, vtype, ref, alt = key
        for seg in aligned:
            allele, bq = get_allele_from_read_at_node_pos(
                seg["off"], seg["seq"], seg["qual"], seg["cigar"],
                pos, node_seq, vtype,
                ref if vtype == "D" else None)
            if allele is None:
                continue
            c = cand[key]
            c["cov"] += 1
            if allele == alt:
                c["alt_count"] += 1
                if bq is not None:
                    c["alt_bqs"].append(bq)
            elif allele == ref or (vtype in ("I", "D") and allele == "REF_STATE_FOR_INDEL"):
                c["ref_count"] += 1
            else:
                c["other"] += 1

    # --- filter + tensor generation ----------------------------------------
    variant_headers = []
    tensors_out    = 0
    out_dir        = os.path.join(worker_base_output_dir, str(node_id))
    os.makedirs(out_dir, exist_ok=True)

    for key, st in cand.items():
        pos, vtype, ref, alt = key
        if st["alt_count"] <= min_alt_count:
            continue
        af = st["alt_count"] / st["cov"] if st["cov"] else 0.0
        if af < min_af:
            continue
        mean_bq = sum(st["alt_bqs"]) / len(st["alt_bqs"]) if st["alt_bqs"] else 0.0
        if mean_bq < min_bq:
            continue
        # (Tensor construction retained from previous version – omitted for brevity)
        # Save minimal header only (you can extend as before):
        variant_headers.append({
            "key": f"{pos}_{vtype}_{ref}_{alt}",
            "alt_count": st["alt_count"],
            "cov": st["cov"],
            "af": round(af,4),
            "mean_alt_bq": round(mean_bq,2)
        })
        tensors_out += 1

    # write summary ---------------------------------------------------------
    if variant_headers:
        with open(os.path.join(out_dir, "variant_summary.json"), "w") as fh:
            json.dump({"node_id": node_id, "variants": variant_headers}, fh, indent=2)

    return node_id, {}, tensors_out

# ────────────────────────────────────────────────────────────────────────────
# Main ----------------------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser("Generate variant‑centred tensors from .dat alignments")
    p.add_argument("dat")
    p.add_argument("idx")
    p.add_argument("output")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--node_id", type=int)
    g.add_argument("--node_id_file")

    p.add_argument("--gfa")
    p.add_argument("--min_af", type=float, default=0.1)
    p.add_argument("--min_variants", type=int, default=2,
                   help="Minimum ALT reads to consider variant")
    p.add_argument("--min_allele_bq", type=float, default=15.0,
                   help="Mean Phred of ALT‑supporting bases at site")
    p.add_argument("--num_workers", type=int, default=os.cpu_count())
    args = p.parse_args()

    # (Index loading, GFA sequence loading, task creation identical to previous
    # version – obtain offset & n_records from idx, sequences from gfa/cache.)
    # For brevity, that scaffolding is omitted here – integrate from the
    # previous full script, ensuring each task tuple now carries args.min_allele_bq.

if __name__ == "__main__":
    main()
