#!/usr/bin/env python3
"""
extract_core_snps_flex.py
-----------------------------------------------
For every record in a vg-deconstructed VCF:

  • keep alleles where len(ALT) == len(REF)  (substitution or MNP)
  • locate every single-base difference between REF and that ALT
  • report the REF base + its node and the ALT base + *the node
    that is unique to the ALT traversal* (if several, choose the
    first one that differs from the REF node at that position).

Output (TAB-separated):
  CHROM  SNP_POS  ID  REF_BASE  REF_NODE  ALT_BASE  ALT_NODE
"""

import argparse, re, sys, pysam
from typing import List

_SPLIT = re.compile(r"[><]+")

def nodes(trav: str) -> List[int]:
    """'>1745>1743>1742' → [1745, 1743, 1742]"""
    return [int(tok) for tok in _SPLIT.split(trav) if tok]

def snps_for_alt(
    ref_seq: str,
    ref_nodes: List[int],
    alt_seq: str,
    alt_nodes: List[int],
):
    """
    Yield tuples   (offset, ref_base, ref_node, alt_base, alt_node)

    • Works even when node lists are shorter/longer than the base string.
    • For k single-base differences, uses the first k node IDs that appear
      *only* on the REF side and the first k node IDs that appear *only*
      on the ALT side.
    """
    if len(ref_seq) != len(alt_seq):          # skip indels
        return

    # 1. Which positions differ at the nucleotide level?
    diff_idx = [i for i, (rb, ab) in enumerate(zip(ref_seq, alt_seq)) if rb != ab]
    if not diff_idx:
        return

    # 2. Node IDs that are unique to each traversal
    ref_only = [n for n in ref_nodes if n not in alt_nodes]
    alt_only = [n for n in alt_nodes if n not in ref_nodes]

    # Pad with last node if we ran out (rare, but keeps script robust)
    while len(ref_only) < len(diff_idx):
        ref_only.append(ref_nodes[-1])
    while len(alt_only) < len(diff_idx):
        alt_only.append(alt_nodes[-1])

    # 3. Emit one row per differing base, pairing node IDs in order
    for k, i in enumerate(diff_idx):
        ref_base = ref_seq[i]
        alt_base = alt_seq[i]
        yield (
            i,                # offset within the record
            ref_base,
            ref_only[k],      # kth node unique to REF
            alt_base,
            alt_only[k],      # kth node unique to ALT
        )

def main(vcf_path):
    vcf = pysam.VariantFile(vcf_path)
    if "AT" not in vcf.header.info:
        sys.exit("ERROR: VCF header lacks AT tag")

    print("CHROM\tSNP_POS\tID\tREF_BASE\tREF_NODE\tALT_BASE\tALT_NODE")
    for rec in vcf:
        ref_seq, ref_trav = rec.ref, rec.info["AT"][0]
        ref_nodes = nodes(ref_trav)

        for alt_idx, alt_seq in enumerate(rec.alts):
            alt_nodes = nodes(rec.info["AT"][alt_idx+1])
            for off, rb, rn, ab, an in snps_for_alt(ref_seq, ref_nodes,
                                                   alt_seq, alt_nodes):
                print(f"{rec.contig}\t{rec.pos+off}\t{rec.id}\t"
                      f"{rb}\t{rn}\t{ab}\t{an}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Extract SNP node IDs without assuming 1-bp-per-node."
    )
    ap.add_argument("vcf", help="vg-deconstructed VCF")
    main(ap.parse_args().vcf)
