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

def snps_for_alt(ref_seq, ref_nodes, alt_seq, alt_nodes):
    """Yield (offset, ref_base, ref_node, alt_base, alt_node) tuples."""
    if len(ref_seq) != len(alt_seq):           # skip indels
        return
    # positions where the bases differ
    diff_idx = [i for i,(rb,ab) in enumerate(zip(ref_seq,alt_seq)) if rb!=ab]
    if not diff_idx:                           # identical allele
        return
    ref_set = set(ref_nodes)
    for i in diff_idx:
        ref_base, alt_base = ref_seq[i], alt_seq[i]
        # heuristic: alt_node = first node in ALT at or after i that’s not in REF
        alt_node = next((n for n in alt_nodes[i:] if n not in ref_set), alt_nodes[i])
        ref_node = ref_nodes[i] if i < len(ref_nodes) else ref_nodes[-1]
        yield i, ref_base, ref_node, alt_base, alt_node

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
