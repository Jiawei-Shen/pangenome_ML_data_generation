#!/usr/bin/env python3
"""
extract_core_snps.py

For every record in a vg-deconstructed VCF:

1.  keep any ALT allele whose sequence length equals the REF length
    (so it is an MNP at most, not an indel);
2.  compare the REF and ALT traversals node-by-node;
3.  for every position where nucleotide and node BOTH differ,
    emit one line:

   CHROM  SNP_POS  ID  REF_BASE  REF_NODE  ALT_BASE  ALT_NODE

   • SNP_POS  = record.POS + within-record offset (1-based)
   • REF_NODE / ALT_NODE = node IDs that actually encode that base
"""

import argparse
import re
import sys
from typing import List

import pysam

# split on orientation symbols '>' or '<'
_SPLIT = re.compile(r"[><]+")


def parse_nodes(trav: str) -> List[int]:
    """Turn a traversal string into a list of (orientation-stripped) node IDs."""
    return [int(tok) for tok in _SPLIT.split(trav) if tok]


def is_substitution_allele(ref_seq: str, alt_seq: str) -> bool:
    """True when ALT is the same length as REF (SNP or MNP, not indel)."""
    return len(ref_seq) == len(alt_seq)


def find_snp_positions(
    ref_seq: str,
    ref_nodes: List[int],
    alt_seq: str,
    alt_nodes: List[int],
) -> List[tuple[int, str, int]]:
    """
    Return a list of tuples (offset, alt_base, alt_node) for every single-base
    difference between REF and this ALT.  A difference is counted only when

        • the nucleotide differs AND
        • the node ID differs
    """
    if (
        len(ref_seq) != len(alt_seq)
        or len(ref_nodes) != len(alt_nodes)
        or len(ref_seq) != len(ref_nodes)  # safety check: 1 node == 1 base
    ):
        return []

    snps = []
    for i, (rb, ab, rn, an) in enumerate(
        zip(ref_seq, alt_seq, ref_nodes, alt_nodes)
    ):
        if rb != ab and rn != an:
            snps.append((i, ab, an))
    return snps


def main(vcf_path: str):
    vcf = pysam.VariantFile(vcf_path)
    if "AT" not in vcf.header.info:
        sys.stderr.write("ERROR: VCF header lacks an AT tag\n")
        sys.exit(1)

    print(
        "CHROM\tSNP_POS\tID\tREF_BASE\tREF_NODE\tALT_BASE\tALT_NODE"
    )

    for rec in vcf:
        ref_seq = rec.ref
        ref_nodes = parse_nodes(rec.info["AT"][0])

        for alt_idx, alt_seq in enumerate(rec.alts):
            if not is_substitution_allele(ref_seq, alt_seq):
                continue  # skip indels

            alt_nodes = parse_nodes(rec.info["AT"][alt_idx + 1])
            for offset, alt_base, alt_node in find_snp_positions(
                ref_seq, ref_nodes, alt_seq, alt_nodes
            ):
                snp_pos = rec.pos + offset  # 1-based VCF coordinate
                ref_base = ref_seq[offset]
                ref_node = ref_nodes[offset]

                print(
                    f"{rec.contig}\t{snp_pos}\t{rec.id}\t"
                    f"{ref_base}\t{ref_node}\t{alt_base}\t{alt_node}"
                )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=(
            "Extract per-base SNP node IDs from a vg-deconstructed VCF, "
            "including cases where the record is a multi-base substitution."
        )
    )
    ap.add_argument("vcf", help="input .vcf or .vcf.gz file")
    main(ap.parse_args().vcf)
