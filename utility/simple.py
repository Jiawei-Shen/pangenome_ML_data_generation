#!/usr/bin/env python3
"""
collect_nodes_from_vcf.py
-------------------------------------------------
Collect ALL node IDs that appear in a vg-deconstructed VCF:

  • ID field  (e.g.  >1745>1742  ➟ 1745, 1742)
  • AT= tag   (e.g.  >1745>1743>1742  ➟ 1745, 1743, 1742)

Output: one node ID per line, unsorted or (optionally) sorted numerically.

Usage
-----
    python collect_nodes_from_vcf.py  input.vcf[.gz]  nodes.txt          # unsorted
    python collect_nodes_from_vcf.py  -s input.vcf.gz nodes.txt          # sorted

Dependencies:  pysam  (pip install pysam)
"""

import argparse
import re
import sys
import pysam

# matches every run of digits, ignoring orientation symbols >  <
_DIGIT_RE = re.compile(r"\d+")

def extract_ids_from_string(s: str):
    """Yield every integer found in a traversal or ID string."""
    for m in _DIGIT_RE.finditer(s):
        yield int(m.group())

def main(vcf_path: str, out_path: str, do_sort: bool):
    seen = set()
    vcf = pysam.VariantFile(vcf_path)

    # -- iterate through records ------------------------------------------------
    for rec in vcf:
        # 1.  ID field  (might be '.')
        if rec.id and rec.id != ".":
            seen.update(extract_ids_from_string(rec.id))

        # 2.  AT tag: list[str] like ['>1745>1743>1742', ...]
        at_vals = rec.info.get("AT")
        if at_vals:
            for trav in at_vals:
                seen.update(extract_ids_from_string(trav))

    # -- write output -----------------------------------------------------------
    with open(out_path, "w") as fout:
        ids = sorted(seen) if do_sort else seen
        for nid in ids:
            fout.write(f"{nid}\n")

    sys.stderr.write(f"Wrote {len(seen):,} unique node IDs ➟ {out_path}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Collect all node IDs referenced in a vg-deconstructed VCF."
    )
    parser.add_argument("vcf",  help="Input VCF (.vcf or .vcf.gz)")
    parser.add_argument("out",  help="Output text file (one node ID per line)")
    parser.add_argument("-s", "--sort", action="store_true",
                        help="sort node IDs numerically before writing")
    args = parser.parse_args()
    main(args.vcf, args.out, args.sort)
