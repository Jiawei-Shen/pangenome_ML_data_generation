#!/usr/bin/env python3

import argparse
import gzip
from collections import Counter


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1].upper()


def parse_gfa_nodes(gfa_path):
    """
    Parse GFA S-lines.

    Expected GFA format:
    S   node_id   sequence
    """
    nodes = {}

    with open_maybe_gz(gfa_path) as f:
        for line in f:
            if not line.startswith("S\t"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue

            node_id = fields[1]
            seq = fields[2].upper()

            if seq == "*":
                continue

            nodes[node_id] = seq

    return nodes


def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def get_context_from_node(nodes, node_id, pos_1based):
    """
    Get trinucleotide context around a 1-based position.

    If the mutation is at the start of a node, use node_id - 1 as left context.
    If the mutation is at the end of a node, use node_id + 1 as right context.

    This follows your requested heuristic:
    missing context base -> use numeric neighbor node id +/- 1.

    Returns:
        context, status
    """

    node_id = str(node_id)

    if node_id not in nodes:
        return None, "missing_node"

    seq = nodes[node_id]
    idx = pos_1based - 1

    if idx < 0 or idx >= len(seq):
        return None, "position_out_of_range"

    center = seq[idx]

    # Left base
    if idx > 0:
        left = seq[idx - 1]
    else:
        prev_node_id = str(int(node_id) - 1)
        if prev_node_id in nodes and len(nodes[prev_node_id]) > 0:
            left = nodes[prev_node_id][-1]
        else:
            return None, "missing_left_context"

    # Right base
    if idx < len(seq) - 1:
        right = seq[idx + 1]
    else:
        next_node_id = str(int(node_id) + 1)
        if next_node_id in nodes and len(nodes[next_node_id]) > 0:
            right = nodes[next_node_id][0]
        else:
            return None, "missing_right_context"

    context = (left + center + right).upper()

    if any(base not in "ACGT" for base in context):
        return None, "non_acgt_context"

    return context, "ok"


def canonical_sbs96(context, ref, alt):
    """
    Convert mutation to canonical pyrimidine-centered SBS96 notation.

    Input:
        context: 3-mer, center should equal REF
        ref: REF allele
        alt: ALT allele

    Output:
        e.g. A[C>T]G
    """

    context = context.upper()
    ref = ref.upper()
    alt = alt.upper()

    if len(ref) != 1 or len(alt) != 1:
        return None

    if ref not in "ACGT" or alt not in "ACGT":
        return None

    if ref == alt:
        return None

    # Canonical mutation should have C or T as REF.
    if ref in {"C", "T"}:
        left, center, right = context
        return f"{left}[{ref}>{alt}]{right}"

    # If REF is A or G, reverse-complement the whole event.
    rc_context = revcomp(context)
    rc_ref = revcomp(ref)
    rc_alt = revcomp(alt)

    left, center, right = rc_context
    return f"{left}[{rc_ref}>{rc_alt}]{right}"


def all_sbs96_channels():
    bases = ["A", "C", "G", "T"]
    mutation_types = [
        ("C", "A"),
        ("C", "G"),
        ("C", "T"),
        ("T", "A"),
        ("T", "C"),
        ("T", "G"),
    ]

    channels = []
    for ref, alt in mutation_types:
        for left in bases:
            for right in bases:
                channels.append(f"{left}[{ref}>{alt}]{right}")

    return channels


def build_signature(vcf_path, nodes, min_prob=None, min_af=None, require_pass=True):
    counts = Counter()
    skipped = Counter()
    examples = []

    with open_maybe_gz(vcf_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                skipped["malformed_vcf_line"] += 1
                continue

            chrom, pos, vid, ref, alt, qual, filt, info_str = fields[:8]

            if require_pass and filt != "PASS":
                skipped["not_pass"] += 1
                continue

            info = parse_info(info_str)

            # Only SNV/X records
            var_type = info.get("TYPE", None)
            if var_type is not None and var_type != "X":
                skipped["not_snv_type_x"] += 1
                continue

            if len(ref) != 1 or len(alt) != 1:
                skipped["not_snv"] += 1
                continue

            if "," in alt:
                skipped["multi_allelic"] += 1
                continue

            try:
                pos_1based = int(pos)
            except ValueError:
                skipped["bad_pos"] += 1
                continue

            if min_prob is not None:
                prob = float(info.get("PROB", 0))
                if prob < min_prob:
                    skipped["low_prob"] += 1
                    continue

            if min_af is not None:
                af = float(info.get("AF", 0))
                if af < min_af:
                    skipped["low_af"] += 1
                    continue

            node_id = info.get("NID", chrom)

            context, status = get_context_from_node(nodes, node_id, pos_1based)
            if context is None:
                skipped[status] += 1
                continue

            # Optional sanity check: center base should match REF.
            if context[1] != ref.upper():
                skipped["ref_context_mismatch"] += 1

                if len(examples) < 10:
                    examples.append(
                        {
                            "node": node_id,
                            "pos": pos_1based,
                            "vcf_ref": ref,
                            "vcf_alt": alt,
                            "context": context,
                            "center_base": context[1],
                        }
                    )

                # Usually better to skip mismatches.
                continue

            channel = canonical_sbs96(context, ref, alt)
            if channel is None:
                skipped["cannot_canonicalize"] += 1
                continue

            counts[channel] += 1

    return counts, skipped, examples


def write_signature_tsv(counts, out_path):
    channels = all_sbs96_channels()
    total = sum(counts.values())

    with open(out_path, "w") as out:
        out.write("channel\tcount\tfrequency\n")
        for ch in channels:
            count = counts.get(ch, 0)
            freq = count / total if total > 0 else 0
            out.write(f"{ch}\t{count}\t{freq:.8f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Build SBS96 mutation signature from graph-style VCF and GFA."
    )

    parser.add_argument("--vcf", required=True, help="Graph-style VCF file")
    parser.add_argument("--gfa", required=True, help="GFA file containing node sequences")
    parser.add_argument("--out", required=True, help="Output SBS96 TSV")
    parser.add_argument("--min-prob", type=float, default=None, help="Minimum PROB filter")
    parser.add_argument("--min-af", type=float, default=None, help="Minimum AF filter")
    parser.add_argument(
        "--allow-non-pass",
        action="store_true",
        help="Do not require FILTER=PASS"
    )

    args = parser.parse_args()

    print(f"[INFO] Loading GFA nodes: {args.gfa}")
    nodes = parse_gfa_nodes(args.gfa)
    print(f"[INFO] Loaded {len(nodes):,} nodes")

    print(f"[INFO] Reading VCF: {args.vcf}")
    counts, skipped, mismatch_examples = build_signature(
        vcf_path=args.vcf,
        nodes=nodes,
        min_prob=args.min_prob,
        min_af=args.min_af,
        require_pass=not args.allow_non_pass,
    )

    write_signature_tsv(counts, args.out)

    print(f"[INFO] Wrote signature: {args.out}")
    print(f"[INFO] Total counted SNVs: {sum(counts.values()):,}")

    print("[INFO] Skipped records:")
    for k, v in skipped.most_common():
        print(f"  {k}: {v:,}")

    if mismatch_examples:
        print("[WARNING] Example REF/context mismatches:")
        for ex in mismatch_examples:
            print(
                f"  node={ex['node']} pos={ex['pos']} "
                f"VCF={ex['vcf_ref']}>{ex['vcf_alt']} "
                f"context={ex['context']} center={ex['center_base']}"
            )


if __name__ == "__main__":
    main()