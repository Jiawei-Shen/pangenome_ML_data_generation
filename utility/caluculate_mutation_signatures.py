#!/usr/bin/env python3

import argparse
import gzip
import os
from collections import Counter, defaultdict


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1].upper()


def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def parse_gfa_nodes(gfa_path):
    """
    Parse node sequences from GFA S lines.

    GFA format:
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


def all_sbs96_channels():
    """
    COSMIC/SigProfiler-style SBS96 order:

    A[C>A]A
    A[C>A]C
    A[C>A]G
    A[C>A]T
    A[C>G]A
    ...
    T[T>G]T

    Order:
        left base outer loop
        mutation type middle loop
        right base inner loop
    """

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

    for left in bases:
        for ref, alt in mutation_types:
            for right in bases:
                channels.append(f"{left}[{ref}>{alt}]{right}")

    return channels


def get_graph_position(fields):
    """
    Return:
        node_id, pos_1based

    For your graph VCF:
        CHROM is often node id
        POS is node-local 1-based position
        INFO may contain NID and NSO

    If NID/NSO are present:
        NID = node id
        NSO = node sequence offset, likely 0-based

    Your examples:
        POS=6, NSO=5
    so:
        POS = NSO + 1
    """

    chrom = fields[0]
    pos = fields[1]
    info = parse_info(fields[7])

    if "NID" in info and "NSO" in info:
        node_id = str(info["NID"])
        pos_1based = int(info["NSO"]) + 1
        return node_id, pos_1based

    node_id = str(chrom)
    pos_1based = int(pos)
    return node_id, pos_1based


def get_context_same_node_only(nodes, node_id, pos_1based):
    """
    Strict version:
    only use same-node trinucleotide context.

    This avoids unsafe node_id +/- 1 assumptions.

    If mutation is at node boundary, return unresolved.
    """

    node_id = str(node_id)

    if node_id not in nodes:
        return None, "missing_node"

    seq = nodes[node_id]
    idx = pos_1based - 1

    if idx < 0 or idx >= len(seq):
        return None, "position_out_of_range"

    if idx == 0:
        return None, "left_boundary_context_missing"

    if idx == len(seq) - 1:
        return None, "right_boundary_context_missing"

    context = seq[idx - 1: idx + 2].upper()

    if len(context) != 3:
        return None, "bad_context_length"

    if any(base not in "ACGT" for base in context):
        return None, "non_acgt_context"

    return context, "ok"


def get_context_with_node_id_neighbor_rescue(nodes, node_id, pos_1based):
    """
    Debug/relaxed version:
    use node_id - 1 or node_id + 1 only when boundary context is missing.

    WARNING:
    This is not graph-topology-safe unless node IDs are path ordered.
    Use only if you intentionally want this heuristic.
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
        try:
            prev_node_id = str(int(node_id) - 1)
        except ValueError:
            return None, "left_boundary_no_numeric_node_id"

        if prev_node_id in nodes and len(nodes[prev_node_id]) > 0:
            left = nodes[prev_node_id][-1]
        else:
            return None, "left_boundary_neighbor_missing"

    # Right base
    if idx < len(seq) - 1:
        right = seq[idx + 1]
    else:
        try:
            next_node_id = str(int(node_id) + 1)
        except ValueError:
            return None, "right_boundary_no_numeric_node_id"

        if next_node_id in nodes and len(nodes[next_node_id]) > 0:
            right = nodes[next_node_id][0]
        else:
            return None, "right_boundary_neighbor_missing"

    context = (left + center + right).upper()

    if any(base not in "ACGT" for base in context):
        return None, "non_acgt_context"

    return context, "ok"


def canonical_sbs96(context, ref, alt):
    """
    Convert SNV to canonical pyrimidine-centered SBS96 channel.

    Example:
        raw G>A with context TGC
        T[G>A]C

    Reverse complement:
        G[C>T]A

    Output:
        G[C>T]A
    """

    context = context.upper()
    ref = ref.upper()
    alt = alt.upper()

    if len(context) != 3:
        return None

    if len(ref) != 1 or len(alt) != 1:
        return None

    if ref not in "ACGT" or alt not in "ACGT":
        return None

    if ref == alt:
        return None

    # Already canonical
    if ref in {"C", "T"}:
        left, center, right = context
        return f"{left}[{ref}>{alt}]{right}"

    # Reverse-complement to pyrimidine representation
    rc_context = revcomp(context)
    rc_ref = revcomp(ref)
    rc_alt = revcomp(alt)

    left, center, right = rc_context
    return f"{left}[{rc_ref}>{rc_alt}]{right}"


def sample_name_from_path(path):
    base = os.path.basename(path)

    for suffix in [".vcf.gz", ".vcf", ".gz"]:
        if base.endswith(suffix):
            base = base[: -len(suffix)]

    return base


def count_sbs96_for_vcf(
    vcf_path,
    nodes,
    min_prob=None,
    min_af=None,
    require_pass=True,
    allow_node_id_neighbor_rescue=False,
):
    counts = Counter()
    skipped = Counter()
    mismatch_examples = []

    if allow_node_id_neighbor_rescue:
        get_context = get_context_with_node_id_neighbor_rescue
    else:
        get_context = get_context_same_node_only

    with open_maybe_gz(vcf_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                skipped["malformed_vcf_line"] += 1
                continue

            ref = fields[3].upper()
            alt = fields[4].upper()
            filt = fields[6]
            info = parse_info(fields[7])

            if require_pass and filt != "PASS":
                skipped["not_pass"] += 1
                continue

            # Your SNV type is TYPE=X
            vtype = info.get("TYPE")
            if vtype is not None and vtype != "X":
                skipped["not_type_x"] += 1
                continue

            if len(ref) != 1 or len(alt) != 1:
                skipped["not_snv"] += 1
                continue

            if "," in alt:
                skipped["multi_allelic"] += 1
                continue

            if min_prob is not None:
                try:
                    prob = float(info.get("PROB", 0))
                except ValueError:
                    skipped["bad_prob"] += 1
                    continue

                if prob < min_prob:
                    skipped["low_prob"] += 1
                    continue

            if min_af is not None:
                try:
                    af = float(info.get("AF", 0))
                except ValueError:
                    skipped["bad_af"] += 1
                    continue

                if af < min_af:
                    skipped["low_af"] += 1
                    continue

            try:
                node_id, pos_1based = get_graph_position(fields)
            except Exception:
                skipped["bad_graph_position"] += 1
                continue

            context, status = get_context(nodes, node_id, pos_1based)
            if context is None:
                skipped[status] += 1
                continue

            # REF should match center base before canonicalization.
            if context[1] != ref:
                skipped["ref_context_mismatch"] += 1

                if len(mismatch_examples) < 10:
                    mismatch_examples.append(
                        {
                            "node_id": node_id,
                            "pos_1based": pos_1based,
                            "ref": ref,
                            "alt": alt,
                            "context": context,
                            "center": context[1],
                        }
                    )

                continue

            channel = canonical_sbs96(context, ref, alt)

            if channel is None:
                skipped["cannot_canonicalize"] += 1
                continue

            counts[channel] += 1

    return counts, skipped, mismatch_examples


def write_cosmic_matrix(sample_to_counts, out_path):
    """
    Write COSMIC-like SBS96 matrix:

        Mutation Types    Sample_1    Sample_2
        A[C>A]A           10          2
        ...
    """

    channels = all_sbs96_channels()
    sample_names = list(sample_to_counts.keys())

    with open(out_path, "w") as out:
        out.write("Mutation Types\t" + "\t".join(sample_names) + "\n")

        for channel in channels:
            row = [channel]
            for sample in sample_names:
                row.append(str(sample_to_counts[sample].get(channel, 0)))
            out.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert graph VCF(s) + GFA into COSMIC-style SBS96 mutation matrix."
    )

    parser.add_argument(
        "--gfa",
        required=True,
        help="GFA file containing node sequences"
    )

    parser.add_argument(
        "--vcf",
        nargs="+",
        required=True,
        help="One or more graph VCF/VCF.GZ files"
    )

    parser.add_argument(
        "--sample-names",
        nargs="+",
        default=None,
        help="Optional sample names, same number/order as --vcf"
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output COSMIC-style SBS96 matrix TSV"
    )

    parser.add_argument(
        "--min-prob",
        type=float,
        default=None,
        help="Minimum PROB threshold"
    )

    parser.add_argument(
        "--min-af",
        type=float,
        default=None,
        help="Minimum AF threshold"
    )

    parser.add_argument(
        "--allow-non-pass",
        action="store_true",
        help="Do not require FILTER=PASS"
    )

    parser.add_argument(
        "--allow-node-id-neighbor-rescue",
        action="store_true",
        help=(
            "Use node_id-1/node_id+1 to rescue boundary trinucleotide context. "
            "Not recommended for final analysis unless node IDs are path ordered."
        )
    )

    args = parser.parse_args()

    if args.sample_names is not None and len(args.sample_names) != len(args.vcf):
        raise ValueError("--sample-names must have the same number of entries as --vcf")

    print(f"[INFO] Loading GFA nodes: {args.gfa}")
    nodes = parse_gfa_nodes(args.gfa)
    print(f"[INFO] Loaded nodes: {len(nodes):,}")

    sample_to_counts = {}
    all_skipped = {}

    for i, vcf_path in enumerate(args.vcf):
        if args.sample_names is not None:
            sample = args.sample_names[i]
        else:
            sample = sample_name_from_path(vcf_path)

        print("============================================================")
        print(f"[INFO] Sample: {sample}")
        print(f"[INFO] VCF: {vcf_path}")

        counts, skipped, mismatch_examples = count_sbs96_for_vcf(
            vcf_path=vcf_path,
            nodes=nodes,
            min_prob=args.min_prob,
            min_af=args.min_af,
            require_pass=not args.allow_non_pass,
            allow_node_id_neighbor_rescue=args.allow_node_id_neighbor_rescue,
        )

        sample_to_counts[sample] = counts
        all_skipped[sample] = skipped

        print(f"[INFO] Counted SBS96 SNVs: {sum(counts.values()):,}")

        if skipped:
            print("[INFO] Skipped records:")
            for k, v in skipped.most_common():
                print(f"  {k}: {v:,}")

        if mismatch_examples:
            print("[WARNING] First REF/context mismatch examples:")
            for ex in mismatch_examples:
                print(
                    f"  node={ex['node_id']} pos={ex['pos_1based']} "
                    f"REF={ex['ref']} ALT={ex['alt']} "
                    f"context={ex['context']} center={ex['center']}"
                )

    write_cosmic_matrix(sample_to_counts, args.out)

    print("============================================================")
    print(f"[INFO] Wrote COSMIC-style SBS96 matrix: {args.out}")
    print(f"[INFO] Samples: {len(sample_to_counts)}")
    print(f"[INFO] Rows: 96")


if __name__ == "__main__":
    main()