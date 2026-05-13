#!/usr/bin/env python3

import argparse
import gzip
from collections import Counter


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


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


def make_graph_key(fields):
    """
    Build a graph-coordinate variant key.

    Preferred:
        NID, NSO, TYPE, REF, ALT

    Fallback:
        CHROM, POS, REF, ALT
    """

    chrom = fields[0]
    pos = fields[1]
    ref = fields[3]
    alt = fields[4]
    info_str = fields[7] if len(fields) > 7 else ""
    info = parse_info(info_str)

    nid = info.get("NID")
    nso = info.get("NSO")
    vtype = info.get("TYPE", "")

    if nid is not None and nso is not None:
        return ("GRAPH", nid, nso, vtype, ref, alt)

    return ("VCF", chrom, pos, ref, alt)


def get_prob(fields):
    """
    Extract PROB from INFO.

    Returns None if PROB is missing or invalid.
    """

    if len(fields) < 8:
        return None

    info = parse_info(fields[7])

    if "PROB" not in info:
        return None

    try:
        return float(info["PROB"])
    except ValueError:
        return None


def load_keys_from_vcf(vcf_path):
    keys = set()
    n_records = 0
    malformed = 0

    with open_maybe_gz(vcf_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                malformed += 1
                continue

            key = make_graph_key(fields)
            keys.add(key)
            n_records += 1

    return keys, n_records, malformed


def recover_dropped(all_vcf, linear_vcf, out_vcf, min_prob=None):
    linear_keys, n_linear, malformed_linear = load_keys_from_vcf(linear_vcf)

    n_all = 0
    n_all_after_prob_filter = 0
    n_kept_in_linear = 0
    n_dropped = 0
    n_low_prob = 0
    n_missing_prob = 0
    malformed_all = 0

    key_counter_all = Counter()

    with open_maybe_gz(all_vcf) as fin, open(out_vcf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                malformed_all += 1
                continue

            n_all += 1

            if min_prob is not None:
                prob = get_prob(fields)

                if prob is None:
                    n_missing_prob += 1
                    continue

                # only keep PROB > x
                if prob <= min_prob:
                    n_low_prob += 1
                    continue

            n_all_after_prob_filter += 1

            key = make_graph_key(fields)
            key_counter_all[key] += 1

            if key in linear_keys:
                n_kept_in_linear += 1
            else:
                fout.write(line)
                n_dropped += 1

    print("[INFO] Done")
    print(f"[INFO] all records: {n_all:,}")

    if min_prob is not None:
        print(f"[INFO] min PROB filter: PROB > {min_prob}")
        print(f"[INFO] all records after PROB filter: {n_all_after_prob_filter:,}")
        print(f"[INFO] skipped low PROB records: {n_low_prob:,}")
        print(f"[INFO] skipped missing/invalid PROB records: {n_missing_prob:,}")

    print(f"[INFO] linear records: {n_linear:,}")
    print(f"[INFO] recovered dropped records: {n_dropped:,}")
    print(f"[INFO] records from all found in linear by graph key: {n_kept_in_linear:,}")
    print(f"[INFO] malformed all records: {malformed_all:,}")
    print(f"[INFO] malformed linear records: {malformed_linear:,}")
    print(f"[INFO] output: {out_vcf}")

    if n_dropped == 0:
        print("[WARNING] No dropped records recovered.")
        print("[WARNING] This may mean the linear VCF does not preserve NID/NSO graph metadata.")
        print("[WARNING] Check INFO fields in the linear VCF.")

    duplicated_all = sum(1 for k, v in key_counter_all.items() if v > 1)
    if duplicated_all:
        print(f"[WARNING] duplicated keys in filtered all VCF: {duplicated_all:,}")


def main():
    parser = argparse.ArgumentParser(
        description="Recover dropped/unconverted graph-format records from all.vcf.gz minus linear.vcf.gz."
    )

    parser.add_argument("--all", required=True, help="All graph-format VCF")
    parser.add_argument("--linear", required=True, help="Linear converted VCF")
    parser.add_argument("--out", required=True, help="Output dropped graph-format VCF")

    parser.add_argument(
        "--min-prob",
        type=float,
        default=None,
        help="Only keep records with PROB > this value"
    )

    args = parser.parse_args()

    recover_dropped(
        all_vcf=args.all,
        linear_vcf=args.linear,
        out_vcf=args.out,
        min_prob=args.min_prob,
    )


if __name__ == "__main__":
    main()