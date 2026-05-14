#!/usr/bin/env python3

import argparse
import gzip
import json
from collections import Counter


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if not item or item == ".":
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def info_to_str(info):
    return ";".join(
        f"{k}={v}" if v is not True else k
        for k, v in info.items()
    ) or "."


def first_present(d, keys, default=None):
    for k in keys:
        if k in d:
            return d[k]
    return default


def to_int(x):
    try:
        return int(x)
    except Exception:
        return None


def iter_node_map_records(json_path):
    with open(json_path, "r") as f:
        data = json.load(f)

    if isinstance(data, list):
        for x in data:
            if isinstance(x, dict):
                yield x
    elif isinstance(data, dict):
        for _, x in data.items():
            if isinstance(x, dict):
                yield x
    else:
        raise ValueError("Unsupported map_json format")


def load_node_map(json_path, node_start_is_1_based=True):
    """
    Returns:
        anchors[node_id] = (chrom, start0, strand, length)
    """
    anchors = {}

    for rec in iter_node_map_records(json_path):
        node_id = to_int(first_present(rec, ["node_id", "node", "id", "nid"]))
        if node_id is None:
            continue

        chrom = first_present(rec, ["chrom", "chr", "chromosome", "contig"])
        strand = str(first_present(rec, ["strand_in_path", "strand"], "+") or "+")
        length = to_int(first_present(rec, ["length", "len", "node_length"]))

        start0 = first_present(rec, ["start_0based", "pos0"])
        start1 = first_present(rec, ["grch38_position_start", "start_1based", "pos1"])

        if start0 is not None:
            s0 = to_int(start0)
        elif start1 is not None:
            s1 = to_int(start1)
            s0 = s1 - 1 if s1 is not None else None
        else:
            generic = first_present(rec, ["start", "pos", "start_pos"])
            sg = to_int(generic)
            if sg is None:
                s0 = None
            else:
                s0 = sg - 1 if node_start_is_1_based else sg

        if chrom is not None and s0 is not None and length is not None:
            anchors[node_id] = (
                str(chrom),
                int(s0),
                strand if strand in ("+", "-") else "+",
                int(length),
            )

    return anchors


def load_neighbor_map(path):
    """
    TSV format:

    node_id    left_node    right_node

    Missing neighbor can be "." or empty.
    """
    neighbor = {}

    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        header_lc = [h.lower() for h in header]

        def idx(name):
            if name not in header_lc:
                raise ValueError(f"neighbor_map missing required column: {name}")
            return header_lc.index(name)

        i_node = idx("node_id")
        i_left = idx("left_node")
        i_right = idx("right_node")

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            node_id = to_int(fields[i_node])
            if node_id is None:
                continue

            left_raw = fields[i_left].strip() if i_left < len(fields) else "."
            right_raw = fields[i_right].strip() if i_right < len(fields) else "."

            left = to_int(left_raw) if left_raw not in ("", ".") else None
            right = to_int(right_raw) if right_raw not in ("", ".") else None

            neighbor[node_id] = (left, right)

    return neighbor


def make_graph_key(fields, require_graph_key=True):
    chrom = fields[0]
    pos = fields[1]
    ref = fields[3]
    alt = fields[4]
    info = parse_info(fields[7])

    nid = info.get("NID")
    nso = info.get("NSO")
    vtype = info.get("TYPE", "")

    if nid is not None and nso is not None:
        return ("GRAPH", nid, nso, vtype, ref, alt)

    if require_graph_key:
        return None

    return ("VCF", chrom, pos, ref, alt)


def load_keys_from_vcf(vcf_path, require_graph_key=True):
    keys = set()
    n_records = 0
    malformed = 0
    missing_graph_key = 0

    with open_maybe_gz(vcf_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                malformed += 1
                continue

            key = make_graph_key(fields, require_graph_key=require_graph_key)
            if key is None:
                missing_graph_key += 1
                continue

            keys.add(key)
            n_records += 1

    return keys, n_records, malformed, missing_graph_key


def rescue_with_neighbor(
    fields,
    anchors,
    neighbor_map,
    offset_is_0_based=True,
    prefer="left",
    require_same_chrom_when_both=True,
):
    """
    Rescue unanchored graph node using left/right anchored neighbor.

    For left anchor:
        rescued POS = left_start0 + left_length + offset0 + 1

    For right anchor:
        rescued POS = right_start0 + offset0 + 1

    This is approximate and should be marked as rescued.
    """
    chrom_node = fields[0]
    ref = fields[3]
    alt = fields[4]
    info = parse_info(fields[7])

    node_id = to_int(info.get("NID", chrom_node))
    offset = to_int(info.get("NSO"))

    if node_id is None or offset is None:
        return None

    off0 = offset if offset_is_0_based else offset - 1
    if off0 < 0:
        return None

    if node_id not in neighbor_map:
        return None

    left_node, right_node = neighbor_map[node_id]

    left_anchor = anchors.get(left_node) if left_node is not None else None
    right_anchor = anchors.get(right_node) if right_node is not None else None

    if require_same_chrom_when_both and left_anchor and right_anchor:
        if left_anchor[0] != right_anchor[0]:
            return None

    candidates = []

    if left_anchor is not None:
        chrom, start0, strand, length = left_anchor
        if strand == "+":
            pos = start0 + length + off0 + 1
            candidates.append(("left", left_node, chrom, pos))

    if right_anchor is not None:
        chrom, start0, strand, length = right_anchor
        if strand == "+":
            pos = start0 + off0 + 1
            candidates.append(("right", right_node, chrom, pos))

    if not candidates:
        return None

    if prefer == "left":
        chosen = next((x for x in candidates if x[0] == "left"), candidates[0])
    elif prefer == "right":
        chosen = next((x for x in candidates if x[0] == "right"), candidates[0])
    else:
        chosen = candidates[0]

    side, neighbor_node, chrom, pos = chosen

    new_info = info.copy()
    new_info["RESCUED_BY_NEIGHBOR"] = side
    new_info["RESCUE_NEIGHBOR_NODE"] = str(neighbor_node)
    new_info["ORIGINAL_GRAPH_NODE"] = str(node_id)
    new_info["ORIGINAL_GRAPH_OFFSET"] = str(offset)

    new_fields = fields[:]
    new_fields[0] = chrom
    new_fields[1] = str(pos)
    new_fields[2] = "."
    new_fields[3] = ref
    new_fields[4] = alt
    new_fields[6] = "NeighborRescued"
    new_fields[7] = info_to_str(new_info)

    return new_fields


def recover_and_rescue(
    all_vcf,
    linear_vcf,
    out_dropped_graph_vcf,
    out_rescued_linear_vcf,
    map_json,
    neighbor_map_path=None,
    rescue_with_neighbor_nodes=False,
    offset_is_0_based=True,
    node_start_is_1_based=True,
    prefer_neighbor="left",
    require_graph_key=True,
):
    linear_keys, n_linear, malformed_linear, missing_linear_graph_key = load_keys_from_vcf(
        linear_vcf,
        require_graph_key=require_graph_key,
    )

    anchors = load_node_map(
        map_json,
        node_start_is_1_based=node_start_is_1_based,
    )

    neighbor_map = {}
    if rescue_with_neighbor_nodes:
        if not neighbor_map_path:
            raise ValueError("--rescue_with_neighbor_nodes requires --neighbor_map")
        neighbor_map = load_neighbor_map(neighbor_map_path)

    n_all = 0
    n_kept_in_linear = 0
    n_dropped = 0
    n_rescued = 0
    n_unrescued = 0
    malformed_all = 0
    missing_all_graph_key = 0

    key_counter_all = Counter()

    with open_maybe_gz(all_vcf) as fin, \
         open(out_dropped_graph_vcf, "w") as fout_drop, \
         open(out_rescued_linear_vcf, "w") as fout_rescue:

        for line in fin:
            if line.startswith("#"):
                fout_drop.write(line)

                if line.startswith("##"):
                    fout_rescue.write(line)
                elif line.startswith("#CHROM"):
                    fout_rescue.write(
                        '##FILTER=<ID=NeighborRescued,Description="Graph record rescued to approximate linear coordinate using anchored neighbor node">\n'
                    )
                    fout_rescue.write(
                        '##INFO=<ID=RESCUED_BY_NEIGHBOR,Number=1,Type=String,Description="Neighbor side used for rescue: left or right">\n'
                    )
                    fout_rescue.write(
                        '##INFO=<ID=RESCUE_NEIGHBOR_NODE,Number=1,Type=String,Description="Anchored neighbor node used for rescue">\n'
                    )
                    fout_rescue.write(
                        '##INFO=<ID=ORIGINAL_GRAPH_NODE,Number=1,Type=String,Description="Original unanchored graph node">\n'
                    )
                    fout_rescue.write(
                        '##INFO=<ID=ORIGINAL_GRAPH_OFFSET,Number=1,Type=String,Description="Original graph offset before rescue">\n'
                    )
                    fout_rescue.write(line)

                continue

            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                malformed_all += 1
                continue

            key = make_graph_key(fields, require_graph_key=require_graph_key)
            if key is None:
                missing_all_graph_key += 1
                continue

            key_counter_all[key] += 1
            n_all += 1

            if key in linear_keys:
                n_kept_in_linear += 1
                continue

            fout_drop.write(line)
            n_dropped += 1

            if rescue_with_neighbor_nodes:
                rescued_fields = rescue_with_neighbor(
                    fields=fields,
                    anchors=anchors,
                    neighbor_map=neighbor_map,
                    offset_is_0_based=offset_is_0_based,
                    prefer=prefer_neighbor,
                    require_same_chrom_when_both=True,
                )

                if rescued_fields is not None:
                    fout_rescue.write("\t".join(rescued_fields) + "\n")
                    n_rescued += 1
                else:
                    n_unrescued += 1
            else:
                n_unrescued += 1

    print("[INFO] Done")
    print(f"[INFO] rescue_with_neighbor_nodes: {rescue_with_neighbor_nodes}")
    print(f"[INFO] all records: {n_all:,}")
    print(f"[INFO] linear records with graph key: {n_linear:,}")
    print(f"[INFO] records from all found in linear by graph key: {n_kept_in_linear:,}")
    print(f"[INFO] dropped/unconverted records: {n_dropped:,}")
    print(f"[INFO] neighbor-rescued records: {n_rescued:,}")
    print(f"[INFO] still-unrescued records: {n_unrescued:,}")
    print(f"[INFO] anchors loaded: {len(anchors):,}")
    print(f"[INFO] neighbor map records: {len(neighbor_map):,}")
    print(f"[INFO] malformed all records: {malformed_all:,}")
    print(f"[INFO] malformed linear records: {malformed_linear:,}")
    print(f"[INFO] missing graph key in all: {missing_all_graph_key:,}")
    print(f"[INFO] missing graph key in linear: {missing_linear_graph_key:,}")
    print(f"[INFO] dropped graph VCF: {out_dropped_graph_vcf}")
    print(f"[INFO] rescued linear VCF: {out_rescued_linear_vcf}")

    duplicated_all = sum(1 for _, v in key_counter_all.items() if v > 1)
    if duplicated_all:
        print(f"[WARNING] duplicated graph keys in all VCF: {duplicated_all:,}")


def main():
    p = argparse.ArgumentParser(
        description="Recover dropped graph records. Neighbor-node rescue is OFF by default."
    )

    p.add_argument("--all", required=True, help="All graph/node-form VCF")
    p.add_argument("--linear", required=True, help="Linear converted VCF")
    p.add_argument("--map_json", required=True, help="Node map JSON with GRCh38 anchors")

    p.add_argument("--out_dropped_graph", required=True, help="Output dropped graph-format VCF")
    p.add_argument("--out_rescued_linear", required=True, help="Output rescued linear-coordinate VCF")

    p.add_argument(
        "--neighbor_map",
        default=None,
        help="Optional TSV with columns: node_id, left_node, right_node",
    )

    p.add_argument(
        "--rescue_with_neighbor_nodes",
        action="store_true",
        help="Enable rescue of dropped graph records using anchored left/right neighbor nodes. Default: OFF.",
    )

    p.add_argument(
        "--offset_is_0_based",
        type=lambda s: str(s).lower() in ("1", "true", "t", "yes", "y"),
        default=True,
        help="Whether NSO/v_pos offset is 0-based. Default: true.",
    )

    p.add_argument(
        "--node_start_is_1_based",
        type=lambda s: str(s).lower() in ("1", "true", "t", "yes", "y"),
        default=True,
        help="Whether grch38_position_start in map_json is 1-based. Default: true.",
    )

    p.add_argument(
        "--prefer_neighbor",
        choices=["left", "right", "any"],
        default="left",
        help="Which anchored neighbor to prefer when both left and right are available. Default: left.",
    )

    p.add_argument(
        "--allow_vcf_pos_fallback",
        action="store_true",
        help="Allow fallback to CHROM/POS/REF/ALT if NID/NSO missing. Not recommended.",
    )

    args = p.parse_args()

    recover_and_rescue(
        all_vcf=args.all,
        linear_vcf=args.linear,
        out_dropped_graph_vcf=args.out_dropped_graph,
        out_rescued_linear_vcf=args.out_rescued_linear,
        map_json=args.map_json,
        neighbor_map_path=args.neighbor_map,
        rescue_with_neighbor_nodes=args.rescue_with_neighbor_nodes,
        offset_is_0_based=args.offset_is_0_based,
        node_start_is_1_based=args.node_start_is_1_based,
        prefer_neighbor=args.prefer_neighbor,
        require_graph_key=not args.allow_vcf_pos_fallback,
    )


if __name__ == "__main__":
    main()