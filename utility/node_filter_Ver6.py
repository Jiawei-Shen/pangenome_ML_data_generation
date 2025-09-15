#!/usr/bin/env python3
import json
import struct
import argparse
import sys
import time
import copy
import re
import os
import pysam  # Replaces subprocess and shutil for VCF operations

# ─────────────────────────────────────────────────────────────────────────────
# File format constants

# Global header (unchanged)
GLOBAL_MAGIC = b"MYFMT\x01"
GLOBAL_VER_PACK = struct.Struct("<BBI16s")   # major, minor, block_count, reserved[16]

# LATEST per-block header (18 bytes): nid, nrec, flags, R, C
BLOCK_HDR_PACK_LATEST = struct.Struct("<I I H I I")
BLOCK_HDR_SIZE_LATEST = BLOCK_HDR_PACK_LATEST.size  # 18

# OLDER per-block header (14 bytes): nid, nrec, flags, node_length
BLOCK_HDR_PACK_OLD = struct.Struct("<I I H I")
BLOCK_HDR_SIZE_OLD = BLOCK_HDR_PACK_OLD.size        # 14

# LATEST .idx entry (30 bytes): nid, offset, block_size, n_records, flags, R, C
IDX_ENTRY_PACK_LATEST = struct.Struct("<I Q I I H I I")
IDX_ENTRY_SIZE_LATEST = IDX_ENTRY_PACK_LATEST.size  # 30

# "New (older)" .idx entry (26 bytes): nid, offset, block_size, n_records, flags, node_length
IDX_ENTRY_PACK_26 = struct.Struct("<I Q I I H I")
IDX_ENTRY_SIZE_26 = IDX_ENTRY_PACK_26.size          # 26

# Oldest .idx entry (22 bytes): nid, offset, block_size, n_records, flags
IDX_ENTRY_PACK_22 = struct.Struct("<I Q I I H")
IDX_ENTRY_SIZE_22 = IDX_ENTRY_PACK_22.size          # 22


# ─────────────────────────────────────────────────────────────────────────────
# IDX loader with format auto-detect (30 / 26 / 22). Optionally peeks into .dat.
# ─────────────────────────────────────────────────────────────────────────────
def load_index(idx_path, dat_path=None):
    """
    Load .idx (supports latest and older layouts):

      • Latest (30B/entry): <I Q I I H I I> → nid, offset, block_size, n_records, flags, R, C
      • 26B/entry:          <I Q I I H I>   → nid, offset, block_size, n_records, flags, node_length
      • 22B/entry:          <I Q I I H>     → nid, offset, block_size, n_records, flags

    Returns dict[nid] = {
        "start": offset,
        "block_size": block_size,
        "n_records": n_records,
        "flags": flags,
        # present if known (0 if unknown):
        "max_read_len": R,
        "max_cigar_len": C,
        "node_length": node_length
    }

    If dat_path is given and any fields are missing (e.g., old idx),
    we'll peek into .dat block headers to fill them (preferring latest header <R,C>,
    falling back to old header <node_length>).
    """
    node_index = {}

    try:
        with open(idx_path, "rb") as f:
            raw = f.read(4)
            if len(raw) != 4:
                print(f"Error: Could not read block count from {idx_path}", file=sys.stderr)
                return {}
            (count,) = struct.unpack("<I", raw)

            f.seek(0, os.SEEK_END)
            size = f.tell()
            remaining = size - 4
            if count <= 0 or remaining <= 0:
                print(f"Error: {idx_path} is empty or corrupt.", file=sys.stderr)
                return {}

            # auto-detect entry size
            entry_size = remaining // count
            if entry_size not in (IDX_ENTRY_SIZE_LATEST, IDX_ENTRY_SIZE_26, IDX_ENTRY_SIZE_22):
                # If not cleanly divisible, try to guess latest first
                if remaining % IDX_ENTRY_SIZE_LATEST == 0:
                    entry_size = IDX_ENTRY_SIZE_LATEST
                elif remaining % IDX_ENTRY_SIZE_26 == 0:
                    entry_size = IDX_ENTRY_SIZE_26
                elif remaining % IDX_ENTRY_SIZE_22 == 0:
                    entry_size = IDX_ENTRY_SIZE_22
                else:
                    print(f"Warning: Unrecognized idx layout (remaining={remaining}, count={count}). "
                          f"Attempting latest (30B) parsing.", file=sys.stderr)
                    entry_size = IDX_ENTRY_SIZE_LATEST

            f.seek(4)
            for i in range(count):
                data = f.read(entry_size)
                if len(data) != entry_size:
                    print(f"Error: Truncated .idx at entry {i} in {idx_path}", file=sys.stderr)
                    return {}

                if entry_size == IDX_ENTRY_SIZE_LATEST:
                    nid, off, blk_sz, nrec, flg, R, C = IDX_ENTRY_PACK_LATEST.unpack(data)
                    node_index[nid] = {
                        "start": off,
                        "block_size": blk_sz,
                        "n_records": nrec,
                        "flags": flg,
                        "max_read_len": int(R),
                        "max_cigar_len": int(C),
                        "node_length": 0,
                    }
                elif entry_size == IDX_ENTRY_SIZE_26:
                    nid, off, blk_sz, nrec, flg, nlen = IDX_ENTRY_PACK_26.unpack(data)
                    node_index[nid] = {
                        "start": off,
                        "block_size": blk_sz,
                        "n_records": nrec,
                        "flags": flg,
                        "max_read_len": 0,
                        "max_cigar_len": 0,
                        "node_length": int(nlen),
                    }
                else:  # 22B
                    nid, off, blk_sz, nrec, flg = IDX_ENTRY_PACK_22.unpack(data)
                    node_index[nid] = {
                        "start": off,
                        "block_size": blk_sz,
                        "n_records": nrec,
                        "flags": flg,
                        "max_read_len": 0,
                        "max_cigar_len": 0,
                        "node_length": 0,
                    }

    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except struct.error as e:
        print(f"Error: Could not unpack data from IDX file {idx_path}. Details: {e}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return {}

    # Fill missing per-block details from .dat (if provided)
    if dat_path and any((info.get("max_read_len", 0) <= 0 and info.get("node_length", 0) <= 0)
                        for info in node_index.values()):
        try:
            with open(dat_path, "rb") as df:
                magic = df.read(len(GLOBAL_MAGIC))
                if magic != GLOBAL_MAGIC:
                    print(f"Warning: .dat magic mismatch for {dat_path}; header peek skipped.", file=sys.stderr)
                else:
                    _ = GLOBAL_VER_PACK.unpack(df.read(GLOBAL_VER_PACK.size))  # version, block_count (not strictly used)
                    for nid, info in node_index.items():
                        # Skip blocks that already have R/C (latest) or node_length
                        if info.get("max_read_len", 0) > 0 and info.get("max_cigar_len", 0) > 0:
                            continue
                        if info.get("node_length", 0) > 0:
                            continue

                        df.seek(info["start"], os.SEEK_SET)
                        # Try latest 18-byte header first
                        hdr = df.read(BLOCK_HDR_SIZE_LATEST)
                        if len(hdr) == BLOCK_HDR_SIZE_LATEST:
                            nid2, nrec2, flg2, R, C = BLOCK_HDR_PACK_LATEST.unpack(hdr)
                            if nid2 == nid:
                                if R > 0 and C > 0:
                                    info["max_read_len"] = int(R)
                                    info["max_cigar_len"] = int(C)
                                    # Header matched: no need to try old layout for this block
                                    continue
                        # Fall back to old 14-byte header
                        df.seek(info["start"], os.SEEK_SET)
                        hdr_old = df.read(BLOCK_HDR_SIZE_OLD)
                        if len(hdr_old) == BLOCK_HDR_SIZE_OLD:
                            nid2, nrec2, flg2, nlen = BLOCK_HDR_PACK_OLD.unpack(hdr_old)
                            if nid2 == nid:
                                info["node_length"] = int(nlen)
        except FileNotFoundError:
            print(f"Note: .dat not found at {dat_path}; could not fill missing block details.", file=sys.stderr)
        except Exception as e:
            print(f"Note: Could not read .dat to fill block details: {e}", file=sys.stderr)

    return node_index


# ─────────────────────────────────────────────────────────────────────────────
# JSON helper (unchanged)
# ─────────────────────────────────────────────────────────────────────────────
def load_json_data_ids_and_map(json_filepath):
    main_json_data = {}
    node_ids_set = set()
    json_nodes_map = {}
    try:
        with open(json_filepath, 'r') as f_json:
            main_json_data = json.load(f_json)
            if not isinstance(main_json_data, dict):
                print(f"Error: Root JSON content in {json_filepath} is not an object/dictionary.", file=sys.stderr)
                return None, set(), {}
            node_list_from_json = main_json_data.get("nodes")
            if not isinstance(node_list_from_json, list):
                print(f"Warning: 'nodes' key not found or not a list in {json_filepath}.", file=sys.stderr)
                return main_json_data, set(), {}
            for item_index, item in enumerate(node_list_from_json):
                if not isinstance(item, dict): continue
                try:
                    node_id_str = item.get("node_id")
                    if node_id_str is None: continue
                    node_id_int = int(node_id_str)
                    node_ids_set.add(node_id_int)
                    json_nodes_map[node_id_int] = item
                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Node ID format error for item {item_index} in 'nodes' (JSON): '{node_id_str}'. Error: {e}. Skipping.",
                        file=sys.stderr)
            return main_json_data, node_ids_set, json_nodes_map
    except FileNotFoundError:
        print(f"Error: JSON file not found: {json_filepath}", file=sys.stderr)
        return None, set(), {}
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {json_filepath}.", file=sys.stderr)
        return None, set(), {}
    except Exception as e:
        print(f"Unexpected error reading JSON {json_filepath}: {e}", file=sys.stderr)
        return None, set(), {}


# ─────────────────────────────────────────────────────────────────────────────
# Chromosome extractor (unchanged)
# ─────────────────────────────────────────────────────────────────────────────
def extract_chromosome_from_path_pattern(path_pattern):
    if not path_pattern or not isinstance(path_pattern, str):
        return None
    match = re.search(r'(chr([0-9A-Za-z_]+|X|Y|M|MT))', path_pattern)
    if match:
        return match.group(1)
    return None


# ─────────────────────────────────────────────────────────────────────────────
# Main filter (unchanged semantics; calls new load_index)
# ─────────────────────────────────────────────────────────────────────────────
def filter_json_nodes_and_write(json_filepath, idx_filepath, output_json_filepath,
                                vcf_file=None, txt_output_path=None, dat_filepath=None):
    print("Step 1: Loading JSON data and node map...")
    main_json_structure, target_node_ids_from_json, json_nodes_map = load_json_data_ids_and_map(json_filepath)
    if main_json_structure is None:
        print("Critical error: Could not load main JSON. Cannot proceed.")
        return
    print(f"Loaded {len(json_nodes_map)} nodes from JSON for potential processing.")

    chromosome_for_vcf = None
    if vcf_file:
        print(f"Using pysam for VCF queries on: {vcf_file}")
        path_pattern = main_json_structure.get("path_name_input_pattern")
        chromosome_for_vcf = extract_chromosome_from_path_pattern(path_pattern)
        if not chromosome_for_vcf:
            print(
                f"Warning: Could not extract chromosome from 'path_name_input_pattern': '{path_pattern}'. Using 'chr1' as fallback.",
                file=sys.stderr)
            chromosome_for_vcf = "chr1"
        else:
            print(f"Extracted chromosome '{chromosome_for_vcf}' for VCF queries.")

    print("\nStep 2: Loading node index (IDX)...")
    # Pass dat_filepath so we can fill R/C or node_length when missing
    idx_data = load_index(idx_filepath, dat_path=dat_filepath)
    available_idx_node_ids = set(idx_data.keys())
    print(f"Found {len(available_idx_node_ids)} unique node IDs in {idx_filepath}")

    print("\nStep 3: Identifying common node IDs...")
    common_node_ids_int = target_node_ids_from_json.intersection(available_idx_node_ids)
    num_initial_common_nodes = len(common_node_ids_int)
    if num_initial_common_nodes == 0:
        print("No common node IDs found between JSON and IDX. Output 'nodes' list will be empty.")
    else:
        print(f"Found {num_initial_common_nodes} common node IDs between JSON and IDX for potential processing.")

    ultimate_filtered_nodes_list = []
    total_nodes_queried_with_pysam = 0
    total_nodes_filtered_out_by_vcf_prereq = 0
    total_nodes_filtered_out_by_empty_vcf = 0
    total_nodes_kept_with_vcf_results = 0

    if num_initial_common_nodes > 0:
        print(
            f"\nStep 4: Processing {num_initial_common_nodes} common nodes for filtering and VCF query (if applicable)...")

        try:
            with pysam.VariantFile(vcf_file) if vcf_file else open(os.devnull, 'w') as vcf:
                sorted_common_ids = sorted(list(common_node_ids_int))

                batch_node_counter_step4 = 0
                batch_queried_for_vcf_step4 = 0
                batch_kept_with_vcf_results_step4 = 0
                batch_start_time_step4 = time.time()

                for node_idx_step4, node_id_int in enumerate(sorted_common_ids):
                    original_node_item = json_nodes_map.get(node_id_int)
                    if not original_node_item:
                        sys.stderr.write(f"Internal Warning: Common node ID {node_id_int} not in JSON map. Skipping.\n")
                        continue

                    current_node_copy = copy.deepcopy(original_node_item)
                    include_node_in_final_output = True

                    if vcf_file and chromosome_for_vcf:
                        try:
                            start_0based_str = current_node_copy.get("grch38_position_start")
                            length_str = current_node_copy.get("length")

                            if start_0based_str is None or length_str is None:
                                include_node_in_final_output = False
                                total_nodes_filtered_out_by_vcf_prereq += 1
                            else:
                                start_0based = int(start_0based_str)
                                length = int(length_str)
                                if length <= 0:
                                    include_node_in_final_output = False
                                    total_nodes_filtered_out_by_vcf_prereq += 1
                                else:
                                    query_end_0based = start_0based + length
                                    vcf_records_iterator = vcf.fetch(chromosome_for_vcf, start_0based, query_end_0based)
                                    vcf_results = [str(record).strip() for record in vcf_records_iterator]

                                    current_node_copy["vcf_query_results"] = vcf_results
                                    total_nodes_queried_with_pysam += 1
                                    batch_queried_for_vcf_step4 += 1

                                    if not vcf_results:
                                        include_node_in_final_output = False
                                        total_nodes_filtered_out_by_empty_vcf += 1
                                    else:
                                        total_nodes_kept_with_vcf_results += 1
                                        batch_kept_with_vcf_results_step4 += 1
                        except (ValueError, KeyError) as e:
                            include_node_in_final_output = False
                            total_nodes_filtered_out_by_vcf_prereq += 1

                    if include_node_in_final_output:
                        ultimate_filtered_nodes_list.append(current_node_copy)

                    batch_node_counter_step4 += 1
                    if (batch_node_counter_step4 % 100000 == 0 and batch_node_counter_step4 > 0) or \
                            (node_idx_step4 + 1 == num_initial_common_nodes and batch_node_counter_step4 > 0):
                        batch_time_step4 = time.time() - batch_start_time_step4
                        print(
                            f"  Step 4 Batch: Processed {batch_node_counter_step4} nodes ({node_idx_step4 + 1}/{num_initial_common_nodes} total common). "
                            f"VCF queried in batch: {batch_queried_for_vcf_step4}. Kept with VCF results in batch: {batch_kept_with_vcf_results_step4}. "
                            f"Time: {batch_time_step4:.2f}s.")
                        batch_node_counter_step4 = 0
                        batch_queried_for_vcf_step4 = 0
                        batch_kept_with_vcf_results_step4 = 0
                        batch_start_time_step4 = time.time()
        except FileNotFoundError:
            print(f"Error: VCF file not found at '{vcf_file}'. Cannot perform VCF filtering.", file=sys.stderr)
            ultimate_filtered_nodes_list = [json_nodes_map[node_id] for node_id in common_node_ids_int]
        except ValueError as e:
            print(f"Error processing VCF file '{vcf_file}'. Is it a valid bgzipped VCF with a .tbi index? Details: {e}",
                  file=sys.stderr)
            ultimate_filtered_nodes_list = [json_nodes_map[node_id] for node_id in common_node_ids_int]

    output_json_structure = copy.deepcopy(main_json_structure)
    output_json_structure["nodes"] = ultimate_filtered_nodes_list
    num_ultimate_nodes = len(ultimate_filtered_nodes_list)
    print(f"\nStep 5: Final processing complete. {num_ultimate_nodes} nodes will be in output JSON.")

    if vcf_file:
        print(f"VCF Query Summary:")
        print(f"  Total nodes eligible and queried with pysam: {total_nodes_queried_with_pysam}")
        print(
            f"  Nodes filtered out due to VCF prerequisites (missing coords/length, etc.): {total_nodes_filtered_out_by_vcf_prereq}")
        print(f"  Nodes filtered out due to empty VCF query results: {total_nodes_filtered_out_by_empty_vcf}")
        print(f"  Nodes kept that had non-empty VCF results: {total_nodes_kept_with_vcf_results}")

    if txt_output_path:
        print(f"\nStep 6: Writing {num_ultimate_nodes} ultimate filtered node ID(s) to {txt_output_path}...")
        if ultimate_filtered_nodes_list:
            try:
                with open(txt_output_path, 'w') as f_txt:
                    ids_to_write = [str(node.get("node_id")) for node in ultimate_filtered_nodes_list if
                                    node.get("node_id") is not None]
                    try:
                        sorted_ids = sorted(ids_to_write, key=int)
                    except ValueError:
                        sorted_ids = sorted(ids_to_write)
                    for node_id_str in sorted_ids: f_txt.write(f"{node_id_str}\n")
                print(f"Successfully wrote ultimate filtered node IDs to {txt_output_path}")
            except Exception as e:
                print(f"Error writing ultimate filtered node IDs to TXT {txt_output_path}: {e}", file=sys.stderr)
        else:
            print(f"No ultimate filtered node IDs to write to {txt_output_path}.")
            try:
                with open(txt_output_path, 'w') as f_txt:
                    pass
                print(f"Created empty TXT file at {txt_output_path}.")
            except Exception as e:
                print(f"Error creating empty TXT file {txt_output_path}: {e}", file=sys.stderr)

    print(f"\nStep 7: Writing final JSON (with {num_ultimate_nodes} nodes) to {output_json_filepath}...")
    try:
        with open(output_json_filepath, 'w') as f_out_json:
            json.dump(output_json_structure, f_out_json, indent=4)
        print(f"Successfully wrote resulting JSON to {output_json_filepath}")
    except Exception as e:
        print(f"Error writing output JSON {output_json_filepath}: {e}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Filter 'nodes' list in a JSON based on IDX file and VCF content using pysam. "
                    "Understands latest 30B .idx (R,C) and older 26B/22B idx. "
                    "Optionally peeks into .dat to fill missing block details."
    )
    parser.add_argument("json_path", help="Path to the input JSON file.")
    parser.add_argument("idx_path", help="Path to the .idx file.")
    parser.add_argument("output_json_path", help="Path for the output filtered JSON file.")
    parser.add_argument("--vcf_file",
                        help="Optional: bgzipped & tabix-indexed VCF (.vcf.gz). Filters nodes whose region has no variants.",
                        default=None)
    parser.add_argument("--txt", dest="txt_output_path",
                        help="Optional: write final node IDs (one per line).", default=None)
    parser.add_argument("--dat", dest="dat_filepath",
                        help="Optional: matching .dat path (used to fill R/C or node_length when missing in old idx).",
                        default=None)
    args = parser.parse_args()

    filter_json_nodes_and_write(
        args.json_path,
        args.idx_path,
        args.output_json_path,
        args.vcf_file,
        args.txt_output_path,
        args.dat_filepath
    )

if __name__ == "__main__":
    main()
