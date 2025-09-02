#!/usr/bin/env python3
"""
single_node_pileup.py

Build Pansoma-style 6-channel tensors for ONE node from .idx/.dat, using streamed
sequence+AF from merged.json (or merged.json.gz).

• Same encodings and tensor layout as your wave script.
• Channel 6 (AF) policy:
    - Accepts per-base digits "0".."7" (string) OR float array → binned into 8 levels:
      0–1e-6→0, 1e-6–1e-5→1, 1e-5–1e-4→2, 1e-4–1e-3→3, 1e-3–1e-2→4, 1e-2–0.1→5, 0.1–0.5→6, 0.5–1.0→7
    - If AF missing: zeros, except the variant column is set to 3
• Outputs .npy files named: {pos}_{type}_{ref}_{alt}.npy  + variant_summary.json

Usage:
  python single_node_pileup.py DAT.dat IDX.idx merged.json[.gz] --node_id 12345 -o OUT_DIR \
      --variant_type snp --min_af 0.1 --min_variants 3 --min_allele_bq 10 --view 10
"""

import argparse
import gzip
import json
import os
import re
import struct
import sys
import time
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Logging

def ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def log(msg: str) -> None:
    print(f"[{ts()}] {msg}", file=sys.stderr, flush=True)

# ─────────────────────────────────────────────────────────────────────────────
# Formats / constants (kept identical to your multi-wave script)

RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size

BASE_TO_INDEX = {'A': 20, 'C': 30, 'G': 50, 'T': 70, 'N': 10, '*': 90, '_PADDING_': 0}
PADDING_BASE_INDEX = 0

CIGAR_OP_TO_INDEX = {'M': 10, 'N': 20, 'S': 30, 'I': 40, 'D': 50, 'H': 60, 'P': 70, '=': 80, 'X': 90, '_PADDING_': 0}
CIGAR_PADDING_INDEX = 0

INDEX_TO_BASE_FOR_VIEW = {20: 'A', 30: 'C', 50: 'G', 70: 'T', 10: 'N', 90: '*', 0: '0'}

TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200
DEFAULT_QUALITY_PADDING = 0
DEFAULT_MAPPING_QUALITY_PADDING = -1
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = -1

# ─────────────────────────────────────────────────────────────────────────────
# Helpers

def reverse_complement(sequence: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(comp)[::-1]

def calculate_window_start(variant_pos: int, window_size: int) -> int:
    center_index = window_size // 2
    return variant_pos - center_index

def decode_cigar_to_int_ops(cigar_string: str) -> List[Tuple[int, str]]:
    if not cigar_string or cigar_string == '*':
        return []
    try:
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        sys.stderr.write(f"Warning: Could not parse CIGAR '{cigar_string}': {e}\n")
        return []

def _open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

# ─────────────────────────────────────────────────────────────────────────────
# IDX / JSON readers for a single node

def find_node_in_idx(idx_path: str, target_node_id: int) -> Optional[Tuple[int, int]]:
    """
    Return (dat_offset, n_records) for the node_id, or None if not found.
    IDX layout assumed:
      header: uint32 N
      entries: struct '<I Q I I H' → (node_id, dat_offset, _, n_records, _)
    """
    with open(idx_path, 'rb') as f:
        hdr = f.read(4)
        if len(hdr) < 4:
            raise RuntimeError("Index file too small")
        n = struct.unpack('<I', hdr)[0]
        for i in range(n):
            rec = f.read(22)
            if len(rec) < 22:
                break
            node_id, dat_off, _u32, n_records, _u16 = struct.unpack('<I Q I I H', rec)
            if node_id == target_node_id:
                return dat_off, n_records
    return None

def af_float_to_bin(x: float) -> int:
    try:
        x = float(x)
    except Exception:
        return 0
    if x <= 0.0:          return 0
    if x < 1e-6:          return 0
    if x < 1e-5:          return 1
    if x < 1e-4:          return 2
    if x < 1e-3:          return 3
    if x < 1e-2:          return 4
    if x < 0.1:           return 5
    if x < 0.5:           return 6
    return 7

def normalize_af_to_bins_compact(af_field, expected_len: int) -> Optional[np.ndarray]:
    """
    Return np.uint8 array length==expected_len with values 0..7, or None if AF truly missing.
    Accept list[float] or digit string "0".."7" per base.
    """
    if af_field is None:
        return None
    if isinstance(af_field, str):
        if not af_field:
            return None
        a = np.frombuffer(af_field.encode('ascii'), dtype=np.uint8)
        a = np.where((a >= ord('0')) & (a <= ord('7')), a - ord('0'), 0).astype(np.uint8)
        if a.size < expected_len:
            out = np.zeros(expected_len, dtype=np.uint8)
            out[:a.size] = a
            return out
        return a[:expected_len].copy()
    if isinstance(af_field, list):
        out = np.zeros(expected_len, dtype=np.uint8)
        upto = min(expected_len, len(af_field))
        out[:upto] = np.fromiter((af_float_to_bin(v) for v in af_field[:upto]), count=upto, dtype=np.uint8)
        return out
    return None

def load_single_node_seq_and_af(merged_path: str, node_id: int) -> Tuple[Optional[str], Optional[np.ndarray]]:
    """
    Stream merged JSON to find the one node_id. Supports:
      • [ {...}, {...}, ... ]
      • {"nodes":[ {...}, ... ], ...}
    Uses a light brace-tracking to avoid loading whole file.
    """
    with _open_maybe_gzip(merged_path, "rt") as f:
        # scan to first '['
        while True:
            ch = f.read(1)
            if not ch:
                return None, None
            if ch == '[':
                break

        buf = []
        depth = 0
        in_str = False
        esc = False

        def _emit_object(s: str):
            try:
                return json.loads(s)
            except Exception:
                return None

        while True:
            ch = f.read(1)
            if not ch:
                break

            if depth == 0:
                if ch.isspace() or ch == ',':
                    continue
                if ch == ']':
                    break
                if ch != '{':
                    continue
                buf = ['{']
                depth = 1
                in_str = False
                esc = False
                continue

            buf.append(ch)

            if in_str:
                if esc:
                    esc = False
                elif ch == '\\':
                    esc = True
                elif ch == '"':
                    in_str = False
            else:
                if ch == '"':
                    in_str = True
                elif ch == '{':
                    depth += 1
                elif ch == '}':
                    depth -= 1
                    if depth == 0:
                        obj = _emit_object(''.join(buf))
                        buf = []
                        if isinstance(obj, dict):
                            nid_raw = obj.get('node_id')
                            if nid_raw is not None:
                                try:
                                    nid = int(str(nid_raw))
                                except Exception:
                                    nid = None
                                if nid == node_id:
                                    seq = str(obj.get('sequence', '')).upper()
                                    if not seq:
                                        return None, None
                                    af_bins = normalize_af_to_bins_compact(obj.get('genomead_af'), len(seq))
                                    return seq, af_bins
        return None, None

# ─────────────────────────────────────────────────────────────────────────────
# Variant extraction + pileup row builders

def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence, read_quality_values, read_cigar_ops_decoded,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None, expected_ref_allele_for_indel=None):
    current_node_pos = read_offset_on_node
    current_read_pos = 0
    for length, op in read_cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            if current_node_pos <= target_node_pos < current_node_pos + length:
                offset_in_block = target_node_pos - current_node_pos
                read_idx = current_read_pos + offset_in_block
                if read_idx < len(read_sequence):
                    allele = read_sequence[read_idx].upper()
                    quality = read_quality_values[read_idx] if read_idx < len(read_quality_values) else 0
                    if expected_var_type in ('I', 'D'):
                        return "REF_STATE_FOR_INDEL", quality
                    return allele, quality
                return None, None
            current_node_pos += length
            current_read_pos += length
        elif op == 'I':
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                if current_read_pos + length <= len(read_sequence):
                    qualities = read_quality_values[current_read_pos: current_read_pos + length]
                    mean_quality = sum(qualities) / len(qualities) if qualities else 0.0
                    return read_sequence[current_read_pos: current_read_pos + length].upper(), mean_quality
                return None, None
            current_read_pos += length
        elif op == 'D':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I':
                    return "OTHER_FOR_INDEL", None
                if expected_var_type == 'D':
                    if 0 <= current_node_pos < len(node_sequence) and current_node_pos + length <= len(node_sequence):
                        deleted_seq_in_ref_context = node_sequence[current_node_pos: current_node_pos + length]
                        if deleted_seq_in_ref_context == expected_ref_allele_for_indel:
                            return "*", None
                        else:
                            return "OTHER_FOR_INDEL", None
                    return "OTHER_FOR_INDEL", None
                return "*", None
            current_node_pos += length
        elif op == 'S':
            current_read_pos += length
        elif op == 'N':
            current_node_pos += length
        if current_node_pos > target_node_pos + 1 and not (
            expected_var_type == 'I' and (current_node_pos - 1) <= target_node_pos):
            break
    return None, None

def detect_variants_from_cigar(offset_on_node, cigar_ops_decoded, read_sequence, node_sequence):
    variants = []
    node_pos, read_pos = offset_on_node, 0
    node_seq_len, read_seq_len = len(node_sequence), len(read_sequence)
    for length, op in cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            for i in range(length):
                cur_node_p, cur_read_p = node_pos + i, read_pos + i
                if cur_node_p < node_seq_len and cur_read_p < read_seq_len:
                    node_base = node_sequence[cur_node_p].upper()
                    read_base = read_sequence[cur_read_p].upper()
                    if node_base != read_base and op != '=':
                        variants.append((cur_node_p, 'X', read_base, node_base))
                else:
                    break
            node_pos += length
            read_pos += length
        elif op == 'I':
            inserted_sequence = read_sequence[read_pos: read_pos + length].upper()
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0
            ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if 0 <= ref_anchor_pos < node_seq_len else "*"
            variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor))
            read_pos += length
        elif op == 'D':
            deleted_sequence_from_ref = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= node_seq_len else ""
            if deleted_sequence_from_ref:
                variants.append((node_pos, 'D', "*", deleted_sequence_from_ref))
            node_pos += length
        elif op == 'S':
            read_pos += length
        elif op == 'N':
            node_pos += length
    return variants

def get_read_tensor_rows_in_window(segment_cigar_ops, segment_offset_on_node,
                                   segment_read_sequence, segment_quality_values,
                                   mapping_quality,
                                   window_start_node, tensor_win_size, node_len):
    bases = [PADDING_BASE_INDEX] * tensor_win_size
    quals = [DEFAULT_QUALITY_PADDING] * tensor_win_size
    mapqs = [DEFAULT_MAPPING_QUALITY_PADDING] * tensor_win_size
    cigar_ops_indices = [CIGAR_PADDING_INDEX] * tensor_win_size
    node_pos, read_pos = segment_offset_on_node, 0
    read_seq_len = len(segment_read_sequence)
    qual_len = len(segment_quality_values)
    for L, op in segment_cigar_ops:
        op_idx = CIGAR_OP_TO_INDEX.get(op, CIGAR_PADDING_INDEX)
        if node_pos >= window_start_node + tensor_win_size and read_pos > 0:
            break
        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = node_pos + i, read_pos + i
                win_idx = n_aln - window_start_node
                if 0 <= win_idx < tensor_win_size:
                    cigar_ops_indices[win_idx] = op_idx
                    mapqs[win_idx] = mapping_quality
                    if r_aln < read_seq_len:
                        base_char = segment_read_sequence[r_aln].upper()
                        bases[win_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])
                        if r_aln < qual_len:
                            quals[win_idx] = segment_quality_values[r_aln]
            node_pos += L
            read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                win_idx = (node_pos + i) - window_start_node
                if 0 <= win_idx < tensor_win_size:
                    cigar_ops_indices[win_idx] = op_idx
                    mapqs[win_idx] = mapping_quality
                    bases[win_idx] = BASE_TO_INDEX['*']
                    quals[win_idx] = DEFAULT_QUALITY_PADDING
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L
        if read_pos >= read_seq_len:
            break
    return bases, quals, mapqs, cigar_ops_indices

# ─────────────────────────────────────────────────────────────────────────────
# View helper (text)

def display_pileup_header(node_id: int, node_len: int, n_reads: int):
    print(f"\n=== Node {node_id} (len={node_len}, reads={n_reads}) ===")

def display_variant_marker(node_seq: str, v_pos: int, v_type: str):
    window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
    window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)
    ref_chars = []
    for j in range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE):
        ref_chars.append(node_seq[j] if 0 <= j < len(node_seq) else '0')
    marker_pos = v_pos - window_start_pos
    marker_line = [' '] * TENSOR_WINDOW_SIZE
    if 0 <= marker_pos < TENSOR_WINDOW_SIZE:
        marker_line[marker_pos] = '^'
    print(f"\n  Ref    : {''.join(ref_chars)}")
    print(f"  Marker : {''.join(marker_line)}")

# ─────────────────────────────────────────────────────────────────────────────
# Main single-node routine

def process_single_node(dat_path: str,
                        idx_path: str,
                        merged_json_path: str,
                        node_id: int,
                        out_dir: str,
                        variant_type: str,
                        min_af: float,
                        min_variants: int,
                        min_allele_bq: float,
                        view_limit: Optional[int],
                        max_reads_view: int,
                        mapq_min: int,
                        max_reads_guard: int) -> Tuple[int, int]:
    """
    Return (n_variants_emitted, n_reads_consumed)
    """
    # 1) Locate node in IDX
    loc = find_node_in_idx(idx_path, node_id)
    if not loc:
        raise SystemExit(f"Node {node_id} not found in IDX.")
    dat_offset, n_records = loc

    # 2) Get sequence + AF bins
    seq, af_bins = load_single_node_seq_and_af(merged_json_path, node_id)
    if not seq:
        raise SystemExit(f"Node {node_id} not found in merged JSON, or empty sequence.")
    node_seq = seq.upper()
    node_len = len(node_seq)

    # 3) Read aligned segments for this node from .dat
    aligned = []
    with open(dat_path, 'rb') as df:
        # replicate your offset behavior (+10 bytes header per node block)
        df.seek(dat_offset + 10)
        bulk = df.read(n_records * RECORD_SIZE)
        if len(bulk) < n_records * RECORD_SIZE:
            n_records = len(bulk) // RECORD_SIZE
            bulk = bulk[: n_records * RECORD_SIZE]

        for (off_from_file, raw_seq, raw_qual, raw_cigar, mapq_val, strand_byte) in RECORD_STRUCT.iter_unpack(bulk):
            if mapq_val < mapq_min:
                continue
            try:
                rseq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                rqual = list(raw_qual.rstrip(b'\0'))
                cigar = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                continue
            if not rseq or len(rseq) != len(rqual):
                continue
            ops = decode_cigar_to_int_ops(cigar)
            if not ops and cigar != '*':
                continue

            cur_seq = rseq
            cur_qual = rqual
            cur_ops = ops
            cur_off = off_from_file
            if strand == '-':
                cur_seq = reverse_complement(rseq)
                cur_qual = rqual[::-1]
                cur_ops = [op for op in reversed(ops)] if ops else []
                aln_span = len(cur_seq)
                cur_off = node_len - aln_span - off_from_file
                if cur_off < 0:
                    continue

            aligned.append({
                "offset_on_node": cur_off,
                "read_sequence": cur_seq,
                "processed_quality_values": cur_qual,
                "cigar_ops": cur_ops,
                "original_cigar_str": cigar,
                "strand": strand,
                "mapping_quality": mapq_val
            })

    n_reads = len(aligned)
    if n_reads == 0:
        log(f"Node {node_id}: no usable reads.")
        return 0, 0
    if n_reads > max_reads_guard:
        log(f"Node {node_id}: too many reads ({n_reads} > {max_reads_guard}), skipping.")
        return 0, 0

    # 4) Variant discovery
    candidate = defaultdict(int)
    for seg in aligned:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_seq):
            candidate[(v_pos, v_type, v_ref, v_alt)] += 1

    # 5) Per-variant counting + filters + tensor build
    out_node_dir = os.path.join(out_dir, str(node_id))
    os.makedirs(out_node_dir, exist_ok=True)

    emitted = 0
    summary = []

    # Optional view
    if view_limit is not None:
        display_pileup_header(node_id, node_len, n_reads)

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate.items():
        if variant_type == 'snp' and v_type != 'X':
            continue
        if variant_type == 'indel' and v_type not in ('I', 'D'):
            continue

        alt = v_alt_from_cigar
        ref = v_ref_from_cigar
        ref_for_indel_ctx = None
        if v_type == 'D':
            alt = "*"
            if 0 <= v_pos < node_len:
                ref = node_seq[v_pos]
            ref_for_indel_ctx = v_ref_from_cigar
        elif v_type == 'I':
            ref = node_seq[v_pos] if 0 <= v_pos < node_len else "*"
            ref_for_indel_ctx = ref

        alt_cnt = ref_cnt = other_cnt = cov = 0
        alt_bqs = []

        for seg in aligned:
            allele, bq = get_allele_from_read_at_node_pos(
                seg["offset_on_node"], seg["read_sequence"], seg["processed_quality_values"], seg["cigar_ops"],
                v_pos, node_seq, v_type, ref_for_indel_ctx)
            if allele is not None:
                cov += 1
                if allele == alt:
                    alt_cnt += 1
                    if bq is not None:
                        alt_bqs.append(bq)
                elif allele == ref or (v_type in ('I', 'D') and allele == "REF_STATE_FOR_INDEL"):
                    ref_cnt += 1
                else:
                    other_cnt += 1

        if alt_cnt < min_variants:
            continue

        if v_type == 'X':
            tmp_af = alt_cnt / cov if cov > 0 else 0.0
            if tmp_af < min_af:
                continue
            mean_bq = sum(alt_bqs) / len(alt_bqs) if alt_bqs else 0.0
            if mean_bq < min_allele_bq:
                continue
        else:
            mean_bq = sum(alt_bqs) / len(alt_bqs) if alt_bqs else 0.0

        cur_af = alt_cnt / cov if cov > 0 else 0.0

        key = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"
        win_center = v_pos + 1 if v_type == 'I' else v_pos
        win_start = calculate_window_start(win_center, TENSOR_WINDOW_SIZE)

        # Row 0 (reference bases)
        ref_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i, pos in enumerate(range(win_start, win_start + TENSOR_WINDOW_SIZE)):
            if 0 <= pos < node_len:
                ref_row[i] = BASE_TO_INDEX.get(node_seq[pos].upper(), BASE_TO_INDEX['N'])

        # Row 0 for AF bins
        af_row = np.zeros(TENSOR_WINDOW_SIZE, dtype=np.uint8)
        if af_bins is not None:
            for i, pos in enumerate(range(win_start, win_start + TENSOR_WINDOW_SIZE)):
                if 0 <= pos < node_len and pos < af_bins.shape[0]:
                    af_row[i] = af_bins[pos]
        else:
            v_idx = v_pos - win_start
            if 0 <= v_idx < TENSOR_WINDOW_SIZE:
                af_row[v_idx] = 3

        # Allocate channels
        H, W = 1 + TENSOR_MAX_READ_ROWS, TENSOR_WINDOW_SIZE
        ch1 = np.full((H, W), PADDING_BASE_INDEX, dtype=np.int8)
        ch2 = np.full((H, W), DEFAULT_QUALITY_PADDING, dtype=np.int8)
        ch3 = np.full((H, W), MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
        ch4 = np.full((H, W), DEFAULT_MAPPING_QUALITY_PADDING, dtype=np.int8)
        ch5 = np.full((H, W), CIGAR_PADDING_INDEX, dtype=np.int8)
        ch6 = np.zeros((H, W), dtype=np.int8)

        ch1[0, :] = np.asarray(ref_row, dtype=np.int8)
        ch6[0, :] = af_row.astype(np.int8, copy=False)

        reads_added = 0
        v_idx = v_pos - win_start
        for seg in aligned:
            if reads_added >= TENSOR_MAX_READ_ROWS:
                break
            base_row, qual_row, mapq_row, cigar_row = get_read_tensor_rows_in_window(
                seg["cigar_ops"], seg["offset_on_node"],
                seg["read_sequence"], seg["processed_quality_values"],
                max(0, min(int(seg["mapping_quality"]), 127)),
                win_start, TENSOR_WINDOW_SIZE, node_len)

            if any(b != PADDING_BASE_INDEX for b in base_row):
                r = 1 + reads_added
                ch1[r, :] = np.asarray(base_row, dtype=np.int8)
                ch2[r, :] = np.asarray(qual_row, dtype=np.int8)
                ref_vec = ch1[0, :].astype(np.int16)
                read_vec = ch1[r, :].astype(np.int16)
                flags = np.full(W, MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
                mask = (read_vec != PADDING_BASE_INDEX) & (ref_vec != PADDING_BASE_INDEX)
                flags[mask] = (read_vec[mask] != ref_vec[mask]).astype(np.int8)
                if 0 <= v_idx < W and mask[v_idx] and flags[v_idx] == 1:
                    flags[v_idx] = 5
                ch3[r, :] = flags
                ch4[r, :] = np.asarray(mapq_row, dtype=np.int8)
                ch5[r, :] = np.asarray(cigar_row, dtype=np.int8)
                ch6[r, :] = ch6[0, :]
                reads_added += 1

        # Save tensor
        npy_name = f"{key}.npy"
        np.save(os.path.join(out_node_dir, npy_name), np.stack([ch1, ch2, ch3, ch4, ch5, ch6], axis=0))

        summary.append({
            "variant_key": key,
            "tensor_file": npy_name,
            "alt_allele_count": alt_cnt,
            "ref_allele_count_at_locus": ref_cnt,
            "other_allele_count_at_locus": other_cnt,
            "coverage_at_locus": cov,
            "alt_allele_frequency": round(cur_af, 4),
            "mean_alt_allele_base_quality": round(mean_bq, 2)
        })
        emitted += 1

        # Optional view
        if view_limit is not None:
            display_variant_marker(node_seq, v_pos, v_type)
            print(f"  Variant: pos={v_pos}, type={v_type}, ref={ref}, alt={alt}, "
                  f"AC={alt_cnt}, RC={ref_cnt}, OC={other_cnt}, DP={cov}, AF={cur_af:.3f}, meanBQ={mean_bq:.1f}")
            if view_limit > 0:
                view_limit -= 1
                if view_limit == 0:
                    print("\n  ... (view limit reached)")
                    view_limit = None  # stop printing more

    # Write summary
    with open(os.path.join(out_node_dir, "variant_summary.json"), "w") as fh:
        json.dump({"node_id": node_id, "node_length": node_len,
                   "variants_passing_af_filter": summary}, fh, indent=2)

    log(f"Node {node_id}: emitted {emitted} tensor(s) from {n_reads} reads.")
    return emitted, n_reads

# ─────────────────────────────────────────────────────────────────────────────
# CLI

def main():
    ap = argparse.ArgumentParser(
        description="Build tensors for ONE .idx/.dat node; JSON streamed; AF channel policy preserved.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("dat", help=".dat alignment file")
    ap.add_argument("idx", help=".idx index file")
    ap.add_argument("merged_json", help="Merged node JSON with sequence and AF; .gz allowed")
    ap.add_argument("-o", "--output", required=True, help="Base output directory")
    ap.add_argument("--node_id", type=int, required=True, help="Node ID to process")

    ap.add_argument("--variant_type", choices=["snp", "indel", "all"], default="all",
                    help="Which variant types to emit")
    ap.add_argument("--min_af", type=float, default=0.10, help="Min AF for SNP to pass")
    ap.add_argument("--min_variants", type=int, default=3, help="Min ALT count to pass")
    ap.add_argument("--min_allele_bq", type=float, default=10.0, help="Min mean base quality for ALT alleles")

    ap.add_argument("--mapq_min", type=int, default=10, help="Minimum MAPQ to accept a read")
    ap.add_argument("--max_reads_guard", type=int, default=10000, help="Skip node if reads exceed this")

    ap.add_argument("--view", nargs='?', const=10, type=int, metavar="N",
                    help="Print brief text view for up to N variants (omit N to default to 10)")
    ap.add_argument("--max_view_reads", type=int, default=20, help="(reserved) reads per variant for view")

    args = ap.parse_args()

    for p in (args.dat, args.idx, args.merged_json):
        if not os.path.isfile(p):
            sys.exit(f"Error: file not found: {p}")
    os.makedirs(args.output, exist_ok=True)

    process_single_node(
        dat_path=args.dat,
        idx_path=args.idx,
        merged_json_path=args.merged_json,
        node_id=args.node_id,
        out_dir=args.output,
        variant_type=args.variant_type,
        min_af=args.min_af,
        min_variants=args.min_variants,
        min_allele_bq=args.min_allele_bq,
        view_limit=args.view,
        max_reads_view=args.max_view_reads,
        mapq_min=args.mapq_min,
        max_reads_guard=args.max_reads_guard
    )

if __name__ == "__main__":
    main()
