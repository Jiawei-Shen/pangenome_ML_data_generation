#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
使用 LMDB 存储 GAM 文件中按 node_id 分段的脚本
每处理 milestone_step 条 reads 就提交一次事务并 flush
"""
import argparse
import gzip
import pickle
import time
import gc
import lmdb
from collections import defaultdict
import vg_pb2

# ─────────────────────────────────────────────────────────────────────────────
# 读取 varint
def read_varint(stream):
    value = 0
    shift = 0
    while True:
        b = stream.read(1)
        if not b:
            raise EOFError
        byte = b[0]
        value |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            return value
        shift += 7

# 检查文件是否 gzip
def file_is_gzip(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"

# 逐条产出 GAM record
def gam_record_iter(path, tag="GAM"):
    open_func = gzip.open if file_is_gzip(path) else open
    with open_func(path, "rb") as f:
        while True:
            try:
                count = read_varint(f)
            except EOFError:
                break
            if count == 0:
                continue
            try:
                tlen = read_varint(f)
                t = f.read(tlen).decode()
            except (EOFError, UnicodeDecodeError):
                break
            if t != tag:
                for _ in range(count - 1):
                    try:
                        skip = read_varint(f)
                        f.seek(skip, 1)
                    except EOFError:
                        break
                continue
            for _ in range(count - 1):
                try:
                    sz = read_varint(f)
                    yield f.read(sz)
                except EOFError:
                    break

# 处理单条对齐，返回 node_id->segments, 以及读数计数
def process_alignment(raw, wanted_nodes, chrom):
    segs = {}
    aln = vg_pb2.Alignment()
    aln.ParseFromString(raw)
    # 染色体过滤
    if chrom and not any(p.name == chrom for p in aln.refpos):
        return segs, 0
    seq = aln.sequence
    qual = aln.quality
    mapq = aln.mapping_quality
    offset = 0
    reads = 1
    for m in aln.path.mapping:
        nid = m.position.node_id
        if nid not in wanted_nodes:
            for e in m.edit:
                offset += max(e.from_length, len(e.sequence))
            continue
        node_off = m.position.offset
        strand = '-' if m.position.is_reverse else '+'
        parts = []
        bq = bytearray()
        for e in m.edit:
            if e.from_length:
                frag = seq[offset:offset+e.from_length]
                parts.append(frag.lower() if e.sequence else frag)
                bq.extend(qual[offset:offset+e.from_length])
                offset += e.from_length
            elif e.sequence:
                ins = len(e.sequence)
                frag = seq[offset:offset+ins]
                parts.append(frag.lower())
                bq.extend(qual[offset:offset+ins])
                offset += ins
        segs.setdefault(nid, []).append({
            'offset': node_off,
            'seq': ''.join(parts),
            'bq': bytes(bq),
            'rq': mapq,
            'strand': strand
        })
    return segs, reads

# 主流程
def run_pipeline(gam_path, stats_path, prefix, milestone_step, chrom_filter):
    print(f"加载统计文件: {stats_path}")
    with open(stats_path, 'rb') as f:
        stats = pickle.load(f)
    # 筛选 wanted nodes
    wanted = {int(n) for n, s in stats.items()
              if s['not_perfect'] > 1 and s['not_perfect']/(s['perfect']+s['not_perfect']) > 0.10}
    del stats
    gc.collect()
    print(f"筛选出 {len(wanted)} 个节点")

    # 打开 LMDB 环境
    lmdb_path = prefix + '.lmdb'
    env = lmdb.open(
        lmdb_path,
        map_size=600 * 1024**3,
        subdir=False,
        max_dbs=1,
        metasync=False,
        sync=False,
        map_async=True
    )
    db = env.open_db(b'segments', dupsort=True)
    print(f"LMDB 存储路径: {lmdb_path}")

    total = 0
    next_ms = milestone_step
    start = time.perf_counter()

    # 长事务写入
    txn = env.begin(write=True, db=db)
    for raw in gam_record_iter(gam_path):
        segdict, cnt = process_alignment(raw, wanted, chrom_filter)
        total += cnt
        # 写入 LMDB
        for nid, segs in segdict.items():
            key = str(nid).encode()
            for seg in segs:
                txn.put(key, pickle.dumps(seg, protocol=pickle.HIGHEST_PROTOCOL), dupdata=True)
        # milestone 提交
        if total >= next_ms:
            txn.commit()
            elapsed = time.perf_counter() - start
            print(f"{total} reads 完成 | 用时 {elapsed:.1f}s")
            next_ms += milestone_step
            txn = env.begin(write=True, db=db)

    # 最后提交并关闭
    txn.commit()
    env.close()
    elapsed = time.perf_counter() - start
    print("\n最终总结:")
    print(f"  总 reads: {total}\n  LMDB 文件: {lmdb_path}\n  总用时: {elapsed:.2f}s")

# CLI
if __name__ == '__main__':
    p = argparse.ArgumentParser(description='GAM->LMDB 读段提取')
    p.add_argument('gam_path')
    p.add_argument('stats_pickle')
    p.add_argument('output_prefix')
    p.add_argument('--milestone', type=int, default=10_000_000)
    p.add_argument('--chr', default='')
    args = p.parse_args()
    run_pipeline(
        args.gam_path,
        args.stats_pickle,
        args.output_prefix,
        args.milestone,
        args.chr
    )
