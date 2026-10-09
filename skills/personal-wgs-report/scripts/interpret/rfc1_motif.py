#!/usr/bin/env python3
"""RFC1（CANVAS）重复单元组成：从 ExpansionHunter 的 BAMlet 里取分配给 RFC1 的读段，统计五核苷酸单元。

ExpansionHunter 用简并基序 (AARRG)* 计数，分不清良性 AAAAG / AAAGG 扩增和致病 AAGGG / ACAGG 扩增。
这里把读段里的五聚体按“旋转 + 反向互补”归一，看全重复读段（IRR，≥80% 是同一种单元）以哪种单元为主。
用法：rfc1_motif.py <sample>
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import collections
import json
import sys

import pysam

s = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
comp = str.maketrans("ACGT", "TGCA")

def canon(k):
    rc = k.translate(comp)[::-1]
    rots = [k[i:] + k[:i] for i in range(len(k))] + [rc[i:] + rc[:i] for i in range(len(rc))]
    return min(rots)

NAMES = {canon("AAAAG"): "AAAAG（参考，良性）", canon("AAAGG"): "AAAGG（良性多态 / 极长时有争议）", canon("AAGGG"): "AAGGG（致病，经典）",
         canon("ACAGG"): "ACAGG（致病，亚洲 / 大洋洲报道）", canon("AGAGG"): "AGAGG（意义不明）", canon("AAGAG"): "AAGAG（良性）",
         canon("AGGGC"): "AGGGC（致病，罕见）", canon("AAGGC"): "AAGGC（报道致病）"}

bam = pysam.AlignmentFile(f"work/{s}/str/{s}.eh_realigned.sorted.bam")
irr = collections.Counter()
partial = collections.Counter()
n = 0
for r in bam.fetch(until_eof=True):
    tags = dict(r.get_tags())
    if not any("RFC1" in str(v) for v in tags.values()):
        continue
    seq = r.query_sequence or ""
    if len(seq) < 100:
        continue
    n += 1
    c = collections.Counter(canon(seq[i:i + 5]) for i in range(0, len(seq) - 4))
    top, cnt = c.most_common(1)[0]
    frac = cnt / max(1, len(seq) - 4)
    if frac >= 0.8:
        irr[top] += 1
    elif frac >= 0.3:
        partial[top] += 1
res = {"sample": s, "reads_assigned": n,
       "in_repeat_reads": {NAMES.get(k, k): v for k, v in irr.most_common()},
       "partial_repeat_reads": {NAMES.get(k, k): v for k, v in partial.most_common(6)}}
json.dump(res, open(f"work/{s}/str/{s}.rfc1_motif.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(res, ensure_ascii=False, indent=1))
