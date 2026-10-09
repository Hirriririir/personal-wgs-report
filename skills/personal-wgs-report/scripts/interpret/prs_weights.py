#!/usr/bin/env python3
"""PGS Catalog 评分文件（GRCh38 harmonized）→ plink2 --score 用的权重表。只用常染色体 SNV；链方向分不清的回文位点（A/T、C/G）去掉。

  prs_weights.py build <评分文件.txt.gz> <输出前缀>
      有另一等位基因（other_allele / hm_inferOtherAllele）的行 ≥90%：写 <前缀>.weights.tsv（ID=chr:pos:ref:alt，两种排列都写）
      否则（只给了效应等位基因，如部分 TPMI 评分）：写 <前缀>.weights_pos.tsv（ID=chr:pos），之后用 1000 Genomes 补全 REF/ALT
      打印 full 或 pos
  prs_weights.py finish-pos <输出前缀>
      读 <前缀>.used.pvar（1000 Genomes 上实际用到的位点及其 REF/ALT）→ 写 <前缀>.weights.tsv（完整 ID，去掉回文位点）
  prs_weights.py union <输出 sites.tsv> <ID ...>
      所有评分在参考人群上实际用到的位点取并集（按染色体、位置排序），供从 gVCF 取基因型
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import csv
import gzip
import os
import sys

PAL = ({"A", "T"}, {"C", "G"})


def rows_of(path):
    with gzip.open(path, "rt") as fh:
        yield from csv.DictReader((x for x in fh if not x.startswith("#")), delimiter="\t")


def build(path, out):
    full, pos, total, with_oa = [], [], 0, 0
    for r in rows_of(path):
        c, p, ea, w = r.get("hm_chr", ""), r.get("hm_pos", ""), r.get("effect_allele", ""), r.get("effect_weight", "")
        if not c or not p or not w or c in ("X", "Y", "MT", "XY") or len(ea) != 1 or ea not in "ACGT":
            continue
        total += 1
        oa = r.get("hm_inferOtherAllele") or r.get("other_allele") or ""
        if len(oa) == 1 and oa in "ACGT" and "/" not in oa:
            with_oa += 1
            if {ea, oa} not in PAL:
                full.append((c, int(p), ea, oa, w))
        pos.append((c, int(p), ea, w))
    if total and with_oa / total >= 0.9:
        seen = set()
        with open(out + ".weights.tsv", "w") as o:
            o.write("ID\tA1\tW\n")
            for c, p, ea, oa, w in full:
                for ref, alt in ((oa, ea), (ea, oa)):
                    vid = f"{c}:{p}:{ref}:{alt}"
                    if vid not in seen:
                        seen.add(vid)
                        o.write(f"{vid}\t{ea}\t{w}\n")
        print("full")
        print(f"{out}: {len(full)} 个 SNV 权重（共 {total} 行可用）", file=sys.stderr)
    else:
        seen = set()
        with open(out + ".weights_pos.tsv", "w") as o:
            o.write("ID\tA1\tW\n")
            for c, p, ea, w in pos:
                k = f"{c}:{p}"
                if k not in seen:
                    seen.add(k)
                    o.write(f"{k}\t{ea}\t{w}\n")
        print("pos")
        print(f"{out}: 只有效应等位基因，按位置匹配 {len(seen)} 个位点", file=sys.stderr)


def finish_pos(out):
    w = {}
    for line in open(out + ".weights_pos.tsv"):
        if not line.startswith("ID"):
            k, a1, wt = line.split()
            w[k] = (a1, wt)
    n = 0
    with open(out + ".used.pvar") as pv, open(out + ".weights.tsv", "w") as ow:
        ow.write("ID\tA1\tW\n")
        for line in pv:
            if line.startswith("#"):
                continue
            c, p, vid, ref, alt = line.split()[:5]
            a1, wt = w[vid]
            if {ref, alt} in PAL or a1 not in (ref, alt):
                continue
            ow.write(f"{c}:{p}:{ref}:{alt}\t{a1}\t{wt}\n")
            n += 1
    print(f"{out}: 补全 REF/ALT 后 {n} 个位点", file=sys.stderr)


def union(out, ids):
    used = set()
    for i in ids:
        p = f"work/prs/{i}.ref.sscore.vars"
        if os.path.exists(p):
            used |= set(open(p).read().split())
    order = {str(c): k for k, c in enumerate(range(1, 23))}
    rows = sorted((order.get(v.split(":")[0], 99), int(v.split(":")[1]), v) for v in used)
    with open(out, "w") as o:
        for _, p, v in rows:
            c, _, r, a = v.split(":")
            o.write(f"{c}\t{p}\t{v}\t{r}\t{a}\n")
    print(f"并集位点 {len(rows)}", file=sys.stderr)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "build":
        build(sys.argv[2], sys.argv[3])
    elif cmd == "finish-pos":
        finish_pos(sys.argv[2])
    elif cmd == "union":
        union(sys.argv[2], sys.argv[3:])
    else:
        sys.exit(__doc__)
