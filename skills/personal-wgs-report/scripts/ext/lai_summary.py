#!/usr/bin/env python3
"""局部祖源汇总：样本和对照个体全基因组各参考组的比例（按物理长度加权、两条单倍型平均）→ work/ext/lai/summary.<s>.json
另存样本的片段表 segments.<s>.tsv（画染色体涂色图用）。用法：lai_summary.py <sample>"""
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402
import pandas as pd  # noqa: E402

s = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
D = "work/ext/lai"
groups = {}  # 编号 → 组名（读 msp 第一行 “#Subpopulation order/codes: North=0	South=1”）
msp = []
for p in sorted(glob.glob(f"{D}/chr/rfmix.{s}.chr*.msp.tsv")):
    with open(p) as f:
        groups = {k: g for g, k in re.findall(r"(\S+?)=(\d+)", f.readline())}
        cols = f.readline().lstrip("#").rstrip("\n").split("\t")
    msp.append(pd.read_csv(p, sep="\t", skiprows=2, header=None, names=cols))
if not msp:
    sys.exit("没有 RFMix 结果")
msp = pd.concat(msp, ignore_index=True).copy()
msp["len"] = msp["epos"] - msp["spos"]
haps = [c for c in msp.columns if re.search(r"\.[01]$", c)]
samples = sorted({h.rsplit(".", 1)[0] for h in haps})
pop = dict(line.split() for line in open(f"{D}/query_1kg.tsv"))
frac = {}
for x in samples:
    tot = 2 * msp["len"].sum()
    frac[x] = {groups[k]: round(float(((msp[f"{x}.0"] == int(k)) * msp["len"]).sum() + ((msp[f"{x}.1"] == int(k)) * msp["len"]).sum()) / tot, 4)
               for k in groups}
ctrl = {}
for x, f in frac.items():
    if x != s:
        ctrl.setdefault(pop.get(x, "?"), []).append(f)
summary = {"sample": s, "fractions": frac[s], "groups": groups, "n_chrom": int(msp["chm"].nunique()),
           "controls": {p: {g: [round(v[g], 4) for v in vs] for g in groups.values()} for p, vs in ctrl.items()}}
json.dump(summary, open(f"{D}/summary.{s}.json", "w"), ensure_ascii=False, indent=1)
msp[["chm", "spos", "epos", "sgpos", "egpos", f"{s}.0", f"{s}.1"]].rename(columns={f"{s}.0": "hap1", f"{s}.1": "hap2"}).to_csv(
    f"{D}/segments.{s}.tsv", sep="\t", index=False)
print(json.dumps({k: v for k, v in summary.items() if k != "controls"}, ensure_ascii=False))
for p, d in summary["controls"].items():
    print(p, {g: round(sum(v) / len(v), 3) for g, v in d.items()})
