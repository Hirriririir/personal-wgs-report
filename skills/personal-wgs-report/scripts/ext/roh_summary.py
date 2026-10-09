#!/usr/bin/env python3
"""ROH 汇总：每人 ≥1 Mb / ≥5 Mb 纯合片段的总长、段数、最长一段，F_ROH；样本在参考人群各亚群中的百分位 → work/ext/roh/summary.json
参照：一级表亲婚配的后代平均约 180 Mb（基因组的 1/16），二级表亲约 45 Mb。"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import popcfg  # noqa: E402
import wgs  # noqa: E402
import pandas as pd  # noqa: E402

D = "work/ext/roh"
cols = ["RG", "sample", "chrom", "start", "end", "len", "nmark", "qual"]
seg = pd.concat([pd.read_csv(f"{D}/pop.roh.RG.txt", sep="\t", names=cols)] +
                [pd.read_csv(f"{D}/{s}.roh.RG.txt", sep="\t", names=cols) for s in wgs.all_samples()
                 if os.path.getsize(f"{D}/{s}.roh.RG.txt") > 0])
seg = seg[seg["nmark"] >= 50]
AUTO = 2875e6
pop = {i: p for i, (_, p) in popcfg.kg_samples().items()}
ids = [x.strip() for x in open("work/ext/panel/pop_unrel.ids")] + wgs.all_samples()
rows = []
for s in ids:
    x = seg[seg["sample"] == s]
    rows.append(dict(sample=s, pop=pop.get(s, "sample"), sum1=round(x.loc[x["len"] >= 1e6, "len"].sum() / 1e6, 2),
                     n1=int((x["len"] >= 1e6).sum()), sum5=round(x.loc[x["len"] >= 5e6, "len"].sum() / 1e6, 2),
                     n5=int((x["len"] >= 5e6).sum()), longest=round(x["len"].max() / 1e6, 2) if len(x) else 0))
t = pd.DataFrame(rows)
t["froh"] = (t["sum1"] * 1e6 / AUTO).round(4)
t.to_csv(f"{D}/per_individual.tsv", sep="\t", index=False)
ref = t[t["pop"] != "sample"]
out = {"reference_median": ref.groupby("pop")[["sum1", "n1", "sum5", "n5", "longest"]].median().round(2).to_dict("index"), "samples": {}}
for s in wgs.all_samples():
    z = t[t["sample"] == s].iloc[0].to_dict()
    z["percentile_by_pop"] = {p: round(float((ref.loc[ref["pop"] == p, "sum1"] < z["sum1"]).mean() * 100), 1) for p in sorted(ref["pop"].unique())}
    out["samples"][s] = z
json.dump(out, open(f"{D}/summary.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(out["samples"], ensure_ascii=False, indent=1))
