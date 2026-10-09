#!/usr/bin/env python3
"""outgroup f3(样本, X; 外群)：样本和每个古今人群共享多少遗传漂移（越大 = 共同的遗传历史越多）。
ADMIXTOOLS qp3Pop，数据为 aadr_subset.py 写出的 work/ext/aadr/sub.*；外群默认姆布蒂人（Mbuti）。
输出 work/ext/aadr/f3.<s>.tsv（每个人群一行，带年代、经纬度、样本数）。用法：aadr_f3.py <sample>"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import popcfg  # noqa: E402
import wgs  # noqa: E402
import pandas as pd  # noqa: E402

s = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
cfg = popcfg.load()["aadr"]
O = "work/ext/aadr"
og = cfg.get("f3_outgroup", "Mbuti")
labels = sorted({x.split()[2] for x in open(f"{O}/sub.ind")})
samples = set(wgs.all_samples())
pops = [x for x in labels if x not in samples and x not in (og, "Chimp")]
open(f"{O}/f3.{s}.pops", "w").write("".join(f"{s} {x} {og}\n" for x in pops))
open(f"{O}/f3.{s}.par", "w").write(f"genotypename: {O}/sub.geno\nsnpname: {O}/sub.snp\nindivname: {O}/sub.ind\n"
                                   f"popfilename: {O}/f3.{s}.pops\nhashcheck: NO\ninbreed: NO\nnumchrom: 22\n")
img = os.environ.get("IMG_ADMIXTOOLS")
with open(f"{O}/f3.{s}.out", "w") as out:
    subprocess.run(["scripts/dr.sh", img, "qp3Pop", "-p", f"{O}/f3.{s}.par"], stdout=out, stderr=subprocess.STDOUT, check=True)
rows = []
for line in open(f"{O}/f3.{s}.out"):
    if line.startswith(" result:"):
        f = line.split()
        rows.append(dict(group=f[2], f3=float(f[4]), se=float(f[5]), z=float(f[6]), nsnp=int(f[7])))
f3 = pd.DataFrame(rows)
meta = pd.read_csv(f"{O}/sub.meta.tsv", sep="\t", low_memory=False)
for c in ("lat", "lon"):
    meta[c] = pd.to_numeric(meta[c], errors="coerce")
g = meta.groupby("label").agg(n=("gid", "size"), bp=("bp", "median"), lat=("lat", "median"), lon=("lon", "median"),
                              country=("country", "first"), locality=("locality", "first")).reset_index()
f3 = f3.merge(g, left_on="group", right_on="label", how="left").drop(columns=["label"])
f3.sort_values("f3", ascending=False).to_csv(f"{O}/f3.{s}.tsv", sep="\t", index=False)
anc = f3[(f3["bp"] > 0) & (f3["nsnp"] >= 30000)].sort_values("f3", ascending=False)
print("古代人群（≥3 万位点）前 10：")
print(anc.head(10)[["group", "bp", "n", "f3", "se"]].to_string(index=False))
