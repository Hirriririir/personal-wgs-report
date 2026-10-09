#!/usr/bin/env python3
"""从 AADR 1240K 里挑出参考地区（assets/population/<POP>.json 的 aadr.countries）的古今个体 + 外群，加上样本，写成小的 TGENO 数据集。
质控：ASSESSMENT 为 PASS 类；1240K 常染色体覆盖位点 ≥ 20,000；同一个体多条记录只留覆盖最多的；去掉 Ignore_ / 离群（-o）/ Dup 等标签。
输出 work/ext/aadr/sub.{geno,snp,ind} + sub.meta.tsv（群体、年代、经纬度、国家、覆盖位点、Y / mt 单倍群）
用法：aadr_subset.py [sample ...]（默认工作目录里所有已编码的样本）"""
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import popcfg  # noqa: E402
import wgs  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from aadr_lib import TGeno, pack, write_tgeno  # noqa: E402

cfg = popcfg.load()["aadr"]
A = "ref/aadr/v66.p1_1240K.aadr"
O = "work/ext/aadr"
t = TGeno(A + ".patch.PUB")
an = pd.read_csv(A + ".PUB.anno", sep="\t", dtype=str, low_memory=False)
c = list(an.columns)
an = an.rename(columns={c[0]: "gid", c[2]: "iid", c[10]: "bp", c[14]: "group", c[15]: "locality", c[16]: "country",
                        c[17]: "lat", c[18]: "lon", c[21]: "dtype", c[26]: "snps", c[34]: "y_term", c[35]: "y_isogg",
                        c[38]: "mt", c[47]: "assess", c[5]: "pub"})
an["bp"] = pd.to_numeric(an["bp"], errors="coerce")
an["snps"] = pd.to_numeric(an["snps"], errors="coerce").fillna(0)
ind = pd.DataFrame(t.ind, columns=["gid", "sex", "label"])
ind["row"] = range(len(ind))
df = ind.merge(an[["gid", "iid", "bp", "group", "locality", "country", "lat", "lon", "dtype", "snps", "y_term", "y_isogg", "mt",
                   "assess", "pub"]], on="gid", how="left")
pas = df["assess"].fillna("").str.contains("PASS", case=False)
lab = df["label"]
bad = lab.str.contains("Ignore|_o$|-o$|-o[A-Z]|Dup|_dup|contam|Outlier|outlier|_rel|DontUse", regex=True)
keep = (df["country"].isin(cfg["countries"]) | lab.isin(cfg["outgroups"])) & pas & (df["snps"] >= 20000) & ~bad
d = df[keep].copy().sort_values("snps", ascending=False).drop_duplicates("iid").sort_values("row")
print("保留个体", len(d), "人群", d["label"].nunique())
samples = sys.argv[1:] or sorted(os.path.basename(p)[:-len(".codes.npy")] for p in glob.glob(f"{O}/*.codes.npy"))
rows = [np.asarray(t.mm[r]) for r in d["row"]] + [pack(np.load(f"{O}/{s}.codes.npy")) for s in samples]
ind_lines = [(g, s_, l) for g, s_, l in zip(d["gid"], d["sex"], d["label"])] + \
            [(s, "M" if wgs.sex_of(s) == "male" else "F", s) for s in samples]
write_tgeno(f"{O}/sub", rows, ind_lines, A + ".patch.PUB.snp")
d.drop(columns=["row"]).to_csv(f"{O}/sub.meta.tsv", sep="\t", index=False)
print("写出", f"{O}/sub.geno", len(rows), "个个体（含样本", samples, "）")
