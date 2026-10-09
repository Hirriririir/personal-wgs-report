#!/usr/bin/env python3
"""几个“拷贝数会变”的基因（读段深度法）→ work/ext/cnv_genes.<s>.json：拷贝数 = 2 × 区域平均深度 ÷ 常染色体平均深度（同一过滤条件）。
多拷贝/高同源基因（AMY1、C4、LPA KIV-2）用 MAPQ≥0 并把参考里的各拷贝加总：总拷贝数 = 2 × Σ各参考拷贝深度 ÷ D；
单拷贝、可能整段缺失的基因（UGT2B17、GSTM1、CYP2A6、RHD）用 MAPQ≥20 的唯一比对读段。
深度法看不出拷贝在哪条染色体上、也分不清基因转换；临床意义的判断要用 MLPA / 长读长等方法确认。
用法：cnv_genes.py <sample>"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import wgs  # noqa: E402
import depth_lib as dl  # noqa: E402
import numpy as np  # noqa: E402
import pysam  # noqa: E402

s = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
bam = dl.open_bam(s)
fa = pysam.FastaFile(wgs.REF)


def depth(chrom, start, end, minq):
    return dl.depth(bam, chrom, start, end, minq, 0)


D = {q: dl.autosomal_depth(bam, fa, minq=q, baseq=0) for q in (0, 20)}
genes = {}
for line in open("ref/annot/gencode/pc_genes.bed"):
    c, a, b, g = line.split()
    genes.setdefault(g, (c, int(a), int(b)))
res = {"sample": s, "auto_depth_q0": round(D[0], 2), "auto_depth_q20": round(D[20], 2)}


def cn_single(g, q=20):
    c, a, b = genes[g]
    return round(2 * float(np.mean(depth(c, a, b, q))) / D[q], 2)


def cn_multi(gs, q=0):
    tot = sum(float(np.mean(depth(*genes[g], q))) for g in gs)
    return round(2 * tot / D[q], 2)


res["AMY1_total"] = cn_multi(["AMY1A", "AMY1B", "AMY1C"])
res["AMY2A"] = cn_single("AMY2A", 0)
res["AMY2B"] = cn_single("AMY2B", 0)
res["C4_total"] = cn_multi(["C4A", "C4B"])
for g in ("UGT2B17", "UGT2B15", "GSTM1", "CYP2A6", "CYP2A7", "RHD", "RHCE", "FCGR3A", "FCGR3B"):
    res[g] = cn_single(g, 20)
# LPA KIV-2：参考里 6 个单元（chr6:160,612,481–160,647,481）；对照用 LPA 其余部分
kiv = depth("chr6", 160612481, 160647481, 0)
lpa_rest = np.concatenate([depth("chr6", 160531481, 160605000, 0), depth("chr6", 160650000, 160664276, 0)])
res["LPA_KIV2_total_copies"] = round(2 * 6 * float(np.mean(kiv)) / D[0], 1)
res["LPA_nonrepeat_CN"] = round(2 * float(np.mean(lpa_rest)) / D[0], 2)
# 深度剖面（1 kb 窗口）供画图
prof = {}
for name, (c, a, b) in {"LPA": ("chr6", 160520000, 160675000), "AMY": ("chr1", 103540000, 103770000),
                        "UGT2B17": ("chr4", 68480000, 68680000), "GSTM1": ("chr1", 109660000, 109720000),
                        "CYP2A6": ("chr19", 40830000, 40890000)}.items():
    q = 0 if name in ("LPA", "AMY") else 20
    d = depth(c, a, b, q)
    prof[name] = dict(chrom=c, start=a, q=q, win=1000,
                      cn=[round(2 * float(np.mean(d[i:i + 1000])) / D[q], 2) for i in range(0, len(d), 1000)])
os.makedirs("work/ext", exist_ok=True)
json.dump(dict(res, profiles=prof), open(f"work/ext/cnv_genes.{s}.json", "w"), indent=1)
print(json.dumps(res, ensure_ascii=False, indent=1))
