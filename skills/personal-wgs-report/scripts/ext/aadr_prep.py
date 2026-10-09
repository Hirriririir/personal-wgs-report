#!/usr/bin/env python3
"""古 DNA 准备：AADR 1240K 位点（hg19）→ GRCh38 坐标，并把样本在这些位点上的基因型编码成 AADR 的 2-bit 格式。
输出 work/ext/aadr/：
  sites38.tsv              chr pos idx hg38REF other a1is（喂给 gvcf_genotype_stream_par.py）
  <s>.1240k.vcf.gz         样本在这些位点的基因型（x02_aadr.sh 生成）
  <s>.codes.npy            长度 = AADR 位点数；值 = AADR 第 5 列等位基因的拷贝数（0/1/2），3 = 缺失
用法：aadr_prep.py lift | encode <sample>
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401
import numpy as np  # noqa: E402
import pysam  # noqa: E402

A = "ref/aadr/v66.p1_1240K.aadr.patch.PUB"
O = "work/ext/aadr"
os.makedirs(O, exist_ok=True)
COMP = str.maketrans("ACGT", "TGCA")


def read_snp():
    snps = []
    with open(A + ".snp") as f:
        for line in f:
            sid, c, gp, pos, a1, a2 = line.split()
            snps.append((sid, c, int(pos), a1, a2))
    return snps


def lift():
    snps = read_snp()
    bed = f"{O}/aadr19.bed"
    with open(bed, "w") as o:
        for k, (sid, c, pos, a1, a2) in enumerate(snps):
            if c.isdigit() and 1 <= int(c) <= 22:
                o.write(f"chr{c}\t{pos - 1}\t{pos}\t{k}\n")
    crossmap = os.path.join(os.path.dirname(sys.executable), "CrossMap")
    subprocess.run([crossmap, "bed", "ref/liftover/hg19ToHg38.over.chain.gz", bed, f"{O}/aadr38.bed"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    fa = pysam.FastaFile(wgs.REF)
    out, seen, stats = [], set(), dict(ok=0, swap=0, flip=0, bad=0, chr=0, dup=0)
    for line in open(f"{O}/aadr38.bed"):
        c38, s, e, k = line.split()[:4]
        k = int(k)
        sid, c19, pos19, a1, a2 = snps[k]
        if c38 != "chr" + c19:
            stats["chr"] += 1
            continue
        pos = int(e)
        if (c38, pos) in seen:
            stats["dup"] += 1
            continue
        seen.add((c38, pos))
        b = fa.fetch(c38, pos - 1, pos).upper()
        pal = {a1, a2} in ({"A", "T"}, {"C", "G"})
        if b == a1:
            out.append((c38, pos, k, a1, a2, "R")); stats["ok"] += 1
        elif b == a2:
            out.append((c38, pos, k, a2, a1, "A")); stats["swap"] += 1
        elif b == a1.translate(COMP) and not pal:   # 链方向翻转：记成互补碱基
            out.append((c38, pos, k, a1.translate(COMP), a2.translate(COMP), "R")); stats["flip"] += 1
        elif b == a2.translate(COMP) and not pal:
            out.append((c38, pos, k, a2.translate(COMP), a1.translate(COMP), "A")); stats["flip"] += 1
        else:
            stats["bad"] += 1
    order = {f"chr{i}": i for i in range(1, 23)}
    out.sort(key=lambda r: (order[r[0]], r[1]))
    with open(f"{O}/sites38.tsv", "w") as o:
        for c, p, k, r, a, a1is in out:
            o.write(f"{c}\t{p}\t{k}\t{r}\t{a}\t{a1is}\n")
    print("AADR 位点", len(snps), "转到 GRCh38 的常染色体位点", len(out), stats)


def encode(s):
    snps = read_snp()
    a1is = {}
    for line in open(f"{O}/sites38.tsv"):
        f = line.split()
        a1is[int(f[2])] = f[5]
    codes = np.full(len(snps), 3, dtype=np.uint8)
    n = dict(called=0, miss=0)
    for rec in pysam.VariantFile(f"{O}/{s}.1240k.vcf.gz"):
        k = int(rec.id)
        gt = rec.samples[0]["GT"]
        if gt is None or None in gt:
            n["miss"] += 1
            continue
        want = 0 if a1is[k] == "R" else 1   # AADR 第 5 列等位基因对应 hg38 的 REF（R）还是 ALT（A）
        codes[k] = sum(g == want for g in gt)
        n["called"] += 1
    np.save(f"{O}/{s}.codes.npy", codes)
    print(s, n, "有基因型的位点", int((codes != 3).sum()), "/", len(codes))


if __name__ == "__main__":
    if sys.argv[1] == "lift":
        lift()
    elif sys.argv[1] == "encode":
        encode(sys.argv[2] if len(sys.argv) > 2 else wgs.SAMPLE)
