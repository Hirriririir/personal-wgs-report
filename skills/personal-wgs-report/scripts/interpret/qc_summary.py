#!/usr/bin/env python3
"""汇总每个样本的数据质量指标 → work/findings/qc.json（报告“数据质量”一节和技术附录用）。

来源（有哪个用哪个）：
  fastp JSON（原始读段数、碱基数、Q30、GC、fastp 估计的重复率）、读段名里的仪器 / 泳道
  samtools stats（比对率、错误率、插入片段）、Parabricks 的 Picard 重复指标或 samtools markdup 统计
  mosdepth（平均深度、≥N× 覆盖比例）、性别推断、VerifyBamID2 污染、联合分型的变异计数、somalier 亲缘
用法：qc_summary.py [样本 ...]（默认工作目录里所有样本）
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import csv
import glob
import gzip
import json
import os
import re
import sys

# 读段名开头的仪器编号 → 机型（Illumina 惯例；只作提示）
INSTRUMENT = (("LH", "NovaSeq X 系列"), ("A0", "NovaSeq 6000"), ("E0", "HiSeq X"), ("K00", "HiSeq 4000"),
              ("J00", "HiSeq 3000"), ("D00", "HiSeq 2500"), ("NB", "NextSeq 500/550"), ("VH", "NextSeq 2000"),
              ("M0", "MiSeq"), ("FS", "iSeq"))


def read_picard(path, key):
    """Picard 风格指标文件里第一张表（表头含 key 的那一行开始）。"""
    if not os.path.exists(path):
        return []
    lines = [x.rstrip("\n") for x in open(path) if x.strip() and not x.startswith("#")]
    for i, line in enumerate(lines):
        if key in line.split("\t"):
            hdr, rows = line.split("\t"), []
            for x in lines[i + 1:]:
                p = x.split("\t")
                if len(p) != len(hdr):
                    break
                rows.append(dict(zip(hdr, p)))
            return rows
    return []


def fastp(s):
    q, js = {}, sorted(glob.glob(f"qc/fastp/{s}.*.json"))
    if not js:
        return q
    reads = bases = q30 = gc = dup = 0.0
    lens = []
    for f in js:
        d = json.load(open(f))
        b = d["summary"]["before_filtering"]
        reads += b["total_reads"]; bases += b["total_bases"]
        q30 += b["q30_rate"] * b["total_bases"]; gc += b["gc_content"] * b["total_bases"]
        dup += d.get("duplication", {}).get("rate", 0) * b["total_reads"]
        lens.append(b.get("read1_mean_length"))
    q.update(raw_read_pairs=int(reads / 2), raw_bases_G=round(bases / 1e9, 1), Q30=round(q30 / bases * 100, 2),
             GC=round(gc / bases * 100, 1), fastp_dup_rate=round(dup / reads * 100, 1), read_length=max(x or 0 for x in lens),
             fastq_pairs=len(js))
    return q


def lanes(s):
    p = f"qc/{s}.lanes.txt"
    if not os.path.exists(p):
        return {}
    rows, mgi = [], False
    for line in open(p):
        n, key = line.split()
        f = key.lstrip("@").split(":")
        if len(f) >= 4:                       # Illumina：仪器:运行号:芯片:泳道
            rows.append((int(n), f[0], f"{f[2]}:{f[3]}"))
        elif re.match(r"[A-Z]*\d+L\d+$", f[0]):   # 华大：芯片号 + L 泳道
            mgi = True
            fc, ln = re.match(r"([A-Z]*\d+)L(\d+)$", f[0]).groups()
            rows.append((int(n), fc, f"{fc}:{ln}"))
    if not rows:
        return {}
    tot = sum(r[0] for r in rows)
    inst = sorted({r[1] for r in rows})
    model = sorted({m for i in inst for pre, m in INSTRUMENT if i.startswith(pre)}) + (["DNBSEQ（华大 / MGI）"] if mgi else [])
    return {"instruments": len(inst) if not mgi else None, "lanes": len(rows), "platform_guess": "、".join(model) or "未知",
            "top_lane_pct": round(max(r[0] for r in rows) / tot * 100, 1)}


def samtools_stats(s):
    q, p = {}, f"work/{s}/qc/{s}.samtools_stats.txt"
    if not os.path.exists(p):
        return q
    sn, ins = {}, []
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if f[0] == "SN":
            sn[f[1].rstrip(":")] = f[2]
        elif f[0] == "IS":
            ins.append((int(f[1]), int(f[2])))
    tot = float(sn.get("raw total sequences", 0) or 0)
    if tot:
        q["pct_aligned"] = round(float(sn["reads mapped"]) / tot * 100, 2)
        q["pct_unmapped"] = round(100 - q["pct_aligned"], 2)
        q["unmapped_reads"] = int(tot - float(sn["reads mapped"]))
        q["pct_properly_paired"] = round(float(sn.get("reads properly paired", 0)) / tot * 100, 2)
        q["pct_dup_samtools"] = round(float(sn.get("reads duplicated", 0)) / tot * 100, 1)
        q["pct_mapq0"] = round(float(sn.get("reads MQ0", 0)) / tot * 100, 2)
    if "error rate" in sn:
        q["error_rate"] = float(sn["error rate"])
    if "insert size average" in sn:
        q["insert_mean"] = round(float(sn["insert size average"]), 1)
        q["insert_sd"] = round(float(sn["insert size standard deviation"]), 1)
    n = sum(c for _, c in ins)
    if n:
        acc = 0
        for size, c in ins:
            acc += c
            if acc >= n / 2:
                q["insert_median"] = size
                break
    return q


def duplicates(s):
    q = {}
    dup = read_picard(f"work/{s}/{s}.dup_metrics.txt", "LIBRARY")   # GPU：Parabricks（Picard 格式）
    if dup:
        d = dup[0]
        pairs, dpairs, opt = int(d["READ_PAIRS_EXAMINED"]), int(d["READ_PAIR_DUPLICATES"]), int(d["READ_PAIR_OPTICAL_DUPLICATES"])
        q.update(dup_pct=round(float(d["PERCENT_DUPLICATION"]) * 100, 1), optical_dup_pct_of_pairs=round(opt / pairs * 100, 1),
                 pcr_dup_pct_of_pairs=round((dpairs - opt) / pairs * 100, 1), est_library_size=int(d["ESTIMATED_LIBRARY_SIZE"] or 0))
        return q
    p = f"work/{s}/{s}.markdup_stats.txt"                              # CPU：samtools markdup -f
    if os.path.exists(p):
        m = {}
        for line in open(p):
            if ":" in line:
                k, v = line.split(":", 1)
                m[k.strip()] = v.strip()
        num = lambda k: float(m.get(k) or 0)  # noqa: E731
        paired, single = num("PAIRED"), num("SINGLE")
        if paired:
            q.update(dup_pct=round((num("DUPLICATE PAIR") + num("DUPLICATE SINGLE")) / (paired + single) * 100, 1),
                     optical_dup_pct_of_pairs=round(num("DUPLICATE PAIR OPTICAL") / paired * 100, 1),
                     pcr_dup_pct_of_pairs=round((num("DUPLICATE PAIR") - num("DUPLICATE PAIR OPTICAL")) / paired * 100, 1),
                     est_library_size=int(num("ESTIMATED_LIBRARY_SIZE")))
    return q


def coverage(s):
    q, p = {}, f"work/{s}/qc/{s}.mosdepth.summary.txt"
    if os.path.exists(p):
        md = {r["chrom"]: float(r["mean"]) for r in csv.DictReader(open(p), delimiter="\t")}
        auto = [md.get(f"chr{i}", 0) for i in range(1, 23)]
        q.update(mean_cov_autosomal=round(sum(auto) / 22, 1), mean_cov_chrX=round(md.get("chrX", 0), 1),
                 mean_cov_chrY=round(md.get("chrY", 0), 1), mean_cov_chrM=round(md.get("chrM", 0)))
        if q["mean_cov_autosomal"]:
            # 线粒体拷贝数 ≈ 2 × chrM 深度 / 常染色体深度（每个细胞的大致线粒体 DNA 份数）
            q["mtDNA_copies_per_cell"] = round(2 * md.get("chrM", 0) / q["mean_cov_autosomal"])
    t = f"work/{s}/qc/{s}.thresholds.bed.gz"
    if os.path.exists(t):
        tot, c = 0, [0] * 5
        with gzip.open(t, "rt") as f:
            next(f)
            for line in f:
                p_ = line.split("\t")
                if not p_[0][3:].isdigit():
                    continue
                tot += int(p_[2]) - int(p_[1])
                for k in range(5):
                    c[k] += int(p_[4 + k])
        for k, th in enumerate((1, 10, 15, 20, 30)):
            q[f"pct_autosome_ge{th}x"] = round(c[k] / tot * 100, 1)
    return q


def misc(s):
    q = {}
    sx = f"work/{s}/qc/{s}.sex.txt"
    if os.path.exists(sx):
        f = open(sx).read().split("\t")
        q["sex_inferred"] = f[0]
        q["sex_detail"] = " ".join(f[1:]).strip()
    vb = f"work/{s}/qc/{s}.verifybamid2.selfSM"
    if os.path.exists(vb):
        r = list(csv.DictReader(open(vb), delimiter="\t"))[0]
        q["freemix_contamination_pct"] = round(float(r["FREEMIX"]) * 100, 2)
    return q


samples = sys.argv[1:] or wgs.all_samples() or [wgs.SAMPLE]
out = {}
for s in samples:
    q = {}
    for fn in (fastp, lanes, samtools_stats, duplicates, coverage, misc):
        q.update(fn(s))
    out[s] = q

# 变异计数（联合分型后、拆多等位）：PSC id sample nRefHom nNonRefHom nHets nTransitions nTransversions nIndels ...
st = "work/joint/cohort.dv.norm.stats.txt"
if os.path.exists(st):
    for line in open(st):
        if line.startswith("PSC\t"):
            p = line.rstrip().split("\t")
            if p[2] in out:
                out[p[2]].update({"n_hom_alt": int(p[4]), "n_het": int(p[5]), "n_snv": int(p[6]) + int(p[7]), "n_indel": int(p[8]),
                                  "het_hom_ratio": round(int(p[5]) / max(1, int(p[4])), 2), "ti_tv": round(int(p[6]) / max(1, int(p[7])), 3)})
rel = "work/ancestry/relate.pairs.tsv"
if os.path.exists(rel):
    keep = ("#sample_a", "sample_b", "relatedness", "ibs0", "ibs2", "hom_concordance", "n")
    out["_relatedness"] = [{k: r[k] for k in r if k in keep} for r in csv.DictReader(open(rel), delimiter="\t")]
os.makedirs("work/findings", exist_ok=True)
json.dump(out, open("work/findings/qc.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(out, ensure_ascii=False, indent=1))
