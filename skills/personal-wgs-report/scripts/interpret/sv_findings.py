#!/usr/bin/env python3
"""SV / CNV 解读：合并 Manta、Delly、CNVpytor 的 AnnotSV 结果 → work/findings/<s>.sv.tsv

健康成年人的大片段致病 CNV 很少，SV 软件的假阳性却很多，所以规则偏严：
  基本过滤：DEL / DUP；≤5 Mb；人群频率（gnomAD-SV / DGV / 1000G 等，AnnotSV B_*_AFmax）<1%；
            断点不在缺口 / ENCODE 黑名单；两端不同时落在片段重复区；不在 PAR；
            Manta / Delly 要 PASS 且样本基因型非参考
  证据：≥2 种方法支持（重叠 ≥50%）；≥5 kb 的事件必须有 CNVpytor 读深支持
  临床相关：AnnotSV ACMG 4/5；或缺失覆盖隐性 / X 连锁基因编码区（→ 携带）；或缺失覆盖 ClinGen HI=3 基因编码区；
            或重复完整覆盖 ClinGen TS=3 基因
report=yes：全部满足；report=review：罕见 + 覆盖编码区，但证据或相关性不足（附在技术附录）
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import collections
import csv
import os

GT = {r["gene"]: r for r in csv.DictReader(open("ref/genelists/gene_table.tsv"), delimiter="\t")}
SAMPLES = [wgs.SAMPLE]
O = f"work/{wgs.SAMPLE}/sv/annot"
PAR = {"X": [(10_001, 2_781_479), (155_701_383, 156_030_895)], "Y": [(10_001, 2_781_479), (56_887_903, 57_217_415)]}

def rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

def gt_of(r, s):
    v = r.get(s, "")
    return v.split(":")[0] if v else ""

def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

HI, TS = {}, {}
for line in open("ref/genelists/clingen_dosage_GRCh38.tsv"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    HI[f[0]] = f[4]
    TS[f[0]] = f[12] if len(f) > 12 else ""

def moi(g):
    return set((GT.get(g, {}).get("moi_strong") or "").split(",")) - {""}

calls = collections.defaultdict(list)
for caller, path in (("Manta", f"{O}/manta.annotsv.tsv"), ("Delly", f"{O}/delly.annotsv.tsv"), ("CNVpytor", f"{O}/cnvpytor.annotsv.tsv")):
    allrows = rows(path)
    # split 行按基因给出编码区重叠；汇总到每个 SV（AnnotSV_ID）
    cds_genes = collections.defaultdict(set)
    for r in allrows:
        if r.get("Annotation_mode") == "split" and (fnum(r.get("Overlapped_CDS_length")) or 0) > 0 and r.get("Gene_name"):
            cds_genes[r["AnnotSV_ID"]].add(r["Gene_name"])
    for r in allrows:
        if r.get("Annotation_mode") != "full":
            continue
        r["_cds_genes"] = sorted(cds_genes.get(r["AnnotSV_ID"], set()))
        svt = r["SV_type"]
        if svt not in ("DEL", "DUP"):
            continue
        chrom, st, en = r["SV_chrom"], int(r["SV_start"]), int(r["SV_end"])
        if caller == "CNVpytor":
            who = [r.get("Samples_ID", "")]
        else:
            if r.get("FILTER") not in ("PASS", "."):
                continue
            who = [s for s in SAMPLES if gt_of(r, s) not in ("", "0/0", "./.", ".", "0|0")]
        for s in who:
            if s in SAMPLES:
                calls[s].append((caller, chrom, st, en, svt, r))

def overlap(a, b):
    if a[1] != b[1] or a[4] != b[4]:
        return False
    lo, hi = max(a[2], b[2]), min(a[3], b[3])
    if hi <= lo:
        return False
    return (hi - lo) >= 0.5 * max(1, a[3] - a[2]) and (hi - lo) >= 0.5 * max(1, b[3] - b[2])

def in_par(chrom, st, en):
    return any(not (en < a or st > b) for a, b in PAR.get(chrom, []))

os.makedirs("work/findings", exist_ok=True)
for s in SAMPLES:
    out = []
    for c in calls[s]:
        caller, chrom, st, en, svt, r = c
        size = en - st
        support = {caller} | {d[0] for d in calls[s] if d is not c and overlap(c, d)}
        af = max([x for x in (fnum(r.get("B_loss_AFmax")), fnum(r.get("B_gain_AFmax"))) if x is not None] or [0])
        genes = [g for g in (r.get("Gene_name") or "").split(";") if g]
        cdsg = r.get("_cds_genes", [])
        cds = bool(cdsg)
        acmg = r.get("ACMG_class", "")
        bad_bp = any(r.get(k) for k in ("Gap_left", "Gap_right", "ENCODE_blacklist_left", "ENCODE_blacklist_right")) \
            or (r.get("SegDup_left") and r.get("SegDup_right"))
        reasons = []
        if size > 5_000_000:
            reasons.append(">5Mb")
        if af >= 0.01:
            reasons.append(f"人群频率 {af:.2g}")
        if bad_bp:
            reasons.append("断点在缺口/黑名单/片段重复")
        if in_par(chrom, st, en):
            reasons.append("PAR")
        if len(support) < 2:
            reasons.append("单一方法")
        if size >= 5000 and "CNVpytor" not in support:
            reasons.append("≥5kb 无读深支持")
        ar = [g for g in cdsg if moi(g) & {"AR", "XL"}]
        hi3 = [g for g in cdsg if GT.get(g, {}).get("acmg_sf")]
        hi_genes = [g for g in cdsg if HI.get(g) == "3"] if svt == "DEL" else []
        ts_genes = [g for g in genes if TS.get(g) == "3"] if svt == "DUP" else []
        relevant = acmg in ("4", "5") or (svt == "DEL" and cds and (ar or hi_genes or hi3)) or (svt == "DUP" and ts_genes)
        if not relevant and not (cds and af < 0.01):
            continue
        if not reasons and relevant:
            report = "yes"
        elif af < 0.01 and cds and size <= 5_000_000 and not in_par(chrom, st, en):
            report = "review"
        else:
            continue
        note = []
        if ar and svt == "DEL":
            note.append("覆盖隐性/X连锁基因编码区：" + ",".join(ar[:5]))
        if hi_genes:
            note.append("ClinGen HI=3：" + ",".join(hi_genes[:5]))
        if acmg in ("4", "5"):
            note.append(f"AnnotSV ACMG {acmg}")
        if reasons:
            note.append("证据不足：" + "、".join(reasons))
        out.append({"sample": s, "type": svt, "region": f"chr{chrom}:{st:,}-{en:,}", "size": f"{size:,}",
                    "genes": ",".join(cdsg[:20]) or ",".join(genes[:8]), "acmg_class": acmg, "pop_af": f"{af:.3g}", "callers": "+".join(sorted(support)),
                    "gt": gt_of(r, s) if caller != "CNVpytor" else "", "location": r.get("Location", ""), "report": report,
                    "note": "；".join(note), "caller": caller})
    seen, uniq = set(), []
    for o in sorted(out, key=lambda x: (x["report"] != "yes", -len(x["callers"].split("+")), x["region"])):
        k = (o["type"], o["region"].split(":")[0], o["genes"][:60])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(o)
    keys = ["sample", "type", "region", "size", "genes", "acmg_class", "pop_af", "callers", "gt", "location", "report", "note", "caller"]
    with open(f"work/findings/{s}.sv.tsv", "w") as f:
        w = csv.DictWriter(f, keys, delimiter="\t")
        w.writeheader()
        for o in uniq:
            w.writerow(o)
    print(s, "SV calls", len(calls[s]), "kept", len(uniq), "yes", sum(o["report"] == "yes" for o in uniq),
          "review", sum(o["report"] == "review" for o in uniq))
    for o in uniq:
        if o["report"] == "yes":
            print("  YES", o["type"], o["region"], o["size"], o["genes"][:60], o["callers"], o["note"])
