#!/usr/bin/env python3
"""小变异解读：从 VEP 注释后的 VCF（work/annot/cohort.*）里挑出要报告 / 要人工复核的变异。

输出 work/findings/：
  <s>.acmg_sf.tsv      A. ACMG SF v3.3 次要发现（按基因规则）
  <s>.carrier.tsv      B. 隐性 / X 连锁携带（P/LP 或高置信功能缺失）
  <s>.biallelic.tsv    同一隐性基因两个致病等位（纯合 / 可能复合杂合）——需要关注
  <s>.nmd.tsv          D. 神经肌肉病基因罕见变异（含 VUS，给专业判断）
  <s>.risk_alleles.tsv 风险等位 / 低外显（ClinVar risk factor、APOE 等）
  pair.shared_carrier.tsv    （只在工作目录里恰好有两个人、比如一对伴侣时）两人在同一隐性基因都携带
  all_candidates.tsv   以上所有候选的完整证据（复核用）

分级（简化的 ACMG 思路，不替代临床判读）：
  P/LP(ClinVar)      ClinVar 致病 / 可能致病，≥1 星且无冲突
  LP(预测功能缺失)    罕见 + 无义 / 移码 / 经典剪接 / 起始密码子丢失，且不在逃逸 NMD 的末端（NMD 插件），
                     基因的致病机制是功能缺失（隐性病基因，或 LOEUF<0.6 的显性病基因）
  VUS-倾向致病        罕见错义 REVEL≥0.773 或 AlphaMissense likely_pathogenic，或 SpliceAI≥0.5
  VUS                 其余罕见编码 / 剪接区变异（只在神经肌肉面板里列出）
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import collections
import csv
import json
import os
import re
import sys

from cyvcf2 import VCF

VCF_PATH = sys.argv[1] if len(sys.argv) > 1 else wgs.annot_vcf()
OUT = os.environ.get("FINDINGS_OUT", "work/findings")
os.makedirs(OUT, exist_ok=True)

# ---------- 基因表 ----------
GT = {}
for r in csv.DictReader(open("ref/genelists/gene_table.tsv"), delimiter="\t"):
    GT[r["gene"]] = r

def moi(g):
    """遗传方式：只认 GenCC Definitive / Strong（Moderate 单独一条不足以据此报告携带 / 致病）。"""
    return set((GT.get(g, {}).get("moi_strong") or "").split(",")) - {""}

def loeuf(g):
    try:
        return float(GT.get(g, {}).get("loeuf") or "nan")
    except ValueError:
        return float("nan")

SEX = {}  # 读 VCF 后按样本填（见下）

# ClinGen 剂量敏感性：HI=3 表示“单倍剂量不足致病”证据充分 → 显性病基因的功能缺失才算致病机制
HI = {}
for line in open("ref/genelists/clingen_dosage_GRCh38.tsv"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    HI[f[0]] = f[4]

PAR = {"chrX": [(10_001, 2_781_479), (155_701_383, 156_030_895)], "chrY": [(10_001, 2_781_479), (56_887_903, 57_217_415)]}

def in_par(chrom, pos):
    return any(a <= pos <= b for a, b in PAR.get(chrom, []))

def canonical_splice(hgvsc):
    """VEP 的 HGVS 默认 3' 平移；剪接受体 / 供体后果里，只有内含子偏移 ≤2 的才算经典剪接位点被破坏。
    重复序列里的缺失左对齐时会压到剪接位点上，3' 平移后其实在内含子深处（如 PEX5 c.147+77_147+121del）。"""
    offs = [int(x) for x in re.findall(r"[+-](\d+)", hgvsc.split(":")[-1])]
    return (not offs) or min(offs) <= 2

def qc_ok(chrom, pos, zyg, vaf, dp, gq, sex):
    """读段层面的可信度：杂合 VAF≥0.30 且 ≥4 条支持读段；纯合 VAF≥0.8；男性 X/Y 非 PAR 区必须像半合子。"""
    alt = vaf * dp
    if gq != "" and gq is not None and int(gq) < 20:
        return False, "GQ<20"
    if sex == "male" and chrom in ("chrX", "chrY") and not in_par(chrom, pos):
        return (vaf >= 0.8, "男性 X/Y 非 PAR 区非半合子" if vaf < 0.8 else "")
    if zyg == "het" and (vaf < 0.30 or alt < 4):
        return False, f"杂合但 VAF={vaf:.2f}、支持读段≈{alt:.0f}"
    if zyg == "hom" and vaf < 0.8:
        return False, f"纯合但 VAF={vaf:.2f}"
    return True, ""

LOF = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
       "start_lost", "transcript_ablation", "stop_lost"}
CODING = LOF | {"missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant",
                "splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant",
                "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant"}
STARS = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "criteria_provided,_single_submitter": 1, "criteria_provided,_conflicting_classifications": 1,
         "criteria_provided,_conflicting_interpretations": 1}

# ---------- VCF ----------
vcf = VCF(VCF_PATH)
samples = vcf.samples
for s in samples:
    try:
        SEX[s] = open(f"work/{s}/qc/{s}.sex.txt").read().split("\t")[0]
    except FileNotFoundError:
        SEX[s] = "unknown"
csq_fmt = re.search(r"Format: ([^\"]+)", vcf.get_header_type("CSQ")["Description"]).group(1).split("|")
IDX = {k: i for i, k in enumerate(csq_fmt)}

def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

def pick_per_gene(csqs):
    """每个基因选一个转录本：MANE Select > Ensembl canonical > 第一个；同时找对应的 RefSeq NM 记录拿 HGVS。"""
    by_gene = collections.defaultdict(list)
    for c in csqs:
        if c["SYMBOL"]:
            by_gene[c["SYMBOL"]].append(c)
    out = {}
    for g, cs in by_gene.items():
        ens = [c for c in cs if c.get("SOURCE") == "Ensembl"] or cs
        mane = [c for c in ens if c.get("MANE_SELECT")]
        canon = [c for c in ens if c.get("CANONICAL") == "YES"]
        best = (mane or canon or ens)[0]
        nm = best.get("MANE_SELECT", "")  # Ensembl 记录的 MANE_SELECT 字段 = 对应的 RefSeq NM 号
        ref = [c for c in cs if nm and c["Feature"] == nm]
        best = dict(best)
        if ref:
            best["HGVSc_NM"] = ref[0].get("HGVSc", "")
            best["HGVSp_NM"] = ref[0].get("HGVSp", "")
        out[g] = (best, cs)
    return out

def clinvar(c):
    sig = c.get("ClinVar_CLNSIG", "") or ""
    rev = c.get("ClinVar_CLNREVSTAT", "") or ""
    stars = STARS.get(rev.split("&")[0], 0) if rev else 0
    s = sig.lower()
    if not sig:
        cls = None
    elif "conflicting" in s:
        cls = "conflicting"
        conf = c.get("ClinVar_CLNSIGCONF", "") or ""
        if "pathogenic" in conf.lower():
            cls = "conflicting(P/LP 在内)"
    elif "pathogenic" in s and "benign" not in s:
        cls = "P" if s.startswith("pathogenic") and "likely" not in s.split("/")[0] else "LP"
        if "low_penetrance" in s:
            cls += "(低外显)"
    elif "uncertain" in s:
        cls = "VUS"
    elif "benign" in s:
        cls = "B/LB"
    elif "risk" in s or "drug_response" in s or "association" in s or "protective" in s:
        cls = sig
    else:
        cls = sig
    return cls, stars, sig, rev, c.get("ClinVar_CLNDN", "") or ""

def spliceai_max(info):
    v = info.get("SpliceAI")
    if not v:
        return None, ""
    best, txt = 0.0, ""
    for e in str(v).split(","):
        p = e.split("|")
        if len(p) >= 6:
            ds = [fnum(x) or 0 for x in p[2:6]]
            m = max(ds)
            if m > best:
                best, txt = m, e
    return best, txt

rows = []
n = 0
for v in vcf:
    n += 1
    if v.FILTER not in (None, "PASS"):
        continue
    gts = v.genotypes
    carriers = [samples[i] for i, g in enumerate(gts) if 1 in g[:2]]
    if not carriers:
        continue
    raw = v.INFO.get("CSQ")
    if not raw:
        continue
    csqs = [dict(zip(csq_fmt, e.split("|"))) for e in raw.split(",")]
    sai, sai_txt = spliceai_max(v.INFO)
    for g, (c, allc) in pick_per_gene(csqs).items():
        cons = set(c["Consequence"].split("&"))
        cls, stars, sig, rev, dn = clinvar(c)
        max_af = fnum(c.get("MAX_AF")) or 0.0
        eas = max(fnum(c.get("gnomADe_EAS_AF")) or 0, fnum(c.get("gnomADg_EAS_AF")) or 0)
        gnomad = max(fnum(c.get("gnomADe_AF")) or 0, fnum(c.get("gnomADg_AF")) or 0)
        revel = fnum(c.get("REVEL"))
        am = fnum(c.get("am_pathogenicity"))
        amc = c.get("am_class", "")
        nmd_flag = c.get("NMD", "")
        in_scope = g in GT
        clinvar_plp = cls in ("P", "LP", "P(低外显)", "LP(低外显)") and stars >= 1
        lof = bool(cons & LOF)
        hgvsc_main = c.get("HGVSc_NM") or c.get("HGVSc") or ""
        # 功能缺失可信度：不在 NMD 逃逸区（VEP NMD 插件），剪接类必须真的落在经典 ±1/2 位
        splice_only = lof and not (cons & (LOF - {"splice_acceptor_variant", "splice_donor_variant"}))
        lof_ok = lof and "NMD_escaping" not in nmd_flag and (not splice_only or canonical_splice(hgvsc_main))
        # 致病机制是功能缺失：隐性 / X 连锁基因（携带意义），或 ClinGen 单倍剂量不足评分 3 的显性基因
        mechanism_lof = ("AR" in moi(g) or "XL" in moi(g)) or HI.get(g) == "3"
        coding = bool(cons & CODING) or (sai is not None and sai >= 0.2)
        if not (clinvar_plp or (cls or "").startswith("conflicting") or ((coding or lof) and max_af < 0.02)
                or (cls and cls not in ("B/LB", "VUS") and "pathogenic" not in (cls or "").lower())):
            continue
        if lof_ok and max_af < 0.005 and mechanism_lof and not clinvar_plp and cls != "B/LB":
            tier = "LP(预测功能缺失)"
        elif clinvar_plp:
            tier = "P/LP(ClinVar)"
        elif ((revel is not None and revel >= 0.773) or amc == "likely_pathogenic" or (sai or 0) >= 0.5) and max_af < 0.01 and cls != "B/LB":
            tier = "VUS-倾向致病"
        elif cls == "B/LB":
            tier = "B/LB(ClinVar)"
        elif coding and max_af < 0.01:
            tier = "VUS"
        else:
            tier = "其他"
        for s in carriers:
            i = samples.index(s)
            gt = gts[i]
            zyg = "hom" if gt[0] == gt[1] == 1 else "het"
            if v.CHROM in ("chrX", "chrY") and SEX.get(s) == "male" and zyg == "hom":
                zyg = "hemi"
            ad = v.format("AD")[i] if v.format("AD") is not None else [0, 0]
            dp = int(sum(x for x in ad if x >= 0)) if ad is not None else 0
            vaf = (ad[1] / dp) if dp else 0
            gq = int(v.format("GQ")[i][0]) if v.format("GQ") is not None else ""
            ok, why = qc_ok(v.CHROM, v.POS, zyg, vaf, dp, gq, SEX.get(s))
            s_tier = tier if ok or tier in ("其他", "B/LB(ClinVar)") else "低质量(疑似假阳性)"
            rows.append({
                "sample": s, "gene": g, "chrom": v.CHROM, "pos": v.POS, "ref": v.REF, "alt": v.ALT[0],
                "zygosity": zyg, "DP": dp, "VAF": round(vaf, 2), "GQ": gq, "qc_note": why,
                "consequence": c["Consequence"], "impact": c["IMPACT"],
                "transcript": c.get("MANE_SELECT") or c["Feature"], "HGVSc": c.get("HGVSc_NM") or c["HGVSc"],
                "HGVSp": c.get("HGVSp_NM") or c["HGVSp"], "exon": c.get("EXON", ""), "rsid": ",".join(x for x in c.get("Existing_variation", "").split("&") if x.startswith("rs")),
                "gnomAD_AF": f"{gnomad:.2e}" if gnomad else "0", "gnomAD_EAS_AF": f"{eas:.2e}" if eas else "0", "max_AF": f"{max_af:.2e}" if max_af else "0",
                "ClinVar": cls or "", "ClinVar_stars": stars if sig else "", "ClinVar_disease": dn.replace("_", " ")[:200],
                "REVEL": revel if revel is not None else "", "AlphaMissense": f"{am:.3f}({amc})" if am is not None else "",
                "SpliceAI": round(sai, 2) if sai is not None else "", "NMD": nmd_flag, "LoF": "Y" if lof else "",
                "tier": s_tier, "tier_before_qc": tier, "clingen_hi": HI.get(g, ""), "gencc_moi": GT.get(g, {}).get("gencc_moi", ""), "gencc_diseases": GT.get(g, {}).get("gencc_diseases", "")[:200],
                "acmg_sf": GT.get(g, {}).get("acmg_sf", ""), "nmd_myopathy": GT.get(g, {}).get("nmd_myopathy", ""),
                "nmd_superpanel": GT.get(g, {}).get("nmd_superpanel", ""), "loeuf": GT.get(g, {}).get("loeuf", ""),
            })
print(f"scanned {n} records; {len(rows)} sample-variant-gene candidate rows", file=sys.stderr)

FIELDS = list(rows[0].keys()) if rows else []
def write(path, rs):
    with open(path, "w") as f:
        w = csv.DictWriter(f, FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rs:
            w.writerow(r)

write(f"{OUT}/all_candidates.tsv", rows)
PATHO = ("P/LP(ClinVar)", "LP(预测功能缺失)")
summary = {}
for s in samples:
    R = [r for r in rows if r["sample"] == s]
    # A. ACMG SF
    sf = []
    by_gene = collections.defaultdict(list)
    for r in R:
        if r["acmg_sf"] and r["tier"] in PATHO:
            by_gene[r["gene"]].append(r)
    for g, rs in by_gene.items():
        rule = GT[g]["acmg_sf"]
        n_alleles = sum(2 if r["zygosity"] == "hom" else 1 for r in rs)
        if rule == "all":
            sf += rs
        elif rule == "biallelic" and n_alleles >= 2:
            sf += rs
        elif rule == "HFE" and any(r["zygosity"] == "hom" and "Cys282Tyr" in r["HGVSp"] for r in rs):
            sf += rs
        elif rule == "TTN" and any(r["LoF"] for r in rs):
            sf += [r for r in rs if r["LoF"]]
    write(f"{OUT}/{s}.acmg_sf.tsv", sf)
    # B. 携带者：隐性 / X 连锁基因上的致病等位（女性 X 连锁杂合 = 携带；男性半合子单列）
    car = [r for r in R if r["tier"] in PATHO and (("AR" in moi(r["gene"])) or ("XL" in moi(r["gene"])))]
    write(f"{OUT}/{s}.carrier.tsv", car)
    # 同一隐性基因 ≥2 个致病等位
    cnt = collections.Counter()
    for r in car:
        if "AR" in moi(r["gene"]):
            cnt[r["gene"]] += 2 if r["zygosity"] == "hom" else 1
    bi = [r for r in car if cnt.get(r["gene"], 0) >= 2 or (r["zygosity"] == "hemi")]
    write(f"{OUT}/{s}.biallelic.tsv", bi)
    # D. 神经肌肉病基因
    nmd = [r for r in R if (r["nmd_myopathy"] in ("2", "3") or r["nmd_superpanel"] == "3")
           and r["tier"] in PATHO + ("VUS-倾向致病", "VUS") and r["tier"] != "其他"]
    write(f"{OUT}/{s}.nmd.tsv", nmd)
    # 风险等位 / 低外显
    # 低外显致病（如 CHEK2、APC I1307K 一类）或 ≥2 星的风险因子；药物反应交给 PharmCAT
    risk = [r for r in R if r["qc_note"] == "" and r["ClinVar"] and ("低外显" in r["ClinVar"] or
            ("risk" in r["ClinVar"].lower() and str(r["ClinVar_stars"]) not in ("", "0", "1")))]
    write(f"{OUT}/{s}.risk_alleles.tsv", risk)
    summary[s] = {"sex": SEX.get(s), "acmg_sf": len(sf), "carrier": len(car), "carrier_genes": sorted({r["gene"] for r in car}),
                  "biallelic_genes": sorted({r["gene"] for r in bi}), "nmd_rows": len(nmd), "risk_rows": len(risk)}

# 两人（如一对伴侣）：同一隐性基因两人都携带致病等位 —— 只有工作目录里恰好两个样本时才做
car = {s: {r["gene"] for r in rows if r["sample"] == s and r["tier"] in PATHO and "AR" in moi(r["gene"])} for s in samples}
shared = sorted(car[samples[0]] & car[samples[1]]) if len(samples) == 2 else []
if len(samples) == 2:
    write(f"{OUT}/pair.shared_carrier.tsv", [r for r in rows if r["gene"] in shared and r["tier"] in PATHO])
# 一方致病 / 可能致病、另一方在同一隐性基因上有罕见意义未明变异（错义 / 框内 / 剪接区）→ 单列给遗传咨询
vus_tiers = ("VUS-倾向致病", "VUS")
partial = []
if len(samples) == 2:
    a, b = samples
    for x, y in ((a, b), (b, a)):
        px = {r["gene"] for r in rows if r["sample"] == x and r["tier"] in PATHO and "AR" in moi(r["gene"])}
        for r in rows:
            if r["sample"] == y and r["gene"] in px and r["gene"] not in shared and r["tier"] in vus_tiers:
                partial.append(r)
        partial += [r for r in rows if r["sample"] == x and r["gene"] in {q["gene"] for q in partial if q["sample"] == y} and r["tier"] in PATHO]
seen_p = set(); partial_u = []
for r in partial:
    k = (r["sample"], r["chrom"], r["pos"], r["alt"], r["gene"])
    if k not in seen_p:
        seen_p.add(k); partial_u.append(r)
if len(samples) == 2:
    write(f"{OUT}/pair.partial_overlap.tsv", partial_u)
# 已共同携带的基因里，另一方额外的 VUS（可能影响本人或孩子）
extra = [r for r in rows if r["gene"] in shared and r["tier"] in vus_tiers]
if len(samples) == 2:
    write(f"{OUT}/pair.shared_gene_extra_vus.tsv", extra)
# X 连锁：女方携带 → 儿子 50%
xl = {s: sorted({r["gene"] for r in rows if r["sample"] == s and r["tier"] in PATHO and "XL" in moi(r["gene"])}) for s in samples}
if len(samples) == 2:
    summary["pair"] = {"shared_AR_carrier_genes": shared, "XL_by_sample": xl,
                       "partial_overlap_genes": sorted({r["gene"] for r in partial_u}),
                       "shared_gene_extra_vus": sorted({(r["sample"], r["gene"], r["HGVSp"].split(":")[-1] or r["HGVSc"].split(":")[-1]) for r in extra})}
json.dump(summary, open(f"{OUT}/summary.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(summary, ensure_ascii=False, indent=1))
