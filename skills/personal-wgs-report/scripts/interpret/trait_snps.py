#!/usr/bin/env python3
"""性状 / 单位点风险等位：从 DeepVariant gVCF 取基因型（参考区块 GQ≥20 记纯合参考）。

位点坐标首次运行时用 Ensembl REST 查 GRCh38 位置并缓存到 ref/genelists/trait_snps.tsv。
输出 work/findings/traits.tsv（每人每位点一行）和 traits.json。
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import csv
import json
import os
import sys
import time
import urllib.request

import pysam

SNPS = [
    # rsid, 基因, 主题, 效应等位基因, 解释（效应等位基因的作用）
    ("rs671", "ALDH2", "饮酒脸红 / 乙醛代谢", "A", "ALDH2*2（Glu504Lys），携带者乙醛清除慢：喝酒脸红、心悸；常饮酒者食管鳞癌风险明显升高"),
    ("rs1229984", "ADH1B", "乙醇代谢速度", "T", "ADH1B*2（His48Arg），乙醇→乙醛快；东亚常见，与 ALDH2*2 叠加时脸红更明显"),
    ("rs1815739", "ACTN3", "快肌纤维 α-辅肌动蛋白-3", "T", "R577X；XX 纯合者快肌纤维缺 α-actinin-3，偏耐力型、爆发力稍弱（普通人影响很小）"),
    ("rs4988235", "LCT/MCM6", "成年乳糖耐受", "A", "-13910*T（正链 A）；带 A 者成年后多能消化乳糖。东亚人几乎都是 G/G（乳糖不耐受很常见）"),
    ("rs17822931", "ABCC11", "耳垢类型 / 体味", "T", "T/T 为干耳垢、腋臭轻（东亚多数）；C 等位基因为湿耳垢"),
    ("rs3827760", "EDAR", "头发粗细 / 铲形门齿", "G", "V370A（G）；东亚高频，头发更粗、铲形门齿、汗腺更多"),
    ("rs12913832", "HERC2/OCA2", "眼睛颜色", "G", "G/G 多为浅色眼睛（东亚几乎都是 A/A 深色）"),
    ("rs1426654", "SLC24A5", "肤色", "A", "A 等位基因与较浅肤色有关（欧洲高频，东亚低频）"),
    ("rs713598", "TAS2R38", "苦味敏感（PTC）", "G", "G（Ala49Pro 中的 Pro）与尝得出 PTC 苦味有关；完整判断需三个位点单倍型"),
    ("rs762551", "CYP1A2", "咖啡因代谢", "A", "A/A 为 CYP1A2*1F 快代谢型（证据中等，仅供参考）"),
    ("rs429358", "APOE", "APOE ε4 判定位点 1", "C", "与 rs7412 一起决定 APOE ε2/ε3/ε4"),
    ("rs7412", "APOE", "APOE ε2 判定位点 2", "T", "与 rs429358 一起决定 APOE ε2/ε3/ε4"),
    ("rs8176719", "ABO", "ABO 血型（O 等位基因，c.261delG）", "-", "“-”= c.261delG 移码 → O 等位（GRCh38 参考本身就是 O 型序列，检出插入表示一个非 O 等位）"),
    ("rs8176746", "ABO", "ABO 血型（B 等位基因标记 c.796C>A）", "T", "B 等位基因的特征位点之一（正链 T）"),
    ("rs1801133", "MTHFR", "MTHFR C677T", "A", "常被商业检测夸大；单独意义很小，不需要据此补充叶酸剂量以外的干预"),
    ("rs6025", "F5", "凝血因子 V Leiden", "T", "东亚极罕见；携带者静脉血栓风险升高"),
    ("rs1799963", "F2", "凝血酶原 G20210A", "A", "东亚极罕见；携带者静脉血栓风险升高"),
    ("rs10455872", "LPA", "脂蛋白(a) 升高", "G", "欧洲常见的 Lp(a) 升高位点；东亚少见"),
    ("rs2814778", "ACKR1", "Duffy 血型阴性", "C", "非洲高频；东亚极少"),
]

CACHE = "ref/genelists/trait_snps.tsv"

def ensembl_pos(rsid):
    url = f"https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                d = json.load(r)
            for m in d.get("mappings", []):
                if m.get("assembly_name") == "GRCh38" and m.get("seq_region_name") in [str(i) for i in range(1, 23)] + ["X", "Y"]:
                    return "chr" + m["seq_region_name"], int(m["start"]), int(m["end"]), m["allele_string"]
        except Exception:  # noqa: BLE001
            time.sleep(2)
    return None

if not os.path.exists(CACHE):
    with open(CACHE, "w") as o:
        o.write("rsid\tchrom\tstart\tend\tallele_string\n")
        for rs, *_ in SNPS:
            p = ensembl_pos(rs)
            print(rs, p, file=sys.stderr)
            if p:
                o.write(f"{rs}\t{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\n")
POS = {r["rsid"]: r for r in csv.DictReader(open(CACHE), delimiter="\t")}

FA = pysam.FastaFile("ref/GRCh38/GRCh38.fa")

def genotype(vf, chrom, start, end):
    """返回该位置的两个等位基因（正链），纯合参考区块给出参考碱基；插入 / 缺失型位点单独处理。"""
    ref_base = FA.fetch(chrom, start - 1, end).upper() if start <= end else "-"
    for rec in vf.fetch(chrom, start - 2, end):
        smp = rec.samples[0]
        alts = rec.alts or ()
        if alts in ((), ("<*>",)):
            if rec.start <= start - 1 < (rec.stop or rec.start + 1):
                return ((ref_base, ref_base), smp.get("GQ"), "ref-block")
            continue
        gt = smp.get("GT")
        if gt is None or None in gt:
            continue
        # SNV 正好在该位置
        if rec.pos == start and len(rec.ref) == 1 and start == end:
            al = [rec.alleles[i] for i in gt]
            return (tuple(al), smp.get("GQ"), "variant")
        # 缺失（VCF 左锚定在前一位）
        if rec.pos == start - 1 and len(rec.ref) > 1:
            al = ["-" if len(rec.alleles[i]) < len(rec.ref) else ref_base for i in gt]
            return (tuple(al), smp.get("GQ"), "variant")
        # 插入：Ensembl 用 start=end+1 表示两碱基之间（GRCh38 参考本身是缺失型，例如 ABO c.261 的 O 等位）
        if start > end and rec.pos == end and any(len(a) > len(rec.ref) for a in rec.alleles[1:]):
            al = [rec.alleles[i][len(rec.ref):] if len(rec.alleles[i]) > len(rec.ref) else "-" for i in gt]
            return (tuple(al), smp.get("GQ"), "variant")
    if start > end:
        # 插入位点没有变异记录：看是否在高质量参考区块里 → 两条都是参考（缺失型）
        for rec in vf.fetch(chrom, end - 1, end):
            if (rec.alts or ()) in ((), ("<*>",)) and (rec.samples[0].get("GQ") or 0) >= 20:
                return (("-", "-"), rec.samples[0].get("GQ"), "ref-block")
    return (None, None, "no-call")

out = []
SAMPLES = sys.argv[1:] or wgs.all_samples() or [wgs.SAMPLE]
for s in SAMPLES:
    vf = pysam.VariantFile(f"work/{s}/{s}.dv.g.vcf.gz")
    for rs, gene, topic, eff, note in SNPS:
        p = POS.get(rs)
        if not p:
            out.append({"sample": s, "rsid": rs, "gene": gene, "topic": topic, "genotype": "坐标未知"})
            continue
        al, gq, how = genotype(vf, p["chrom"], int(p["start"]), int(p["end"]))
        gt = "/".join(al) if al else "未检出"
        n_eff = sum(a == eff for a in al) if al else ""
        out.append({"sample": s, "rsid": rs, "gene": gene, "topic": topic, "pos": f'{p["chrom"]}:{p["start"]}',
                    "ensembl_alleles": p["allele_string"], "genotype": gt, "effect_allele": eff, "effect_allele_count": n_eff,
                    "GQ": gq, "source": how, "note": note})

# APOE：rs429358(C) + rs7412(C) → ε4；T/T → ε3 / ε2 判定
def apoe(rows, s):
    g = {r["rsid"]: r["genotype"] for r in rows if r["sample"] == s}
    a, b = g.get("rs429358", ""), g.get("rs7412", "")
    if "/" not in a or "/" not in b:
        return "无法判定"
    a1, a2 = a.split("/"); b1, b2 = b.split("/")
    # 未定相：按常见组合推断（ε1 极罕见忽略）
    hap = []
    c429 = [a1, a2].count("C"); t412 = [b1, b2].count("T")
    e4, e2 = c429, t412
    e3 = 2 - e4 - e2
    if e3 < 0:
        return f"罕见组合（rs429358 {a}, rs7412 {b}）"
    hap = ["ε2"] * e2 + ["ε3"] * e3 + ["ε4"] * e4
    return "/".join(sorted(hap))

def abo(rows, s):
    """粗略推断：O 等位 = c.261delG 个数，B 等位 = c.796A（rs8176746 正链 T）个数，其余为 A。
    未定相、未覆盖少见的 O2 / 亚型，只作参考，以血型检测为准。"""
    g = {r["rsid"]: r for r in rows if r["sample"] == s}
    o, b = g.get("rs8176719", {}), g.get("rs8176746", {})
    if o.get("effect_allele_count") in ("", None) or b.get("effect_allele_count") in ("", None):
        return "无法判定"
    n_o, n_b = int(o["effect_allele_count"]), int(b["effect_allele_count"])
    n_a = 2 - n_o - n_b
    if n_a < 0:
        return f"无法判定（O 等位 {n_o}、B 标记 {n_b}）"
    alle = ["A"] * n_a + ["B"] * n_b + ["O"] * n_o
    pheno = {"AA": "A", "AO": "A", "BB": "B", "BO": "B", "AB": "AB", "OO": "O"}["".join(sorted(alle))]
    return f"{pheno} 型（推测基因型 {'/'.join(sorted(alle))}）"

summary = {s: {"APOE": apoe(out, s), "ABO": abo(out, s)} for s in SAMPLES}
os.makedirs("work/findings", exist_ok=True)
with open("work/findings/traits.tsv", "w") as f:
    keys = ["sample", "rsid", "gene", "topic", "pos", "ensembl_alleles", "genotype", "effect_allele", "effect_allele_count", "GQ", "source", "note"]
    w = csv.DictWriter(f, keys, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in out:
        w.writerow(r)
json.dump(summary, open("work/findings/traits.json", "w"), ensure_ascii=False, indent=1)
for r in out:
    print(r["sample"], r["rsid"], r["gene"], r.get("genotype"), r.get("source"))
print(summary)
