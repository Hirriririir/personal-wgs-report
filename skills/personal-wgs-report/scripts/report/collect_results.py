#!/usr/bin/env python3
"""把各步骤的结果收拢成一份结构化摘要 → results/<样本>/summary.json + summary.md
summary.md 是给写报告的人（或 AI）看的“数据底稿”：每一节的关键数字、表格和对应的原始文件，不含任何解读措辞。
写报告时只用这里的数字；需要细节时按“来源”回到原始文件查。哪个分析没跑，对应一节就写“未运行”。
用法：collect_results.py [sample]
"""
import csv
import datetime
import glob
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402

S = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
R = f"results/{S}"
os.makedirs(R, exist_ok=True)
OUT = {"sample": S, "generated": datetime.datetime.now().isoformat(timespec="seconds"), "pop": wgs.POP, "sex": wgs.sex_of(S)}
MD = []
PATHO = ("P/LP(ClinVar)", "LP(预测功能缺失)")


def tsv(p, skip=0):
    if not os.path.exists(p):
        return None
    with open(p) as f:
        for _ in range(skip):
            f.readline()
        return list(csv.DictReader(f, delimiter="\t"))


def js(p):
    return json.load(open(p)) if os.path.exists(p) else None


def section(title, body, src=()):
    MD.append(f"## {title}\n")
    MD.append(body.rstrip() + "\n")
    if src:
        MD.append("来源：" + "、".join(f"`{x}`" for x in src) + "\n")


def table(rows, cols, head=None):
    if not rows:
        return "（无）\n"
    head = head or cols
    s = "| " + " | ".join(head) + " |\n|" + "---|" * len(cols) + "\n"
    for r in rows:
        s += "| " + " | ".join(str(r.get(c, "")).replace("|", "/").replace("\n", " ") for c in cols) + " |\n"
    return s


def short_variant(r):
    p = (r.get("HGVSp") or "").split(":")[-1]
    c = (r.get("HGVSc") or "").split(":")[-1]
    return p if p and p != "." else c


# ---------------------------------------------------------------- 0. 版本
env = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib", "env.sh")).read()
OUT["tools"] = dict(re.findall(r"^IMG_(\w+)=\S*?([^/\s]+?)(?:\}|$)", env, re.M))
try:
    OUT["pipeline_commit"] = subprocess.run(["git", "-C", os.path.dirname(os.path.realpath(__file__)), "rev-parse", "--short", "HEAD"],
                                            capture_output=True, text=True).stdout.strip()
except OSError:
    pass

# ---------------------------------------------------------------- 1. 数据质量
qc = (js("work/findings/qc.json") or {}).get(S)
OUT["qc"] = qc
if qc:
    keys = [("raw_bases_G", "原始数据量（Gb）"), ("raw_read_pairs", "读段对数"), ("read_length", "读长"), ("Q30", "Q30（%）"),
            ("GC", "GC（%）"), ("platform_guess", "测序平台（按读段名推测）"), ("lanes", "泳道数"),
            ("pct_aligned", "比对率（%）"), ("unmapped_reads", "未比对读段"), ("dup_pct", "重复率（%）"),
            ("optical_dup_pct_of_pairs", "光学重复（%）"), ("insert_median", "插入片段中位数（bp）"),
            ("mean_cov_autosomal", "常染色体平均深度（×）"), ("pct_autosome_ge10x", "≥10× 比例（%）"), ("pct_autosome_ge20x", "≥20× 比例（%）"),
            ("mean_cov_chrX", "X 深度"), ("mean_cov_chrY", "Y 深度"), ("sex_inferred", "推断性别"),
            ("freemix_contamination_pct", "污染估计 FREEMIX（%）"), ("error_rate", "测序错误率"), ("mtDNA_copies_per_cell", "线粒体拷贝数 / 细胞"),
            ("n_snv", "SNV 数"), ("n_indel", "插入缺失数"), ("n_het", "杂合"), ("n_hom_alt", "纯合变异"), ("het_hom_ratio", "杂合 / 纯合比"),
            ("ti_tv", "Ti/Tv")]
    section("1. 数据质量", table([{"k": lab, "v": qc.get(k, "")} for k, lab in keys if k in qc], ["k", "v"], ["指标", "数值"]),
            ["work/findings/qc.json"])
else:
    section("1. 数据质量", "未运行（scripts/interpret/qc_summary.py）")

# ---------------------------------------------------------------- 2. 小变异
F = "work/findings"
fs = (js(f"{F}/summary.json") or {}).get(S)
OUT["findings_summary"] = fs
cols = ["gene", "HGVSp_short", "zygosity", "tier", "ClinVar", "ClinVar_stars", "gnomAD_AF", "max_AF", "gencc_moi", "gencc_diseases", "qc_note"]
head = ["基因", "变异", "合子", "分级", "ClinVar", "星级", "gnomAD AF", "最高人群 AF", "遗传方式", "相关疾病", "读段复核"]
body = ""
for key, title in (("acmg_sf", "A. ACMG SF v3.3 次要发现（可干预的显性病基因等）"), ("carrier", "B. 隐性 / X 连锁携带"),
                   ("biallelic", "同一隐性基因两个致病等位（需要关注）"), ("risk_alleles", "风险等位 / 低外显"),
                   ("nmd", "D. 神经肌肉病基因罕见变异（含意义未明）")):
    rows = tsv(f"{F}/{S}.{key}.tsv")
    if rows is None:
        continue
    for r in rows:
        r["HGVSp_short"] = short_variant(r)
        r["gencc_diseases"] = (r.get("gencc_diseases") or "")[:60]
    OUT.setdefault("findings", {})[key] = rows
    body += f"### {title}：{len(rows)} 条\n\n" + table(rows[:40], cols, head) + ("\n（只列前 40 条）\n" if len(rows) > 40 else "") + "\n"
for p in sorted(glob.glob(f"{F}/pair.*.tsv")):
    rows = tsv(p) or []
    for r in rows:
        r["HGVSp_short"] = short_variant(r)
    OUT.setdefault("pair", {})[os.path.basename(p)] = rows
    body += f"### 两人联合：{os.path.basename(p)}：{len(rows)} 条\n\n" + table(rows[:30], ["sample"] + cols[:6], ["样本"] + head[:6]) + "\n"
ko = f"{F}/{S}.knockouts.txt"
if os.path.exists(ko):
    OUT["knockouts"] = open(ko).read().strip()
    body += "### 纯合功能缺失基因（“人类敲除”）\n\n```\n" + OUT["knockouts"][:3000] + "\n```\n"
section("2. 小变异解读（SNV / 插入缺失）", body or "未运行（11_findings.sh）",
        [f"{F}/{S}.*.tsv", f"{F}/all_candidates.tsv", f"results/{S}/evidence_igv.html"])

# ---------------------------------------------------------------- 3. 结构变异
sv = tsv(f"{F}/{S}.sv.tsv")
if sv is not None:
    yes = [r for r in sv if r["report"] == "yes"]
    rev = [r for r in sv if r["report"] == "review"]
    OUT["sv"] = {"report": yes, "review_count": len(rev), "review": rev[:50]}
    section("3. 结构变异 / CNV", f"报告级 {len(yes)} 条；复核级 {len(rev)} 条（证据或临床相关性不足，只放附录）\n\n" +
            table(yes, ["type", "region", "size", "genes", "acmg_class", "pop_af", "callers", "note"],
                  ["类型", "区域", "大小", "基因", "AnnotSV ACMG", "人群频率", "支持方法", "说明"]), [f"{F}/{S}.sv.tsv"])
else:
    section("3. 结构变异 / CNV", "未运行（15_sv.sh，可选）")

# ---------------------------------------------------------------- 4. 重复扩增
str_rows = tsv(f"work/{S}/str/{S}.str_summary.tsv")
if str_rows is not None:
    OUT["str"] = str_rows
    flag = [r for r in str_rows if any(k in (r.get("status") or "").lower() for k in ("pre_mutation", "full_mutation"))]
    rf = js(f"work/{S}/str/{S}.rfc1_motif.json")
    OUT["rfc1"] = rf
    section("4. 重复扩增（ExpansionHunter + stranger）",
            f"共 {len(str_rows)} 个位点；进入前突变 / 致病范围的：{len(flag)} 个\n\n" + table(flag or str_rows[:0], list(str_rows[0].keys()) if str_rows else []) +
            ("\nRFC1 重复单元：" + json.dumps(rf, ensure_ascii=False)[:600] + "\n" if rf else ""),
            [f"work/{S}/str/{S}.str_summary.tsv", f"work/{S}/str/reviewer/"])
else:
    section("4. 重复扩增", "未运行（05_str.sh）")

# ---------------------------------------------------------------- 5. 难测基因
sp = {}
for key, p in (("SMN1_SMN2", f"work/{S}/special/{S}.smn.tsv"), ("CYP2D6", f"work/{S}/special/{S}.cyp2d6.tsv"),
               ("GBA1", f"work/{S}/special/{S}.gba.tsv")):
    rows = tsv(p)
    if rows:
        sp[key] = rows[0]
hba = f"work/{S}/special/{S}.hba.txt"
if os.path.exists(hba):
    sp["HBA_alpha_thalassemia_depth"] = open(hba).read().strip()
OUT["special_loci"] = sp
body = ""
if "SMN1_SMN2" in sp:
    r = sp["SMN1_SMN2"]
    body += f"- SMN1 拷贝数 {r.get('SMN1_CN')}，SMN2 拷贝数 {r.get('SMN2_CN')}；SMA 携带者 {r.get('isCarrier')}，SMA 患者 {r.get('isSMA')}\n"
if "CYP2D6" in sp:
    body += f"- CYP2D6（Cyrius）：{sp['CYP2D6'].get('Genotype')}（{sp['CYP2D6'].get('Filter')}）\n"
if "GBA1" in sp:
    r = sp["GBA1"]
    body += f"- GBA1（Gauchian）：GBA+GBAP1 拷贝数 {r.get('CN(GBA+GBAP1)')}，携带 {r.get('is_carrier(GBAP1-like_variant_exon9-11)')}，其他变异 {r.get('other_unphased_variants')}\n"
if "HBA_alpha_thalassemia_depth" in sp:
    body += "- α-地贫缺失（读深法）：\n\n```\n" + sp["HBA_alpha_thalassemia_depth"] + "\n```\n"
section("5. 普通流程测不准的基因", body or "未运行（06_special_loci.sh）", [f"work/{S}/special/"])

# ---------------------------------------------------------------- 6. HLA / KIR
def t1k(p):
    out = []
    if not os.path.exists(p):
        return None
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        a1 = f[2] if f[2] != "." else ""
        a2 = f[5] if len(f) > 5 and f[5] != "." else ""
        out.append(dict(gene=f[0], allele1=a1.split(",")[0], allele2=a2.split(",")[0], q1=f[4], q2=f[7] if len(f) > 7 else ""))
    return out


hla, kir = t1k(f"work/{S}/hla/{S}_genotype.tsv"), t1k(f"work/{S}/kir/{S}_genotype.tsv")
OUT["hla"], OUT["kir"] = hla, kir
body = ""
if hla:
    body += "### HLA（T1K，IPD-IMGT/HLA）\n\n" + table([r for r in hla if r["allele1"]], ["gene", "allele1", "allele2"], ["基因", "等位基因 1", "等位基因 2"])
if kir:
    pres = [r["gene"] for r in kir if r["allele1"]]
    body += "\n### KIR（T1K，IPD-KIR）\n\n有：" + "、".join(pres) + "\n\n无：" + "、".join(r["gene"] for r in kir if not r["allele1"]) + "\n"
section("6. HLA / KIR", body or "未运行（07_hla_kir.sh）", [f"work/{S}/hla/{S}_genotype.tsv", f"work/{S}/kir/{S}_genotype.tsv"])

# ---------------------------------------------------------------- 7. 线粒体 / Y
mt = {}
hg = f"work/{S}/mito/{S}.haplogrep3.txt"
if os.path.exists(hg):
    rows = list(csv.DictReader(open(hg), delimiter="\t"))
    if rows:
        r = rows[0]
        mt = {"haplogroup": r.get("Haplogroup"), "quality": r.get("Quality"), "top_hits": [x.get("Haplogroup") for x in rows[:3]]}
het = []
mv = f"work/{S}/mito/{S}.mutserve.vcf.gz"
# 线粒体里比对 / 测序伪影集中的位置：poly-C 区（D-loop 302–316、513–525、16180–16195）、3107 占位 N；插入缺失也不可靠
MT_NOISY = [(302, 316), (513, 525), (3105, 3109), (16180, 16195)]
if os.path.exists(mv):
    try:
        import pysam
        pysam.set_verbosity(0)
        for rec in pysam.VariantFile(mv):
            af = rec.samples[0].get("AF")
            af = af[0] if isinstance(af, tuple) else af
            alt = (rec.alts or ("",))[0]
            if af is None or not 0.03 <= float(af) < 0.9 or len(rec.ref) != 1 or alt not in ("A", "C", "G", "T"):
                continue
            if any(a <= rec.pos <= b for a, b in MT_NOISY):
                continue
            het.append(dict(pos=rec.pos, ref=rec.ref, alt=alt, level=round(float(af), 3)))
    except Exception as e:  # noqa: BLE001
        het = [{"error": str(e)}]
mt["heteroplasmies"] = het
OUT["mito"] = mt
y = {}
for p in sorted(glob.glob(f"work/{S}/ychr/yleaf/hg_prediction*.hg")):
    rows = tsv(p)
    if rows:
        y[os.path.basename(p)] = rows[0]
OUT["ychr"] = y
body = ""
if mt.get("haplogroup"):
    body += f"- 母系（线粒体）单倍群：{mt['haplogroup']}（HaploGrep3 质量 {mt['quality']}；候选 {', '.join(mt['top_hits'])}）\n"
    body += f"- 异质性 SNV（水平 3%–90%，mutserve；已排除 poly-C 等伪影区和插入缺失）：{len(het)} 个 " + json.dumps(het[:10], ensure_ascii=False) + "\n"
if y:
    body += "- 父系（Y）单倍群（Yleaf）：" + json.dumps(y.get("hg_prediction_combined.hg") or next(iter(y.values())), ensure_ascii=False) + "\n"
elif OUT["sex"] != "male":
    body += "- Y 染色体：非男性样本，不适用\n"
section("7. 母系 / 父系单倍群", body or "未运行（08_mito.sh / 09_ychr.sh）", [hg, f"work/{S}/ychr/yleaf/"])

# ---------------------------------------------------------------- 8. 药物基因组
pg = tsv(f"work/{S}/pgx/{S}.report.tsv", skip=1)
if pg is not None:
    keep = ["Gene", "Source Diplotype", "Phenotype", "Activity Score", "Outside Call", "Missing positions"]
    OUT["pgx"] = [{k: r.get(k, "") for k in keep} for r in pg]
    called = [r for r in OUT["pgx"] if r["Source Diplotype"] and "No Result" not in r["Source Diplotype"]]
    section("8. 药物基因组（PharmCAT）", f"有结果的基因 {len(called)} 个\n\n" +
            table(called, keep, ["基因", "双倍型", "表型", "活性评分", "外部分型", "缺失位点"]) +
            "\n用药建议见 PharmCAT 报告（按 CPIC / DPWG / FDA 指南）。\n", [f"work/{S}/pgx/{S}.report.html", f"work/{S}/pgx/{S}.report.json"])
else:
    section("8. 药物基因组", "未运行（12_pgx.sh）")

# ---------------------------------------------------------------- 9. 性状位点
tr = [r for r in (tsv(f"{F}/traits.tsv") or []) if r["sample"] == S]
trj = (js(f"{F}/traits.json") or {}).get(S)
OUT["traits"], OUT["traits_summary"] = tr, trj
section("9. 性状与单位点", (f"APOE：{trj.get('APOE')}；ABO（粗略）：{trj.get('ABO')}\n\n" if trj else "") +
        table(tr, ["rsid", "gene", "topic", "genotype", "effect_allele", "effect_allele_count", "source", "note"],
              ["rsID", "基因", "主题", "基因型", "效应等位", "效应等位个数", "来源", "说明"]), [f"{F}/traits.tsv"])

# ---------------------------------------------------------------- 10. 祖源
anc = {}
sa = tsv("work/ancestry/ancestry.somalier-ancestry.tsv")
if sa:
    me = [r for r in sa if r["#sample_id"] == S]
    if me:
        anc["somalier"] = {k: me[0][k] for k in ("predicted_ancestry", "AFR_prob", "EUR_prob", "EAS_prob", "AMR_prob", "SAS_prob")}
pc = js("work/ancestry/pca_summary.json")
if pc:
    anc["pca"] = {k: v["samples"].get(S) for k, v in pc.items()}
rel = tsv("work/ancestry/relate.pairs.tsv")
if rel:
    anc["relatedness"] = [r for r in rel if S in (r.get("#sample_a"), r.get("sample_b"))]
OUT["ancestry"] = anc
body = ""
if "somalier" in anc:
    body += f"- 大洲人群（somalier，对照 1000 Genomes）：{json.dumps(anc['somalier'], ensure_ascii=False)}\n"
for k, v in (anc.get("pca") or {}).items():
    if v:
        body += f"- PCA {k}：最近的人群中心 {v['nearest']}；距离 {v['dist']}\n"
for r in anc.get("relatedness") or []:
    body += f"- 亲缘：{r.get('#sample_a')}–{r.get('sample_b')} relatedness={r.get('relatedness')} IBS0={r.get('ibs0')}\n"
section("10. 祖源（大洲 / 人群）", body or "未运行（13_ancestry.sh）", ["work/ancestry/"])

# ---------------------------------------------------------------- 11. PRS
prs = js("work/prs/summary.json")
if prs:
    rows = [dict(trait=r["trait"], pgs=r["pgs_id"], pct=r["pct"].get(S), z=r["z"].get(S), n_used=r["n_used"], n_weights=r["n_weights"],
                 note=r["note"], sex_specific=r.get("sex_specific") or "") for r in prs if S in r["pct"]]
    rows = [r for r in rows if r["sex_specific"] in ("", OUT["sex"])]
    OUT["prs"] = rows
    section(f"11. 多基因风险评分（对照 1000 Genomes {wgs.POP}）", table(sorted(rows, key=lambda r: -r["pct"]),
            ["trait", "pct", "z", "pgs", "n_used", "n_weights", "note"], ["性状", "百分位", "Z", "PGS", "用到位点", "评分位点", "评分来源"]),
            ["work/prs/summary.json"])
else:
    section("11. 多基因风险评分", "未运行（14_prs.sh）")

# ---------------------------------------------------------------- 12. 扩展分析
E = "work/ext"
ext = {}
lai = js(f"{E}/lai/summary.{S}.json")
if lai:
    ext["lai"] = {"fractions": lai["fractions"], "controls_mean": {p: {g: round(sum(v) / len(v), 3) for g, v in d.items()} for p, d in lai["controls"].items()}}
f3 = tsv(f"{E}/aadr/f3.{S}.tsv")
if f3:
    a = [r for r in f3 if float(r["bp"] or 0) > 0 and int(r["nsnp"]) >= 30000][:15]
    m = [r for r in f3 if float(r["bp"] or 0) == 0 and int(r["nsnp"]) >= 300000][:10]
    ext["f3_top_ancient"] = [{k: r[k] for k in ("group", "f3", "se", "nsnp", "n", "bp", "country", "locality")} for r in a]
    ext["f3_top_modern"] = [{k: r[k] for k in ("group", "f3", "se", "n", "country")} for r in m]
qp = tsv(f"{E}/aadr/qpadm/results.tsv")
if qp:
    ext["qpadm"] = [r for r in qp if r["target"] == S]
    ext["qpadm_compare"] = [r for r in qp if r["target"] != S and r["right"] == "base"]
qpi = tsv(f"{E}/aadr/qpadm/results_indiv.tsv")
if qpi:
    ext["qpadm_indiv"] = qpi
arc = js(f"{E}/archaic/summary.json")
if arc and S in arc["samples"]:
    ext["archaic"] = arc["samples"][S]
roh = js(f"{E}/roh/summary.json")
if roh and S in roh["samples"]:
    ext["roh"] = roh["samples"][S]
    ext["roh_reference_median"] = roh["reference_median"]
for k, p in (("aging", f"{E}/aging/aging.{S}.json"), ("bloodgroups", f"{E}/bloodgroup/{S}.json"), ("cnv_genes", f"{E}/cnv_genes.{S}.json"),
             ("genome_stats", f"{E}/genome_stats.{S}.json"), ("virome", f"{E}/virome/{S}.summary.json")):
    v = js(p)
    if v is not None:
        if k == "cnv_genes":
            v = {kk: vv for kk, vv in v.items() if kk != "profiles"}
        ext[k] = v
OUT["ext"] = ext
body = ""
if "lai" in ext:
    body += f"### 局部祖源\n\n样本：{ext['lai']['fractions']}\n\n对照人群平均：{json.dumps(ext['lai']['controls_mean'], ensure_ascii=False)}\n\n"
if "f3_top_ancient" in ext:
    body += "### outgroup f3 最近的古代人群\n\n" + table(ext["f3_top_ancient"], ["group", "f3", "se", "bp", "n", "country", "locality"]) + "\n"
if "qpadm" in ext:
    body += "### qpAdm（样本）\n\n" + table(ext["qpadm"], [k for k in ext["qpadm"][0] if k not in ("pop",)]) + "\n"
if "qpadm_indiv" in ext:
    body += "### qpAdm 个体级对照（单人当目标时的离散程度）\n\n" + table(ext["qpadm_indiv"], list(ext["qpadm_indiv"][0].keys())) + "\n"
if "archaic" in ext:
    x = ext["archaic"]
    body += (f"### 古人类片段\n\n共 {x['n_segments']} 段、{x['total_Mb']} Mb；尼安德特型 {x['nea_Mb']} Mb，丹尼索瓦型 {x['den_Mb']} Mb；"
             f"对照百分位 {x.get('percentile_vs_controls')}\n\n著名渗入基因上的片段：{json.dumps(x['known_gene_hits'], ensure_ascii=False)}\n\n")
if "roh" in ext:
    x = ext["roh"]
    body += f"### 纯合片段\n\n≥1 Mb 总长 {x['sum1']} Mb（{x['n1']} 段），≥5 Mb 总长 {x['sum5']} Mb，最长 {x['longest']} Mb；在参考各人群中的百分位 {x['percentile_by_pop']}\n\n"
for k, title in (("aging", "衰老指标"), ("bloodgroups", "红细胞血型"), ("cnv_genes", "拷贝数会变的基因"), ("genome_stats", "基因组概览"),
                 ("virome", "未比对读段分类")):
    if k in ext:
        v = ext[k]
        if k == "aging":
            v = {kk: vv for kk, vv in v.items() if kk != "chip_candidates"} | {"chip_candidates_n": len(v.get("chip_candidates", [])),
                                                                               "chip_candidates": v.get("chip_candidates", [])[:10]}
        body += f"### {title}\n\n```json\n{json.dumps(v, ensure_ascii=False, indent=1)[:2500]}\n```\n\n"
section("12. 扩展分析（祖源细节、古 DNA、古人类、ROH、衰老、血型、CNV、病毒）", body or "未运行（scripts/ext/）", [f"{E}/"])

# ---------------------------------------------------------------- 13. 图
figs = js(f"{R}/figs/index.json")
OUT["figures"] = figs
section("13. 图", "\n".join(f"- {k}：" + "、".join(f"`{p}`" for p in v) for k, v in (figs or {}).items()) or "未生成（scripts/report/make_figures.py）")

json.dump(OUT, open(f"{R}/summary.json", "w"), ensure_ascii=False, indent=1, default=str)
open(f"{R}/summary.md", "w").write(f"# {S} 结果底稿（{OUT['generated'][:10]}）\n\n参考人群 {wgs.POP}；性别 {OUT['sex']}\n\n" + "\n".join(MD))
print(f"写出 {R}/summary.json、{R}/summary.md")
