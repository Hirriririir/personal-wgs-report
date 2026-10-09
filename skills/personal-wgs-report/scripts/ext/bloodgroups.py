#!/usr/bin/env python3
"""红细胞血型扩展分型：关键 SNP → 氨基酸 → 抗原（ABO 在 interpret/trait_snps.py 里）→ work/ext/bloodgroup/<s>.json
每个位点的 GRCh38 坐标与“等位基因→氨基酸”由 Ensembl REST（VEP）给出（只发 rsID，不发任何个人数据；缓存于 ref/genelists/bg_cache.json），
抗原规则只写“哪个氨基酸 = 哪个抗原”，从而避免链方向 / 参考等位基因写错。只作参考，输血配型以血库检测为准。
用法：bloodgroups.py <sample>   （需要能访问 rest.ensembl.org）"""
import json
import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402
import pysam  # noqa: E402

s = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
CACHE = "ref/genelists/bg_cache.json"
cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}

# 系统, 基因, rsID, 密码子位置（蛋白）, {氨基酸: 抗原}
RULES = [
    ("Duffy", "ACKR1", "rs12075", 42, {"G": "Fy(a)", "D": "Fy(b)"}),
    ("Duffy", "ACKR1", "rs2814778", None, None),  # GATA 盒 -67T>C：C 纯合 = 红细胞不表达 Duffy（非洲常见）
    ("Kidd", "SLC14A1", "rs1058396", 280, {"D": "Jk(a)", "N": "Jk(b)"}),
    ("Kell", "KEL", "rs8176058", 193, {"T": "k", "M": "K"}),
    ("Kell", "KEL", "rs8176059", 281, {"R": "Kp(b)", "W": "Kp(a)"}),
    ("MNS", "GYPB", "rs7683365", 48, {"M": "S", "T": "s"}),
    ("Diego", "SLC4A1", "rs2285644", 854, {"P": "Di(b)", "L": "Di(a)"}),
    ("Diego", "SLC4A1", "rs75731670", 658, {"E": "Wr(b)", "K": "Wr(a)"}),
    ("Rh", "RHCE", "rs609320", 226, {"A": "e", "P": "E"}),
    ("Dombrock", "ART4", "rs11276", 265, {"N": "Do(a)", "D": "Do(b)"}),
    ("Colton", "AQP1", "rs28362692", 45, {"A": "Co(a)", "V": "Co(b)"}),
    ("Yt", "ACHE", "rs1799805", 353, {"H": "Yt(a)", "N": "Yt(b)"}),
    ("Lutheran", "BCAM", "rs28399653", 77, {"R": "Lu(b)", "H": "Lu(a)"}),
    ("Junior", "ABCG2", "rs72552713", 126, {"Q": "Jr(a) 正常", "*": "Jr 无效等位（Q126*）"}),
    ("Lewis", "FUT3", "rs28362459", 20, {"L": "FUT3 正常", "R": "FUT3 失活（L20R）"}),
    ("Lewis", "FUT3", "rs3745635", 170, {"G": "FUT3 正常", "S": "FUT3 失活（G170S）"}),
    ("Lewis", "FUT3", "rs3894326", 356, {"I": "FUT3 正常", "K": "FUT3 失活（I356K）"}),
    ("Secretor", "FUT2", "rs1047781", 129, {"I": "FUT2 正常", "F": "FUT2 弱化（I129F，东亚弱分泌型）"}),
    ("Secretor", "FUT2", "rs601338", 143, {"W": "FUT2 正常", "*": "FUT2 失活（W143*）"}),
    ("RhD-DEL", "RHD", "rs549616139", None, None),  # RHD c.1227G>A 亚洲型 DEL
]


def rest(url):
    for _ in range(4):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.load(r)
        except Exception:  # noqa: BLE001
            time.sleep(3)
    return None


def info(rs, gene):
    if rs in cache:
        return cache[rs]
    v = rest(f"https://rest.ensembl.org/variation/human/{rs}?content-type=application/json")
    m = next((m for m in (v or {}).get("mappings", []) if m.get("assembly_name") == "GRCh38" and m["seq_region_name"].isalnum()
              and len(m["seq_region_name"]) <= 2), None)
    vep = rest(f"https://rest.ensembl.org/vep/human/id/{rs}?content-type=application/json&canonical=1&hgvs=1")
    aa = {}
    hg = ""
    for rec in vep or []:
        for tc in rec.get("transcript_consequences", []):
            if tc.get("gene_symbol") == gene and tc.get("canonical") and tc.get("amino_acids"):
                aa[tc["variant_allele"]] = tc["amino_acids"]
                hg = tc.get("hgvsp", "") or hg
    cache[rs] = {"chrom": "chr" + m["seq_region_name"] if m else None, "pos": m["start"] if m else None,
                 "end": m["end"] if m else None, "alleles": m["allele_string"] if m else None, "aa": aa, "hgvsp": hg}
    return cache[rs]


vf = pysam.VariantFile(f"work/{s}/{s}.dv.g.vcf.gz")
fa = pysam.FastaFile(wgs.REF)
out = []
for system, gene, rs, codon, amap in RULES:
    d = info(rs, gene)
    if not d["chrom"]:
        out.append(dict(system=system, gene=gene, rs=rs, genotype="坐标未知"))
        continue
    c, p = d["chrom"], d["pos"]
    ref = fa.fetch(c, p - 1, p).upper()
    gt = None
    q = None
    for rec in vf.fetch(c, p - 1, p):
        alts = rec.alts or ()
        smp = rec.samples[0]
        if alts in ((), ("<*>",)):
            if rec.start <= p - 1 < rec.stop:
                gt, q = (ref, ref), smp.get("GQ")
            continue
        if rec.pos == p and smp.get("GT") and None not in smp["GT"]:
            al = rec.alleles
            gt, q = tuple(al[i] for i in smp["GT"]), smp.get("GQ")
            break
    row = dict(system=system, gene=gene, rs=rs, pos=f"{c}:{p}", ref=ref, genotype="/".join(gt) if gt else "未覆盖", GQ=q,
               hgvsp=d.get("hgvsp", ""))
    if gt and amap:
        # 参考碱基对应的氨基酸 = aa 字段 "X/Y" 的左边（Ensembl 写法：ref/alt）
        any_aa = next(iter(d["aa"].values()), None)
        if any_aa:
            ref_aa = any_aa.split("/")[0]
            ants = []
            for b in gt:
                if b == ref:
                    ants.append(amap.get(ref_aa, f"?({ref_aa})"))
                else:
                    alt_aa = d["aa"].get(b, "/?").split("/")[-1]
                    ants.append(amap.get(alt_aa, f"?({alt_aa})"))
            row["antigens"] = " / ".join(ants)
    out.append(row)
json.dump(cache, open(CACHE, "w"), indent=1)
os.makedirs("work/ext/bloodgroup", exist_ok=True)
json.dump(out, open(f"work/ext/bloodgroup/{s}.json", "w"), ensure_ascii=False, indent=1)
for r in out:
    print(r)
