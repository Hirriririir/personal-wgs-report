#!/usr/bin/env python3
"""α-珠蛋白基因簇（chr16p13.3）读深分析，筛查中国常见的 α-地贫缺失。

坐标为 GRCh38（= hg19 坐标减 49,999，此段两版本平移一致）：
  HBA2 chr16:172,876-173,710   HBA1 chr16:176,680-177,522   HBQ1 chr16:180,459-181,302
  --SEA  ≈ chr16:165,400-184,700（19.3 kb，HBA2+HBA1+HBQ1 全缺）
  -α3.7  ≈ chr16:173,300-177,100（右向缺失，两个 Z 盒之间）
  -α4.2  ≈ chr16:169,820-174,030（左向缺失，含 HBA2）
HBA1/HBA2 高度同源，比对质量常为 0，所以这里统计所有比对上的读段（含 MAPQ=0，去掉重复 / 次要比对），
按两侧单拷贝区（150–162 kb、190–220 kb）的中位深度归一化。二倍体正常 ≈1.0，杂合缺失 ≈0.5。
深度法只能提示，阳性需 gap-PCR / MLPA 等临床方法确认。

用法：hba_depth.py sample.bam out_prefix
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import sys
import statistics
import pysam

bam_path, out = sys.argv[1:3]
CHROM, START, END = "chr16", 140_000, 230_000
depth = [0] * (END - START)
bam = pysam.AlignmentFile(bam_path)
for r in bam.fetch(CHROM, START, END):
    if r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_unmapped:
        continue
    for a, b in r.get_blocks():
        for p in range(max(a, START), min(b, END)):
            depth[p - START] += 1

def mean(a, b):
    return statistics.fmean(depth[a - START:b - START])

flank = statistics.median(depth[150_000 - START:162_000 - START] + depth[190_000 - START:220_000 - START])
segments = [
    ("SEA 5' 段（4.2/3.7 之外）", 165_500, 169_700),
    ("-α4.2 段（含 HBA2 上游）", 169_900, 172_800),
    ("HBA2 基因", 172_876, 173_710),
    ("HBA2–HBA1 间（-α3.7 缺失段）", 173_800, 176_600),
    ("HBA1 基因", 176_680, 177_522),
    ("SEA 3' 段（含 HBQ1）", 177_700, 184_600),
    ("对照：SEA 外 3'", 185_000, 189_000),
]
if flank < 5:
    open(out + ".hba.txt", "w").write(f"flank_median_depth\t{flank:.1f}\ncall\t深度不足，无法判断\n")
    sys.exit("flank depth too low")
lines = [f"flank_median_depth\t{flank:.1f}"]
for name, a, b in segments:
    m = mean(a, b)
    lines.append(f"{name}\tchr16:{a}-{b}\tmean_depth={m:.1f}\tratio={m / flank:.2f}")
r = {name: mean(a, b) / flank for name, a, b in segments}
calls = []
if r["SEA 5' 段（4.2/3.7 之外）"] < 0.7 and r["SEA 3' 段（含 HBQ1）"] < 0.7:
    calls.append("疑似 --SEA 或其他大缺失（两侧段都降到约一半）")
if r["HBA2–HBA1 间（-α3.7 缺失段）"] < 0.7 and r["SEA 5' 段（4.2/3.7 之外）"] > 0.8:
    calls.append("疑似 -α3.7（两基因之间降低）")
if r["-α4.2 段（含 HBA2 上游）"] < 0.7 and r["SEA 5' 段（4.2/3.7 之外）"] > 0.8:
    calls.append("疑似 -α4.2（HBA2 上游降低）")
if max(r.values()) > 1.35:
    calls.append("某段读深升高，可能有 αααanti3.7 / αααanti4.2 三联体")
lines.append("call\t" + ("；".join(calls) if calls else "未见常见 α-地贫缺失的读深信号"))
open(out + ".hba.txt", "w").write("\n".join(lines) + "\n")
with open(out + ".hba_depth_500bp.tsv", "w") as f:
    f.write("pos\tdepth\tratio\n")
    for p in range(START, END, 500):
        m = mean(p, p + 500)
        f.write(f"{p}\t{m:.1f}\t{m / flank:.3f}\n")
print("\n".join(lines))
