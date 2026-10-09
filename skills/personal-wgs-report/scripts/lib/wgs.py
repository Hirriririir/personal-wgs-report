"""Python 脚本的公共入口：进入工作目录、读取配置（WGS_CONFIG，默认 config.env；KEY=VALUE），给出样本名等常用设置。

用法（脚本在 scripts/ 的任一子目录里）：
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
    import wgs
    s = wgs.SAMPLE
"""
import os
import shlex

WORKDIR = os.environ.get("WGS_WORKDIR", os.getcwd())
os.chdir(WORKDIR)

CONFIG_FILE = os.environ.get("WGS_CONFIG", "config.env")
CONFIG = {}
if os.path.exists(CONFIG_FILE):
    for line in open(CONFIG_FILE, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        v = v.strip()
        try:
            v = shlex.split(v, comments=True)
            v = v[0] if v else ""
        except ValueError:
            pass
        CONFIG[k.strip()] = v

SAMPLE = CONFIG.get("SAMPLE", "")
POP = CONFIG.get("POP", "EAS")
THREADS = int(CONFIG.get("THREADS") or os.cpu_count() or 4)
REF = "ref/GRCh38/GRCh38.fa"


def sex():
    """config 里写了 male/female 就用；否则读 03_bam_qc 推断的结果（ambiguous 返回 female）。"""
    s = CONFIG.get("SEX", "auto")
    if s in ("male", "female"):
        return s
    try:
        s = open(f"work/{SAMPLE}/qc/{SAMPLE}.sex.txt").read().split()[0]
    except (OSError, IndexError):
        return "female"
    return "male" if s == "male" else "female"


def annot_vcf():
    """VEP（+SpliceAI）注释后的 VCF：有 SpliceAI 版本就用它。"""
    p = "work/annot/cohort.vep.spliceai.vcf.gz"
    return p if os.path.exists(p + ".tbi") else "work/annot/cohort.vep.vcf.gz"


def all_samples():
    """工作目录里已有变异检测结果（work/<s>/<s>.dv.g.vcf.gz）的所有样本。"""
    import glob
    return sorted(os.path.basename(p)[:-len(".dv.g.vcf.gz")] for p in glob.glob("work/*/*.dv.g.vcf.gz"))


def sex_of(s):
    """任一样本的性别（当前配置的样本按 config 的 SEX；其他样本读推断结果）。"""
    if s == SAMPLE:
        return sex()
    try:
        x = open(f"work/{s}/qc/{s}.sex.txt").read().split()[0]
    except (OSError, IndexError):
        return "female"
    return "male" if x == "male" else "female"
