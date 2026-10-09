#!/usr/bin/env bash
# 古人类渗入（尼安德特人 / 丹尼索瓦人片段）：hmmix（Skov 2018）二倍体流程
#   观测：样本的 PASS SNP，去掉外群（1000G+HGDP 非洲人）里出现过的、只留衍生等位基因 → 训练 HMM → 解码片段
#   片段归类：与 4 个高覆盖古人类基因组（Altai、Vindija、Chagyrskaya、Denisova）共享的衍生变异数
#   对照：参考人群各亚群随机若干人（assets/population/<POP>.json 的 controls）跑完全相同的流程
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
O=work/ext/archaic; K=$O/kg; H=ref/hmmix
mkdir -p $O $K tools/bin-docker
# hmmix 读古人类 BCF 要调用 bcftools：用容器包一层放进 PATH
printf '#!/usr/bin/env bash\nexec "%s/scripts/dr.sh" %s bcftools "$@"\n' "$(pwd -P)" "$IMG_BCFTOOLS" > tools/bin-docker/bcftools
chmod +x tools/bin-docker/bcftools
export PATH=$(pwd -P)/tools/bin-docker:$PATH
HM=tools/venv-ext/bin/hmmix
ARGS="-weights=$H/hg38_strick_callability_mask.bed -mutrates=$H/hg38_mutationrate.bed"
decode() {  # decode <obs> <前缀>
  [ -s $2.diploid.txt ] && return 0
  $HM train -obs=$1 $ARGS -out=$2.trained.json > $2.train.log 2>&1
  $HM decode -obs=$1 $ARGS -param=$2.trained.json -admixpop="$H/archaic/individuals_highcov.*.bcf" -out=$2 > $2.decode.log 2>&1
}
export -f decode; export HM ARGS H
# 样本
for s in $(all_samples); do
  [ -s $O/obs.$s.txt ] || $PYB scripts/ext/hmmix_obs.py work/$s/$s.dv.vcf.gz $O/obs.$s.txt 2> $O/obs.$s.log
  decode $O/obs.$s.txt $O/decoded.$s
  log "hmmix $s 完成"
done
# 对照
if [ ! -s $K/ids.tsv ]; then
  $PYB - <<'PY'
import random, sys
sys.path.insert(0, "scripts/ext")
import popcfg
cfg = popcfg.load()["controls"]
by = {}
for i, (sp, p) in popcfg.kg_samples().items():
    by.setdefault(p, []).append(i)
random.seed(7)
out = []
for p in cfg["pops"]:
    ids = sorted(by.get(p, [])); random.shuffle(ids); out += [(i, p) for i in ids[:cfg.get("per_pop", 8)]]
open("work/ext/archaic/kg/ids.tsv", "w").write("".join(f"{i}\t{p}\n" for i, p in out))
open("work/ext/archaic/kg/ids.txt", "w").write("".join(f"{i}\n" for i, _ in out))
print(len(out), "个对照")
PY
fi
if [ ! -s $K/kg.vcf.gz.csi ]; then
  $D $IMG_PLINK2 plink2 --pfile ref/pca/all_hg38 vzs --keep $K/ids.txt --chr 1-22 --snps-only just-acgt --mac 1 \
      --export vcf bgz id-paste=iid --output-chr chrM --out $K/kg --threads $THREADS > $K/kg.log 2>&1
  $D $IMG_BCFTOOLS bcftools index -f $K/kg.vcf.gz
fi
[ -s $K/obs.done ] || { $PYB scripts/ext/hmmix_obs.py $K/kg.vcf.gz "$K/obs.{sample}.txt" 2> $K/obs.log && touch $K/obs.done; }
cut -f1 $K/ids.tsv | xargs -P ${JOBS:-4} -I{} bash -c 'decode work/ext/archaic/kg/obs.{}.txt work/ext/archaic/kg/decoded.{}'
$PYB scripts/ext/archaic_summary.py > $O/summary.log 2>&1 && cat $O/summary.log
log "ARCHAIC DONE"
