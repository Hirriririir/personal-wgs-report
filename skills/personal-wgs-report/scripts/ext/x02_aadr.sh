#!/usr/bin/env bash
# 古 DNA 对照（AADR v66.p1，1240K 位点；参数在 assets/population/<POP>.json 的 aadr 段）：
#   1) 位点转 GRCh38、样本基因型编码（所有样本）；2) 挑参考地区的古今个体 + 外群
#   3) smartpca：现代人群定轴，古人和样本投影；4) outgroup f3：和哪些古代人群最近；5) qpAdm 混合模型（参数里有模型时）
# 用法：scripts/ext/x02_aadr.sh [prep|pca|f3|qpadm|all]
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
O=work/ext/aadr; mkdir -p $O/pca; STEP=${1:-all}
export IMG_ADMIXTOOLS IMG_EIGENSOFT
case $STEP in prep|all)
  [ -s $O/sites38.tsv ] || $PYB scripts/ext/aadr_prep.py lift
  for s in $(all_samples); do
    if [ ! -s $O/$s.codes.npy ]; then
      $PY scripts/interpret/gvcf_genotype_stream_par.py $O/sites38.tsv $O/$s.1240k.vcf.gz $s=work/$s/$s.dv.g.vcf.gz > $O/$s.1240k.log
      $PYB scripts/ext/aadr_prep.py encode $s
    fi
  done
  $PYB scripts/ext/aadr_subset.py ;;
esac
case $STEP in pca|all)
  $PYB - <<'PY'
import sys
sys.path.insert(0, "scripts/ext")
import popcfg, wgs
cfg = popcfg.load()["aadr"]
have = {x.split()[2] for x in open("work/ext/aadr/sub.ind")}
pops = [p for p in cfg["modern_pops"] if p in have]
open("work/ext/aadr/pca/modern.pops", "w").write("\n".join(pops) + "\n")
print("PCA 定轴人群", len(pops))
PY
  cat > $O/pca/pca.par <<PAR
genotypename: $O/sub.geno
snpname: $O/sub.snp
indivname: $O/sub.ind
evecoutname: $O/pca/pca.evec
evaloutname: $O/pca/pca.eval
poplistname: $O/pca/modern.pops
lsqproject: YES
numoutevec: 4
numoutlieriter: 0
hashcheck: NO
numthreads: $THREADS
PAR
  $D $IMG_EIGENSOFT smartpca -p $O/pca/pca.par > $O/pca/smartpca.log 2>&1 ;;
esac
case $STEP in f3|all) for s in $(all_samples); do $PYB scripts/ext/aadr_f3.py $s > $O/f3.$s.log 2>&1; tail -12 $O/f3.$s.log; done ;; esac
case $STEP in qpadm|all)
  $PYB scripts/ext/aadr_qpadm.py > $O/qpadm.log 2>&1 && $PYB scripts/ext/aadr_qpadm.py --indiv > $O/qpadm_indiv.log 2>&1 \
    || log "qpAdm 失败（见 $O/qpadm*.log）" ;;
esac
log "AADR $STEP DONE"
