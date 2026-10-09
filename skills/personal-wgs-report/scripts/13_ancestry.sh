#!/usr/bin/env bash
# 样本身份 / 亲缘 / 大洲祖源（多人步骤：工作目录里所有样本一起算）
#   1) somalier relate：两人以上时算亲缘系数（无血缘 ≈0，亲子 / 同胞 ≈0.5）、核对性别
#   2) somalier ancestry：对照 1000 Genomes 2504 人判断大洲人群
#   3) plink2 PCA：1000 Genomes 高深度 3202 人、LD 修剪后的常见 SNP，全球 + 参考人群（config 的 POP）两套主成分，把样本投影上去
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
SAMPLES=$(all_samples); pop=$(echo "$POP" | tr 'A-Z' 'a-z')
O=work/ancestry; P=ref/pca; mkdir -p $O
SOMS=$(for s in $SAMPLES; do echo work/somalier/$s.somalier; done)
# ---- 1) relate
if [ "$(echo $SAMPLES | wc -w)" -ge 2 ]; then
  { printf "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\n"; for s in $SAMPLES; do printf "fam\t%s\t0\t0\t0\t-9\n" $s; done; } > $O/cohort.ped
  $D $IMG_SOMALIER somalier relate --ped $O/cohort.ped -o $O/relate $SOMS > $O/relate.log 2>&1
fi
# ---- 2) ancestry
$D $IMG_SOMALIER sh -c "somalier ancestry --labels ref/qc/somalier/ancestry-labels-1kg.tsv ref/qc/somalier/1kg-somalier/*.somalier ++ $SOMS -o $O/ancestry" > $O/ancestry.log 2>&1
# ---- 3) PCA 参考面板（只做一次）
PL="$D $IMG_PLINK2 plink2"
if [ ! -s $P/kg_pruned.pgen ]; then
  # 常染色体、双等位 SNP、MAF≥5%，去掉二级以内亲属；ID 统一成 chr:pos:ref:alt
  $PL --pfile $P/all_hg38 vzs --remove $P/deg2_hg38.king.cutoff.out.id --chr 1-22 --snps-only just-acgt \
     --max-alleles 2 --maf 0.05 --geno 0.02 --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1 truncate \
     --rm-dup exclude-all --exclude range $P/high_ld_hg38.txt --make-pgen --out $P/kg_qc --threads $THREADS > $O/plink_qc.log 2>&1
  $PL --pfile $P/kg_qc --indep-pairwise 500kb 0.1 --out $P/kg_prune --threads $THREADS >> $O/plink_qc.log 2>&1
  $PL --pfile $P/kg_qc --extract $P/kg_prune.prune.in --make-pgen --out $P/kg_pruned --threads $THREADS >> $O/plink_qc.log 2>&1
fi
[ -s $P/$pop.ids ] || awk -F'\t' -v p="$POP" 'NR==1{for(i=1;i<=NF;i++) if($i=="SuperPop") c=i} NR>1 && $c==p{print $1}' $P/all_hg38.psam > $P/$pop.ids
if [ ! -s $P/kg_${pop}_pruned.pgen ]; then   # 参考人群内部重新按 MAF / LD 筛
  $PL --pfile $P/kg_qc --keep $P/$pop.ids --maf 0.05 --indep-pairwise 500kb 0.1 --out $P/kg_${pop}_prune --threads $THREADS >> $O/plink_qc.log 2>&1
  $PL --pfile $P/kg_qc --keep $P/$pop.ids --extract $P/kg_${pop}_prune.prune.in --make-pgen --out $P/kg_${pop}_pruned --threads $THREADS >> $O/plink_qc.log 2>&1
fi
for set in kg_pruned kg_${pop}_pruned; do
  [ -s $P/$set.pca.eigenvec.allele ] || $PL --pfile $P/$set --nonfounders --freq counts --pca 10 allele-wts --out $P/$set.pca --threads $THREADS >> $O/plink_pca.log 2>&1
  # 参考样本自己也走一遍投影，用来把投影分数校准到 eigenvec 的尺度（pca_summary.py 用）
  [ -s $P/$set.selfproj.sscore ] || $PL --pfile $P/$set --read-freq $P/$set.pca.acount \
     --score $P/$set.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
     --score-col-nums 6-15 --out $P/$set.selfproj --threads $THREADS >> $O/plink_pca.log 2>&1
done
# ---- 样本投影：从 gVCF 取这些位点的基因型（参考区块 GQ≥20 记 0/0），REF/ALT 与参考面板一致
GV=$(for s in $SAMPLES; do echo "$s=work/$s/$s.dv.g.vcf.gz"; done)
for set in kg_pruned kg_${pop}_pruned; do
  grep -v "^#" $P/$set.pvar | cut -f1-5 > $O/$set.sites.tsv
  $PY scripts/interpret/gvcf_genotype_sites.py $O/$set.sites.tsv $O/cohort.$set.vcf.gz $GV > $O/cohort.$set.genotype.log
  $PL --vcf $O/cohort.$set.vcf.gz --make-pgen --out $O/cohort.$set >> $O/plink_proj.log 2>&1
  $PL --pfile $O/cohort.$set --read-freq $P/$set.pca.acount \
     --score $P/$set.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
     --score-col-nums 6-15 --out $O/cohort.$set.proj >> $O/plink_proj.log 2>&1
done
$PY scripts/interpret/pca_summary.py > $O/pca_summary.log 2>&1
log "ANCESTRY DONE"
