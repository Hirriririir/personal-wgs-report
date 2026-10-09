#!/usr/bin/env bash
# 多基因风险评分（PRS，多人步骤）：ref/pgs/selected.tsv 里的 PGS Catalog 评分（GRCh38 harmonized）
#   1) 每个评分 → 权重表 → 在 1000 Genomes 参考人群（config 的 POP）上打分：得到参照分布和实际用到的位点
#   2) 所有评分用到的位点取并集，从每个人的 gVCF 取基因型（参考区块 GQ≥20 记 0/0）
#   3) 样本按同一套位点打分 → prs_summary.py 算百分位
# 评分清单：默认 assets/pgs_selected_<POP>.tsv（仓库自带东亚 16 个）；其他人群先用 scripts/interpret/pgs_search.py 挑
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
O=work/prs; P=ref/pca; G=ref/pgs; mkdir -p $O $G
pop=$(echo "$POP" | tr 'A-Z' 'a-z')
LIST=${PGS_LIST:-$G/selected.tsv}
if [ ! -s $LIST ]; then
  A=$(readlink -f scripts)/../assets/pgs_selected_$POP.tsv
  [ -s "$A" ] || { echo "没有 $LIST：先运行 $PY scripts/interpret/pgs_search.py 挑评分，按 assets/pgs_selected_EAS.tsv 的格式写好" >&2; exit 1; }
  cp "$A" $LIST
fi
[ -s $P/$pop.ids ] || awk -F'\t' -v p="$POP" 'NR==1{for(i=1;i<=NF;i++) if($i=="SuperPop") c=i} NR>1 && $c==p{print $1}' $P/all_hg38.psam > $P/$pop.ids
PL="$D $IMG_PLINK2 plink2"
IDS=$(grep -v "^#" $LIST | cut -f2)
for id in $IDS; do
  f=$G/${id}_hmPOS_GRCh38.txt.gz
  python3 scripts/setup/fetch.py url "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/$id/ScoringFiles/Harmonized/${id}_hmPOS_GRCh38.txt.gz" $f > /dev/null \
    || { log "下载失败：$id，跳过"; continue; }
  [ -s $O/$id.ref.sscore ] && continue
  kind=$($PY scripts/interpret/prs_weights.py build $f $O/$id 2> $O/$id.weights.log)
  if [ "$kind" = pos ]; then
    # 只有效应等位基因：先按 chr:pos 在 1000G 双等位 SNP 上匹配 → 取用到位点的 REF/ALT → 去掉回文位点后重新按完整 ID 打分
    $PL --pfile $P/all_hg38 vzs --keep $P/$pop.ids --snps-only just-acgt --max-alleles 2 --set-all-var-ids '@:#' --rm-dup exclude-all \
       --score $O/$id.weights_pos.tsv 1 2 3 header no-mean-imputation cols=+scoresums list-variants --out $O/$id.pos --threads $THREADS > $O/$id.plink.log 2>&1
    $PL --pfile $P/all_hg38 vzs --snps-only just-acgt --max-alleles 2 --set-all-var-ids '@:#' --rm-dup exclude-all \
       --extract $O/$id.pos.sscore.vars --make-just-pvar --out $O/$id.used --threads $THREADS >> $O/$id.plink.log 2>&1
    $PY scripts/interpret/prs_weights.py finish-pos $O/$id 2>> $O/$id.weights.log
  fi
  $PL --pfile $P/all_hg38 vzs --keep $P/$pop.ids --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1 truncate \
     --score $O/$id.weights.tsv 1 2 3 header no-mean-imputation cols=+scoresums list-variants --out $O/$id.ref --threads $THREADS >> $O/$id.plink.log 2>&1 \
     || log "参考人群打分失败：$id（见 $O/$id.plink.log）"
  log "$id 参考人群打分完成（$kind）"
done
# 样本：所有评分用到的位点并集，一次取基因型
$PY scripts/interpret/prs_weights.py union $O/union.sites.tsv $IDS
GV=$(for s in $(all_samples); do echo "$s=work/$s/$s.dv.g.vcf.gz"; done)
$PY scripts/interpret/gvcf_genotype_stream_par.py $O/union.sites.tsv $O/union.cohort.vcf.gz $GV > $O/union.genotype.log
$PL --vcf $O/union.cohort.vcf.gz --make-pgen --out $O/union.cohort --threads $THREADS > $O/union.plink.log 2>&1
for id in $IDS; do
  [ -s $O/$id.weights.tsv ] || continue
  $PL --pfile $O/union.cohort --score $O/$id.weights.tsv 1 2 3 header no-mean-imputation cols=+scoresums \
    --out $O/$id.cohort --threads 4 >> $O/$id.plink.log 2>&1 || log "样本打分失败：$id"
done
$PY scripts/interpret/prs_summary.py > $O/prs_summary.log 2>&1
log "PRS DONE"
