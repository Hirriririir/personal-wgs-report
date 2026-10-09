#!/usr/bin/env bash
# 小变异注释：只注释解读用得到的范围（GENCODE 外显子 ±50 bp + ClinVar 收录位点，约占全部变异的 5%），
#   VEP 116（Ensembl + RefSeq merged 缓存，含 gnomAD v4.1 人群频率、MANE 转录本）
#   + AlphaMissense + REVEL（下载了才用）+ NMD 插件 + ClinVar（下载当周的版本）
# 输入 work/joint/cohort.dv.norm.vcf.gz（04_cohort_genotype.sh）→ 输出 work/annot/cohort.vep.vcf.gz
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
O=work/annot; mkdir -p $O logs
W=$(pwd -P)
$D $IMG_BCFTOOLS sh -c "bcftools view -R ref/annot/regions/exons50_clinvar.bed work/joint/cohort.dv.norm.vcf.gz -Ou \
  | bcftools sort -Oz -o $O/cohort.regions.vcf.gz && bcftools index -f -t $O/cohort.regions.vcf.gz"
plugins=(--plugin "AlphaMissense,file=$W/ref/annot/alphamissense/AlphaMissense_hg38.tsv.gz" --plugin NMD)
REVEL=ref/annot/revel/new_tabbed_revel_grch38.tsv.gz
if [ -s $REVEL.tbi ]; then plugins+=(--plugin "REVEL,file=$W/$REVEL,no_match=1")
else log "没有 REVEL 文件，跳过 REVEL（错义变异仍有 AlphaMissense 打分）"; fi
log "VEP：$($D $IMG_BCFTOOLS sh -c "bcftools view -H $O/cohort.regions.vcf.gz | wc -l") 个位点"
$D $IMG_VEP vep --input_file $O/cohort.regions.vcf.gz --output_file $O/cohort.vep.vcf.gz --vcf --compress_output bgzip --force_overwrite \
  --offline --cache --dir_cache ref/vep --cache_version 116 --assembly GRCh38 --merged \
  --fasta $REF --fork $(( THREADS > 16 ? 16 : THREADS )) --buffer_size 20000 \
  --everything --flag_pick_allele_gene --mane --no_stats \
  --dir_plugins /plugins "${plugins[@]}" \
  --custom file=$W/ref/annot/clinvar/clinvar.chr.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF%ALLELEID%GENEINFO%CLNVC \
  > logs/vep.log 2>&1 || { echo "VEP 失败，见 logs/vep.log" >&2; exit 1; }
$D $IMG_BCFTOOLS tabix -f -p vcf $O/cohort.vep.vcf.gz
log "ANNOTATION DONE"
