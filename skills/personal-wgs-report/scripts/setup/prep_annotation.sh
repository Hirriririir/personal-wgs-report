#!/usr/bin/env bash
# 注释前的预处理：ClinVar 改 chr 命名；AlphaMissense / REVEL 建索引；注释范围 BED；基因 - 疾病表。
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
A=ref/annot
# ClinVar：1..22 X Y MT → chr…，去掉 NT_ 片段
if [ ! -s $A/clinvar/clinvar.chr.vcf.gz.tbi ]; then
  ( for c in $(seq 1 22) X Y; do echo -e "$c\tchr$c"; done; echo -e "MT\tchrM" ) > $A/clinvar/chr_rename.txt
  $D $IMG_BCFTOOLS sh -c "bcftools annotate --rename-chrs $A/clinvar/chr_rename.txt $A/clinvar/clinvar.vcf.gz \
    | bcftools view -t \$(cut -f2 $A/clinvar/chr_rename.txt | paste -sd,) -Oz -o $A/clinvar/clinvar.chr.vcf.gz && tabix -p vcf $A/clinvar/clinvar.chr.vcf.gz"
fi
# AlphaMissense：VEP 插件需要 tabix 索引
[ -s $A/alphamissense/AlphaMissense_hg38.tsv.gz.tbi ] || \
  $D $IMG_BCFTOOLS tabix -s 1 -b 2 -e 2 -f -S 1 $A/alphamissense/AlphaMissense_hg38.tsv.gz
# REVEL v1.3（可选）：官方 CSV → 制表符 → 取 GRCh38 坐标（第 3 列）排序 → bgzip + tabix
if [ ! -s $A/revel/new_tabbed_revel_grch38.tsv.gz.tbi ] && [ -s $A/revel/revel-v1.3_all_chromosomes.zip ]; then
  ( cd $A/revel
    python3 -c "import zipfile; zipfile.ZipFile('revel-v1.3_all_chromosomes.zip').extractall('.')"
    f=$(ls | grep -i -E "^revel.*(csv|ids)$" | head -1)
    head -1 "$f" | tr "," "\t" | sed "1s/^/#/" > h.tsv
    tail -n +2 "$f" | tr "," "\t" | awk -F"\t" '$3 != "."' | LC_ALL=C sort -S 8G --parallel=8 -k1,1 -k3,3n | cat h.tsv - > revel_grch38.tsv
    rm -f "$f" h.tsv )
  $D $IMG_BCFTOOLS sh -c "bgzip -@ 8 -c $A/revel/revel_grch38.tsv > $A/revel/new_tabbed_revel_grch38.tsv.gz && tabix -f -s 1 -b 3 -e 3 $A/revel/new_tabbed_revel_grch38.tsv.gz"
  rm -f $A/revel/revel_grch38.tsv
fi
[ -s $A/regions/exons50_clinvar.bed ] || python3 scripts/setup/build_regions.py
# 蛋白编码基因区间（扩展分析标注基因用）；克隆性造血（CHIP）常见 25 个基因的 MANE 编码区
G=$A/gencode/gencode.v50.annotation.gtf.gz
[ -s $A/gencode/pc_genes.bed ] || zcat $G | awk -F'\t' '$3=="gene" && $9 ~ /gene_type "protein_coding"/ {
  match($9, /gene_name "[^"]+"/); print $1"\t"$4-1"\t"$5"\t"substr($9, RSTART+11, RLENGTH-12)}' | sort -k1,1V -k2,2n > $A/gencode/pc_genes.bed
CHIP="ASXL1 BCOR BCORL1 CALR CBL CHEK2 DNMT3A EZH2 GNAS GNB1 IDH1 IDH2 JAK2 KRAS MPL NRAS PPM1D RUNX1 SF3B1 SRSF2 STAG2 TET2 TP53 U2AF1 ZRSR2"
[ -s $A/gencode/chip_cds.bed ] || zcat $G | awk -F'\t' -v genes="$CHIP" 'BEGIN{n=split(genes,a," "); for(i=1;i<=n;i++) keep[a[i]]=1}
  $3=="CDS" && $9 ~ /tag "MANE_Select"/ { match($9, /gene_name "[^"]+"/); g=substr($9, RSTART+11, RLENGTH-12); if (g in keep) print $1"\t"$4-1"\t"$5"\t"g }' \
  | sort -k1,1V -k2,2n -u > $A/gencode/chip_cds.bed
$PY scripts/interpret/build_gene_tables.py
echo PREP DONE
