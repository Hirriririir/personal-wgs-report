#!/usr/bin/env bash
# 核心资源（约 40 GB）：质控、重复扩增目录、注释数据库、基因 - 疾病表、HLA / KIR 分型索引、线粒体参考。
# 可选：REVEL 需要手动下载（见下方提示）；PanelApp 神经肌肉面板拿不到也不影响主流程。
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
F="python3 scripts/setup/fetch.py url"
A=$(readlink -f scripts)/../assets
mkdir -p ref/{qc/somalier,qc/verifybamid2,str,annot/clinvar,annot/gnomad,annot/alphamissense,annot/revel,annot/gencode,genelists,vep,hla,kir,mt}
# --- 样本身份 / 污染
$F https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz ref/qc/somalier/sites.hg38.vcf.gz
$F https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv ref/qc/somalier/ancestry-labels-1kg.tsv
if [ ! -d ref/qc/somalier/1kg-somalier ]; then
  $F https://zenodo.org/records/3479773/files/1kg.somalier.tar.gz ref/qc/somalier/1kg.somalier.tar.gz
  tar -xzf ref/qc/somalier/1kg.somalier.tar.gz -C ref/qc/somalier
fi
for x in UD mu bed V; do
  $F https://raw.githubusercontent.com/Griffan/VerifyBamID/master/resource/1000g.phase3.100k.b38.vcf.gz.dat.$x ref/qc/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.$x
done
# --- 重复扩增位点目录
$F https://raw.githubusercontent.com/Clinical-Genomics/stranger/main/stranger/resources/variant_catalog_grch38.json ref/str/stranger_variant_catalog_grch38.json
$F https://raw.githubusercontent.com/Illumina/ExpansionHunter/master/variant_catalog/grch38/variant_catalog.json ref/str/illumina_variant_catalog_grch38.json
# --- 注释
$F https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz ref/annot/clinvar/clinvar.vcf.gz
$F https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi ref/annot/clinvar/clinvar.vcf.gz.tbi
$F https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv ref/annot/gnomad/gnomad.v4.1.constraint_metrics.tsv
$F https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz ref/annot/alphamissense/AlphaMissense_hg38.tsv.gz
$F https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz ref/annot/gencode/gencode.v50.annotation.gtf.gz
if [ ! -s ref/annot/revel/revel-v1.3_all_chromosomes.zip ] && [ ! -s ref/annot/revel/new_tabbed_revel_grch38.tsv.gz ]; then
  curl -fsSL --retry 3 -o ref/annot/revel/revel-v1.3_all_chromosomes.zip https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip \
    || { rm -f ref/annot/revel/revel-v1.3_all_chromosomes.zip
         echo "[可选] REVEL 自动下载失败：浏览器打开 https://sites.google.com/site/revelgenomics/downloads ，"
         echo "       下载 revel-v1.3_all_chromosomes.zip 放到 ref/annot/revel/ 再运行 prep_annotation.sh；不放也能跑，只是少一个错义变异打分。"; }
fi
# --- 基因 - 疾病 - 遗传方式
$F "https://search.thegencc.org/download/action/submissions-export-tsv" ref/genelists/gencc_submissions.tsv
$F https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv ref/genelists/clingen_dosage_GRCh38.tsv
cp $A/genelists/acmg_sf_v3.3.tsv $A/genelists/trait_snps.tsv ref/genelists/
for p in 3101 4092; do   # PanelApp Australia：肌病 / 神经肌肉超级面板（可选）
  [ -s ref/genelists/panelapp_aus_$p.json ] || curl -fsSL -A "Mozilla/5.0" -o ref/genelists/panelapp_aus_$p.json "https://panelapp-aus.org/api/v1/panels/$p/" \
    || { rm -f ref/genelists/panelapp_aus_$p.json; echo "[可选] PanelApp 面板 $p 下载失败（部分地区会被拦），跳过"; }
done
date > ref/genelists/DOWNLOADED_AT
# --- VEP 缓存（Ensembl + RefSeq merged，约 28 GB）
if [ ! -d ref/vep/homo_sapiens_merged/116_GRCh38 ]; then
  $F https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_merged_vep_116_GRCh38.tar.gz ref/vep/homo_sapiens_merged_vep_116_GRCh38.tar.gz
  tar -xzf ref/vep/homo_sapiens_merged_vep_116_GRCh38.tar.gz -C ref/vep && rm ref/vep/homo_sapiens_merged_vep_116_GRCh38.tar.gz*
fi
# --- HLA / KIR 分型索引（T1K；坐标用 GENCODE v50）
[ -s ref/annot/gencode/gencode.v50.gtf ] || gzip -dc ref/annot/gencode/gencode.v50.annotation.gtf.gz > ref/annot/gencode/gencode.v50.gtf
if [ ! -s ref/hla/hlaidx/hlaidx_dna_seq.fa ]; then
  $F https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla.dat.zip ref/hla/hla.dat.zip
  python3 -c "import zipfile; zipfile.ZipFile('ref/hla/hla.dat.zip').extract('hla.dat', 'ref/hla')"
  $D $IMG_T1K t1k-build.pl -o ref/hla/hlaidx -d ref/hla/hla.dat -g ref/annot/gencode/gencode.v50.gtf > ref/hla/t1k-build.log 2>&1
fi
if [ ! -s ref/kir/kiridx/kiridx_dna_seq.fa ]; then
  $F https://raw.githubusercontent.com/ANHIG/IPDKIR/Latest/kir.dat ref/kir/kir.dat
  $D $IMG_T1K t1k-build.pl -o ref/kir/kiridx -d ref/kir/kir.dat -g ref/annot/gencode/gencode.v50.gtf > ref/kir/t1k-build.log 2>&1
fi
# --- 线粒体参考（rCRS，就是 GRCh38 的 chrM）
[ -s ref/mt/chrM.fa ] || $D $IMG_SAMTOOLS samtools faidx $REF chrM -o ref/mt/chrM.fa
echo RESOURCES DONE
