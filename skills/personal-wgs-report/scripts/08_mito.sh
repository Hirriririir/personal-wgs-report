#!/usr/bin/env bash
# 线粒体：mutserve（同质 / 异质性变异，rCRS）+ GATK Mutect2 线粒体模式交叉验证 + HaploGrep3 单倍群
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
s=$SAMPLE
O=work/$s/mito; mkdir -p $O
$D $IMG_MUTSERVE mutserve call --reference ref/mt/chrM.fa --contig-name chrM --deletions --insertions \
  --output $O/$s.mutserve.vcf.gz --threads 4 --level 0.01 work/$s/$s.bam > $O/mutserve.log 2>&1
$D $IMG_GATK gatk --java-options "-Xmx8g" Mutect2 -R $REF -L chrM --mitochondria-mode \
  -I work/$s/$s.bam -O $O/$s.chrM.mutect2.raw.vcf.gz > $O/mutect2.log 2>&1
$D $IMG_GATK gatk FilterMutectCalls -R $REF --mitochondria-mode \
  -V $O/$s.chrM.mutect2.raw.vcf.gz -O $O/$s.chrM.mutect2.vcf.gz >> $O/mutect2.log 2>&1
$D $IMG_HAPLOGREP haplogrep3 classify --in $O/$s.mutserve.vcf.gz \
  --tree phylotree-fu-rcrs@1.2 --out $O/$s.haplogrep3.txt --extend-report --hits 3 > $O/haplogrep3.log 2>&1
echo "MITO DONE $s"
