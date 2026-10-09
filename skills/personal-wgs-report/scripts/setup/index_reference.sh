#!/usr/bin/env bash
# 比对索引与序列字典：
#   GPU（Parabricks fq2bam）用经典 BWA 索引（.bwt/.sa，约 1 小时），放在 ref/GRCh38/bwa/（Parabricks 会解析软链接，所以用硬链接）
#   CPU 用 bwa-mem2 索引（约 1–1.5 小时，内存约 70–90 GB）
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
W=$PWD
cd ref/GRCh38
M=GRCh38_no_alt_plus_hs38d1_GRCmasked
[ -s $M.dict ] || "$W/$D" $IMG_SAMTOOLS samtools dict -a GRCh38 -s "Homo sapiens" -o $M.dict $M.fa
ln -sfn $M.dict GRCh38.dict
if use_gpu; then
  mkdir -p bwa && for x in .fa .fa.fai .dict; do rm -f bwa/GRCh38$x; ln $M$x bwa/GRCh38$x; done
  [ -s bwa/GRCh38.fa.sa ] || ( cd bwa && "$W/$D" $IMG_BWA bwa index -a bwtsw GRCh38.fa )
else
  [ -s GRCh38.fa.bwt.2bit.64 ] || "$W/tools/bin/bwa-mem2" index GRCh38.fa
fi
echo INDEX DONE
