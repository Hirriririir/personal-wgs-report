#!/usr/bin/env bash
# Y 染色体（父系）单倍群：只取 chrY 读段建小 BAM（Yleaf 的 -p 模式对整个基因组 BAM 会生成几十 GB 的 pileup），
# 再用 Yleaf 同时对照 YFull、FTDNA、ISOGG 三棵树。女性跳过。
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
s=$SAMPLE
[ "$(sample_sex)" = male ] || { log "非男性样本，跳过 Y 染色体"; exit 0; }
O=work/$s/ychr; mkdir -p $O
[ -s $O/$s.chrY.bam.bai ] || $D $IMG_SAMTOOLS sh -c "samtools view -b -F 0x400 work/$s/$s.bam chrY -o $O/$s.chrY.bam && samtools index $O/$s.chrY.bam"
W=$(pwd -P)
( cd $O && "$W/tools/venv-yleaf/bin/Yleaf" -bam "$W/$O/$s.chrY.bam" -rg hg38 --ref-fasta "$W/$REF" -o yleaf -t 8 -r 3 -q 20 -b 90 \
    -tree yfull ftdna isogg -p --no-update-check --report-json "$W/$O/report.json" > yleaf.log 2>&1 )
log "Y DONE $s：$(head -3 $O/yleaf/hg_prediction*.hg 2>/dev/null | tail -1)"
