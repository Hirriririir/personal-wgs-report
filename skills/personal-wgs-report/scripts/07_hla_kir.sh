#!/usr/bin/env bash
# HLA 与 KIR 分型（T1K）：HLA 用 IPD-IMGT/HLA（含用药风险等位基因 HLA-B*15:02 / B*58:01 / A*31:01 等），KIR 用 IPD-KIR
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
s=$SAMPLE
O=work/$s/hla; mkdir -p $O
$D $IMG_T1K run-t1k -b work/$s/$s.bam -c ref/hla/hlaidx/hlaidx_dna_coord.fa -f ref/hla/hlaidx/hlaidx_dna_seq.fa \
  --preset hla-wgs --abnormalUnmapFlag -t 8 -o $s --od $O > $O/t1k.log 2>&1
O=work/$s/kir; mkdir -p $O
$D $IMG_T1K run-t1k -b work/$s/$s.bam -c ref/kir/kiridx/kiridx_dna_coord.fa -f ref/kir/kiridx/kiridx_dna_seq.fa \
  --preset kir-wgs --abnormalUnmapFlag -t 8 -o $s --od $O > $O/t1k.log 2>&1
log "HLA/KIR DONE $s"
