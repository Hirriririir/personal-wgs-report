#!/usr/bin/env bash
# 参考基因组：NCBI GRCh38 no-alt + hs38d1 decoy analysis set（UCSC 风格 chr 命名；chrY PAR 已屏蔽），
# 再按 GRC_exclusions.bed 硬屏蔽假重复区（同 GIAB v3 的 maskedGRC 做法），坐标不变。
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
W=$PWD
mkdir -p ref/GRCh38 && cd ref/GRCh38
U=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
N=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
M=GRCh38_no_alt_plus_hs38d1_GRCmasked.fa
[ -s $M ] && [ -L GRCh38.fa ] && { echo "参考已就绪"; exit 0; }
curl -fSL -C - --retry 5 -o $N.gz $U/$N.gz
curl -fSL --retry 5 -o $N.fai $U/$N.fai
curl -fSL --retry 5 -o md5checksums.txt $U/md5checksums.txt
grep -E " \./$N\.(gz|fai)$" md5checksums.txt | sed "s| \./| |" | md5sum -c -
gzip -dc $N.gz > $N
curl -fSL --retry 5 -o GRCh38_GRC_exclusions.bed $U/GCA_000001405.15_GRCh38_GRC_exclusions.bed
python3 "$W/scripts/setup/mask_reference.py" $N GRCh38_GRC_exclusions.bed $M
cp $N.fai $M.fai && rm $N
ln -sfn $M GRCh38.fa; ln -sfn $M.fai GRCh38.fa.fai
echo REF DONE
