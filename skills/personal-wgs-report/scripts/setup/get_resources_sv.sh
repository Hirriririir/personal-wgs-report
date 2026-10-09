#!/usr/bin/env bash
# 结构变异（15_sv.sh，可选）用的资源，约 20 GB：AnnotSV 3.5 人类注释包（gnomAD-SV / DGV / ClinGen / OMIM / ACMG CNV 打分规则等）、
# Delly 的 GRCh38 排除区（着丝粒、端粒等容易出假阳性的区域）。
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
F="python3 scripts/setup/fetch.py url"
mkdir -p ref/sv/annotsv
$F https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv ref/sv/delly_human.hg38.excl.tsv
if [ ! -d ref/sv/annotsv/Annotations_Human ]; then
  $F https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz ref/sv/annotsv/Annotations_Human_3.5.tar.gz
  tar -xzf ref/sv/annotsv/Annotations_Human_3.5.tar.gz -C ref/sv/annotsv && rm -f ref/sv/annotsv/Annotations_Human_3.5.tar.gz*
fi
echo SV RESOURCES DONE
