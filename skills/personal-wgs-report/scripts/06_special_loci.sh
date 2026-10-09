#!/usr/bin/env bash
# 普通流程测不准的同源基因：SMN1/SMN2（SMA 携带）、CYP2D6（Cyrius）、GBA1（Gauchian）
# α-地贫缺失（HBA1/HBA2）在 scripts/interpret/hba_depth.py 里按深度判断
# 三个工具互相独立，一个失败不影响其他
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
set +e
s=$SAMPLE
O=work/$s/special; mkdir -p $O; echo "$(pwd -P)/work/$s/$s.bam" > $O/manifest.txt
$D $IMG_SMN smn_caller.py --manifest $O/manifest.txt --genome 38 \
  --prefix $s.smn --outDir $O --threads 8 > $O/smn.log 2>&1 || echo "SMNCopyNumberCaller failed (see $O/smn.log)"
tools/venv-illumina/bin/python tools/Cyrius/star_caller.py --manifest $O/manifest.txt --genome 38 --prefix $s.cyp2d6 \
  --outDir $O --threads 8 > $O/cyrius.log 2>&1 || echo "Cyrius failed (see $O/cyrius.log)"
PYTHONPATH=tools/Gauchian tools/venv-illumina/bin/python -c "import sys; from gauchian.gauchian import run; sys.argv[0]='gauchian'; run()" \
  -m $O/manifest.txt -g 38 -o $O -p $s.gba -t 8 > $O/gauchian.log 2>&1 || echo "Gauchian failed (see $O/gauchian.log)"
$PY scripts/interpret/hba_depth.py work/$s/$s.bam $O/$s > $O/hba.log 2>&1 || echo "HBA depth failed (see $O/hba.log)"
echo "SPECIAL DONE $s"
