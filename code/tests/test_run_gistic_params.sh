#!/usr/bin/env bash
# The GISTIC parameters are pre-registered in the spec. A silent drift here (e.g. a
# stale -maxseg 100000, or an hg19 refgene) invalidates every peak downstream, and
# neither failure is visible in the output. Assert them literally.
set -euo pipefail
S=code/run_gistic.sh

for p in "-brlen 0.7" "-conf 0.99" "-ta 0.3" "-td 0.3" "-cap 1.5" \
         "-smallmem 0" "-v 30" "-js 4" "-genegistic 1" "-broad 1" \
         "-armpeel 1" "-savegene 1" "-gcm extreme" "-rx 0"; do
  grep -q -- "$p" "$S" || { echo "MISSING PARAM: $p"; exit 1; }
done

grep -qE 'MAXSEG:?=\{?MAXSEG:-46000\}?|MAXSEG="\$\{MAXSEG:-46000\}"' "$S" \
  || { echo "MISSING: maxseg default 46000"; exit 1; }

grep -q "hg38" "$S" || { echo "MISSING: hg38 refgene reference"; exit 1; }

for g in mmrp_high mmrp_low mmrd_low mmrd_high; do
  grep -q "$g" "$S" || { echo "MISSING GROUP: $g"; exit 1; }
done
grep -qi "leave-one-out\|loo" "$S" || { echo "MISSING: LOO mode"; exit 1; }

bash -n "$S" || { echo "SYNTAX ERROR"; exit 1; }

echo "run_gistic.sh params + groups + LOO: OK"
