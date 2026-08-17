#!/usr/bin/env bash
set -euo pipefail
CAPD_SHA=731079217a9254ea2948d742df2b170895effe7f
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
python3 analytic_prechecks.py
if [ ! -d CAPD/.git ]; then
  git clone https://github.com/CAPDGroup/CAPD.git CAPD
fi
git -C CAPD fetch origin master
git -C CAPD checkout "$CAPD_SHA"
cmake -S CAPD -B CAPD/build -DCMAKE_BUILD_TYPE=Release
cmake --build CAPD/build --parallel 2
make CAPD_CONFIG="$ROOT/CAPD/build/bin/capd-config"
set +e
./si3bp_capd 2>&1 | tee capd-proof.log
code=${PIPESTATUS[0]}
set -e
grep -F "C1: PASS" capd-proof.log
grep -F "C4 (Burns--Weiss crossing): PASS" capd-proof.log
grep -F "ALL ENTROPY CERTIFICATES PASSED." capd-proof.log
exit "$code"
