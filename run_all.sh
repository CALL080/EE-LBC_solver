#!/bin/bash
# Regenerate outputs/ (run from this directory).  Needs Python 3.8+ (standard library only) and gcc with OpenMP.
#   ./run_all.sh          checks, exact searches, hull sums and figure data (a few minutes)
#   ./run_all.sh --full   also the reduced-round experiments of Table 3 (about one hour on 20 cores)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
mkdir -p bin outputs/figs
gcc -O3 -march=native -o bin/matsui matsui.c -lm
gcc -O3 -o bin/hull hull.c -lm
gcc -O3 -march=native -fopenmp -o bin/experiments experiments.c -lm
gcc -O3 -march=native -fopenmp -o bin/fixedhull fixedhull.c -lm

$PY eelbc.py > outputs/eelbc_selftest.log
$PY check_spec.py > outputs/check_spec.log

# exact minimum weights (Theorem 1, Table 2, Fig. 2); the last argument is a time limit per round in seconds
bin/matsui eelbc lin 34 600 > outputs/matsui_eelbc_lin_34.log
bin/matsui eelbc diff 34 600 > outputs/matsui_eelbc_diff_34.log
bin/matsui eelbc aslin 31 600 > outputs/matsui_eelbc_aslin.log
bin/matsui eelbc asdiff 31 600 > outputs/matsui_eelbc_asdiff.log
bin/matsui present lin 19 600 > outputs/matsui_present_lin.log
bin/matsui present diff 31 600 > outputs/matsui_present_diff.log
bin/matsui present aslin 31 600 > outputs/matsui_present_aslin.log
bin/matsui present asdiff 31 600 > outputs/matsui_present_asdiff.log

# linear hull and differential restricted to at most two active S-boxes per round (Sect. 3.1)
bin/hull eelbc lin 21 21 31 > outputs/hull_eelbc_lin_21.log
bin/hull eelbc diff 0 0 31 > outputs/hull_eelbc_diff_0.log
bin/hull present lin 21 21 31 > outputs/hull_present_lin_21.log

# Sect. 4.3: fixed-key hull for the keys of corr_16keys.log, and its distribution over 2000 keys
{ for L in 8 10 12; do bin/fixedhull check $L 201 216; done; } > outputs/fixedhull_check.log
{ for L in 8 10 12 13 20 28 30; do bin/fixedhull stats $L 2000; done; } > outputs/fixedhull_stats.log

if [ "${1:-}" = "--full" ]; then
  # Table 3: linear r = 4..13 (four keys with the key schedule, two with independent round keys), r = 14 (two keys);
  # differential r = 3..8 (four keys).  r rounds cover L = r - 1 trail rounds; the last round has no key.
  {
    for r in 4 5 6 7 8 9 10 11 12 13; do lg=$((2*(r-1)+8))
      for s in 1 2 3 4; do bin/experiments lin $r $lg $s; done
      for s in 5 6; do bin/experiments lin $r $lg $s indep; done
    done
    for s in 1 2; do bin/experiments lin 14 34 $s; done
    for r in 3 4 5 6 7 8; do lg=$((3*(r-1)+8)); for s in 1 2 3 4; do bin/experiments diff $r $lg $s; done; done
  } > outputs/experiments.log
  # Sect. 3.2: spread of the correlation over 16 keys (L = 8, 10, 12)
  { for r in 9 11 13; do lg=$((2*(r-1)+8)); for s in $(seq 201 216); do bin/experiments lin $r $lg $s; done; done; } > outputs/corr_16keys.log
  # Table 5: key recovery with a 20-bit guess on 12 rounds (N = 2^5, 2^6, 2^7 / c^2) and on 14 and 16 rounds (2^6 / c^2)
  { for lg in 23 24 25; do for s in $(seq 101 132); do bin/experiments keyrec 12 $lg $s; done; done; } > outputs/keyrec_r12_32seeds.log
  { for a in '14 28' '16 32'; do set -- $a; for s in $(seq 301 316); do bin/experiments keyrec $1 $2 $s; done; done; } > outputs/keyrec_longer.log
fi

$PY complexity.py > outputs/complexity.log

$PY summarize_experiments.py > outputs/summarize_experiments.log
$PY make_figs.py > /dev/null
echo "done; see outputs/"
