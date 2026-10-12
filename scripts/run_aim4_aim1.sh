#!/bin/bash
# Aim 4, step 1: Aim 1 outcome M for the two matched pairs and the three new mixed pairs
# (plan §0 orientation, M = A - B), frozen pipeline pipeline-v1.1. Half-split terms are stored
# for the knockdown bootstrap (as on the development pair). Then the group-mean variant
# (sensitivity) for the two matched pairs. RPE1_vs_K562 already exists and is reused.
cd "$(dirname "$0")/.."
unset PYTHONPATH; source env/venv/bin/activate; export PYTHONPATH=src
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p data/aim1_cache results/aim4/logs
for pair in "K562 Jurkat" "RPE1 HepG2" "K562 HepG2" "Jurkat RPE1" "Jurkat HepG2"; do
  set -- $pair
  python -c "
from cdtest import aim1
aim1.run_pair('$1','$2', save_half_terms='data/aim1_cache/$1_vs_$2_halfD.npy')
" > results/aim4/logs/aim1_$1_vs_$2.log 2>&1 &
done
wait
for pair in "K562 Jurkat" "RPE1 HepG2"; do
  set -- $pair
  python -c "
from cdtest import aim1
aim1.run_pair('$1','$2', generic='group')
" > results/aim4/logs/aim1_$1_vs_$2_groupmean.log 2>&1 &
done
wait
echo AIM4_AIM1_DONE
