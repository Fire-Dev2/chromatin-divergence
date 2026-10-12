#!/bin/bash
# Aim 1 on the development pair: rerun primary (stores half-split terms for the knockdown
# bootstrap; outputs must reproduce the earlier run) and the group-mean variant.
cd "$(dirname "$0")/.."
unset PYTHONPATH; source env/venv/bin/activate; export PYTHONPATH=src
mkdir -p data/aim1_cache
python -c "
from cdtest import aim1
aim1.run_pair('RPE1','K562', save_half_terms='data/aim1_cache/RPE1_vs_K562_halfD.npy')
" > results/aim1/RPE1_vs_K562_rerun.log 2>&1 &
python -c "
from cdtest import aim1
aim1.run_pair('RPE1','K562', generic='group')
aim1.p53_gate('RPE1_vs_K562_groupmean')
" > results/aim1/RPE1_vs_K562_groupmean.log 2>&1 &
wait
echo AIM1_DEV_DONE
