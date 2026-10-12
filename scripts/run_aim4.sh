#!/bin/bash
# Aim 4, steps 2-8 (after scripts/run_aim4_aim1.sh). Frozen pipeline pipeline-v1.1.
cd "$(dirname "$0")/.."
unset PYTHONPATH; source env/venv/bin/activate; export PYTHONPATH=src
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
L=results/aim4/logs; mkdir -p $L
# main test + comparators + sensitivities (v1.3) for the 5 new pairs; v1.2 for the matched pairs
for pair in "K562 Jurkat" "RPE1 HepG2" "K562 HepG2" "Jurkat RPE1" "Jurkat HepG2"; do
  set -- $pair; python scripts/run_aim4.py tests $1 $2 > $L/tests_$1_vs_$2.log 2>&1 &
done
wait
for pair in "K562 Jurkat" "RPE1 HepG2"; do
  set -- $pair
  python scripts/run_aim4.py tests_v12 $1 $2 > $L/tests_v12_$1_vs_$2.log 2>&1 &
  python scripts/run_aim4.py groupmean $1 $2 > $L/groupmean_$1_vs_$2.log 2>&1 &
done
python scripts/run_aim4.py extras > $L/extras.log 2>&1 &
wait
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for pair in "K562 Jurkat" "RPE1 HepG2"; do
  set -- $pair; python scripts/run_aim4.py planted $1 $2 > $L/planted_$1_vs_$2.log 2>&1
done
echo AIM4_RUN_DONE
