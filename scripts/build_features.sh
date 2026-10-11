#!/usr/bin/env bash
# Rebuild accessibility features (ANALYSIS_PLAN §3) for every line in config/feature_lines.tsv whose
# inputs exist, then sanity checks + heatmap. Re-run when Jurkat / local RPE1 outputs appear.
set -euo pipefail
cd "$(dirname "$0")/.."
unset PYTHONPATH
source env/venv/bin/activate
export PYTHONPATH=src
mkdir -p results/features
python -m cdtest.features --config config/feature_lines.tsv --procs "${PROCS:-6}" 2>&1 | tee results/features/build.log
python scripts/features_sanity.py || echo "WARNING: some sanity checks failed (see results/features/sanity.md)"
