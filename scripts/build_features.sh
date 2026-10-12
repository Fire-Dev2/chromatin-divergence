#!/usr/bin/env bash
# Rebuild accessibility features (ANALYSIS_PLAN §3; plan v1.2: insertion CPM primary, FC sensitivity)
# for every line in config/feature_lines.tsv whose inputs exist, then sanity checks + heatmaps.
set -euo pipefail
cd "$(dirname "$0")/.."
unset PYTHONPATH
source env/venv/bin/activate
export PYTHONPATH=src
mkdir -p results/features
python -m cdtest.features --signal both --config config/feature_lines.tsv --procs "${PROCS:-10}" 2>&1 | tee results/features/build.log
python scripts/features_sanity.py || echo "WARNING: some sanity checks failed (see results/features/sanity.md)"
