#!/bin/bash
# Project status at a glance: ./status.sh  (run from the app terminal, Ctrl+`)
cd "$(dirname "$0")"

echo "== Downloads (verified by size + md5)"
unset PYTHONPATH
env/venv/bin/python data/fetch_all.py --status 2>/dev/null | grep -v ' done$'
pgrep -f fetch_all.py >/dev/null && echo "  download manager: running" || echo "  download manager: NOT running (rerun: env/venv/bin/python data/fetch_all.py)"

echo "== Running background processes"
pgrep -f caffeinate >/dev/null && echo "  caffeinate: keeping Mac awake" || echo "  caffeinate: NOT running"
for p in bowtie2 macs3 fastp samtools pytest idr; do pgrep -x "$p" >/dev/null && echo "  $p: running"; done
pgrep -f 'cdtest|simulate' >/dev/null && echo "  simulations: running"

echo "== Last commits"
git log --oneline -5 | sed 's/^/  /'
echo "== Disk free: $(df -h ~ | tail -1 | awk '{print $4}')"
