#!/bin/bash
# Project status at a glance: ./status.sh  (run from the app terminal, Ctrl+`)
cd "$(dirname "$0")"
gb() { awk -v b="$(stat -f%z "$1" 2>/dev/null || echo 0)" 'BEGIN{printf "%6.2f", b/1e9}'; }

echo "== Perturb-seq screens (GB downloaded / total)"
for f in K562_essential_raw_singlecell_01:10.66 rpe1_raw_singlecell_01:8.70 \
         GSE264667_jurkat_raw_singlecell_01:9.37 GSE264667_hepg2_raw_singlecell_01:5.61; do
  n=${f%%:*}; t=${f##*:}; printf "  %-40s %s / %s\n" "$n" "$(gb data/perturbseq/$n.h5ad)" "$t"
done

echo "== ATAC signal + peaks (11 lines)"
printf "  files: %s / 20   size: %s\n" "$(ls data/atac 2>/dev/null | grep -cE '\.(bigWig|bw|bed\.gz|bb)$')" "$(du -sh data/atac 2>/dev/null | cut -f1)"

echo "== Jurkat raw reads (GSE273644, ~7.35 GB)"
printf "  size: %s\n" "$(du -sh data/atac_fastq/jurkat_GSE273644 2>/dev/null | cut -f1)"

echo "== Running background processes"
pgrep -fl 'curl .*(figshare|ncbi|encodeproject|4dn|ebi)' >/dev/null && echo "  downloads: $(pgrep -f 'curl .*(figshare|ncbi|encodeproject|4dn|ebi)' | wc -l | tr -d ' ') curl transfers"
pgrep -f caffeinate >/dev/null && echo "  caffeinate: keeping Mac awake" || echo "  caffeinate: NOT running"
for p in bowtie2 macs3 fastp samtools pytest; do pgrep -x "$p" >/dev/null && echo "  $p: running"; done
pgrep -f 'cdtest|simulate' >/dev/null && echo "  simulations: running"

echo "== Last commits"
git log --oneline -5 | sed 's/^/  /'
echo "== Disk free: $(df -h ~ | tail -1 | awk '{print $4}')"
