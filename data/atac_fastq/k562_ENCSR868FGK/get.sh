#!/usr/bin/env bash
# K562 ENCSR868FGK validation FASTQs for the local ATAC pipeline.
# All 3 bio reps total ~45 GB (rep1 4 runs x 43M pairs; rep2/rep3 4 runs x 29-30M pairs),
# so each replicate is subsampled to its FIRST 40,000,000 read pairs, taken from
# runs in the order listed in files.tsv (rep1: run ENCFF962SRP/ENCFF420XJE only;
# rep2: ENCFF098UCE/ENCFF703BGR then ENCFF903FIK/ENCFF924YGE; rep3: ENCFF354EXH/ENCFF575ZTZ
# then ENCFF827WXQ/ENCFF742FCG). Whole runs are downloaded, md5-checked against ENCODE,
# subsampled with head, and the full runs deleted.
set -uo pipefail
cd "$(dirname "$0")"
N_LINES=160000000   # 40M reads x 4 lines

dl() {  # accession md5
  local a=$1 m=$2 f=$1.fastq.gz
  if [[ -f $f.ok ]]; then return 0; fi
  for i in 1 2 3 4 5 6 7 8; do
    curl -sSL -C - --retry 5 --retry-delay 10 -o "$f" "https://www.encodeproject.org/files/$a/@@download/$a.fastq.gz"
    if [[ "$(md5 -q "$f")" == "$m" ]]; then touch "$f.ok"; echo "$(date +%T) $a ok"; return 0; fi
    echo "$(date +%T) $a md5 mismatch/incomplete (try $i)"; sleep 20
  done
  return 1
}

rep_runs() { awk -v r=$1 '$2==r && $4==1 {print $1"\t"$5}' files.tsv; }

for rep in 1 2 3; do
  out1=rep${rep}_R1.fastq.gz; out2=rep${rep}_R2.fastq.gz
  if [[ -f rep${rep}.done ]]; then continue; fi
  r1s=(); r2s=()
  while IFS=$'\t' read -r a b; do
    dl "$a" "$(awk -v a=$a '$1==a{print $8}' files.tsv)" || exit 1
    dl "$b" "$(awk -v a=$b '$1==a{print $8}' files.tsv)" || exit 1
    r1s+=("$a.fastq.gz"); r2s+=("$b.fastq.gz")
  done < <(rep_runs $rep)
  ( cat "${r1s[@]}" | gzip -dc | head -n $N_LINES | bgzip -@3 -l 3 > $out1 ) &
  ( cat "${r2s[@]}" | gzip -dc | head -n $N_LINES | bgzip -@3 -l 3 > $out2 ) &
  wait
  n1=$(gzip -dc $out1 | wc -l); n2=$(gzip -dc $out2 | wc -l)
  echo "$(date +%T) rep$rep subsample lines R1=$n1 R2=$n2"
  if [[ $n1 -eq $N_LINES && $n2 -eq $N_LINES ]]; then
    rm -f "${r1s[@]}" "${r2s[@]}"; touch rep${rep}.done
  else
    echo "rep$rep subsample wrong length"; exit 1
  fi
done
echo "$(date +%T) all done"
