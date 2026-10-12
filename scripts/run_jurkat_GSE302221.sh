#!/bin/bash
# Jurkat CRISPRi-background sensitivity ATAC set (GSE302221 / PRJNA1289670; dCas9-KRAB,
# non-targeting guide, 3 reps, 2x151 NovaSeq). Download from ENA (resume + md5), take the
# first 40M read pairs per replicate (head of file; ~primary Jurkat depth, 31-40M pairs),
# full-length reads (no 50-bp trim), local pipeline with --threads 8 (shared machine), pool,
# insertion arrays data/atac_ins/JurkatCRISPRi__repK.npz. FASTQs deleted at the end.
set -euo pipefail
cd /Users/tharunnambi/chromatin-divergence
unset PYTHONPATH; source env/venv/bin/activate
SAMPLE=jurkat_GSE302221
F=data/atac_fastq/$SAMPLE; S=$F/sub40M; mkdir -p $S
RUNS=(SRR34472975 SRR34473029 SRR34473025)   # rep1 rep2 rep3
N=$((40000000*4))

# ---- download (resume + md5) -------------------------------------------------------------
dl() {  # url md5 dest bytes
  local url=$1 md5=$2 out=$3 bytes=$4
  [ -f $out.ok ] && return 0
  for t in 1 2 3 4 5 6 7 8; do
    curl -sS -L --retry 5 --retry-delay 10 -C - -o $out "$url" || true
    if [ "$(md5 -q $out 2>/dev/null)" == "$md5" ]; then touch $out.ok; echo "OK $out"; return 0; fi
    echo "md5 mismatch/incomplete $out (try $t)"
    # complete-size file with wrong md5 = corrupt: restart from scratch
    [ -f $out ] && [ $(stat -f %z $out) -ge $bytes ] && { echo "corrupt full-size $out, deleting"; rm -f $out; }
    sleep 20
  done
  echo "FAILED $out"; return 1
}
: > $F/ena_filereport.tsv
for r in "${RUNS[@]}"; do
  meta=$(curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$r&result=read_run&fields=fastq_ftp,fastq_bytes,fastq_md5,read_count" | tail -1)
  echo "$meta" >> $F/ena_filereport.tsv
  urls=$(echo "$meta" | cut -f2); md5s=$(echo "$meta" | cut -f4); sizes=$(echo "$meta" | cut -f3)
  for k in 1 2; do
    u=$(echo $urls | cut -d';' -f$k); m=$(echo $md5s | cut -d';' -f$k); b=$(echo $sizes | cut -d';' -f$k)
    dl "https://$u" $m $F/${r}_$k.fastq.gz $b &
  done
done
# replicates are processed in the order their downloads complete (one at a time, --threads 8)
ready() { [ -f $F/$1_1.fastq.gz.ok ] && [ -f $F/$1_2.fastq.gz.ok ]; }
DONE=""
while [ $(echo $DONE | wc -w) -lt 3 ]; do
  i=0; picked=""
  for r in "${RUNS[@]}"; do
    i=$((i+1)); case " $DONE " in *" rep$i "*) continue;; esac
    if ready $r; then picked="$r rep$i"; break; fi
  done
  if [ -z "$picked" ]; then
    if ! pgrep -f "curl -sS -L --retry 5" > /dev/null; then echo "downloads stopped with replicates missing"; exit 1; fi
    sleep 60; continue
  fi
  set -- $picked; r=$1; rep=$2
# ---- per replicate: first 40M pairs, pipeline -----------------------------------------------
  echo "processing $rep ($r)"
  if [ ! -f data/atac_local/$SAMPLE/$rep/.align.done ]; then
    pids=""
    for k in 1 2; do
      if [ ! -f $S/${rep}_R$k.done ]; then
        ( gzip -dc $F/${r}_$k.fastq.gz | head -n $N > $S/${rep}_R$k.fastq || true; touch $S/${rep}_R$k.done ) &
        pids="$pids $!"
      fi
    done
    for p in $pids; do wait $p; done
  fi
  if [ -f $S/${rep}_R1.fastq ]; then
    l1=$(wc -l < $S/${rep}_R1.fastq); l2=$(wc -l < $S/${rep}_R2.fastq)
    [ $l1 -eq $l2 ] && [ $l1 -eq $N ] || { echo "subsample length problem $rep $l1 $l2"; exit 1; }
  fi
  python scripts/atac_pipeline.py rep --sample $SAMPLE --rep $rep --r1 $S/${rep}_R1.fastq --r2 $S/${rep}_R2.fastq --threads 8
  rm -f $S/${rep}_R1.fastq $S/${rep}_R2.fastq
  DONE="$DONE $rep"
done
wait
python scripts/atac_pipeline.py pool --sample $SAMPLE --reps rep1 rep2 rep3 --threads 8
for k in 1 2 3; do
  python scripts/atac_bam_to_insertions.py --tagalign data/atac_local/$SAMPLE/rep$k/rep$k.tn5.tagAlign.gz JurkatCRISPRi__rep$k
done
echo JURKAT_CRISPRI_DONE
