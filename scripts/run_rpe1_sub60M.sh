#!/bin/bash
# RPE1 (4DN 4DNES48BA2C8): first 60M read pairs per replicate, same as the K562 validation
# subsampling (head of file). Then the standard local pipeline, then insertion arrays.
set -euo pipefail
cd "$(dirname "$0")/.."
unset PYTHONPATH; source env/venv/bin/activate
F=data/atac_fastq/rpe1_4DNES48BA2C8
S=$F/sub60M; mkdir -p $S
N=$((60000000*4))
sub() { [ -f $S/$2.done ] || { gzip -dc $F/$1.fastq.gz | head -n $N > $S/$2.fastq; touch $S/$2.done; }; }
sub 4DNFIDNZZ2A9 rep1_R1 & sub 4DNFI14DT5NE rep1_R2 & sub 4DNFIWGX7K9Y rep2_R1 & sub 4DNFIY8ZSQC1 rep2_R2 & wait
for r in rep1 rep2; do
  [ $(wc -l < $S/${r}_R1.fastq) -eq $(wc -l < $S/${r}_R2.fastq) ] || { echo "pair length mismatch $r"; exit 1; }
  python scripts/atac_pipeline.py rep --sample rpe1_4DNES48BA2C8 --rep $r --r1 $S/${r}_R1.fastq --r2 $S/${r}_R2.fastq --threads 12
  rm -f $S/${r}_R1.fastq $S/${r}_R2.fastq
done
python scripts/atac_pipeline.py pool --sample rpe1_4DNES48BA2C8 --reps rep1 rep2
for r in rep1 rep2; do
  python scripts/atac_bam_to_insertions.py --tagalign data/atac_local/rpe1_4DNES48BA2C8/$r/$r.tn5.tagAlign.gz RPE1__$r
done
echo RPE1_DONE
