#!/bin/bash
# Parallel byte-range download of the 65.8 GB K562 genome-wide h5ad from Figshare.
# The single stream runs at ~1 MB/s; Figshare's S3 backend serves ranges, so N streams
# in parallel. Keeps the bytes already downloaded as the first piece. md5 checked at end.
set -u
cd "$(dirname "$0")/perturbseq"
URL=https://ndownloader.figshare.com/files/35775507
F=K562_gwps_raw_singlecell_01.h5ad
SIZE=65830941948
MD5=887e3e6a8c8df6eadf7a3030a53c9546
N=12
pkill -f "curl.*35775507" ; sleep 2
HAVE=$(stat -f%z $F 2>/dev/null || echo 0)
REM=$((SIZE-HAVE)); CH=$(( (REM + N - 1) / N ))
mkdir -p gwps_parts
for i in $(seq 0 $((N-1))); do
  s=$((HAVE + i*CH)); e=$((s + CH - 1)); [ $e -ge $SIZE ] && e=$((SIZE-1))
  ( p=gwps_parts/part$(printf %02d $i); want=$((e-s+1))
    while :; do
      got=$(stat -f%z $p 2>/dev/null || echo 0)
      [ $got -ge $want ] && break
      curl -sS -L --speed-limit 20000 --speed-time 60 -r $((s+got))-$e "$URL" >> $p || sleep 5
    done; echo "$(date +%T) part $i done" ) &
done
wait
cat gwps_parts/part* >> $F && rm -rf gwps_parts
got=$(md5 -q $F)
if [ "$got" = "$MD5" ]; then echo "perturbseq/$F" >> ../manifest_done.txt; echo GWPS_OK; else echo "GWPS_MD5_FAIL $got"; fi
