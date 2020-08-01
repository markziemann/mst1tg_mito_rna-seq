#!/bin/bash
for FQZ in *.fq.gz ; do
  ../sw/skewer/skewer -q 10 -t 8 $FQZ

  FQ1=$(echo $FQZ | sed 's/.fq.gz/.fq-trimmed.fastq/')
  ../sw/kallisto/kallisto quant \
  --single \
  --fragment-length 150 \
  --sd 30 \
  -i ../ref/gencode.vM24.transcripts.fa.idx \
  -o ${FQ1}_kal -t 16 $FQ1
done

for TSV in */*abundance.tsv ; do
  NAME=$(echo $TSV | cut -d '_' -f1) ; cut -f1,4 $TSV | sed 1d | sed "s/^/${NAME}\t/"
done > 3col.tsv

