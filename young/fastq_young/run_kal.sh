#!/bin/bash
for FQZ1 in *_1.fq.gz ; do
  FQZ2=$( echo $FQZ1 | sed 's/_1.fq.gz/_2.fq.gz/' )
  ../sw/skewer/skewer -q 10 -t 8 $FQZ1 $FQZ2

  FQ1=$(echo $FQZ1 | sed 's/_1.fq.gz/_1.fq-trimmed-pair1.fastq/')
  FQ2=$(echo $FQ1 | sed 's/_1.fq-trimmed-pair1.fastq/_1.fq-trimmed-pair2.fastq/')
  ../sw/kallisto/kallisto quant \
  -i ../ref/gencode.vM24.transcripts.fa.idx \
  -o ${FQ1}_kal -t 16 $FQ1 $FQ2
done


for TSV in */*abundance.tsv ; do
  NAME=$(echo $TSV | cut -d '_' -f1) ; cut -f1,4 $TSV | sed 1d | sed "s/^/${NAME}\t/"
done > 3col.tsv

