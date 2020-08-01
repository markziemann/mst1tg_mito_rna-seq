#!/bin/bash

conv(){
line=$1
  NAME_DESC=`echo $line | cut -d ' ' -f-2`

  GENES=`echo $line | cut -d ' ' -f3- \
  | tr ' ' '\n' | sort -uk 1b,1 \
  | join -1 1 -2 1 - \
  <(cut -f3,4 mart_export.txt \
  | sed 1d | awk '$1!="" && $2!=""' \
  | sort -uk 1b,1) | cut -d ' ' -f2 \
  | sort -u | tr '\n' '\t' \
  | sed 's/\t$/\n/'`

  echo $NAME_DESC $GENES | tr ' ' '\t'
}
export -f conv

for GMT in `ls *gmt | grep -v mouse ` ; do
  NAME=`echo $GMT | sed 's/.gmt/_mouse.gmt/'`
  sed -i 's@ @-@g' $GMT
  parallel -k conv < $GMT > $NAME
done

