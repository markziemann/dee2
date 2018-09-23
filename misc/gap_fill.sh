#!/bin/bash
>gap_fill.txt
for ORG in athaliana celegans dmelanogaster drerio ecoli hsapiens mmusculus rnorvegicus scerevisiae ; do
  MDATA=../sradb/${ORG}_metadata.tsv.cut
  DATA=../data/${ORG}
  for SRP in $(sed 1d $MDATA | cut -f5 | sort -u  ) ; do
    TOT=0 ; DONE=0 ; NOT_DONE=""
    for SRR in $(grep -w $SRP $MDATA | cut -f1 ) ;  do
      TOT=$((TOT+1))
      if [ -d $DATA/$SRR ] ; then
        DONE=$((DONE+1))
      else
        NOT_DONE=${NOT_DONE},$SRR
      fi
    done
    echo $SRP $TOT $DONE $NOT_DONE
  done | awk '$3/$2<1 && $3/$2>.5 {print $1,$2,$3,$3/$2,$4}' | sort -k5g \
  |cut -d ' ' -f5 | tr -d  '\n' |sed 's/$/\n/' | sed "s/^/${ORG} /" | sed 's/,//'
done | tee gap_fill.txt

