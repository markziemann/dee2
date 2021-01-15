#!/bin/bash

set -x
# This script is designed to accommodat requests from users who nominate
# specific accessions to be processed

# because it is cron scheduled it needs to cd
cd "$(dirname "$0")";

TODO=todo.txt

PROG=inprogress

# check if already working
if [ -r $PROG ] ; then exit ; fi

scp ubuntu@118.138.234.94:/var/www/request.txt .

comm -23 <(sort request.txt) <(sort oldrequest.txt) > $TODO

len=$(wc -l < $TODO)

if [ $len -gt 0 ] ; then

  while IFS= read -r line; do echo $line

    ORG=$(echo $line | cut -d ' ' -f1)
    SRP=$(echo $line | cut -d ' ' -f2)
    EMAIL=$(echo $line | cut -d ' ' -f3)
    ACCS=$(grep -w $SRP ../sradb/${ORG}.csv | cut -d ',' -f1 | paste -s -d ',')
    len2=$(echo $ACCS | wc -c)
    N=$(echo $ACCS | tr ',' '\n' | wc -l)

    if [ $len2 -gt 0 ] ; then

      > $PROG
      docker run mziemann/tallyup $ORG $ACCS

      CONTAINER=$(docker ps -aql)
      mkdir data
      docker cp ${CONTAINER}:/dee2/data/$ORG data
      ACCS2=$(grep -w $SRP ../sradb/${ORG}.csv | cut -d ',' -f1 | paste -s -d '|')
      mkdir $SRP

      # README
      cp contents.md $SRP/README.md

      # Gene and tx info
      ORG3=$(echo $ORG | cut -c-3)
      cp ../gene_info/${ORG3}_gene_info.tsv $SRP/GeneInfo.tsv
      cp ../gene_info/${ORG3}_tx_info.tsv $SRP/TxInfo.tsv

      # Prepare genecounts file
      SE=$(find data/$ORG/ | egrep $ACCS2 | grep se.tsv)
      NUMSE=$(echo $SE | wc -w)
      NUMCOL=$(echo $((NUMSE * 2)) )
      COLS=$(seq 2 2 $NUMCOL)
      COLS=$(echo $COLS | tr ' ' ',')
      echo GeneID $ACCS2 | tr '| ' '\t' > $SRP/GeneCountMatrix.tsv
      paste $SE | sed 1d | cut -f1,$COLS >> $SRP/GeneCountMatrix.tsv
      ERROR_CNT=0
      RES_COL=$(head $SRP/GeneCountMatrix.tsv | tail -1 | wc -w )
      RES_COL=$((RES_COL-1))
      if [ $RES_COL -ne $N ] ; then
        ERROR_CNT=$((ERROR_CNT+1))
      fi

      # Prepare txcounts file
      KE=$(find data/$ORG/ | egrep $ACCS2 | grep ke)
      NUMKE=$(echo $KE | wc -w)
      NUMCOL=$(echo $((NUMKE * 5)) )
      COLS=$(seq 4 5 $NUMCOL)
      COLS=$(echo $COLS | tr ' ' ',')
      paste $KE | cut -f$COLS | head -1 | sed 's/^/TranscriptID\t/' \
      | sed 's/_est_counts//g' > $SRP/TxCountMatrix.tsv
      paste $KE | sed 1d | cut -f1,$COLS >> $SRP/TxCountMatrix.tsv
      RES_COL=$(head $SRP/TxCountMatrix.tsv | tail -1 | wc -w )
      RES_COL=$((RES_COL-1))
      if [ $RES_COL -ne $N ] ; then
        ERROR_CNT=$((ERROR_CNT+1))
      fi

      # prepare qc file
      QC=$(find data/$ORG/ | egrep $ACCS2 | grep qc)
      NUMQC=$(echo $QC | wc -w)
      NUMQC=$(echo $((NUMQC * 2)) )
      COLS=$(seq 2 2 $NUMQC)
      COLS=$(echo 1 $COLS | tr ' ' ',')
      echo QCmetric $ACCS2 | tr '| ' '\t' > $SRP/QC_Matrix.tsv
      paste $QC | tr ':' '\t' | head -28 | cut -f$COLS >> $SRP/QC_Matrix.tsv
      RES_COL=$(head $SRP/QC_Matrix.tsv | tail -1 | wc -w )
      RES_COL=$((RES_COL-1))
      if [ $RES_COL -ne $N ] ; then
        ERROR_CNT=$((ERROR_CNT+1))
      fi

      LOG=$(find data/$ORG/ | egrep $ACCS2 | grep log )
      mkdir -p $SRP/log
      cp $LOG $SRP/log
      RES_COL=$(echo $LOG | wc -w )
      if [ $RES_COL -ne $N ] ; then
        ERROR_CNT=$((ERROR_CNT+1))
      fi

      cp contents.md $SRP/README.md

      if [ $ERROR_CNT -gt 0 ] ; then
        echo "There were errors processing ${SRP}" >> $SRP.txt
        scp $SRP.txt ubuntu@118.138.234.94:/dee2_data/requests
        cp -r data $SRP.data
        scp $SRP.data ubuntu@118.138.234.94:/dee2_data/requests
        rm $SRP.data
        zip -r $SRP.zip $SRP
        scp $SRP.zip ubuntu@118.138.234.94:/dee2_data/requests
      else
        zip -r $SRP.zip $SRP
        scp $SRP.zip ubuntu@118.138.234.94:/dee2_data/requests
      fi

      rm -rf data

      rm $PROG

    fi

  done < $TODO

fi

cp request.txt oldrequest.txt





