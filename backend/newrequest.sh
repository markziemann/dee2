#!/bin/bash

# This script is designed to accommodat requests from users who nominate
# specific accessions to be processed

# because it is cron scheduled it needs to cd
cd "$(dirname "$0")";

TODO=todo.txt

PROG=inprogress

# check if already working
if [ -r $PROG ] ; then exit ; fi

scp ubuntu@118.138.234.94:/var/www/html/request.txt .

comm -23 <(sort request.txt) <(sort oldrequest.txt) > $TODO

len=$(wc -l < $TODO)

if [ $len -gt 0 ] ; then

  while IFS= read -r line; do echo $line

    ORG=$(cat $TODO | cut -d ' ' -f1)
    SRP=$(cat $TODO | cut -d ' ' -f2)
    EMAIL=$(cat $TODO | cut -d ' ' -f3)
    ACCS=$(grep -w $SRP ../sradb/${ORG}.csv | cut -d ',' -f1 | paste -s -d ',')
    len2=$(echo $ACCS | wc -c)

    if [ $len2 -gt 0 ] ; then

      > $PROG
      docker run mziemann/tallyup $ORG $ACCS

      CONTAINER=$(docker ps -aq)
      docker cp ${CONTAINER}:/dee2/data/$ORG .
      ACCS2=$(grep -w $SRP ../sradb/${ORG}.csv | cut -d ',' -f1 | paste -s -d '|')
      mkdir $SRP

      # Prepare genecounts file
      SE=$(find data/$ORG/ | egrep $ACCS2 | grep se.tsv)
      NUMSE=$(echo $SE | wc -w)
      NUMCOL=$(echo $((NUMSE * 2)) )
      COLS=$(seq 2 2 $NUMCOL)
      COLS=$(echo $COLS | tr ' ' ',')
      echo GeneID $ACCS2 | tr '| ' '\t' > $SRP/genecounts.tsv
      paste $SE | sed 1d | cut -f1,$COLS >> $SRP/genecounts.tsv

      # Prepare txcounts file
      KE=$(find data/$ORG/ | egrep $ACCS2 | grep ke)
      NUMKE=$(echo $KE | wc -w)
      NUMCOL=$(echo $((NUMKE * 5)) )
      COLS=$(seq 4 5 $NUMCOL)
      COLS=$(echo $COLS | tr ' ' ',')
      paste $KE | cut -f$COLS | head -1 | sed 's/^/TranscriptID\t/' \
      | sed 's/_est_counts//g'  > $SRP/txcounts.tsv
      paste $KE | sed 1d | cut -f1,$COLS >> $SRP/txcounts.tsv

      # prepare qc file
      QC=$(find data/$ORG/ | egrep $ACCS2 | grep qc)
      NUMQC=$(echo $QC | wc -w)
      NUMQC=$(echo $((NUMQC * 2)) )
      COLS=1$(seq 2 2 $NUMQC)
      COLS=$(echo 1 $COLS | tr ' ' ',')
      echo QCmetric $ACCS2 | tr '| ' '\t' > $SRP/qc.tsv
      paste $QC | tr ':' '\t' | head -28 | cut -f$COLS >> $SRP/qc.tsv

      LOG=$(find data/$ORG/ | egrep $ACCS2 | grep log )
      mkdir -p $SRP/log
      cp $LOG $SRP/log

      zip -r $SRP.zip $SRP
      scp $SRP.zip ubuntu@118.138.234.94:/dee2_data/requests

      rm -rf data

      rm $PROG

    fi

  done < $TODO

fi

cp request.txt oldrequest.txt





