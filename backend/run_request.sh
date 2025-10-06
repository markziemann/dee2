#!/bin/bash
cd /home/mdz/dee2/request
set -x

if [ -r LOCK ] ; then
  echo LOCK file present. Exiting now.
  exit
fi

# Start by `ls` then create the SRP folder, do analysis in there
ssh -i ~/.ssh/dee2_2025 ubuntu@dee2.io "ls /usr/lib/cgi-bin/newrequests/*confirmed" \
| sed 's#/usr/lib/cgi-bin/newrequests/##' > CONFIRMED

for FILE  in  $(cat CONFIRMED) ; do
  echo $FILE
  if [ ! -r $FILE ] ; then
    echo file $FILE not found
    FILE_PRESENT=FALSE

    #START ROUTINE
    touch LOCK
    SRP=$(echo $FILE | cut -d '.' -f1)
    scp -i ~/.ssh/dee2_2025 ubuntu@dee2.io:/usr/lib/cgi-bin/newrequests/$FILE .
    mkdir $SRP
    cp $FILE $SRP
    cd $SRP
    RUNS=$(grep -w $SRP $FILE | cut -f22)
    NRUNS=$(echo $RUNS | wc -w)
    ORG=$(grep -w $SRP $FILE | head -1 | cut -f7)
    L1=$(echo $ORG | cut -c1 | tr '[:upper:]' '[:lower:]') ; echo $L1
    W2=$(echo $ORG | cut -d ' ' -f2)
    ORG2=$L1$W2

    # DOWNLOAD
    for SRR in $RUNS ; do

      SRACNT=$(find . -maxdepth 1 | grep -c sra$ )

      while [ $SRACNT -gt 3 ] ; do
        sleep 10
        SRACNT=$(find . -maxdepth 1 | grep -c sra$ )
      done

      while [ -r LOCK1 ] ; do
        sleep 1
      done

      touch LOCK1
      prefetch -X 9999999999999 -o ${ORG}_${SRR}.sra $SRR
      rm LOCK1

      {
        while [ -r LOCK2 ] ; do
          sleep 1
        done

        touch LOCK2
        apptainer run -w -B ${PWD}:/dee2/mnt/ tallyup -s $ORG -t 16 -d
        rm LOCK2
      } &

    done
    wait

    for ZIP in *zip ; do
      unzip $ZIP
    done

    mkdir $SRP
    ORG3=$(echo $ORG2 | cut -c-3)
    # get gene and tx info
    cp /mnt/md0/dee2/gene_info/${ORG3}_gene_info.tsv $SRP/GeneInfo.tsv
    cp /mnt/md0/dee2/gene_info/${ORG3}_tx_info.tsv $SRP/TxInfo.tsv

    # get star gene counts
    SE=$(find . | grep se.tsv)
    NUMSE=$(echo $SE | wc -w)
    NUMCOL=$(echo $((NUMSE * 2)) )
    COLS=$(seq 2 2 $NUMCOL)
    COLS=$(echo $COLS | tr ' ' ',')
    ACCS2=$(find . | grep se | cut -d '/' -f2)
    echo GeneID $ACCS2 | tr '| ' '\t' > $SRP/GeneCountMatrix.tsv
    paste $SE | sed 1d | cut -f1,$COLS >> $SRP/GeneCountMatrix.tsv
    ERROR_CNT=0
    RES_COL=$(head $SRP/GeneCountMatrix.tsv | tail -1 | wc -w )
    RES_COL=$((RES_COL-1))
    if [ $RES_COL -ne $NRUNS ] ; then
      ERROR_CNT=$((ERROR_CNT+1))
    fi

    # get kallisto tx counts
    KE=$(find . | grep ke)
    NUMKE=$(echo $KE | wc -w)
    NUMCOL=$(echo $((NUMKE * 5)) )
    COLS=$(seq 4 5 $NUMCOL)
    COLS=$(echo $COLS | tr ' ' ',')
    paste $KE | cut -f$COLS | head -1 | sed 's/^/TranscriptID\t/' \
    | sed 's/_est_counts//g' > $SRP/TxCountMatrix.tsv
    paste $KE | sed 1d | cut -f1,$COLS >> $SRP/TxCountMatrix.tsv
    RES_COL=$(head $SRP/TxCountMatrix.tsv | tail -1 | wc -w )
    RES_COL=$((RES_COL-1))
    if [ $RES_COL -ne $NRUNS ] ; then
      ERROR_CNT=$((ERROR_CNT+1))
    fi

    # get QC data
    QC=$(find . | grep qc)
    NUMQC=$(echo $QC | wc -w)
    NUMQC=$(echo $((NUMQC * 2)) )
    COLS=$(seq 2 2 $NUMQC)
    COLS=$(echo 1 $COLS | tr ' ' ',')
    ACCS2=$(find . | grep se | cut -d '/' -f2)
    echo QCmetric $ACCS2 | tr '| ' '\t' > $SRP/QC_Matrix.tsv
    paste $QC | tr ':' '\t' | head -28 | cut -f$COLS >> $SRP/QC_Matrix.tsv
    RES_COL=$(head $SRP/QC_Matrix.tsv | tail -1 | wc -w )
    RES_COL=$((RES_COL-1))
    if [ $RES_COL -ne $NRUNS ] ; then
      ERROR_CNT=$((ERROR_CNT+1))
    fi

    # get logs
    LOG=$(find . | grep log )
    mkdir -p $SRP/log
    cp $LOG $SRP/log
    RES_COL=$(echo $LOG | wc -w )
    if [ $RES_COL -ne $NRUNS ] ; then
      ERROR_CNT=$((ERROR_CNT+1))
    fi

    cp ../contents.md $SRP/README.md

    zip -r $SRP.zip $SRP
    scp -i ~/.ssh/dee2_2025 $SRP.zip ubuntu@dee2.io:/dee2_data/requests

    cd ..
    rm LOCK
  fi
done
