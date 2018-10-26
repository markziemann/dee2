#!/bin/bash

## Package large datasets for public use

for ORG in athaliana celegans dmelanogaster drerio ecoli hsapiens mmusculus rnorvegicus scerevisiae ; do

FULL_METADATA=${ORG}_metadata.tsv.cut
DEE2_ACCESSIONS=${ORG}_accessions.tsv

DIR=big_proj/${ORG}
if [ ! -d $DIR ] ; then mkdir -p $DIR ; fi

cut -d ' ' -f2 $DEE2_ACCESSIONS \
| sed 1d | sort | uniq -c | awk '$1>200' \
| head -3 | while read line ; do
  SRP=$(echo $line | cut -d ' ' -f2)
  CNT_DEE=$(echo $line | cut -d ' ' -f1)
  CNT_SRA=$(grep -wc $SRP $FULL_METADATA )
  echo $SRP $CNT_SRA $CNT_DEE
  if [ $CNT_DEE -eq $CNT_SRA ] ; then
    #fetch the geo series number
    GSE=$(grep -wm1 $SRP ${ORG}_metadata.tsv.cut | cut -f6)
    SRP_DIR=$DIR/${SRP}_${GSE}
    if [ ! -r $SRP_DIR.zip ] ; then
      mkdir -p $SRP_DIR

      #package the gene counts
      SRRLIST=$(grep -w $SRP $DEE2_ACCESSIONS | cut -d ' ' -f4 \
      | awk '{print $0"/"$0}' | sed "s@^@../data/${ORG}/@" \
      | sed 's@$@_gene.cnt@'  )
      paste ../data/$ORG/rownames_gene.txt $SRRLIST > $SRP_DIR/GeneCountMatrix.tsv

      #package the transcript counts
      SRRLIST=$(grep -w $SRP $DEE2_ACCESSIONS | cut -d ' ' -f4 \
      | awk '{print $0"/"$0}' | sed "s@^@../data/${ORG}/@" \
      | sed 's@$@_tx.cnt@'  )
      paste ../data/$ORG/rownames_tx.txt $SRRLIST > $SRP_DIR/TxCountMatrix.tsv

      #package the QC metrics
      SRRLIST=$(grep -w $SRP $DEE2_ACCESSIONS | cut -d ' ' -f4 \
      | awk '{print $0"/"$0}' | sed "s@^@../data/${ORG}/@" \
      | sed 's@$@.qcl@'  )

      SRRLIST2=$(grep -w $SRP $DEE2_ACCESSIONS | cut -d ' ' -f4 | tr '\n' ' ')

      paste ../data/$ORG/rownames_qc.txt $SRRLIST \
      | cat <(echo SeqMetric $SRRLIST2 | tr ' ' '\t') - > $SRP_DIR/QC_Matrix.tsv

      #package the logs
      mkdir -p $SRP_DIR/logs

      SRRLIST=$(grep -w $SRP $DEE2_ACCESSIONS | cut -d ' ' -f4 \
      | awk '{print $0"/"$0}' | sed "s@^@../data/${ORG}/@" \
      | sed 's@$@.log@'  )

      cp $SRRLIST $SRP_DIR/logs

      #zip up folder
      zip -r $SRP_DIR.zip $SRP_DIR
      rm -rf $SRP_DIR
    fi
  fi
done

#scp -i ~/.ssh/monash/cloud2.key $DIR ubuntu@118.138.234.131:/dee2_data/bulk
rsync -azvh -e "ssh -i ~/.ssh/monash/cloud2.key" $DIR ubuntu@118.138.234.131:/dee2_data/bulk
done
