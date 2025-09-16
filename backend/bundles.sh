#!/bin/bash

# ulimit issue
#if [ $(ulimit -n) -lt 1000000 ] ; then
#  echo Error: need to increase ULIMIT for this to work properly
#  echo 'sudo sh -c "ulimit -n 1000000 && exec su $LOGNAME"'
#  exit 1
#fi

cd ../sradb

ORG=$1

#ORG in athaliana celegans dmelanogaster drerio ecoli rnorvegicus scerevisiae bdistachyon gmax hvulgare ptrichocarpa sbicolor slycopersicum stuberosum taestivum vvinifera  hsapiens mmusculus ; do

FULL_METADATA=${ORG}_metadata.tsv.cut
DEE2_ACCESSIONS=${ORG}_accessions.tsv
DIR=big_proj/${ORG}
if [ ! -d $DIR ] ; then mkdir -p $DIR ; fi

cut -d ' ' -f2 $DEE2_ACCESSIONS \
| sed 1d | sort | uniq -c | awk '$1>5' \
| while read line ; do
  SRP=$(echo $line | cut -d ' ' -f2)
  CNT_DEE=$(echo $line | cut -d ' ' -f1)
  CNT_SRA=$(grep -wc $SRP $FULL_METADATA )
  echo $SRP $CNT_SRA $CNT_DEE
  if [ $CNT_DEE -eq $CNT_SRA ] ; then

    #fetch the geo series number
    GSE=$(grep -wm1 $SRP ${ORG}_metadata.tsv.cut | cut -f7 | grep GSE | cut -d ';' -f1)
    if [ -z "$GSE" ] ; then
      GSE="NA"
    fi
    echo GSE is $GSE

    SRP_DIR=$DIR/${SRP}_${GSE}
    if [ ! -r $SRP_DIR.zip ] ; then
      mkdir -p $SRP_DIR

      #package the gene and tx info
      PFX=$(echo $ORG | cut -c-3)
      cp ../gene_info/${PFX}_gene_info.tsv  $SRP_DIR/GeneInfo.tsv
      cp ../gene_info/${PFX}_tx_info.tsv $SRP_DIR/TxInfo.tsv

      #package the metadata summary
      head -1 $FULL_METADATA > $SRP_DIR/MetadataSummary.tsv
      grep -w $SRP $DEE2_ACCESSIONS | cut -d ' ' -f4 | grep -wFf - $FULL_METADATA >> $SRP_DIR/MetadataSummary.tsv

      #package the metadata in full
      COMPLETE_METADATA=${ORG}_metadata.tsv
      head -1 $COMPLETE_METADATA > $SRP_DIR/MetadataFull.tsv
      grep -w $SRP $DEE2_ACCESSIONS | cut -d ' ' -f4 | grep -wFf - $COMPLETE_METADATA >> $SRP_DIR/MetadataFull.tsv

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
      cd $DIR
      zip -r ${SRP}_${GSE}.zip ${SRP}_${GSE} ../README.md
      rm -rf ${SRP}_${GSE}
      cd -
    fi
  fi
done

echo rsync to webserver
IP_ADD=$(dig +short dee2.io)
chmod -R 775 $DIR
rsync  --delete -azvh -e "ssh -i ~/.ssh/dee2_2025" $DIR ubuntu@${IP_ADD}:/dee2_data/huge

cd ../code
