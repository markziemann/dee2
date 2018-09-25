#!/bin/bash
set -x

#QUERY_STRING='org=ecoli&x=SRR2044879'

CNT=`echo "$QUERY_STRING" | grep -c '&'`
if [ $CNT -eq "0" ] ; then
  echo "Content-type: text/html"
  echo ""
  echo No datasets selected
  echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;"></FORM>'
  exit
fi

echo "Content-type: text"
echo "Content-Type: application/zip"
echo "Content-Encoding: zip"
echo "Content-disposition: attachment; filename=Data.zip"
echo ""


ID=${RANDOM}_${RANDOM}
USRDIR=/mnt/tmp/$ID/
mkdir -p $USRDIR
LOGDIR=$USRDIR/logs
mkdir -p $LOGDIR
ORG=`echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2 | tr 'A-Z' 'a-z'`
DATA_DIR=/dee2_data/data/${ORG}/
LIST=$DATA_DIR/${ORG}_val_list.txt

cd $DATA_DIR

#Format query string to be compatible with egrep multi-query search
QS=`echo $QUERY_STRING | cut -d '&' -f2- | sed 's/DataSetList=on&//' \
| sed 's/x=/\ /g' | sed 's/ //' | tr -d '&' | sed 's/|//' \
| sed 's/ /_gene.cnt /g' | sed 's/$/_gene.cnt/'`

#Omit entries with no datasets
##Check list age
LIST_TIME=`expr $(date +%s) - $(date +%s -r $LIST)`

##remove SRRs with files absent
QS2=`echo $QS | tr ' ' '\n' | cut -d '_' -f1 | sort - $LIST | uniq -d | tr '\n' ' '`

#################################################
# STAR Gene counts
#################################################
GENECOUNTS=$(echo $QS2 | tr ' ' '\n' | awk '{print $1"/"$1"_gene.cnt"}' | tr '\n' ' ' | sed 's/$/\n/')
ROWNAMES_GENE=rownames_gene.txt
paste $ROWNAMES_GENE $GENECOUNTS > $USRDIR/GeneCountMatrix.tsv

#################################################
# QC matrix
#################################################
QCS=$(echo $QS2 | tr ' ' '\n' | awk '{print $1"/"$1".qcl"}' | tr '\n' ' ' | sed 's/$/\n/')
ROWNAMES_QC=rownames_qc.txt
paste $ROWNAMES_QC $QCS > $USRDIR/QC_Matrix.tsv

#################################################
# Kallisto Transcript counts
#################################################
TXCOUNTS=$(echo $QS2 | tr ' ' '\n' | awk '{print $1"/"$1"_tx.cnt"}' | tr '\n' ' ' | sed 's/$/\n/')
ROWNAMES_TX=rownames_tx.txt
paste $ROWNAMES_TX $TXCOUNTS > $USRDIR/TxCountMatrix.tsv

#################################################
# Collect run logs
#################################################
# Collect reports
LOGS=$(echo $QS2 | tr ' ' '\n' | awk '{print $1"/"$1".log"}' | tr '\n' ' ' | sed 's/$/\n/')
cp $LOGS $LOGDIR

cd $USRDIR
zip -r - GeneCountMatrix.tsv QC_Matrix.tsv TxCountMatrix.tsv logs
find $USRDIR -type d -mmin +60 -maxdepth 1 -exec rm -r {} \;
