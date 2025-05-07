#!/bin/bash
set -x

CNT=`echo "$QUERY_STRING" | grep -c '&'`
if [ $CNT -eq "0" ] ; then
  echo "Content-type: text/html"
  echo ""
  echo No datasets selected
  echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;"></FORM>'
  exit
fi

echo "Content-type: text/html"
echo "Content-Type: application/zip"
echo "Content-Encoding: zip"
echo "Content-disposition: attachment; filename=Data.zip"
echo ""

#SPEC=$(echo $QUERY_STRING  | egrep -c '(:|;|}|\{|\[|\]|\/|\\|\@|\<|\>)')
#SPEC=$(echo $QUERY_STRING  | egrep -c '(:|;|}|\{|\[|\]|\/|\\|\@)')
SPEC=$(echo $QUERY_STRING  | egrep -c '(\!|\@|\#|\$|\%|\^|\*|\:|\;|\}|\{|\[|\]|\/|\\|\@)')
if [ $SPEC -gt 0 ] ; then
  echo $QUERY_STRING
  echo "<br>"
  echo "Avoid special characters"
  exit
fi
QUERY_STRING=$(echo $QUERY_STRING | tr -d ':;{}()[]\/<>' )

#ID=$(< /dev/urandom tr -dc A-Za-z0-9 | head -c${1:-32};echo;)
ID=${RANDOM}_${RANDOM}
USRDIR=/dee2_data/tmp/$ID/
mkdir -p $USRDIR
LOGDIR=$USRDIR/logs
mkdir -p $LOGDIR
ORG=`echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2 | tr 'A-Z' 'a-z'`
DATA_DIR=/dee2_data/data/${ORG}/
ACCESSIONS=/dee2_data/metadata/${ORG}_metadata.tsv.cut
METADATA=/dee2_data/metadata/${ORG}_metadata.tsv

cd $DATA_DIR

#Format query string to be compatible with egrep multi-query search
QS=`echo $QUERY_STRING | cut -d '&' -f2- | sed 's/DataSetList=on&//' \
| sed 's/x=/\ /g' | sed 's/ //' | tr -d '&' | sed 's/|//' \
| sed 's/ /_gene.cnt /g' | sed 's/$/_gene.cnt/'`

##remove SRRs with files absent
QS2=`echo $QS | tr ' ' '\n' | cut -d '_' -f1 | tr '\n' ' '`

##get the 3 digit species prefix
PFX=$(echo $ORG | cut -c-3)

#################################################
# Accession numbers
#################################################
head -1 $ACCESSIONS > $USRDIR/MetadataSummary.tsv
echo $QS2 | tr ' ' '\n' | grep -wFf - $ACCESSIONS | sort >> $USRDIR/MetadataSummary.tsv

#################################################
# Accession numbers
#################################################
head -1 $METADATA > $USRDIR/MetadataFull.tsv
echo $QS2 | tr ' ' '\n' | grep -wFf - $METADATA | sort >> $USRDIR/MetadataFull.tsv

#################################################
# Gene info - names and length
#################################################
cp ${PFX}_gene_info.tsv $USRDIR/GeneInfo.tsv
cp ${PFX}_tx_info.tsv $USRDIR/TxInfo.tsv

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
paste $ROWNAMES_QC $QCS | cat <(echo SeqMetric2 $QS2 | tr ' ' '\t') - > $USRDIR/QC_Matrix.tsv

#################################################
# Kallisto Transcript counts
#################################################
TXCOUNTS=$(echo $QS2 | tr ' ' '\n' | awk '{print $1"/"$1"_tx.cnt"}' | tr '\n' ' ' | sed 's/$/\n/')
ROWNAMES_TX=rownames_tx.txt
paste $ROWNAMES_TX $TXCOUNTS > $USRDIR/TxCountMatrix.tsv

#################################################
# README file based on contents
#################################################
cp /usr/lib/cgi-bin/contents.md $USRDIR/README.md

#################################################
# Collect run logs
#################################################
# Collect reports
LOGS=$(echo $QS2 | tr ' ' '\n' | awk '{print $1"/"$1".log"}' | tr '\n' ' ' | sed 's/$/\n/')
cp $LOGS $LOGDIR

cd $USRDIR
zip -r - README.md MetadataSummary.tsv MetadataFull.tsv GeneCountMatrix.tsv QC_Matrix.tsv TxCountMatrix.tsv GeneInfo.tsv TxInfo.tsv logs
find $USRDIR -type d -mmin +60 -maxdepth 1 -exec rm -r {} \;

