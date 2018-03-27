#!/bin/bash
set -x

##########################
# matrix aggregation script
##########################

mxagg(){
ORG=$1
CODE_DIR=$(pwd)
PROJ_DIR=$(echo $CODE_DIR | rev | cut -d '/' -f2- | rev )
DATA_DIR=$PROJ_DIR/data/$ORG/
MX_DIR=$PROJ_DIR/mx/

SE_MX=$MX_DIR/${ORG}_se.tsv
SE_LIST=$DATA_DIR/${ORG}_se_list.txt
SE_GENES=$DATA_DIR/${ORG}_se_genes.txt

KE_MX=$MX_DIR/${ORG}_ke_counts.tsv
KE_TPM=$MX_DIR/${ORG}_ke_tpm.tsv
KE_LIST=$DATA_DIR/${ORG}_ke_list.txt
KE_GENES=$DATA_DIR/${ORG}_ke_genes.txt

QC_MX=$MX_DIR/${ORG}_qc.tsv
QC_LIST=$DATA_DIR/${ORG}_qc_list.txt
QC_GENES=$DATA_DIR/${ORG}_qc_genes.txt

###################################################
# start with extracting single column data from se qc and ke files
###################################################
grep -v x $SE_LIST > $SE_LIST.tmp ; mv $SE_LIST.tmp $SE_LIST

cd $DATA_DIR

#list of folders and files
LIST=$DATA_DIR/list.txt
if [ ! -r $LIST ] ; then ls > $LIST ; fi
LIST_TIME=$(expr $(date '+%s') - $(stat -c%Y  $LIST) )
if [ "$LIST_TIME" -gt 86400 ] ; then ls > $LIST ; fi

#single column of star gene counts
tsv1(){
TSV=$1
OUT=$(echo $TSV | sed 's/.se.tsv/_gene.cnt/')
if [ ! -r $OUT ] ; then
  awk '{print $NF}' $TSV > $OUT
fi
}
export -f tsv1
find . | grep se.tsv | parallel -j16 tsv1 {}

#rownames genes
echo "GeneID" > rownames_gene.txt
ACC1=$(ls | grep -m1 RR)
cut -f1 $ACC1/$ACC1.se.tsv | sed 1d >> rownames_gene.txt

#single column of qc data
tsv2(){
TSV=$1
OUT=${TSV}l
PFX=$(basename $TSV .qc)
if [ ! -r $OUT ] ; then
  echo $PFX > $OUT
  awk -F: '{print $NF}' $TSV >> $OUT
fi
}
export -f tsv2
find . | grep .qc$ | parallel -j16 tsv2 {}

#rownames qc
echo "SeqMetric" > rownames_qc.txt
ACC1=$(ls | grep -m1 RR)
cut -d ':' -f1 $ACC1/$ACC1.qc >> rownames_qc.txt

#single column of kallisto tx counts
tsv3(){
TSV=$1
OUT=$(echo $TSV | sed 's/.ke.tsv/_tx.cnt/')
PFX=$(basename $TSV .ke.tsv)
if [ ! -r $OUT ] ; then
  echo $PFX > $OUT
  awk '{print $4}' $TSV | sed 1d >> $OUT
fi
}
export -f tsv3
find . | grep ke.tsv | parallel -j16 tsv3 {}

#rownames transcripts
echo "TxID" > rownames_tx.txt
ACC1=$(ls | grep -m1 RR)
cut -f1 $ACC1/$ACC1.ke.tsv | sed 1d >> rownames_tx.txt

#################################################
# Now paste STAR gene count matrix
################################################
grep -v x $SE_LIST > $SE_LIST.tmp ; mv $SE_LIST.tmp $SE_LIST

#divide the list into sets of 1000

rm *split
split -l 500 --additional-suffix=split $SE_LIST

paste1(){
  SPLIT=$1
  cut -d '/' -f3 $SPLIT | sed 's/se.tsv/se/' | tr '\n' '\t' | sed 's#\t$#\n#'  > $SPLIT.tsv
  paste  $(cat $SPLIT) \
  | awk '{for(i=2;i<=NF;i=i+2){printf "%s\t", $i}{printf "%s", RS}}' \
  | sed 1d >> $SPLIT.tsv
}
export -f paste1
parallel -j16 paste1 ::: *split

echo GeneID > $SE_GENES
cut -f1 $(head -1 $SE_LIST ) | sed 1d >> $SE_GENES

paste $SE_GENES *split*.tsv > $SE_MX
head $SE_MX | cut -f-5

#tidy up temporary files
rm $SE_GENES *split *split.tsv

#zip the mx
pbzip2 -kf $SE_MX

###################################################
# Now kallisto counts (ke)
###################################################
grep -v x $KE_LIST > $KE_LIST.tmp ; mv $KE_LIST.tmp $KE_LIST

#divide the list into sets of 1000
cd $DATA_DIR
rm *split
split -l 500 --additional-suffix=split $KE_LIST

paste2(){
  SPLIT=$1
  cut -d '/' -f3 $SPLIT | sed 's/ke.tsv/ke_counts/' | tr '\n' '\t' | sed 's#\t$#\n#'  > $SPLIT.tsv
  paste  $(cat $SPLIT) \
  | awk '{for(i=4;i<=NF;i=i+5){printf "%s\t", $i}{printf "%s", RS}}' \
  | sed 1d >> $SPLIT.tsv
}
export -f paste2
parallel -j16 paste2 ::: *split

echo TranscriptID > $KE_GENES
cut -f1 $(head -1 $KE_LIST ) | sed 1d >> $KE_GENES
head $KE_GENES

paste $KE_GENES *split*.tsv > $KE_MX
head $KE_MX | cut -f-5

#tidy up temporary files
rm $KE_GENES *split *split.tsv

#zip the mx
pbzip2 -kf $KE_MX

###################################################
# Now kallisto tpm (ke)
###################################################
grep -v x $KE_LIST > $KE_LIST.tmp ; mv $KE_LIST.tmp $KE_LIST

#divide the list into sets of 1000
cd $DATA_DIR
rm *split
split -l 500 --additional-suffix=split $KE_LIST

paste3(){
  SPLIT=$1
  cut -d '/' -f3 $SPLIT | sed 's/ke.tsv/ke_tpm/' | tr '\n' '\t' | sed 's#\t$#\n#'  > $SPLIT.tsv
  paste  $(cat $SPLIT) \
  | awk '{for(i=5;i<=NF;i=i+5){printf "%s\t", $i}{printf "%s", RS}}' \
  | sed 1d >> $SPLIT.tsv
}
export -f paste3
parallel -j16 paste3 ::: *split

echo TranscriptID > $KE_GENES
cut -f1 $(head -1 $KE_LIST ) | sed 1d >> $KE_GENES
head $KE_GENES

paste $KE_GENES *split*.tsv > $KE_TPM
head $KE_TPM | cut -f-5

#tidy up temporary files
rm $KE_GENES *split *split.tsv

#zip the mx
pbzip2 -kf $KE_TPM

###################################################
# Now qc metrics
###################################################
grep -v x $QC_LIST > $QC_LIST.tmp ; mv $QC_LIST.tmp $QC_LIST

#divide the list into sets of 1000
cd $DATA_DIR
rm *split
split -l 500 --additional-suffix=split $QC_LIST

paste4(){
  SPLIT=$1
  cut -d '/' -f3 $SPLIT | sed 's/se.tsv/se/' | tr '\n' '\t' | sed 's#\t$#\n#'  > $SPLIT.tsv
  paste  $(cat $SPLIT) \
  | tr ':' '\t' \
  | awk '{for(i=2;i<=NF;i=i+2){printf "%s\t", $i}{printf "%s", RS}}' \
  | sed 1d >> $SPLIT.tsv
}
export -f paste4
parallel -j16 paste4 ::: *split

echo QC_metric > $QC_GENES
cut -d ':' -f1 $(head -1 $QC_LIST ) | sed 1d >> $QC_GENES

paste $QC_GENES *split*.tsv > $QC_MX
head $QC_MX | cut -f-5

#tidy up temporary files
rm $QC_GENES *split *split.tsv

#zip the mx
pbzip2 -kf $QC_MX

#SCP
scp -i ~/.ssh/cloud/cloud2.key $SE_MX.bz2 $KE_MX.bz2 $KE_TPM.bz2 $QC_MX.bz2 ubuntu@118.138.240.228:/mnt/dee2_data/mx/
cd $CODE_DIR
}

export -f mxagg

#for ORG in athaliana celegans dmelanogaster drerio ecoli hsapiens mmusculus rnorvegicus scerevisiae ; do
#  mxagg $ORG
#done

ORG=$1
mxagg $ORG
