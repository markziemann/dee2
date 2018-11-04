#!/bin/bash

CORES=16
ORG=$1
DIR=/mnt/md0/dee2/data/$ORG
FINLIST=$DIR/${ORG}_fin_list.txt
SELIST=$DIR/${ORG}_se_list.txt
KELIST=$DIR/${ORG}_ke_list.txt
QCLIST=$DIR/${ORG}_qc_list.txt
VALLIST=$DIR/${ORG}_val_list.txt
QUEUELIST=$DIR/${ORG}_queue_list.txt
MXDIR=/mnt/md0/dee2/mx
SEMX=$MXDIR/${ORG}_se.tsv
KEMX=$MXDIR/${ORG}_ke.tsv
QCMX=$MXDIR/${ORG}_qc.tsv

cd $DIR

##################
# Run some checks
##################
#remove directories without finished or validated files
for i in $DIR/*RR* ; do
  CNT=$(ls $i | egrep -c '(finished|validated)')
  if [ $CNT -lt 1 ] ; then rm -fr $i ; fi
done

####
echo "set up list of files to be checked and added to the repo"
###
find $DIR | grep finished | rev | cut -d '/' -f1 |rev | cut -d '.' -f1 > $FINLIST

####
echo "remove directories that are incomplete and add the accession number to the queue"
###
for SRR in $(cat $FINLIST ) ; do
  chmod -R +w $DIR/$SRR
  SE_TSV=$DIR/$SRR/$SRR.se.tsv
  KE_TSV=$DIR/$SRR/$SRR.ke.tsv
  QC_TSV=$DIR/$SRR/$SRR.qc

  if [ ! -e  "$SE_TSV" ] || [ ! -e  "$KE_TSV" ] || [ ! -e  "$QC_TSV" ] ; then
    rm -rf $DIR/$SRR
    echo $SRR
  fi
done >> $QUEUELIST

#remove duplicate entries
sort -u -o $QUEUELIST $QUEUELIST

####
echo "update list of files to be checked and added to the repo"
###
find $DIR | grep finished | rev | cut -d '/' -f1 |rev | cut -d '.' -f1 > $FINLIST

####
echo "SE checks"
####
parallel -j$CORES wc -l {}/{}.se.tsv  :::: $FINLIST > $SELIST
NROW=$(awk '{print $1}' $SELIST | sort | uniq -c | sort -k1nr | head -1 | awk '{print $NF}')
awk -v n=$NROW '$1!=n {print $2}' $SELIST | cut -d '/' -f2 | parallel -j$CORES rm -rf {}

awk -v n=$NROW '$1==n' $SELIST \
| rev | cut -d '/' -f1 | rev \
| cut -d '.' -f1 > $SELIST.tmp \
 && mv $SELIST.tmp $SELIST

####
echo "KE checks"
####
parallel -j$CORES wc -l {}/{}.ke.tsv :::: $FINLIST > $KELIST
NROW=$(awk '{print $1}' $KELIST | sort | uniq -c | sort -k1nr | head -1 | awk '{print $NF}')
awk -v n=$NROW '$1!=n {print $2}' $KELIST | cut -d '/' -f2 | parallel -j$CORES rm -rf {}

awk -v n=$NROW '$1==n' $KELIST \
| rev | cut -d '/' -f1 | rev \
| cut -d '.' -f1 > $KELIST.tmp \
 && mv $KELIST.tmp $KELIST

####
echo "QC checks"
####
parallel -j$CORES wc -l {}/{}.qc :::: $FINLIST > $QCLIST
NROW=$(awk '{print $1}' $QCLIST | sort | uniq -c | sort -k1nr | head -1 | awk '{print $NF}')
awk -v n=$NROW '$1!=n {print $2}' $QCLIST | cut -d '/' -f2 | parallel -j$CORES rm -rf {}

awk -v n=$NROW '$1==n' $QCLIST \
| rev | cut -d '/' -f1 | rev \
| cut -d '.' -f1 > $QCLIST.tmp \
 && mv $QCLIST.tmp $QCLIST

##################
echo "Run aggregation"
##################

cat $SELIST $KELIST $QCLIST | sort | uniq -c | awk '$1==3 {print $2}' > $VALLIST

##################
echo "validate"
##################
fin_agg(){
DIRPATH=$1
chmod -R +w $DIRPATH
rename -f 's/.finished/.validated/' $DIRPATH/*.finished
}
export -f fin_agg
parallel -j$CORES fin_agg :::: $VALLIST

find $DIR | grep validated | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | sort -u > $VALLIST

##################
echo "Upload validated files"
##################
while read line ; do
  rsync -avzh -e "ssh -i  ~/.ssh/monash/cloud2.key " $line ubuntu@118.138.234.131:/dee2_data/data/${ORG}/
done <  $VALLIST

####
echo "se agg"
####
se_agg(){
ACC=$1
chmod +w -R $ACC
awk '{print $NF}' $ACC/$ACC.se.tsv > $ACC/${ACC}_gene.cnt
sed 1d $ACC/$ACC.se.tsv | sed "s/^/${ACC}\t/"
}
export -f se_agg
parallel -j$CORES se_agg :::: $VALLIST | pbzip2 -c -j$CORES > $SEMX.bz2

####
echo "ke_agg"
####
ke_agg(){
ACC=$1
awk '{print $4}' $ACC/$ACC.ke.tsv | sed 's/_est_counts//'> $ACC/${ACC}_tx.cnt
sed 1d $ACC/$ACC.ke.tsv | cut -f1,4 | sed "s/^/${ACC}\t/"
}
export -f ke_agg
parallel -j$CORES ke_agg :::: $VALLIST | pbzip2 -c -j$CORES > $KEMX.bz2

####
echo "qc agg"
####
qc_agg(){
ACC=$1
cut -d ':' -f2 $ACC/$ACC.qc > $ACC/$ACC.qcl
sed 's/:/\t/' $ACC/$ACC.qc | sed "s/^/${ACC}\t/"
}
export -f qc_agg
parallel -j$CORES qc_agg :::: $VALLIST | pbzip2 -c -j$CORES > $QCMX.bz2

scp -i ~/.ssh/monash/cloud2.key $SEMX.bz2 $KEMX.bz2 $QCMX.bz2  ubuntu@118.138.234.131:/dee2_data/mx

fin_agg(){
DIRPATH=$1
chmod 0544 $DIRPATH
}
export -f fin_agg
parallel -j$CORES fin_agg :::: $VALLIST

comm -23 <(sort $FINLIST) <(sort $VALLIST) > $QUEUELIST

#dat share $MXDIR &
#cd $MXDIR
#for BZ2 in *bz2 ; do rclone copy $BZ2 drive:Public/DEE ; done

