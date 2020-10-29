#!/bin/bash

CORES=16
ORG=$1
DIR=/mnt/md0/dee2/data/$ORG
MDDIR=/mnt/md0/dee2/sradb
MXDIR=/mnt/md0/dee2/mx
SEMX=$MXDIR/${ORG}_se.tsv
KEMX=$MXDIR/${ORG}_ke.tsv
QCMX=$MXDIR/${ORG}_qc.tsv
MD=$MDDIR/${ORG}_metadata.tsv.cut
MDC=$MDDIR/${ORG}_accessions.tsv.bz2
cd $DIR

# prepa a list of validated SRA runs
cut -f1 $MD | sed 1d > $MD.tmp


####
echo "se agg"
####
se_agg(){
ACC=$1
chmod +w -R $ACC
zcat $ACC/$ACC.se.tsv.gz | sed 1d | sed "s/^/${ACC}\t/"
}
export -f se_agg
parallel -j$CORES se_agg :::: $MD.tmp | pbzip2 -c -j$CORES > $SEMX.bz2
md5sum $SEMX.bz2 > $SEMX.bz2.md5

####
echo "ke_agg"
####
ke_agg(){
ACC=$1
zcat $ACC/$ACC.ke.tsv.gz | sed 1d | cut -f1,4 | sed "s/^/${ACC}\t/"
}
export -f ke_agg
parallel -j$CORES ke_agg :::: $MD.tmp | pbzip2 -c -j$CORES > $KEMX.bz2
md5sum $KEMX.bz2 > $KEMX.bz2.md5

####
echo "qc agg"
####
qc_agg(){
ACC=$1
sed -i -n '/:/p'  $ACC/$ACC.qc
cut -d ':' -f2 $ACC/$ACC.qc > $ACC/$ACC.qcl
sed 's/:/\t/' $ACC/$ACC.qc | sed "s/^/${ACC}\t/"
}
export -f qc_agg
parallel -j$CORES qc_agg :::: $MD.tmp | pbzip2 -c -j$CORES > $QCMX.bz2
md5sum $QCMX.bz2 > $QCMX.bz2.md5

####
echo "fix permissions"
####
perm775(){
ACC=$1
chmod 775 $ACC
}
export -f perm775
parallel -j$CORES perm775 :::: $MD.tmp

#remove the tmp file
rm $MD.tmp

# metadata
pbzip2 -c $MD > $MDC

# checksums
rm  $MXDIR/checksums.md5
cat $MXDIR/*md5 | sed 's#/mnt/md0/dee2/mx/##' > $MXDIR/checksums.md5
# transfer to webserver
scp -i ~/.ssh/monash/cloud2.key $SEMX.bz2 $KEMX.bz2 $QCMX.bz2 $MDC $MXDIR/checksums.md5 \
ubuntu@118.138.234.94:/dee2_data/mx

