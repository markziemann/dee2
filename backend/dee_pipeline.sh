#!/bin/bash

CORES=16
ORG=$1
DIR=/mnt/hdd1/dee2/data/$ORG
MDDIR=/mnt/hdd1/dee2/sradb
MXDIR=/mnt/hdd1/dee2/mx
SEH5=$MXDIR/${ORG}_se.h5
KEH5=$MXDIR/${ORG}_ke.h5
QCH5=$MXDIR/${ORG}_qc.h5
MD=$MDDIR/${ORG}_metadata.tsv.cut
MDC=$MDDIR/${ORG}_accessions.tsv.bz2
cd $DIR

# prepa a list of validated SRA runs
cut -f1 $MD | sed 1d > $MD.tmp

# get checksums
md5sum $SEH5 > $SEH5.md5 &
md5sum $KEH5 > $KEH5.md5 &
md5sum $QCH5 > $QCH5.md5 &
wait

####
echo "fix permissions"
####
perm775(){
ACC=$1
chmod 775 $ACC 2> /dev/null
}
export -f perm775
parallel -j$CORES perm775 :::: $MD.tmp

#remove the tmp file
rm $MD.tmp

# metadata
pbzip2 -p$CORES -c $MD > $MDC

# checksums
rm  $MXDIR/checksums.md5
wait
cat $MXDIR/*md5 | sed 's#/mnt/hdd1/dee2/mx/##' > $MXDIR/checksums.md5
# transfer to webserver
scp -i ~/.ssh/dee2_2025 $SEH5 $KEH5 $QCH5 $MDC $MXDIR/checksums.md5 \
ubuntu@118.138.235.221:/dee2_data/bulk/

