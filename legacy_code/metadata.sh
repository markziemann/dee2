#!/bin/bash
ORG=$1
CODEDIR=/data/projects/mziemann/geo2mx_project/v1/code
METAWD=/data/projects/mziemann/geo2mx_project/v1/metadata/metadata_${ORG}
DATADIR=/data/projects/mziemann/geo2mx_project/v1/data/${ORG}
BLACKLIST=${METAWD}/blacklist.txt
ACCESSIONS=${METAWD}/SRXaccessions.txt



cp mineSRAdb.R $METAWD
cd $METAWD
ln $DATADIR/queue/SRR_complete.txt .
cp $CODEDIR/mineSRAdb.R .
Rscript mineSRAdb.R

mine() {
ACC=$1
SRR=`echo $ACC | cut -d _ -f6`
SRA=`echo $ACC | cut -d _ -f3`
SRP=`echo $ACC | cut -d _ -f4`
SRS=`echo $ACC | cut -d _ -f5`
SRX=`echo $ACC | cut -d _ -f2`
grep -w $SRX experiment.txt > ${SRR}.tmp
GSEX=`grep -w GEO $SRR.tmp | cut -f10 | tr -d '"'`
GSE=`echo $GSEX NA | cut -d ' ' -f1`
GSMX=`grep -w GEO $SRR.tmp | cut -f9 | cut -d ':' -f1 | tr -d '"'`
GSM=`echo $GSMX NA | cut -d ' ' -f1`
rm ${SRR}.tmp
( echo $SRR $SRX $SRS $SRP $SRA $GSE $GSM | tr ' ' '\t'
grep -w $SRR run.txt
grep -w $SRX experiment.txt
grep -w $SRS sample.txt
grep -w $SRP study.txt ) | tr '\n' '\t' | tr -s ' \t' | sed 's/$/\n/'
}
export -f mine
sed 1d $ACCESSIONS | grep -vwFf $BLACKLIST  \
| tr '\t' '_' | parallel mine {} > ${ORG}_metadata.txt
