#!/bin/bash
ORG=$1
METAWD=/data/projects/mziemann/geo2mx_project/v1/metadata/metadata_${ORG}
cd $METAWD

CODE=/data/projects/mziemann/geo2mx_project/v1/code
DATADIR=/data/projects/mziemann/geo2mx_project/v1/data/${ORG}
BLACKLIST=/data/projects/mziemann/geo2mx_project/v1/metadata/metadata_${ORG}/blacklist.txt
HEADER=/data/projects/mziemann/geo2mx_project/v1/metadata/metadata_header.txt
QC_REPORTS=/data/projects/mziemann/geo2mx_project/v1/reports/${ORG}/

cp $CODE/mineSRAdb.R $METAWD
ACCESSIONS=$METAWD/SRXaccessions.txt
Rscript mineSRAdb.R

mine() {
ACC=$1
ORG=$2
DATADIR=/data/projects/mziemann/geo2mx_project/v1/data/${ORG}
SRR=`echo $ACC | cut -d _ -f2`
QC=`tail -1 $DATADIR/$SRR/$SRR.qc | cut -d ':' -f2-`
if [ -z "$QC" ]; then QC=INCOMPLETE_UNAVAILABLE ; fi
SRA=`echo $ACC | cut -d _ -f3`
SRP=`echo $ACC | cut -d _ -f4`
SRS=`echo $ACC | cut -d _ -f5`
SRX=`echo $ACC | cut -d _ -f6`
grep -w $SRX experiment.txt > ${SRR}.tmp
GSEX=`grep -w GEO $SRR.tmp | cut -f10 | tr -d '"'`
GSE=`echo $GSEX NA | cut -d ' ' -f1`
GSMX=`grep -w GEO $SRR.tmp | cut -f9 | cut -d ':' -f1 | tr -d '"'`
GSM=`echo $GSMX NA | cut -d ' ' -f1`
rm ${SRR}.tmp
( echo $SRR $QC $SRX $SRS $SRP $SRA $GSE $GSM | tr ' ' '\t'
grep -w $SRR run.txt
grep -w $SRX experiment.txt
grep -w $SRS sample.txt
grep -w $SRP study.txt ) | tr '\n' '\t' | tr -s ' \t' | sed 's/$/\n/'
}
export -f mine
sed 1d $ACCESSIONS | grep -vwFf $BLACKLIST  \
| tr '\t' '_' | parallel mine {} $ORG | cat $HEADER - > ${ORG}_metadata.txt
