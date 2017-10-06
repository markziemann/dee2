#!/bin/bash
set -x

ORG=$1
#Set main directory structure
cd /data/projects/mziemann/geo2mx_project/v1/code
WD=/data/projects/mziemann/geo2mx_project/v1

CODEDIR=$WD/code
PRLL=$CODEDIR/prll_jan18.sh
MXAGG=$CODEDIR/mxagg.sh
TMX=$CODEDIR/tmx.awk
MDSH=$CODEDIR/metadata2.sh

CONFIGDIR=$WD/config
CONFIG=$CONFIGDIR/${ORG}.cfg
DATADIR=$WD/data/${ORG}
INPROGRESS=$WD/data/${ORG}/inprogress.tmp
AGGDIR=$DATADIR/agg
INDV_FILES=$DATADIR/agg/individual_files
UPDATE=$WD/firefox_update/${ORG}_sra_result.csv

REFGENOME=$WD/refgenomes/refgenome_${ORG}
GENENAMES=$REFGENOME/${ORG}_genenames.txt
MXDIR=$WD/mx

METADATA=$WD/metadata/metadata_${ORG}

#below variable caused bug in SRR handling
#SRRCOMPLETEDLIST=$METADATA/SRRdone.txt

ANALYSISREPORTS=$WD/reports/reports_${ORG}
QCREPORTS=$WD/reports/qc_${ORG}


#############################
# Touch progress file
#############################
touch $INPROGRESS

#############################
echo Queue SRRs to run
#############################
#Specify files for queuing
QDIR=$METADATA/queue
SRX=$QDIR/SRX.txt
URLSALL=$QDIR/URLs_all.txt
SRRALL=$QDIR/SRR_all.txt
SRRCOMPLETE=$QDIR/SRR_complete.txt
SRRTODO=$QDIR/SRR_todo.txt
URLSTODO=$QDIR/URLS_todo.txt
ACC=$METADATA/SRXaccessions.txt
BASEURL=ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra

mkdir -p $QDIR

cd $METADATA

cut -d '"' -f2 $UPDATE | sed 1d > $SRX

cp $CODEDIR/queueSRAdb.R .
cp $SRX .
Rscript queueSRAdb.R

find $DATADIR | grep _gene.cnt.summary | rev | cut -d '/' -f1 | rev \
| cut -d '_' -f1 | sort > $SRRCOMPLETE

sed 1d $ACC | cut -f6 | sort > $SRRALL

comm -23 $SRRALL $SRRCOMPLETE > $SRRTODO

for SRR in `cat $SRRTODO` ; do
PFX1=`echo $SRR | cut -c-3`
PFX2=`echo $SRR | cut -c-6`
URL=${BASEURL}/${PFX1}/${PFX2}/${SRR}/${SRR}.sra
echo $URL
done > $URLSTODO

cp $URLSTODO $DATADIR
CNT=`wc -l < $URLSTODO`
if [ $CNT -eq 0 ] ; then
echo No new data sets to analyse, quitting
exit
fi
echo There are $CNT datasets in the queue

#############################
echo Actually run the analysis
#############################
cd $DATADIR
cat $CONFIG
nohup bash $PRLL $CONFIG

#############################
echo Aggregate the data
#############################
mkdir -p $AGGDIR
cd $AGGDIR
cp $MXAGG $TMX .
bash $MXAGG $ORG $GENENAMES
cp ${ORG}*mx $MXDIR

#############################
echo Collect reports
#############################
mkdir -p $ANALYSISREPORTS
cd $ANALYSISREPORTS
find $DATADIR | grep RR | egrep '(log$)' | parallel cp {} .
find $DATADIR | grep RR | egrep '(_gene.cnt.summary$)' | parallel cp {} .

find . | grep log | parallel "sed -i 's#/data/projects/mziemann/#/path/#g;s#/data/home/mziemann/#/path/#' {}"

mkdir -p $QCREPORTS
cd $QCREPORTS
find $DATADIR | grep RR | egrep '(.qc$)' | parallel cp {} .
for QC in *qc ; do
  HEADER=`echo $QC | sed 's/.qc/v1/'`
  echo $HEADER > $QC.t
  cut -d ':' -f2- $QC >> $QC.t
done

echo Pipeline finished

#############################
echo Fetch metadata
#############################
mkdir -p $METADATA
cd $METADATA
cut -f1 $AGGDIR/${ORG}_gene.mx | sed 1d | cut -d '.' -f1 > $SRRCOMPLETE
bash $MDSH $ORG

#############################
# Remove progress file
#############################
rm $INPROGRESS

#############################
echo Send data to webserver
#############################
cd $DATADIR

#Transfer individual count files to webserver
scp -r $INDV_FILES mziemann@172.16.0.53:/var/www/data/$ORG
touch ${ORG}.finished
scp ${ORG}.finished mziemann@172.16.0.53:/var/www/data/

#Transfer full reports to webserver
scp -r /data/projects/mziemann/geo2mx_project/v1/reports/reports_${ORG} mziemann@172.16.0.53:/var/www/metadata/reports_${ORG}2

#Transfer qc reports to webserver
scp -r /data/projects/mziemann/geo2mx_project/v1/reports/qc_${ORG}  mziemann@172.16.0.53:/var/www/metadata/qc_${ORG}2

#Transfer metadata to webserver
scp /data/projects/mziemann/geo2mx_project/v1/metadata/metadata_${ORG}/${ORG}_metadata.txt mziemann@172.16.0.53:/var/www/metadata/${ORG}_metadata.txt2

#Transfer finished file
scp ${ORG}.finished mziemann@172.16.0.53:/var/www/metadata/${ORG}.finished

#Bulk data
rm ${ORG}.finished
cd /data/projects/mziemann/geo2mx_project/v1/bulk/${ORG}_bulk
ln /data/projects/mziemann/geo2mx_project/v1/mx/${ORG}_gene_n.mx .
ln -s /data/projects/mziemann/geo2mx_project/v1/reports/reports_${ORG} .
ln /data/projects/mziemann/geo2mx_project/v1/metadata/metadata_${ORG}/${ORG}_metadata.txt .
cd ..
rm ${ORG}_bulk.zip
zip -qr ${ORG}_bulk.zip ${ORG}_bulk
scp ${ORG}_bulk.zip mziemann@172.16.0.53:/var/www/html/bulkdata/${ORG}_bulk.zip2
touch ${ORG}_bulk.zip.finished
scp ${ORG}_bulk.zip.finished mziemann@172.16.0.53:/var/www/html/bulkdata/


