
#!/bin/bash

echo "Content-type: text"
echo "Content-Type: application/zip"
echo "Content-Encoding: zip"
echo "Content-disposition: attachment; filename=Data.zip"
echo ""

#use ssh to run a command on another server
#ssh -i ~/.ssh/monash/id_rsa mziemann@118.138.246.227 'progress'
#QUERY_STRING='org=Celegans&DataSetList=on&x=SRR1175711&x=SRR1175712&x=SRR1175713&x=SRR1175714'
CNT=`echo "$QUERY_STRING" | grep -c '&'`
if [ $CNT -eq "0" ] ; then
  exit
fi

ID=`< /dev/urandom tr -dc A-Za-z0-9 | head -c${1:-32};echo;`
TMPDIR=/tmp
USRDIR=$TMPDIR/$ID
LOGS=$TMPDIR/$ID/logs
mkdir -p $LOGS
#echo $QUERY_STRING
ORG=`echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2 | tr 'A-Z' 'a-z'`
#echo $ORG
#DATA=/var/www/data/${ORG}_gene_n.mx

#Format query string to be compatible with egrep multi-query search
QS=`echo $QUERY_STRING | cut -d '&' -f2- | sed 's/DataSetList=on&//' \
| sed 's/x=/\ /g' | sed 's/ //' | tr -d '&' | sed 's/|//' \
| sed 's/ /_gene.cnt /g' | sed 's/$/_gene.cnt/'`
#echo "$QS"

#echo 1 >> /var/log/apache2/error.log

#exit

cd /var/www/data/${ORG}_data/

#Omit entries with no datasets
##Check list age
LIST_TIME=`expr $(date +%s) - $(date +%s -r list.txt)`

##Generate list if older than 1 day
if [ ! -r list.txt ] ; then ls > list.txt ; fi
if [ "$LIST_TIME" -gt 86400 ] ; then ls > list.txt ; fi

##remove SRRs with files absent
QS2=`echo $QS | tr ' ' '\n' | sort - list.txt | uniq -d | tr '\n' ' '`

paste rownames.txt $QS2 > $USRDIR/CountMatrix.xls
#echo $USRDIR/CountMatrix.xls

# Collect qc reports and makeqc matrix
QCDIR=/var/www/metadata/qc_$ORG
cd $QCDIR

##Generate list if not existing
if [ ! -r list.txt ] ; then ls > list.txt ; fi
##Check ages of lists
LIST_TIME=`expr $(date +%s) - $(date +%s -r list.txt)`
if [ "$LIST_TIME" -gt 86400 ] ; then ls *.qc.t > list.txt ; fi

#omit absent files
FILES1=`echo $QS | sed 's/_gene.cnt/.qc.t/g' | tr ' ' '\n' \
| sort - list.txt | uniq -d | tr '\n' ' '`

paste ../rownames.txt $FILES1 > $USRDIR/QC_Matrix.xls

# Collect reports
REPORTS=/var/www/metadata/reports_${ORG}/
cd $REPORTS

##Generate list if not existing
if [ ! -r list.log.txt ] ; then ls > list.log.txt ; fi
if [ ! -r list.sum.txt ] ; then ls > list.sum.txt ; fi

##Check ages of lists
LIST_TIME=`expr $(date +%s) - $(date +%s -r list.log.txt)`
if [ "$LIST_TIME" -gt 86400 ] ; then ls *log > list.log.txt ; fi
LIST_TIME=`expr $(date +%s) - $(date +%s -r list.sum.txt)`
if [ ! -r list.sum.txt ] ; then ls > list.sum.txt ; fi

##Omit absent files
#QS3=`echo $QS | sed 's/_gene.cnt//g'`
FILES2=`echo $QS | sed 's/_gene.cnt//g' | tr ' ' '\n' \
| sed 's/$/.log/' | sort - list.log.txt | uniq -d | tr '\n' ' '`
FILES3=`echo $QS |  sed 's/_gene.cnt//g' | tr ' ' '\n' \
| sed 's/$/_gene.cnt.summary/' | sort - list.sum.txt | uniq -d | tr '\n' ' '`

##Copy logs
cp $FILES2 $FILES3 $LOGS

cd $USRDIR
zip -r - CountMatrix.xls QC_Matrix.xls logs
cd ..
