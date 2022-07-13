#!/bin/bash

#run like this
#bash downloader_request.sh hsapiens sanjid.txt

ORG=$1
SRRLIST=$2
TXOUT=TxCountMatrix.tsv
GOUT=GeneCountMatrix.tsv

if [ -z $ORG ] ; then

  echo "Error: ORG not set!"
  exit

fi

for SRR in $(cat $SRRLIST) ; do

    echo $SRR

    prefetch -X 9999999999999 -o ${ORG}_${SRR}.sra $SRR

done

docker run -v $(pwd):/dee2/mnt mziemann/tallyup hsapiens -d

for ZIP in $(ls *${ORG}.zip) ; do

  echo $ZIP
  scp -i ~/.ssh/dee2 $ZIP ubuntu@dee2.io:~/upload

done

for ZIP in *zip ; do unzip $ZIP ; done

rm $(find . | grep tmp)

for KE in $(find . | grep ke) ; do
  cut -f4 $KE | sed 's/_est_counts//' > $KE.tmp
done

FILE1=$(find . | grep ke.tsv$ | head -1)
echo Tx_ID > $TXOUT
cut -f1 $FILE1 > $TXOUT
paste $TXOUT $(find . | grep ke.tsv.tmp) > tmp && mv tmp $TXOUT

for SE in $(find . | grep se) ; do
  cut -f2 $SE > $SE.tmp
done

FILE1=$(find . | grep se.tsv$ | head -1)
echo Gene_ID > $GOUT
cut -f1 $FILE1 > $GOUT
paste $GOUT $(find . | grep se.tsv.tmp) > tmp && mv tmp $GOUT

mkdir data
cp -r $SRRLIST SRR* GeneCountMatrix.tsv TxCountMatrix.tsv data
rm -rf $(find data | grep zip)
zip -r data.zip data
