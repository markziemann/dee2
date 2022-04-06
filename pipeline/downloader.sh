#!/bin/bash

ORG=$1

if [ -z $ORG ] ; then

  echo "Error: ORG not set!"
  exit

fi

while true ; do

  FILECNT=$(ls ${ORG}_*.sra | wc -l)

  if [ $FILECNT -lt 100 ] ; then

    SRR=$(curl "https://dee2.io/cgi-bin/acc.sh?ORG=${ORG}&submit" \
    | grep ACCESSION \
    | cut -d '=' -f2 )

    echo $SRR

    prefetch -X 9999999999999 -o ${ORG}_${SRR}.sra $SRR

    FILECNT=$(ls ${ORG}_*.sra | wc -l)

  else

    date
    sleep 10

  fi

  ZIPCNT=$(ls *${ORG}.zip | grep -c zip$ )

  if [ $ZIPCNT -gt 0 ] ; then

    for ZIP in $(ls *${ORG}.zip) ; do

      SIZE=$(du -s $ZIP | awk '$1')

      if [ $SIZE -gt 5000 ] ; then
        rm $ZIP
      else
        scp -i ~/.ssh/dee2 $ZIP ubuntu@dee2.io:~/upload && rm -f $ZIP
      fi
    done
  fi
done
