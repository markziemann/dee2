#!/bin/bash

ORG=$1

if [ -z $ORG ] ; then

  echo "Error: ORG not set!"
  exit

fi

while true ; do

  FILECNT=$(ls ${ORG}_*.sra | wc -l)


  if [ $FILECNT -lt 100 ] ; then

    SRR=$(curl "http://dee2.io/cgi-bin/acc.sh?ORG=${ORG}&submit" \
    | grep ACCESSION \
    | cut -d '=' -f2 )

    echo $SRR

    prefetch -X 9999999999999 -o ${ORG}_${SRR}.sra $SRR

    FILECNT=$(ls ${ORG}_*.sra | wc -l)

  else

    sleep 10

  fi

  ZIPCNT=$(ls *${ORG}.zip | grep -c zip$ )

  if [ $ZIPCNT -gt 0 ] ; then

    SFTP_URL=$(curl dee2.io/ip)

    for ZIP in $(ls *${ORG}.zip) ; do

      sftp -i ~/.ssh/guestuser guestuser@$SFTP_URL << EOF
put $ZIP
EOF

      rm -f $ZIP

    done
  fi
done

