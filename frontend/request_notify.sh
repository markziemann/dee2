#!/bin/bash

cd "$(dirname "$0")";

REQ=/var/www/html/request.txt

while read line ; do

  SRP=$(echo $line | awk '{print $2}')

  ADDRESS=(echo $line | awk '{print $2}')

  FILENAME=/dee2_data/requests/$SRP.zip

  if [ -r $FILENAME ] ; then

    AGE=$(( `date +%s` - `stat -L --format %Y $FILENAME` ))

    if [ $AGE -gt 60 ] ; then

      if [ ! -z $ADDRESS ] ; then

        sed "s/SRP_ACCESSION/$SRP/" mail1_template.txt > mail.txt

        mail -s "$SRP is ready at dee2.io !" \
        -r "mdz@dee2.io" "$ADDRESS" < mail.txt

        grep -wv $SRP $REQ > tmp && mv tmp $REQ

      if

    fi

  else

    echo $SRP does not exist

  fi

done < $REQ

