#!/bin/bash

set -x

cd "$(dirname "$0")";

REQ=/var/www/request.txt

while read line ; do

  SRP=$(echo $line | awk '{print $2}')

  ADDRESS=$(echo $line | awk '{print $3}')

  FILENAME=/dee2_data/requests/$SRP.zip

  ERR=/dee2_data/requests/$SRP.txt

  if [ -r $FILENAME ] ; then

    AGE=$(( `date +%s` - `stat -L --format %Y $FILENAME` ))

    if [ $AGE -gt 60 ] ; then

      if [ ! -z $ADDRESS ] ; then

        sed "s/SRP_ACCESSION/$SRP/" mail1_template.txt > mail.txt

        mail -s "$SRP is ready at dee2.io !" \
        -r "mdz@dee2.io" "$ADDRESS" < mail.txt

        grep -wv $SRP $REQ > tmp && mv tmp $REQ
        sudo chown www-data:www-data $REQ

      fi

    fi

  elif [ -r $ERR ] ; then

    if [ ! -z $ADDRESS ] ; then

      sed "s/SRP_ACCESSION/$SRP/" mail2_template.txt > mail.txt

      mail -s "Problems with dataset ${SRP}" \
      -r "mdz@dee2.io" "$ADDRESS" < mail.txt

      grep -wv $SRP $REQ > tmp && mv tmp $REQ
      sudo chown www-data:www-data $REQ

    fi

  else

    echo $SRP does not exist

  fi

done < $REQ

