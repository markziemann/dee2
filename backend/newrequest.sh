#!/bin/bash

# This script is designed to accommodat requests from users who nominate
# specific accessions to be processed

TODO=todo.txt

PROG=inprogress

# check if already working
if [ -r $PROG ] ; then exit ; fi

scp ubuntu@118.138.234.94:/var/www/html/request.txt .

comm -23 <(sort request.txt) <(sort oldrequest.txt) > $TODO

len=$(wc -l < $TODO)

if [ $len -gt 0 ] ; then

  while IFS= read -r line; do echo $line

    ORG=$(cat $TODO | cut -d ' ' -f1)
    SRP=$(cat $TODO | cut -d ' ' -f2)
    EMAIL=$(cat $TODO | cut -d ' ' -f3)
    ACCS=$(grep -w $SRP ../sradb/${ORG}.csv | cut -d ',' -f1 | paste -s -d ',')
    len2=$(echo $ACCS | wc -c)

    if [ $len2 -gt 0 ] ; then

      > $PROG
      docker run mziemann/tallyup $ORG $ACCS
      rm $PROG

    fi

  done < $TODO

fi

cp request.txt oldrequest.txt
