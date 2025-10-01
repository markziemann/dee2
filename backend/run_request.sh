#!/bin/bash

if [ -r LOCK ] ; then
  echo LOCK file present. Exiting now.
  exit
fi

# Start by `ls` then create the SRP folder, do analysis in there
ssh -i ~/.ssh/dee2_2025 ubuntu@dee2.io "ls /usr/lib/cgi-bin/newrequests/*confirmed" \
| sed 's#/usr/lib/cgi-bin/newrequests/##' > CONFIRMED

for FILE  in  $(cat CONFIRMED) ; do
  echo $FILE
  if [ ! -r $FILE ] ; then
    echo file $FILE not found
    FILE_PRESENT=FALSE

    #START ROUTINE
    touch LOCK
    SRP=$(echo $FILE | cut -d '.' -f1)
    scp -i ~/.ssh/dee2_2025 ubuntu@dee2.io:/usr/lib/cgi-bin/newrequests/$FILE .
    mkdir $SRP
    cp $FILE $SRP
    cd $SRP
    RUNS=$(grep -w $SRP $SRP | cut -f22)
    ORG=$(grep -w $SRP $SRP | head -1 | cut -f7)
    L1=$(echo $ORG | cut -c1 | tr '[:upper:]' '[:lower:]') ; echo $L1
    W2=$(echo $ORG | cut -d ' ' -f2)
    ORG2=$L1$W2

    # DOWNLOAD
    for SRR in $RUNS ; do
      prefetch -X 9999999999999 -o ${ORG2}_${SRR}.sra $SRR
    done

    # RUN
    ORG3=$echo $ORG2 | cut -c -3)
    ln -s ../tallyup .
    apptainer run -w -B ${PWD}:/dee2/mnt/ tallyup -s $ORG2 -t 8 -d
    cd ..
    rm LOCK
  fi
done
