#!/bin/bash

set -x

ORG=$1

if [ -z $ORG ] ; then

  echo "Error: ORG not set!"
  exit

fi

while true ; do

  FILECNT=$(ssh mz3382@gadi.nci.org.au -t 'find ~/dee2/nocontainer/mmu2/tfr/ | grep -c sra$ ')

  FILECNT=$(echo $FILECNT | dos2unix)

  JOBCNT=$(ssh mz3382@gadi.nci.org.au -t 'qstat | grep -c mmu2' )

  JOBCNT=$(echo $JOBCNT | dos2unix)

  if [ $JOBCNT -eq 0 ] ; then

    if [ $FILECNT -ge 20 ] ; then

      ssh mz3382@gadi.nci.org.au -t 'qsub /home/565/mz3382/dee2/nocontainer/mmu2/mmu2/dee2_pbs.sh'

    fi

  fi

  if [ $FILECNT -lt 20 ] ; then

    SRR=$(curl "https://dee2.io/cgi-bin/acc.sh?ORG=${ORG}&submit" \
    | grep ACCESSION \
    | cut -d '=' -f2 )

    echo $SRR

    prefetch -X 9999999999999 -o ${ORG}_${SRR}.sra $SRR \
    && scp ${ORG}_${SRR}.sra mz3382@gadi-dm.nci.org.au:~/dee2/nocontainer/mmu2/tfr \
    && rm ${ORG}_${SRR}.sra

  else

    date
    sleep 60

  fi

  sleep 10

done
