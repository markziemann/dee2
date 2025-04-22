#!/bin/bash

#PBS -N mmu2
#PBS -l ncpus=12
#PBS -l mem=50GB
#PBS -l jobfs=300GB
#PBS -q normal
#PBS -P rv38
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/rv38+scratch/rv38
#PBS -l wd

set -x

JOBCNT=$(qstat | grep -c mmu2)

if [ $JOBCNT -eq 1 ] ; then

  cd ~/dee2/nocontainer/mmu2/mmu2

  FILECNT=$(find dee2/mnt/ | grep -c sra$ )

  if [ $FILECNT -ge 10 ] ; then

    bash dee2/code/volunteer_pipeline.sh -s mmusculus -t 12 -d -v

  fi

  FILECNT2=$(find ../tfr/ -cmin +10 | grep -c sra$ )

  if [ $FILECNT2 -ge 10 ] ; then

    mv $(find ../tfr/ -cmin +10 | grep sra$ | head -10 ) dl

    bash dee2/code/volunteer_pipeline.sh -s mmusculus -t 12 -d -v

  fi

fi
