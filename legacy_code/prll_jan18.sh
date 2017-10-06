#!/bin/bash
#sra2mx by mark ziemann mark.ziemann@gmail.com
#2015-01-18

#Specify config file, pipeline script, genome, parallel pipelines
CFG=$1
PIPELINE=`grep ^PIPELINE= $CFG | cut -d '=' -f2`
URLS=`grep ^URLS= $CFG | cut -d '=' -f2`
GENOME=`grep ^GENOME= $CFG | cut -d '=' -f2`
PARALLEL=`grep ^PARALLEL= $CFG | cut -d '=' -f2`
DISKLIM=`grep ^DISKLIM= $CFG | cut -d '=' -f2`
DISK=`df . | awk 'END{print$4}'`
#Calc mem ensure we don't run too many alignments and max out the ram
MEM=`free | awk '/buffers\/cache/{print $4}'`
SAINDEX=`du -c $GENOME | tail -1 | cut -f1`
ALNLIM=4
#`grep ^ALNLIM= $CFG | cut -d '=' -f2`
MEMALNLIM=4
#`echo $MEM $SAINDEX | awk '{printf "%d\n", ($1/$2)-1}'`
grep -v MEMALNLIM $CFG > $CFG.tmp ; mv $CFG.tmp $CFG
echo MEMALNLIM=$MEMALNLIM >> $CFG

echo "Starting sra2mx pipeline
Parameters:
  Pipeline script = $PIPELINE
  Config file = $CFG
  URLS file = $URLS
  minimum disk space = $DISKLIM
  current disk space = $DISK
  free memory now = $MEM
  genome index size = $SAINDEX
  max concurrent alignments specified = $ALNLIM
  max allowed alignments with free mem = $MEMALNLIM "

#Load genome
STAR --genomeLoad LoadAndExit --genomeDir $GENOME/

#Run it
parallel -u -j $PARALLEL $PIPELINE $CFG < $URLS
#parallel -u -j $PARALLEL -S omics3,omics4,omics5 $PIPELINE $CFG < $URLS

#unload genome from memory
STAR --genomeLoad Remove --genomeDir $GENOME/



