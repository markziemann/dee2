#!/bin/bash

ORG=$1

while true ; do

  SRR=$(curl "http://dee2.io/cgi-bin/acc.sh?ORG=${ORG}&submit" \
  | grep ACCESSION \
  | cut -d '=' -f2 )

  echo $SRR

  prefetch -X 9999999999999 -o ${ORG}_${SRR}.sra $SRR

done
