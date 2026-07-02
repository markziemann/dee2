#!/bin/bash
source .venv/bin/activate

ORG=$1

if [ $ORG = "athaliana" ] ; then ORGANISM='Arabidopsis thaliana' ; fi
if [ $ORG = "celegans" ] ; then ORGANISM='Caenorhabditis elegans' ; fi
if [ $ORG = "dmelanogaster" ] ; then ORGANISM='Drosophila melanogaster' ; fi
if [ $ORG = "drerio" ] ; then ORGANISM='Danio rerio' ; fi
if [ $ORG = "ecoli" ] ; then ORGANISM='Escherichia coli' ; fi
if [ $ORG = "hsapiens" ] ; then ORGANISM='Homo sapiens' ; fi
if [ $ORG = "mmusculus" ] ; then ORGANISM='Mus musculus' ; fi
if [ $ORG = "scerevisiae" ] ; then ORGANISM='Saccharomyces cerevisiae' ; fi
if [ $ORG = "rnorvegicus" ] ; then ORGANISM='Rattus norvegicus' ; fi
if [ $ORG = "osativa" ] ; then ORGANISM='Oryza sativa' ; fi
if [ $ORG = "zmays" ] ; then ORGANISM='Zea mays' ; fi
if [ $ORG = "bdistachyon" ] ; then ORGANISM='Brachypodium distachyon' ; fi
if [ $ORG = "gmax" ] ; then ORGANISM='Glycine max' ; fi
if [ $ORG = "hvulgare" ] ; then ORGANISM='Hordeum vulgare' ; fi
if [ $ORG = "ptrichocarpa" ] ; then ORGANISM='Populus trichocarpa' ; fi
if [ $ORG = "sbicolor" ] ; then ORGANISM='Sorghum bicolor' ; fi
if [ $ORG = "slycopersicum" ] ; then ORGANISM='Solanum lycopersicum' ; fi
if [ $ORG = "stuberosum" ] ; then ORGANISM='Solanum tuberosum' ; fi
if [ $ORG = "taestivum" ] ; then ORGANISM='Triticum aestivum' ; fi
if [ $ORG = "vvinifera" ] ; then ORGANISM='Vitis vinifera' ; fi


# SRA metadata

if [ !  -d ../sradb/${ORG} ] ; then
  mkdir ../sradb/${ORG}
fi

start="2007-01-01"
end=$(date --date="2 days ago" +"%Y-%m-%d")

while [[ $start < $end ]] ;  do
  echo "$start"
  YEAR=$(echo $start | cut -d '-' -f1)
  MM=$(echo $start | cut -d '-' -f2)
  DD=$(echo $start | cut -d '-' -f3)
  OUTFILE=../sradb/${ORG}/$YEAR/${YEAR}-${MM}-${DD}.csv

  if [ ! -d ../sradb/${ORG}/$YEAR ] ; then
    mkdir ../sradb/${ORG}/$YEAR
  fi

  if [ ! -r ${OUTFILE} ] ; then
    sleep 3
    python -m pysradb.cli search --organism="${ORGANISM}" \
      --publication-date ${DD}-${MM}-${YEAR}:${DD}-${MM}-${YEAR} \
      --source="transcriptomic" --max=999000 --query="Public[Access]" --saveto $OUTFILE
  fi

  # set date
  start=$(date -d "$start + 1 day" +"%Y-%m-%d")
done

# GEO metadata

if [ !  -d ../sradb/${ORG}_geo ] ; then
  mkdir ../sradb/${ORG}_geo
fi

start="2007-01-01"
end=$(date --date="2 days ago" +"%Y-%m-%d")

while [[ $start < $end ]] ;  do
  echo "$start"
  YEAR=$(echo $start | cut -d '-' -f1)
  MM=$(echo $start | cut -d '-' -f2)
  DD=$(echo $start | cut -d '-' -f3)
  OUTFILE=../sradb/${ORG}_geo/$YEAR/${YEAR}-${MM}-${DD}.csv

  if [ ! -d ../sradb/${ORG}/$YEAR ] ; then
    mkdir ../sradb/${ORG}_geo/$YEAR
  fi

  if [ ! -r ${OUTFILE} ] ; then
    sleep 3
    python -m pysradb.cli search -d geo --organism="${ORGANISM}" \
      --publication-date ${DD}-${MM}-${YEAR}:${DD}-${MM}-${YEAR} \
      -C="TRANSCRIPTOMIC" --max=999000 --query="Public[Access]" --saveto $OUTFILE
  fi

  # set date
  start=$(date -d "$start + 1 day" +"%Y-%m-%d")
done
