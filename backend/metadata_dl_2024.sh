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

> $ORG.csv
#for YEAR in $(seq 2007 2022) ; do
for YEAR in 2024 ; do
  for DD in $(seq -w 01 31) ; do
    for MM in $(seq -w 01 12) ; do

      >tmp.csv
      pysradb search --organism="${ORGANISM}" \
        --publication-date ${DD}-${MM}-${YEAR}:${DD}-${MM}-${YEAR} \
        --source="transcriptomic" --max=999000 --query="Public[Access]" --saveto tmp.csv
      cat tmp.csv >> ../sradb/${ORG}_2024.csv
    done
  done
done
rm tmp.csv

#  for MONTH in $(seq -w 01 12) ; do
#    >tmp.csv
#    pysradb search --organism="${ORGANISM}" \
#      --publication-date 01-${MONTH}-${YEAR}:31-${MONTH}-${YEAR} \
#      --source="transcriptomic" --max=999000 --query="Public[Access]" --saveto tmp.csv
#    cat tmp.csv >> ../sradb/${ORG}_2024.csv
#  done
#done
#rm tmp.csv
