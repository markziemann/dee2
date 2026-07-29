#!/bin/bash

set -e
set -x

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

PAUSE=3
BATCH_SIZE=5000
API_KEY=$(cat /mnt/hdd1/dee2/sradb/API_KEY)
ESEARCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH_DIR="sra_batches_${DATE_SAFE}"

# SRA metadata

if [ !  -d ../sradb_new/${ORG} ] ; then
  mkdir -p ../sradb_new/${ORG}
fi

start="2026-01-01"
end=$(date --date="2 days ago" +"%Y-%m-%d")

while [[ $start < $end ]] ;  do
  echo "$start"
  YEAR=$(echo $start | cut -d '-' -f1)
  MM=$(echo $start | cut -d '-' -f2)
  DD=$(echo $start | cut -d '-' -f3)
  DATE_SAFE=${YEAR}-${MM}-${DD}
  OUTFILE=../sradb_new/${ORG}/$YEAR/${DATE_SAFE}.csv

  if [ ! -d ../sradb_new/${ORG}/$YEAR ] ; then
    mkdir ../sradb_new/${ORG}/$YEAR
  fi

  if [ ! -r ${OUTFILE} ] ; then
    sleep $PAUSE
#    python -m pysradb.cli search --organism="${ORGANISM}" \
#      --publication-date ${DD}-${MM}-${YEAR}:${DD}-${MM}-${YEAR} \
#      --source="transcriptomic" --max=999000 --query="Public[Access]" --saveto $OUTFILE

    TERM="${ORGANISM}[Organism] AND transcriptomic[Source] AND ${YEAR}/${MM}/${DD}[PDAT] : ${YEAR}/${MM}/${DD}[PDAT]"
    SEARCH_XML=$(curl -s -G "$ESEARCH_URL" \
      --data-urlencode "db=sra" \
      --data-urlencode "term=${TERM}" \
      --data-urlencode "usehistory=y" \
      --data-urlencode "retmax=0" \
      ${API_KEY:+--data-urlencode "api_key=${API_KEY}"} )


    WEBENV=$(grep -oP '(?<=<WebEnv>)[^<]+' <<< "$SEARCH_XML")
    QUERY_KEY=$(grep -oP '(?<=<QueryKey>)[^<]+' <<< "$SEARCH_XML")
    COUNT=$(grep -oP '(?<=<Count>)[^<]+' <<< "$SEARCH_XML")
    COUNT=$(echo $COUNT | cut -d ' ' -f1)

    BATCH_DIR="sra_batches_${DATE_SAFE}"

    if [[ "$COUNT" -eq 0 ]]; then
      touch $OUTFILE
    else
      mkdir $BATCH_DIR
      BATCH_FILES=()
      for ((RETSTART=0; RETSTART<COUNT; RETSTART+=BATCH_SIZE)); do
        BATCH_FILE="${BATCH_DIR}/batch_${RETSTART}.xml"

        echo "Fetching records ${RETSTART}-$((RETSTART+BATCH_SIZE-1)) of ${COUNT}..." >&2

        curl -s -G "$EFETCH_URL" \
          --data-urlencode "db=sra" \
          --data-urlencode "query_key=${QUERY_KEY}" \
          --data-urlencode "WebEnv=${WEBENV}" \
          --data-urlencode "retstart=${RETSTART}" \
          --data-urlencode "retmax=${BATCH_SIZE}" \
          --data-urlencode "rettype=full" \
          --data-urlencode "retmode=xml" \
          -o "$BATCH_FILE" \
          ${API_KEY:+--data-urlencode "api_key=${API_KEY}"}

          BATCH_FILES+=("$BATCH_FILE")
          sleep $PAUSE
      done

      echo "Parsing ${#BATCH_FILES[@]} batch file(s) into CSV..." >&2
      python3 parse_sra_xml.py "${BATCH_FILES[@]}" "$OUTFILE"
      rm -rf $BATCH_DIR
    fi

  fi

  # set date
  start=$(date -d "$start + 1 day" +"%Y-%m-%d")
done
