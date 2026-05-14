#!/bin/bash

ORG="kphaffii"
TAX_ID="460519"

BASE_DIR="/home/user/x-alife/user/dee2"
LIST_FILE="${BASE_DIR}/${ORG}_accessions.txt"
RESULTS_DIR="${BASE_DIR}/results"
LOGS_DIR="${BASE_DIR}/logs"

# OPTIMIZED FOR n1-standard-96 (96 vCPUs, 360 GB RAM)
MAX_JOBS=24

mkdir -p "$RESULTS_DIR"
mkdir -p "$LOGS_DIR"
mkdir -p "${BASE_DIR}/ref"

# 1. Download the list of SRA runs (if missing or empty)
if [ ! -f "$LIST_FILE" ] || [ ! -s "$LIST_FILE" ]; then
    echo "Fetching transcriptomic SRA accessions for $ORG from EBI ENA API..."
    curl -s "https://www.ebi.ac.uk/ena/portal/api/search?query=tax_tree(${TAX_ID})&result=read_run&fields=run_accession,library_source&format=tsv&limit=0" \
        | awk -F'\t' '$2 ~ /TRANSCRIPTOMIC/ {print $1}' > "$LIST_FILE"
    echo "Found $(wc -l < "$LIST_FILE") RNA-Seq runs to process."
fi

# 2. Define the core processing worker function
process_srr() {
    local SRR="$1"
    local ORG="$2"
    local RESULTS_DIR="$3"
    local BASE_DIR="$4"
    local LOGS_DIR="$5"

    # Safely skip completely processed, non-empty files
    if [ -s "$RESULTS_DIR/${SRR}.${ORG}.zip" ]; then
        echo "Skipping $SRR (already successfully processed)."
        return 0
    fi

    # Staggered start up to 30 seconds to protect network I/O from getting throttled
    local DELAY=$(( (RANDOM % 30) + 1 ))
    sleep $DELAY

    echo "=== Starting background job for $SRR (stagger delay: ${DELAY}s) ==="

    # Route all internal container logs to individual files inside the logs/ folder
    docker run --name "dee2_${SRR}" \
        -v "${BASE_DIR}/pipeline/volunteer_pipeline.sh:/dee2/code/volunteer_pipeline.sh" \
        -v "${BASE_DIR}/ref:/dee2/ref" \
        mziemann/tallyup -s "$ORG" -a "$SRR" > "${LOGS_DIR}/${SRR}.log" 2>&1

    # Extract the resulting zip file back to the host machine
    docker cp "dee2_${SRR}:/dee2/${SRR}.${ORG}.zip" "$RESULTS_DIR/" 2>/dev/null
    
    if [ $? -eq 0 ]; then
        echo "Successfully saved $SRR"
        # Delete successful logs to save storage space, keeping failed ones for troubleshooting
        rm -f "${LOGS_DIR}/${SRR}.log"
    else
        echo "Error processing $SRR (Check logs/${SRR}.log for details)."
    fi

    # Cleanup container storage layer
    docker rm -v "dee2_${SRR}" >/dev/null 2>&1
}

# Export environment variables so the sub-shells spawned by xargs can access them
export -f process_srr
export ORG RESULTS_DIR BASE_DIR LOGS_DIR

echo "------------------------------------------------"
echo "Launching parallel bulk processing queue..."
echo "Concurrency limit: $MAX_JOBS workers | Total runs: $(wc -l < "$LIST_FILE")"
echo "------------------------------------------------"

# Wipe out any dead or lingering 'dee2_' containers to prevent name collisions
echo "Cleaning up any old, lingering containers..."
docker rm -f $(docker ps -a -q --filter name=dee2_) >/dev/null 2>&1

# 3. Stream the accessions into xargs to scale out across your CPU cores
cat "$LIST_FILE" | xargs -P "$MAX_JOBS" -I {} bash -c 'process_srr "{}" "$ORG" "$RESULTS_DIR" "$BASE_DIR" "$LOGS_DIR"'

echo "Bulk processing complete! Parallel execution queue has finished processing."
