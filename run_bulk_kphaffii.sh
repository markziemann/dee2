#!/bin/bash

ORG="kphaffii"
TAX_ID="460519"
LIST_FILE="${ORG}_accessions.txt"
RESULTS_DIR="results"

mkdir -p "$RESULTS_DIR"
mkdir -p ref

# 1. Download the list of SRA runs
# We use the EBI ENA API (which mirrors NCBI SRA perfectly) because it provides a much faster and 
# cleaner text-based API for fetching accessions compared to NCBI E-utilities.
if [ ! -f "$LIST_FILE" ]; then
    echo "Fetching transcriptomic SRA accessions for $ORG from EBI ENA API..."
    
    # Query for all runs matching the Taxon ID, returning the Run Accession and Library Source
    # We filter only for "TRANSCRIPTOMIC" (RNA-Seq) data, matching DEE2's standard behavior.
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=tax_eq(${TAX_ID})&result=read_run&fields=run_accession,library_source&format=tsv" \
        | awk '$2=="TRANSCRIPTOMIC" {print $1}' > "$LIST_FILE"
        
    echo "Found $(wc -l < "$LIST_FILE") RNA-Seq runs to process."
fi

# 2. Process each accession
for SRR in $(cat "$LIST_FILE"); do
    # Skip if we already have the completed zip file (makes the script resumable)
    if [ -f "$RESULTS_DIR/${SRR}.${ORG}.zip" ]; then
        echo "Skipping $SRR (already processed)."
        continue
    fi

    echo "------------------------------------------------"
    echo "Processing $SRR..."
    echo "------------------------------------------------"

    # We name the container so we can reliably extract data from it before deleting it.
    # Note: No --rm flag, otherwise it deletes itself before we can copy the ZIP file out!
    docker run --name "dee2_${SRR}" \
        -v $(pwd)/pipeline/volunteer_pipeline.sh:/dee2/code/volunteer_pipeline.sh \
        -v $(pwd)/ref:/dee2/ref \
        mziemann/tallyup -s $ORG -a $SRR

    # Retrieve the processed ZIP file from the container's working directory (/dee2)
    docker cp "dee2_${SRR}:/dee2/${SRR}.${ORG}.zip" "$RESULTS_DIR/" 2>/dev/null
    
    if [ $? -eq 0 ]; then
        echo "Successfully saved $RESULTS_DIR/${SRR}.${ORG}.zip"
    else
        echo "Error processing $SRR. It may have failed QC or SRA download."
    fi

    # Cleanup the finished container to free up Docker storage space
    docker rm -v "dee2_${SRR}" >/dev/null
done

echo "Bulk processing complete! All processed data is in the $RESULTS_DIR/ directory."
