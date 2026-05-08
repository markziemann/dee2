#!/bin/bash

RESULTS_DIR="results"
AGG_DIR="aggregated_data"

if [ ! -d "$RESULTS_DIR" ] || [ -z "$(ls -A $RESULTS_DIR/*.zip 2>/dev/null)" ]; then
    echo "Error: No .zip files found in $RESULTS_DIR/. Run the bulk processor first!"
    exit 1
fi

mkdir -p "$AGG_DIR"
TMP_DIR=$(mktemp -d)

echo "Unzipping results..."
for ZIP in $RESULTS_DIR/*.zip; do
    # -j ignores directories inside the zip, dumping files flat
    # -q is quiet mode
    unzip -q -j "$ZIP" -d "$TMP_DIR"
done

cd "$TMP_DIR"

#####################################################
# 1. Aggregate Gene Counts (from STAR .se.tsv files)
#####################################################
echo "Aggregating Gene Counts..."
SE_FILES=(*.se.tsv)
if [ -e "${SE_FILES[0]}" ]; then
    # Grab the gene names from the first file. 
    # STAR's ReadsPerGene.out.tab (se.tsv) has 4 lines of metadata at the top, so we use tail -n +5
    echo "Gene_ID" > "../$AGG_DIR/GeneCountMatrix.tsv"
    tail -n +5 "${SE_FILES[0]}" | cut -f1 >> "../$AGG_DIR/GeneCountMatrix.tsv"

    for SE in *.se.tsv; do
        SRR=$(basename "$SE" .se.tsv)
        # Add the SRR header, then grab column 2 (unstranded counts)
        echo "$SRR" > "$SE.tmp"
        tail -n +5 "$SE" | cut -f2 >> "$SE.tmp"
    done
    
    # Paste the Gene IDs and all the SRR count columns together
    paste "../$AGG_DIR/GeneCountMatrix.tsv" *.se.tsv.tmp > "../$AGG_DIR/tmp_gene"
    mv "../$AGG_DIR/tmp_gene" "../$AGG_DIR/GeneCountMatrix.tsv"
fi

#####################################################
# 2. Aggregate Transcript Counts (from Kallisto .ke.tsv files)
#####################################################
echo "Aggregating Transcript Counts..."
KE_FILES=(*.ke.tsv)
if [ -e "${KE_FILES[0]}" ]; then
    # Grab the transcript IDs from the first file.
    # Kallisto's abundance.tsv (ke.tsv) has 1 header line, so we use tail -n +2
    echo "Tx_ID" > "../$AGG_DIR/TxCountMatrix.tsv"
    tail -n +2 "${KE_FILES[0]}" | cut -f1 >> "../$AGG_DIR/TxCountMatrix.tsv"

    for KE in *.ke.tsv; do
        SRR=$(basename "$KE" .ke.tsv)
        # Add the SRR header, then grab column 4 (estimated counts)
        echo "$SRR" > "$KE.tmp"
        tail -n +2 "$KE" | cut -f4 >> "$KE.tmp"
    done
    
    # Paste Transcript IDs and all the SRR count columns together
    paste "../$AGG_DIR/TxCountMatrix.tsv" *.ke.tsv.tmp > "../$AGG_DIR/tmp_tx"
    mv "../$AGG_DIR/tmp_tx" "../$AGG_DIR/TxCountMatrix.tsv"
fi

#####################################################
# 3. Aggregate QC Metrics (from .log files)
#####################################################
echo "Aggregating QC Metadata..."
LOG_FILES=(*.log)
if [ -e "${LOG_FILES[0]}" ]; then
    # The DEE2 pipeline outputs a key:value summary block at the end of each .log file.
    # We will extract the Keys from the first file to form our header.
    grep -E '^[A-Za-z0-9_]+:' "${LOG_FILES[0]}" | cut -d ':' -f1 | awk 'BEGIN{printf "SRR_Accession\t"} {printf "%s\t", $0} END{print ""}' > "../$AGG_DIR/QC_Matrix.tsv"

    for LOG in *.log; do
        SRR=$(basename "$LOG" .log)
        # Extract the values matching the keys
        grep -E '^[A-Za-z0-9_]+:' "$LOG" | cut -d ':' -f2 | awk -v srr="$SRR" 'BEGIN{printf "%s\t", srr} {printf "%s\t", $0} END{print ""}' >> "../$AGG_DIR/QC_Matrix.tsv"
    done
fi

cd ..
rm -rf "$TMP_DIR"

echo "------------------------------------------------"
echo "Aggregation complete! Final DEE2 matrices saved in the $AGG_DIR/ directory:"
ls -lh $AGG_DIR/
echo "------------------------------------------------"
