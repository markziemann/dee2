# Validation of DEE2 data

In this folder are scripts related to the validation of DEE2 pipeline. There are two approaches to validation; simulation
and comparison to published datasets.

## Simulation study

The approach I took was very similar to the one used in my previous papers, particularly the one on short DNA seq mappers
(https://www.biorxiv.org/content/early/2016/05/16/053686). The general approach is to generate realistic RNA-seq reads 
using the complete Ensembl cDNA as the template and then process the simulated reads with DEE2 and then find out whether the 
RPM values from DEE2 match the simulated data.

The script *dee2_quant.sh* does most of the work by generating the simulated reads, processing them with DEE2 pipeline and 
tabulating the raw expression counts on a gene wise and transcript wise level. The *dee2_quent.R* script just calculates the 
RPM values, determines the Spearman correlation and plots the expected vs measured RPM values for each gene/transcript.

The simulated reads were generated with ART software (version art_bin_MountRainier) with a HiSeq2500 error profile and seed 
of 1540165885 (https://doi.org/10.1093/bioinformatics/btr708). The correct version of ART is a requirement that the script 
runs properly.

## Comparison study

In the comparison study, I searched for GEO series that had three or more replicates and summarised expression data on Ensembl
gene names. Data were fetched from DEE2 and GEO and the corresponding datasets underwent Spearman correlation analysis. In 
addition, contrasts were performed on DEE2 and GEO datasets and the results were merged side-by-side. The spearman correlation
of differential expression sign(logFC)/log(p-val) was calculated. 

The GEO series used in this study are GSE53078, GSE46344, GSE43180, GSE80768, GSE80251, GSE63776, GSE59970, GSE65715 and 
GSE76444.
