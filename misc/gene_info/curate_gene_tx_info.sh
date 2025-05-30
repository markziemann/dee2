#!/bin/bash

#which gtftools > /dev/null || echo "GTFtools_0.6.5 is required"
#which gtftools > /dev/null || exit 1
#INSTALLED_VERSION=$(gtftools -v 2>&1 | tr -d ' ' )
#REQUIRED_VERSION="GTFtoolsversion:0.8.5"
#if [ $REQUIRED_VERSION != $INSTALLED_VERSION ] ; then echo "GTFtools_0.8.5 is required" ; exit 1 ; fi

ATH_GTFURL="ftp://ftp.ensemblgenomes.org/pub/release-36/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.36.gtf.gz"
ATH_CDNAURL="ftp://ftp.ensemblgenomes.org/pub/release-36/plants/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
ATH_GTF="Arabidopsis_thaliana.TAIR10.36.gtf"
ATH_CDNA="Arabidopsis_thaliana.TAIR10.cdna.all.fa"
ATH_GENEINFO="ath_gene_info.tsv"
ATH_TXINFO="ath_tx_info.tsv"

CEL_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.90.gtf.gz"
CEL_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz"
CEL_GTF="Caenorhabditis_elegans.WBcel235.90.gtf"
CEL_CDNA="Caenorhabditis_elegans.WBcel235.cdna.all.fa"
CEL_GENEINFO="cel_gene_info.tsv"
CEL_TXINFO="cel_tx_info.tsv"

DME_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.90.gtf.gz"
DME_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz"
DME_GTF="Drosophila_melanogaster.BDGP6.90.gtf"
DME_CDNA="Drosophila_melanogaster.BDGP6.cdna.all.fa"
DME_GENEINFO="dme_gene_info.tsv"
DME_TXINFO="dme_tx_info.tsv"

DRE_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/danio_rerio/Danio_rerio.GRCz10.90.gtf.gz"
DRE_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/danio_rerio/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz"
DRE_GTF="Danio_rerio.GRCz10.90.gtf"
DRE_CDNA="Danio_rerio.GRCz10.cdna.all.fa"
DRE_GENEINFO="dre_gene_info.tsv"
DRE_TXINFO="dre_tx_info.tsv"

ECO_GTFURL="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.36.gtf.gz"
ECO_CDNAURL="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cdna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa.gz"
ECO_GTF="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.36.gtf"
ECO_CDNA="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa"
ECO_GENEINFO="eco_gene_info.tsv"
ECO_TXINFO="eco_tx_info.tsv"

HSA_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz"
HSA_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
HSA_GTF="Homo_sapiens.GRCh38.90.gtf"
HSA_CDNA="Homo_sapiens.GRCh38.cdna.all.fa"
HSA_GENEINFO="hsa_gene_info.tsv"
HSA_TXINFO="hsa_tx_info.tsv"

MMU_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz"
MMU_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
MMU_GTF="Mus_musculus.GRCm38.90.gtf"
MMU_CDNA="Mus_musculus.GRCm38.cdna.all.fa"
MMU_GENEINFO="mmu_gene_info.tsv"
MMU_TXINFO="mmu_tx_info.tsv"

RNO_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf.gz"
RNO_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz"
RNO_GTF="Rattus_norvegicus.Rnor_6.0.90.gtf"
RNO_CDNA="Rattus_norvegicus.Rnor_6.0.cdna.all.fa"
RNO_GENEINFO="rno_gene_info.tsv"
RNO_TXINFO="rno_tx_info.tsv"

SCE_GTFURL="ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.36.gtf.gz"
SCE_CDNAURL="ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
SCE_GTF="Saccharomyces_cerevisiae.R64-1-1.36.gtf"
SCE_CDNA="Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa"
SCE_GENEINFO="sce_gene_info.tsv"
SCE_TXINFO="sce_tx_info.tsv"

OSA_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.59.gtf.gz"
OSA_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/oryza_sativa/cdna/Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz"
OSA_GTF="Oryza_sativa.IRGSP-1.0.59.gtf"
OSA_CDNA="Oryza_sativa.IRGSP-1.0.cdna.all.fa"
OSA_GENEINFO="osa_gene_info.tsv"
OSA_TXINFO="osa_tx_info.tsv"

ZMA_GTFURL="https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-59/plants/gtf/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.gtf.gz"
ZMA_CDNAURL="https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-59/plants/fasta/zea_mays/cdna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz"
ZMA_GTF="Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.gtf"
ZMA_CDNA="Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa"
ZMA_GENEINFO="zma_gene_info.tsv"
ZMA_TXINFO="zma_tx_info.tsv"

TAE_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/triticum_aestivum/Triticum_aestivum.IWGSC.59.gtf.gz"
TAE_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/triticum_aestivum/cdna/Triticum_aestivum.IWGSC.cdna.all.fa.gz"
TAE_GTF="Triticum_aestivum.IWGSC.59.gtf"
TAE_CDNA="Triticum_aestivum.IWGSC.cdna.all.fa"
TAE_GENEINFO="tae_gene_info.tsv"
TAE_TXINFO="tae_tx_info.tsv"

SLY_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/solanum_lycopersicum/Solanum_lycopersicum.SL3.0.59.gtf.gz"
SLY_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/solanum_lycopersicum/cdna/Solanum_lycopersicum.SL3.0.cdna.all.fa.gz"
SLY_GTF="Solanum_lycopersicum.SL3.0.59.gtf"
SLY_CDNA="Solanum_lycopersicum.SL3.0.cdna.all.fa"
SLY_GENEINFO="sly_gene_info.tsv"
SLY_TXINFO="sly_tx_info.tsv"

SBI_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.59.gtf.gz"
SBI_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/sorghum_bicolor/cdna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cdna.all.fa.gz"
SBI_GTF="Sorghum_bicolor.Sorghum_bicolor_NCBIv3.59.gtf"
SBI_CDNA="Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cdna.all.fa"
SBI_GENEINFO="sbi_gene_info.tsv"
SBI_TXINFO="sbi_tx_info.tsv"

GMA_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/glycine_max/Glycine_max.Glycine_max_v2.1.59.gtf.gz"
GMA_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/glycine_max/cdna/Glycine_max.Glycine_max_v2.1.cdna.all.fa.gz"
GMA_GTF="Glycine_max.Glycine_max_v2.1.59.gtf"
GMA_CDNA="Glycine_max.Glycine_max_v2.1.cdna.all.fa"
GMA_GENEINFO="gma_gene_info.tsv"
GMA_TXINFO="gma_tx_info.tsv"

PTR_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/populus_trichocarpa/Populus_trichocarpa.Pop_tri_v4.59.gtf.gz"
PTR_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/populus_trichocarpa/cdna/Populus_trichocarpa.Pop_tri_v4.cdna.all.fa.gz"
PTR_GTF="Populus_trichocarpa.Pop_tri_v4.59.gtf"
PTR_CDNA="Populus_trichocarpa.Pop_tri_v4.cdna.all.fa"
PTR_GENEINFO="ptr_gene_info.tsv"
PTR_TXINFO="ptr_tx_info.tsv"

VVI_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/vitis_vinifera/Vitis_vinifera.PN40024.v4.59.gtf.gz"
VVI_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/vitis_vinifera/cdna/Vitis_vinifera.PN40024.v4.cdna.all.fa.gz"
VVI_GTF="Vitis_vinifera.PN40024.v4.59.gtf"
VVI_CDNA="Vitis_vinifera.PN40024.v4.cdna.all.fa"
VVI_GENEINFO="vvi_gene_info.tsv"
VVI_TXINFO="vvi_tx_info.tsv"

HVU_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/hordeum_vulgare/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.59.gtf.gz"
HVU_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/hordeum_vulgare/cdna/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.cdna.all.fa.gz"
HVU_GTF="Hordeum_vulgare.MorexV3_pseudomolecules_assembly.59.gtf"
HVU_CDNA="Hordeum_vulgare.MorexV3_pseudomolecules_assembly.cdna.all.fa"
HVU_GENEINFO="hvu_gene_info.tsv"
HVU_TXINFO="hvu_tx_info.tsv"

STU_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/solanum_tuberosum/Solanum_tuberosum.SolTub_3.0.59.gtf.gz"
STU_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/solanum_tuberosum/cdna/Solanum_tuberosum.SolTub_3.0.cdna.all.fa.gz"
STU_GTF="Solanum_tuberosum.SolTub_3.0.59.gtf"
STU_CDNA="Solanum_tuberosum.SolTub_3.0.cdna.all.fa"
STU_GENEINFO="stu_gene_info.tsv"
STU_TXINFO="stu_tx_info.tsv"

BDI_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/brachypodium_distachyon/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.59.gtf.gz"
BDI_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/brachypodium_distachyon/cdna/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cdna.all.fa.gz"
BDI_GTF="Brachypodium_distachyon.Brachypodium_distachyon_v3.0.59.gtf"
BDI_CDNA="Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cdna.all.fa"
BDI_GENEINFO="bdi_gene_info.tsv"
BDI_TXINFO="bdi_tx_info.tsv"

###########################################################
# ATH
###########################################################
wget -N $ATH_GTFURL && gunzip -kf $ATH_GTF.gz
wget -N $ATH_CDNAURL && gunzip -kf $ATH_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $ATH_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $ATH_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $ATH_TXINFO

# prep the gene lengths
grep '#' $ATH_GTF > $ATH_GTF.tmp
grep -v '#' $ATH_GTF | awk '{OFS="\t"} $1=1' >> $ATH_GTF.tmp
gtftools  -l $ATH_GTF.genelength $ATH_GTF.tmp
rm $ATH_GTF.tmp

# prep gene names
grep -w gene $ATH_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $ATH_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $ATH_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $ATH_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $ATH_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $ATH_GENEINFO

###########################################################
# CEL
###########################################################
wget -N $CEL_GTFURL && gunzip -kf $CEL_GTF.gz
wget -N $CEL_CDNAURL && gunzip -kf $CEL_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $CEL_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $CEL_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $CEL_TXINFO

# prep the gene lengths
grep '#' $CEL_GTF > $CEL_GTF.tmp
grep -v '#' $CEL_GTF | awk '{OFS="\t"} $1=1' >> $CEL_GTF.tmp
gtftools  -l $CEL_GTF.genelength  $CEL_GTF.tmp
rm $CEL_GTF.tmp

# prep gene names
grep -w gene $CEL_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $CEL_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $CEL_GENEINFO
awk '{print $0,NR}' $CEL_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $CEL_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $CEL_GENEINFO


###########################################################
# DME
###########################################################
wget -N $DME_GTFURL && gunzip -kf $DME_GTF.gz
wget -N $DME_CDNAURL && gunzip -kf $DME_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $DME_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $DME_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $DME_TXINFO

# prep the gene lengths
grep '#' $DME_GTF > $DME_GTF.tmp
grep -v '#' $DME_GTF | awk '{OFS="\t"} $1=1' >> $DME_GTF.tmp
gtftools  -l $DME_GTF.genelength  $DME_GTF.tmp
rm $DME_GTF.tmp

# prep gene names
grep -w gene $DME_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $DME_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $DME_GENEINFO
awk '{print $0,NR}' $DME_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $DME_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $DME_GENEINFO


###########################################################
# DRE
###########################################################
wget -N $DRE_GTFURL && gunzip -kf $DRE_GTF.gz
wget -N $DRE_CDNAURL && gunzip -kf $DRE_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $DRE_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $DRE_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $DRE_TXINFO

# prep the gene lengths
grep '#' $DRE_GTF > $DRE_GTF.tmp
grep -v '#' $DRE_GTF | awk '{OFS="\t"} $1=1' >> $DRE_GTF.tmp
gtftools  -l $DRE_GTF.genelength  $DRE_GTF.tmp
rm $DRE_GTF.tmp

# prep gene names
grep -w gene $DRE_GTF | cut -d '"' -f2,6 | tr '"' '\t' | cut -d ' ' -f1  > $DRE_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $DRE_GENEINFO
awk '{print $0,NR}' $DRE_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $DRE_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $DRE_GENEINFO


###########################################################
# ECO
###########################################################
wget -N $ECO_GTFURL && gunzip -kf $ECO_GTF.gz
wget -N $ECO_CDNAURL && gunzip -kf $ECO_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $ECO_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $ECO_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $ECO_TXINFO

# prep the gene lengths
grep '#' $ECO_GTF > $ECO_GTF.tmp
grep -v '#' $ECO_GTF | awk '{OFS="\t"} $1=1' >> $ECO_GTF.tmp
gtftools  -l $ECO_GTF.genelength  $ECO_GTF.tmp
rm $ECO_GTF.tmp

# prep gene names
grep -w gene $ECO_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $ECO_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $ECO_GENEINFO
awk '{print $0,NR}' $ECO_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $ECO_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $ECO_GENEINFO


###########################################################
# HSA
###########################################################
wget -N $HSA_GTFURL && gunzip -kf $HSA_GTF.gz
wget -N $HSA_CDNAURL && gunzip -kf $HSA_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $HSA_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $HSA_CDNA \
| sed 1d | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $HSA_TXINFO

# prep the gene lengths
grep '#' $HSA_GTF > $HSA_GTF.tmp
grep -v '#' $HSA_GTF | awk '{OFS="\t"} $1=1' >> $HSA_GTF.tmp
gtftools  -l $HSA_GTF.genelength  $HSA_GTF.tmp
rm $HSA_GTF.tmp

# prep gene names
grep -w gene $HSA_GTF | cut -d '"' -f2,6 | tr '"' '\t' > $HSA_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $HSA_GENEINFO
awk '{print $0,NR}' $HSA_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $HSA_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $HSA_GENEINFO


###########################################################
# MMU
###########################################################
wget -N $MMU_GTFURL && gunzip -kf $MMU_GTF.gz
wget -N $MMU_CDNAURL && gunzip -kf $MMU_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $MMU_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $MMU_CDNA \
| sed 1d | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $MMU_TXINFO

# prep the gene lengths
grep '#' $MMU_GTF > $MMU_GTF.tmp
grep -v '#' $MMU_GTF | awk '{OFS="\t"} $1=1' >> $MMU_GTF.tmp
gtftools  -l $MMU_GTF.genelength  $MMU_GTF.tmp
rm $MMU_GTF.tmp

# prep gene names
grep -w gene $MMU_GTF | cut -d '"' -f2,6 | tr '"' '\t' > $MMU_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $MMU_GENEINFO
awk '{print $0,NR}' $MMU_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $MMU_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $MMU_GENEINFO


###########################################################
# RNO
###########################################################
wget -N $RNO_GTFURL && gunzip -kf $RNO_GTF.gz
wget -N $RNO_CDNAURL && gunzip -kf $RNO_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $RNO_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $RNO_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $RNO_TXINFO

# prep the gene lengths
grep '#' $RNO_GTF > $RNO_GTF.tmp
grep -v '#' $RNO_GTF | awk '{OFS="\t"} $1=1' >> $RNO_GTF.tmp
gtftools  -l $RNO_GTF.genelength  $RNO_GTF.tmp
rm $RNO_GTF.tmp

# prep gene names
grep -w gene $RNO_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $RNO_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $RNO_GENEINFO
awk '{print $0,NR}' $RNO_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $RNO_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $RNO_GENEINFO


###########################################################
# SCE
###########################################################
wget -N $SCE_GTFURL && gunzip -kf $SCE_GTF.gz
wget -N $SCE_CDNAURL && gunzip -kf $SCE_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $SCE_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $SCE_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $SCE_TXINFO

# prep the gene lengths
grep '#' $SCE_GTF > $SCE_GTF.tmp
grep -v '#' $SCE_GTF | awk '{OFS="\t"} $1=1' >> $SCE_GTF.tmp
gtftools  -l $SCE_GTF.genelength  $SCE_GTF.tmp
rm $SCE_GTF.tmp

# prep gene names
grep -w gene $SCE_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $SCE_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $SCE_GENEINFO
awk '{print $0,NR}' $SCE_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $SCE_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $SCE_GENEINFO

###########################################################
# ZMA
###########################################################
wget -N $ZMA_GTFURL && gunzip -kf $ZMA_GTF.gz
wget -N $ZMA_CDNAURL && gunzip -kf $ZMA_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $ZMA_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $ZMA_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $ZMA_TXINFO

# prep the gene lengths
grep '#' $ZMA_GTF > $ZMA_GTF.tmp
grep -v '#' $ZMA_GTF | awk '{OFS="\t"} $1=1' >> $ZMA_GTF.tmp
gtftools  -l $ZMA_GTF.genelength  $ZMA_GTF.tmp
rm $ZMA_GTF.tmp

# prep gene names
grep -w gene $ZMA_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $ZMA_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $ZMA_GENEINFO
awk '{print $0,NR}' $ZMA_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $ZMA_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $ZMA_GENEINFO

###########################################################
# OSA
###########################################################
wget -N $OSA_GTFURL && gunzip -kf $OSA_GTF.gz
wget -N $OSA_CDNAURL && gunzip -kf $OSA_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $OSA_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $OSA_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $OSA_TXINFO

# prep the gene lengths
grep '#' $OSA_GTF > $OSA_GTF.tmp
grep -v '#' $OSA_GTF | awk '{OFS="\t"} $1=1' >> $OSA_GTF.tmp
gtftools  -l $OSA_GTF.genelength  $OSA_GTF.tmp
rm $OSA_GTF.tmp

# prep gene names
grep -w gene $OSA_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $OSA_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $OSA_GENEINFO
awk '{print $0,NR}' $OSA_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $OSA_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $OSA_GENEINFO

###########################################################
# TAE
###########################################################
wget -N $TAE_GTFURL && gunzip -kf $TAE_GTF.gz
wget -N $TAE_CDNAURL && gunzip -kf $TAE_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $TAE_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $TAE_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $TAE_TXINFO

# prep the gene lengths
grep '#' $TAE_GTF > $TAE_GTF.tmp
grep -v '#' $TAE_GTF | awk '{OFS="\t"} $1=1' >> $TAE_GTF.tmp
gtftools  -l $TAE_GTF.genelength  $TAE_GTF.tmp
rm $TAE_GTF.tmp

# prep gene names
grep -w gene $TAE_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $TAE_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $TAE_GENEINFO
awk '{print $0,NR}' $TAE_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $TAE_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $TAE_GENEINFO


###########################################################
# SLY
###########################################################
wget -N $SLY_GTFURL && gunzip -kf $SLY_GTF.gz
wget -N $SLY_CDNAURL && gunzip -kf $SLY_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $SLY_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $SLY_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $SLY_TXINFO

# prep the gene lengths
grep '#' $SLY_GTF > $SLY_GTF.tmp
grep -v '#' $SLY_GTF | awk '{OFS="\t"} $1=1' >> $SLY_GTF.tmp
gtftools  -l $SLY_GTF.genelength  $SLY_GTF.tmp
rm $SLY_GTF.tmp

# prep gene names
grep -w gene $SLY_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $SLY_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $SLY_GENEINFO
awk '{print $0,NR}' $SLY_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $SLY_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $SLY_GENEINFO

###########################################################
# SBI
###########################################################
wget -N $SBI_GTFURL && gunzip -kf $SBI_GTF.gz
wget -N $SBI_CDNAURL && gunzip -kf $SBI_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $SBI_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $SBI_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $SBI_TXINFO

# prep the gene lengths
grep '#' $SBI_GTF > $SBI_GTF.tmp
grep -v '#' $SBI_GTF | awk '{OFS="\t"} $1=1' >> $SBI_GTF.tmp
gtftools  -l $SBI_GTF.genelength  $SBI_GTF.tmp
rm $SBI_GTF.tmp

# prep gene names
grep -w gene $SBI_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $SBI_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $SBI_GENEINFO
awk '{print $0,NR}' $SBI_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $SBI_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $SBI_GENEINFO

###########################################################
# GMA
###########################################################
wget -N $GMA_GTFURL && gunzip -kf $GMA_GTF.gz
wget -N $GMA_CDNAURL && gunzip -kf $GMA_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $GMA_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $GMA_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $GMA_TXINFO

# prep the gene lengths
grep '#' $GMA_GTF > $GMA_GTF.tmp
grep -v '#' $GMA_GTF | awk '{OFS="\t"} $1=1' >> $GMA_GTF.tmp
gtftools  -l $GMA_GTF.genelength  $GMA_GTF.tmp
rm $GMA_GTF.tmp

# prep gene names
grep -w gene $GMA_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $GMA_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $GMA_GENEINFO
awk '{print $0,NR}' $GMA_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $GMA_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $GMA_GENEINFO

###########################################################
# PTR
###########################################################
wget -N $PTR_GTFURL && gunzip -kf $PTR_GTF.gz
wget -N $PTR_CDNAURL && gunzip -kf $PTR_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $PTR_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $PTR_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $PTR_TXINFO

# prep the gene lengths
grep '#' $PTR_GTF > $PTR_GTF.tmp
grep -v '#' $PTR_GTF | awk '{OFS="\t"} $1=1' >> $PTR_GTF.tmp
gtftools  -l $PTR_GTF.genelength  $PTR_GTF.tmp
rm $PTR_GTF.tmp

# prep gene names
grep -w gene $PTR_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $PTR_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $PTR_GENEINFO
awk '{print $0,NR}' $PTR_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $PTR_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $PTR_GENEINFO

###########################################################
# VVI
###########################################################
wget -N $VVI_GTFURL && gunzip -kf $VVI_GTF.gz
wget -N $VVI_CDNAURL && gunzip -kf $VVI_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $VVI_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $VVI_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $VVI_TXINFO

# prep the gene lengths
grep '#' $VVI_GTF > $VVI_GTF.tmp
grep -v '#' $VVI_GTF | awk '{OFS="\t"} $1=1' >> $VVI_GTF.tmp
gtftools  -l $VVI_GTF.genelength  $VVI_GTF.tmp
rm $VVI_GTF.tmp

# prep gene names
grep -w gene $VVI_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $VVI_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $VVI_GENEINFO
awk '{print $0,NR}' $VVI_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $VVI_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $VVI_GENEINFO

###########################################################
# HVU
###########################################################
wget -N $HVU_GTFURL && gunzip -kf $HVU_GTF.gz
wget -N $HVU_CDNAURL && gunzip -kf $HVU_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $HVU_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $HVU_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $HVU_TXINFO

# prep the gene lengths
grep '#' $HVU_GTF > $HVU_GTF.tmp
grep -v '#' $HVU_GTF | awk '{OFS="\t"} $1=1' >> $HVU_GTF.tmp
gtftools  -l $HVU_GTF.genelength  $HVU_GTF.tmp
rm $HVU_GTF.tmp

# prep gene names
grep -w gene $HVU_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $HVU_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $HVU_GENEINFO
awk '{print $0,NR}' $HVU_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $HVU_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $HVU_GENEINFO

###########################################################
# STU
###########################################################
wget -N $STU_GTFURL && gunzip -kf $STU_GTF.gz
wget -N $STU_CDNAURL && gunzip -kf $STU_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $STU_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $STU_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $STU_TXINFO

# prep the gene lengths
grep '#' $STU_GTF > $STU_GTF.tmp
grep -v '#' $STU_GTF | awk '{OFS="\t"} $1=1' >> $STU_GTF.tmp
gtftools  -l $STU_GTF.genelength  $STU_GTF.tmp
rm $STU_GTF.tmp

# prep gene names
grep -w gene $STU_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $STU_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $STU_GENEINFO
awk '{print $0,NR}' $STU_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $STU_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $STU_GENEINFO

###########################################################
# BDI
###########################################################
wget -N $BDI_GTFURL && gunzip -kf $BDI_GTF.gz
wget -N $BDI_CDNAURL && gunzip -kf $BDI_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $BDI_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $BDI_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $BDI_TXINFO

# prep the gene lengths
grep '#' $BDI_GTF > $BDI_GTF.tmp
grep -v '#' $BDI_GTF | awk '{OFS="\t"} $1=1' >> $BDI_GTF.tmp
gtftools  -l $BDI_GTF.genelength  $BDI_GTF.tmp
rm $BDI_GTF.tmp

# prep gene names
grep -w gene $BDI_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $BDI_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $BDI_GENEINFO
awk '{print $0,NR}' $BDI_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $BDI_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $BDI_GENEINFO
