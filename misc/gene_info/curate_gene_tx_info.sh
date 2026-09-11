#!/bin/bash

echo GTFtools version 0.9.0 from pip

#which gtftools.py > /dev/null || echo "GTFtools_0.6.5 is required"
#which gtftools.py > /dev/null || exit 1
#INSTALLED_VERSION=$(gtftools.py -v 2>&1 | tr -d ' ' )
#REQUIRED_VERSION="GTFtoolsversion:0.6.5"
#if [ $REQUIRED_VERSION != $INSTALLED_VERSION ] ; then echo "GTFtools_0.6.5 is required" ; exit 1 ; fi

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

OSA_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.59.gtf.gz"
OSA_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/oryza_sativa/cdna/Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz"
OSA_GTF="Oryza_sativa.IRGSP-1.0.59.gtf"
OSA_CDNA="Oryza_sativa.IRGSP-1.0.cdna.all.fa"
OSA_GENEINFO="osa_gene_info.tsv"
OSA_TXINFO="osa_tx_info.tsv"

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

ZMA_GTFURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/plants/gtf/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.gtf.gz"
ZMA_CDNAURL="ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/zea_mays/cdna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz"
ZMA_GTF="Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.gtf"
ZMA_CDNA="Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa"
ZMA_GENEINFO="zma_gene_info.tsv"
ZMA_TXINFO="zma_tx_info.tsv"

#mmulatta
MMUL_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/003/339/765/3/ensembl/2019_12/geneset/genes.gtf.gz"
MMUL_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/003/339/765/3/ensembl/2019_12/geneset/cdna.fa.bgz"
MMUL_GTF="mmul_genes.gtf"
MMUL_CDNA="mmul_cdna.fa"
MMUL_GENEINFO="mmul_gene_info.tsv"
MMUL_TXINFO="mmul_tx_info.tsv"

#ggallus
GGA_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/016/699/485/1/ensembl/2022_01/geneset/genes.gtf.bgz"
GGA_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/016/699/485/1/ensembl/2022_01/geneset/cdna.fa.bgz"
GGA_GTF="gga_genes.gtf"
GGA_CDNA="gga_cdna.fa"
GGA_GENEINFO="gga_gene_info.tsv"
GGA_TXINFO="gga_tx_info.tsv"

#sscrofa
SSC_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/003/025/6/ensembl/2022_02/geneset/genes.gtf.bgz"
SSC_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/003/025/6/ensembl/2022_02/geneset/cdna.fa.bgz"
SSC_GTF="ssc_genes.gtf"
SSC_CDNA="ssc_cdna.fa"
SSC_GENEINFO="ssc_gene_info.tsv"
SSC_TXINFO="ssc_tx_info.tsv"

#btaurus
BTA_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/002/263/795/4/ensembl/2024_11/geneset/genes.gtf.bgz"
BTA_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/002/263/795/4/ensembl/2024_11/geneset/cdna.fa.bgz"
BTA_GTF="bta_genes.gtf"
BTA_CDNA="bta_cdna.fa"
BTA_GENEINFO="bta_gene_info.tsv"
BTA_TXINFO="bta_tx_info.tsv"

#oaries
OAR_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/016/772/045/2/ensembl/2024_12/geneset/genes.gtf.bgz"
OAR_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/016/772/045/2/ensembl/2024_12/geneset/cdna.fa.bgz"
OAR_GTF="oar_genes.gtf"
OAR_CDNA="oar_cdna.fa"
OAR_GENEINFO="oar_gene_info.tsv"
OAR_TXINFO="oar_tx_info.tsv"

#mfascicularis
MFA_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/011/100/615/1/ensembl/2020_08/geneset/genes.gtf.bgz"
MFA_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/011/100/615/1/ensembl/2020_08/geneset/cdna.fa.bgz"
MFA_GTF="mfa_genes.gtf"
MFA_CDNA="mfa_cdna.fa"
MFA_GENEINFO="mfa_gene_info.tsv"
MFA_TXINFO="mfa_tx_info.tsv"

#pfalciparum
PFA_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/002/765/2/community/2017_10/geneset/genes.gtf.bgz"
PFA_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/002/765/2/community/2017_10/geneset/cdna.fa.bgz"
PFA_GTF="pfa_genes.gtf"
PFA_CDNA="pfa_cdna.fa"
PFA_GENEINFO="pfa_gene_info.tsv"
PFA_TXINFO="pfa_tx_info.tsv"

#pvivax
PVI_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/900/093/535/1/ensembl/2022_12/geneset/genes.gtf.bgz"
PVI_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/900/093/535/1/ensembl/2022_12/geneset/cdna.fa.bgz"
PVI_GTF="pvi_genes.gtf"
PVI_CDNA="pvi_cdna.fa"
PVI_GENEINFO="pvi_gene_info.tsv"
PVI_TXINFO="pvi_tx_info.tsv"

#aaegypti
AAE_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/002/204/515/1/veupathdb/2020_06/geneset/genes.gtf.gz"
AAE_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/002/204/515/1/veupathdb/2020_06/geneset/cdna.fa.bgz"
AAE_GTF="aae_genes.gtf"
AAE_CDNA="aae_cdna.fa"
AAE_GENEINFO="aae_gene_info.tsv"
AAE_TXINFO="aae_tx_info.tsv"

#aalbopictus
AAL_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCF/035/046/485/1/veupathdb/2025_06/geneset/genes.gtf.gz"
AAL_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCF/035/046/485/1/veupathdb/2025_06/geneset/cdna.fa.bgz"
AAL_GTF="aal_genes.gtf"
AAL_CDNA="aal_cdna.fa"
AAL_GENEINFO="aal_gene_info.tsv"
AAL_TXINFO="aal_tx_info.tsv"

#agambiae
AGA_GTFURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/005/575/1/veupathdb/2015_10/geneset/genes.gtf.gz"
AGA_CDNAURL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/005/575/1/veupathdb/2015_10/geneset/cdna.fa.bgz"
AGA_GTF="aga_genes.gtf"
AGA_CDNA="aga_cdna.fa"
AGA_GENEINFO="aga_gene_info.tsv"
AGA_TXINFO="aga_tx_info.tsv"

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
gtftools -l $CEL_GTF.genelength  $CEL_GTF.tmp
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
gtftools -l $DME_GTF.genelength  $DME_GTF.tmp
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
gtftools -l $DRE_GTF.genelength  $DRE_GTF.tmp
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
gtftools -l $ECO_GTF.genelength  $ECO_GTF.tmp
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
gtftools -l $HSA_GTF.genelength  $HSA_GTF.tmp
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
gtftools -l $MMU_GTF.genelength  $MMU_GTF.tmp
rm $MMU_GTF.tmp

# prep gene names
grep -w gene $MMU_GTF | cut -d '"' -f2,6 | tr '"' '\t' > $MMU_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $MMU_GENEINFO
awk '{print $0,NR}' $MMU_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $MMU_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $MMU_GENEINFO


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
gtftools  -l $OSA_GTF.genelength $OSA_GTF.tmp
rm $OSA_GTF.tmp

# prep gene names
grep -w gene $OSA_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $OSA_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $OSA_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $OSA_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $OSA_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $OSA_GENEINFO

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
gtftools -l $RNO_GTF.genelength  $RNO_GTF.tmp
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
gtftools -l $SCE_GTF.genelength  $SCE_GTF.tmp
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
gtftools  -l $ZMA_GTF.genelength $ZMA_GTF.tmp
rm $ZMA_GTF.tmp

# prep gene names
grep -w gene $ZMA_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $ZMA_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $ZMA_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $ZMA_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $ZMA_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $ZMA_GENEINFO

###########################################################
# M. mulatta
###########################################################

wget -N $MMUL_GTFURL && mv genes.gtf.gz $MMUL_GTF.gz && gunzip -kf $MMUL_GTF.gz
wget -N $MMUL_CDNAURL && mv cdna.fa.bgz $MMUL_CDNA.gz && gunzip -kf $MMUL_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $MMUL_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $MMUL_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $MMUL_TXINFO

# prep the gene lengths
grep '#' $MMUL_GTF > $MMUL_GTF.tmp
grep -v '#' $MMUL_GTF | awk '{OFS="\t"} $1=1' >> $MMUL_GTF.tmp
gtftools  -l $MMUL_GTF.genelength $MMUL_GTF.tmp
rm $MMUL_GTF.tmp

# prep gene names
grep -w gene $MMUL_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $MMUL_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $MMUL_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $MMUL_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $MMUL_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $MMUL_GENEINFO

###########################################################
# G. gallus
###########################################################

wget -N $GGA_GTFURL && mv genes.gtf.bgz $GGA_GTF.gz && gunzip -kf $GGA_GTF.gz
wget -N $GGA_CDNAURL && mv cdna.fa.bgz $GGA_CDNA.gz && gunzip -kf $GGA_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $GGA_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $GGA_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $GGA_TXINFO

# prep the gene lengths
grep '#' $GGA_GTF > $GGA_GTF.tmp
grep -v '#' $GGA_GTF | awk '{OFS="\t"} $1=1' >> $GGA_GTF.tmp
gtftools  -l $GGA_GTF.genelength $GGA_GTF.tmp
rm $GGA_GTF.tmp

# prep gene names
grep -w gene $GGA_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $GGA_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $GGA_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $GGA_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $GGA_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $GGA_GENEINFO

###########################################################
# S. scrofa
###########################################################

wget -N $SSC_GTFURL && mv genes.gtf.bgz $SSC_GTF.gz && gunzip -kf $SSC_GTF.gz
wget -N $SSC_CDNAURL && mv cdna.fa.bgz $SSC_CDNA.gz && gunzip -kf $SSC_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $SSC_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $SSC_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $SSC_TXINFO

# prep the gene lengths
grep '#' $SSC_GTF > $SSC_GTF.tmp
grep -v '#' $SSC_GTF | awk '{OFS="\t"} $1=1' >> $SSC_GTF.tmp
gtftools  -l $SSC_GTF.genelength $SSC_GTF.tmp
rm $SSC_GTF.tmp

# prep gene names
grep -w gene $SSC_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $SSC_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $SSC_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $SSC_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $SSC_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $SSC_GENEINFO

###########################################################
# B. taurus
###########################################################

wget -N $BTA_GTFURL && mv genes.gtf.bgz $BTA_GTF.gz && gunzip -kf $BTA_GTF.gz
wget -N $BTA_CDNAURL && mv cdna.fa.bgz $BTA_CDNA.gz && gunzip -kf $BTA_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $BTA_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $BTA_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $BTA_TXINFO

# prep the gene lengths
grep '#' $BTA_GTF > $BTA_GTF.tmp
grep -v '#' $BTA_GTF | awk '{OFS="\t"} $1=1' >> $BTA_GTF.tmp
gtftools  -l $BTA_GTF.genelength $BTA_GTF.tmp
rm $BTA_GTF.tmp

# prep gene names
grep -w gene $BTA_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $BTA_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $BTA_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $BTA_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $BTA_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $BTA_GENEINFO

###########################################################
# O. aries
###########################################################

wget -N $OAR_GTFURL && mv genes.gtf.bgz $OAR_GTF.gz && gunzip -kf $OAR_GTF.gz
wget -N $OAR_CDNAURL && mv cdna.fa.bgz $OAR_CDNA.gz && gunzip -kf $OAR_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $OAR_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $OAR_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $OAR_TXINFO

# prep the gene lengths
grep '#' $OAR_GTF > $OAR_GTF.tmp
grep -v '#' $OAR_GTF | awk '{OFS="\t"} $1=1' >> $OAR_GTF.tmp
gtftools  -l $OAR_GTF.genelength $OAR_GTF.tmp
rm $OAR_GTF.tmp

# prep gene names
grep -w gene $OAR_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $OAR_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $OAR_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $OAR_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $OAR_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $OAR_GENEINFO

###########################################################
# M. fascicularis
###########################################################

wget -N $MFA_GTFURL && mv genes.gtf.bgz $MFA_GTF.gz && gunzip -kf $MFA_GTF.gz
wget -N $MFA_CDNAURL && mv cdna.fa.bgz $MFA_CDNA.gz && gunzip -kf $MFA_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $MFA_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $MFA_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $MFA_TXINFO

# prep the gene lengths
grep '#' $MFA_GTF > $MFA_GTF.tmp
grep -v '#' $MFA_GTF | awk '{OFS="\t"} $1=1' >> $MFA_GTF.tmp
gtftools  -l $MFA_GTF.genelength $MFA_GTF.tmp
rm $MFA_GTF.tmp

# prep gene names
grep -w gene $MFA_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $MFA_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $MFA_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $MFA_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $MFA_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $MFA_GENEINFO

###########################################################
# P. falciparum
###########################################################

wget -N $PFA_GTFURL && mv genes.gtf.bgz $PFA_GTF.gz && gunzip -kf $PFA_GTF.gz
wget -N $PFA_CDNAURL && mv cdna.fa.bgz $PFA_CDNA.gz && gunzip -kf $PFA_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $PFA_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $PFA_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $PFA_TXINFO

# prep the gene lengths
grep '#' $PFA_GTF > $PFA_GTF.tmp
grep -v '#' $PFA_GTF | awk '{OFS="\t"} $1=1' >> $PFA_GTF.tmp
gtftools  -l $PFA_GTF.genelength $PFA_GTF.tmp
rm $PFA_GTF.tmp

# prep gene names
grep -w gene $PFA_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $PFA_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $PFA_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $PFA_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $PFA_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $PFA_GENEINFO

###########################################################
# P. vivax
###########################################################

wget -N $PVI_GTFURL && mv genes.gtf.bgz $PVI_GTF.gz && gunzip -kf $PVI_GTF.gz
wget -N $PVI_CDNAURL && mv cdna.fa.bgz $PVI_CDNA.gz && gunzip -kf $PVI_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $PVI_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $PVI_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $PVI_TXINFO

# prep the gene lengths
grep '#' $PVI_GTF > $PVI_GTF.tmp
grep -v '#' $PVI_GTF | awk '{OFS="\t"} $1=1' >> $PVI_GTF.tmp
gtftools  -l $PVI_GTF.genelength $PVI_GTF.tmp
rm $PVI_GTF.tmp

# prep gene names
grep -w gene $PVI_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $PVI_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $PVI_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $PVI_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $PVI_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $PVI_GENEINFO

###########################################################
# A. aegypti
###########################################################

wget -N $AAE_GTFURL && mv genes.gtf.gz $AAE_GTF.gz && gunzip -kf $AAE_GTF.gz
wget -N $AAE_CDNAURL && mv cdna.fa.bgz $AAE_CDNA.gz && gunzip -kf $AAE_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $AAE_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $AAE_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $AAE_TXINFO

# prep the gene lengths
grep '#' $AAE_GTF > $AAE_GTF.tmp
grep -v '#' $AAE_GTF | awk '{OFS="\t"} $1=1' >> $AAE_GTF.tmp
gtftools  -l $AAE_GTF.genelength $AAE_GTF.tmp
rm $AAE_GTF.tmp

# prep gene names
grep -w gene $AAE_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $AAE_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $AAE_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $AAE_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $AAE_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $AAE_GENEINFO

###########################################################
# A. albopictus
###########################################################

wget -N $AAL_GTFURL && mv genes.gtf.gz $AAL_GTF.gz && gunzip -kf $AAL_GTF.gz
wget -N $AAL_CDNAURL && mv cdna.fa.bgz $AAL_CDNA.gz && gunzip -kf $AAL_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $AAL_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $AAL_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $AAL_TXINFO

# prep the gene lengths
grep '#' $AAL_GTF > $AAL_GTF.tmp
grep -v '#' $AAL_GTF | awk '{OFS="\t"} $1=1' >> $AAL_GTF.tmp
gtftools  -l $AAL_GTF.genelength $AAL_GTF.tmp
rm $AAL_GTF.tmp

# prep gene names
grep -w gene $AAL_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $AAL_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $AAL_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $AAL_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $AAL_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $AAL_GENEINFO

###########################################################
# A. gambiae
###########################################################

wget -N $AGA_GTFURL && mv genes.gtf.gz $AGA_GTF.gz && gunzip -kf $AGA_GTF.gz
wget -N $AGA_CDNAURL && mv cdna.fa.bgz $AGA_CDNA.gz && gunzip -kf $AGA_CDNA.gz

# prep the cDNA lengths
echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $AGA_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $AGA_CDNA \
| sed 1d | sed '/>/s/$/ gene_symbol:NA/' | paste - - -d '!' | tr -d '>' \
| sed 's/ gene:/\n/' | sed 's/ gene_symbol:/\n/' | sed 's/!/\n/' \
| cut -d ' ' -f1  | paste - - - - \
| awk '{OFS="\t"} {print $1,$2,$3,length($4)}' >> $AGA_TXINFO

# prep the gene lengths
grep '#' $AGA_GTF > $AGA_GTF.tmp
grep -v '#' $AGA_GTF | awk '{OFS="\t"} $1=1' >> $AGA_GTF.tmp
gtftools  -l $AGA_GTF.genelength $AGA_GTF.tmp
rm $AGA_GTF.tmp

# prep gene names
grep -w gene $AGA_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $AGA_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $AGA_GENEINFO
awk '{OFS="\t"}{print $0,NR}' $AGA_GTF.genenames | sed 's/ /_/g' | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $AGA_GTF.genelength) \
| sort -k3g | sort -k3g | cut -d ' ' -f-2,4- | tr ' ' '\t' >> $AGA_GENEINFO
