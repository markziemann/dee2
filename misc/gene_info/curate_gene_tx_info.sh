#!/bin/bash
set -x

which gtftools.py > /dev/null || echo "GTFtools_0.6.5 is required"
which gtftools.py > /dev/null || exit 1
INSTALLED_VERSION=$(gtftools.py -v 2>&1 | tr -d ' ' )
REQUIRED_VERSION="GTFtoolsversion:0.6.5"
if [ $REQUIRED_VERSION != $INSTALLED_VERSION ] ; then echo "GTFtools_0.6.5 is required" ; exit 1 ; fi

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
gtftools.py  -l $ATH_GTF.genelength $ATH_GTF.tmp
rm $ATH_GTF.tmp

# prep gene names
grep -w gene $ATH_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
| sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - - > $ATH_GTF.genenames
# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $ATH_GENEINFO
awk '{print $0,NR}' $ATH_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $ATH_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $ATH_GENEINFO


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
gtftools.py  -l $CEL_GTF.genelength  $CEL_GTF.tmp
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
gtftools.py  -l $DME_GTF.genelength  $DME_GTF.tmp
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
gtftools.py  -l $DRE_GTF.genelength  $DRE_GTF.tmp
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
gtftools.py  -l $ECO_GTF.genelength  $ECO_GTF.tmp
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
gtftools.py  -l $HSA_GTF.genelength  $HSA_GTF.tmp
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
gtftools.py  -l $MMU_GTF.genelength  $MMU_GTF.tmp
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
gtftools.py  -l $RNO_GTF.genelength  $RNO_GTF.tmp
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
gtftools.py  -l $SCE_GTF.genelength  $SCE_GTF.tmp
rm $SCE_GTF.tmp

# prep gene names
grep -w gene $SCE_GTF | sed 's/$/gene_name "NA"/' | sed 's/gene_id "/\ngene_id "/' \
|  sed 's/gene_name "/\ngene_name "/' | grep ^gene | cut -d '"' -f2 | paste - -  > $SCE_GTF.genenames

# merge gene names
echo "GeneID GeneSymbol mean median longest_isoform merged" | tr ' ' '\t' > $SCE_GENEINFO
awk '{print $0,NR}' $SCE_GTF.genenames | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 $SCE_GTF.genelength) \
| sort -k3g | tr ' '  '\t' | cut -f-2,4- >> $SCE_GENEINFO
