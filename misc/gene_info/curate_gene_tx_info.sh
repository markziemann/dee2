#!/bin/bash
set -x

ATH_GTFURL="ftp://ftp.ensemblgenomes.org/pub/release-36/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.36.gtf.gz"
ATH_CDNAURL="ftp://ftp.ensemblgenomes.org/pub/release-36/plants/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
CEL_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.90.gtf.gz"
CEL_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz"
DME_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.90.gtf.gz"
DME_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz"
DRE_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/danio_rerio/Danio_rerio.GRCz10.90.gtf.gz"
DRE_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/danio_rerio/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz"
ECO_GTFURL="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.36.gtf.gz"
ECO_CDNAURL="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cdna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa.gz"

HSA_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz"
HSA_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
HSA_GTF="Homo_sapiens.GRCh38.90.gtf"
HSA_CDNA="Homo_sapiens.GRCh38.cdna.all.fa"
HSA_GENEINFO=hsa_gene_info.tsv
HSA_TXINFO=hsa_tx_info.tsv

MMU_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz"
MMU_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
RNO_GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf.gz"
RNO_CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz"
SCE_GTFURL="ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.36.gtf.gz"
SCE_CDNAURL="ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"



# HSA Gene length
wget -N $HSA_GTFURL && gunzip -kf $HSA_GTF.gz
wget -N $HSA_CDNAURL && gunzip -kf $HSA_CDNA.gz

grep -w gene $HSA_GTF | cut -d '"' -f2,6 | tr '"' '\t' > $HSA_GENEINFO

echo "GeneID GeneSymbol GeneLength" > $HSA_GENEINFO.tmp

length(){
 ACCESSION=$1
 GENESYMBOL=$2
 HSA_GTF=$3
 GENE_LENGTH=$(grep -w $ACCESSION  $HSA_GTF | cut -f1,4,5 | bedtools sort | bedtools merge | awk '{print sqrt( ($3-$2)^2) }')
 echo $ACCESSION $GENESYMBOL $GENE_LENGTH | tr ' ' '\t'
}
export -f length
cat $HSA_GENEINFO | parallel --colsep '\t' length {1} {2} $HSA_GTF >> $HSA_GENEINFO.tmp && mv $HSA_GENEINFO.tmp $HSA_GENEINFO

echo "TxID GeneID GeneSymbol TxLength" | tr ' ' '\t' > $HSA_TXINFO
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $HSA_CDNA \
| sed 1d | paste - - -d '@'  | while read line ; do
  TX_ACCESSION=$(echo $line | cut -d ' ' -f1 | tr -d '>' )
  GENE_ACCESSION=$(echo $line | tr ' ' '\n' | grep gene: | cut -d ':' -f2 | cut -d '.' -f1)
  GENE_SYMBOL=$(echo $line | tr ' ' '\n' | grep gene_symbol: | cut -d ':' -f2 )
  LENGTH=$(echo $line  | rev | cut -d '@' -f1 | awk '{print length($1)}')
  echo $TX_ACCESSION $GENE_ACCESSION $GENE_SYMBOL $LENGTH
done | tr ' ' '\t' >> $HSA_TXINFO

