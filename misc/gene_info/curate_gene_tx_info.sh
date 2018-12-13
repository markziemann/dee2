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
SCE_GTFURL="ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.36.gtf.gz"
SCE_CDNAURL="ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"






# MMU
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





















# HSA
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





