#!/bin/bash

agg(){
awk  '{array[$1]+=$2} END { for (i in array) {print i, array[i]}}'
}
export -f agg

################################
echo Run some simulations
################################

REF=../../genomes/Homo_sapiens.GRCh38.cdna.all.fa
ART=../../software/art_bin_MountRainier/art_illumina
PFX=hsa_art
SEED=1540165885
for LEN in 50 100 ; do
  $ART --rndSeed $SEED -sam -i $REF -l $LEN -f 2 -o ${PFX}_se_${LEN} &
done

for LEN in 50 100 ; do
  $ART --rndSeed $SEED -sam -i $REF -l $LEN -p -m 300 -s 30 -f 2 -o ${PFX}_pe_${LEN} &
done

wait

################################
echo establish simulated counts
################################
grep '>' ../../genomes/Homo_sapiens.GRCh38.cdna.all.fa | sed 's/gene:/\n/' \
| cut -d ' ' -f1 | paste - - | sed 's/>//' | cut -d '.' -f-2 > tx2gene.tsv

#tx counts
zcat hsa_art_se_50.fq.gz | sed -n '1~4p' | cut -d '-' -f1 \
| sed 's/@//' | sort | uniq -c | awk '{print $2"\t"$1}' \
| sort -k 1b,1 > hsa_art_se_50_tx.tsv

zcat hsa_art_se_100.fq.gz | sed -n '1~4p' | cut -d '-' -f1 \
| sed 's/@//' | sort | uniq -c | awk '{print $2"\t"$1}' \
| sort -k 1b,1 > hsa_art_se_100_tx.tsv

zcat hsa_art_pe_50_1.fq.gz | sed -n '1~4p' | cut -d '-' -f1 \
| sed 's/@//' | sort | uniq -c | awk '{print $2"\t"$1}' \
| sort -k 1b,1 > hsa_art_pe_50_tx.tsv

zcat hsa_art_pe_100_1.fq.gz | sed -n '1~4p' | cut -d '-' -f1 \
| sed 's/@//' | sort | uniq -c | awk '{print $2"\t"$1}' \
| sort -k 1b,1 > hsa_art_pe_100_tx.tsv

#gene counts
join -1 1 -2 1 <(sort -k 1b,1 tx2gene.tsv)  hsa_art_se_50.tsv | sort -k 2b,2 \
| cut -d ' ' -f2- | sort -k 1b,1 | agg | sort -k 1b,1 > hsa_art_se_50_g.tsv

join -1 1 -2 1 <(sort -k 1b,1 tx2gene.tsv)  hsa_art_se_100.tsv | sort -k 2b,2 \
| cut -d ' ' -f2- | sort -k 1b,1 | agg | sort -k 1b,1 > hsa_art_se_100_g.tsv

join -1 1 -2 1 <(sort -k 1b,1 tx2gene.tsv)  hsa_art_pe_50_tx.tsv | sort -k 2b,2 \
| cut -d ' ' -f2- | sort -k 1b,1 | agg | sort -k 1b,1 > hsa_art_pe_50_g.tsv

join -1 1 -2 1 <(sort -k 1b,1 tx2gene.tsv)  hsa_art_pe_100_tx.tsv | sort -k 2b,2 \
| cut -d ' ' -f2- | sort -k 1b,1 | agg | sort -k 1b,1 > hsa_art_pe_100_g.tsv


################################
echo Run mapping
################################
docker run -v /mnt/data/rna-sim/hsapiens:/dee2/mnt mziemann/tallyup_hsa hsapiens \
-f hsa_art_se_100.fq.gz,hsa_art_se_100.fq.gz
docker cp $(docker ps -alq):/dee2/data/ se

docker run -v $(pwd):/mnt mziemann/tallyup_hsa hsapiens \
-f hsa_art_pe_100_1.fq.gz,hsa_art_pe_50_1.fq.gz  hsa_art_pe_100_2.fq.gz,hsa_art_pe_50_2.fq.gz
docker cp $(docker ps -alq):/dee2/data/ pe


################################
echo join simulated and processed data
################################
#SINGLE END
#50 bp SE reads
#sim vs star
join -1 1 -2 1 hsa_art_se_50_g.tsv <(sort -k 1b,1 se/hsapiens/hsa_art_se_50/hsa_art_se_50.se.tsv) \
> hsa_art_se_50_star.tsv

#sim vs kal transcripts
cut -f1,4 se/hsapiens/hsa_art_se_50/hsa_art_se_50.ke.tsv | sed 1d \
| sort -k 1b,1 | join -1 1 -2 1 hsa_art_se_50.tsv - | tr ' ' '\t' > hsa_art_se_50_kalt.tsv

#sim vs kal genes
cut -f1,4 se/hsapiens/hsa_art_se_50/hsa_art_se_50.ke.tsv | sed 1d | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 tx2gene.tsv ) | awk '{print $3"\t"$2}' \
| sort | agg | sort -k 1b,1 \
| join -1 1 -2 1 hsa_art_se_50_g.tsv - > hsa_art_se_50_kalg.tsv

#100 bp SE reads
#sim vs star
join -1 1 -2 1 hsa_art_se_100_g.tsv <(sort -k 1b,1 se/hsapiens/hsa_art_se_100/hsa_art_se_100.se.tsv) \
> hsa_art_se_100_star.tsv

#sim vs kal transcripts
cut -f1,4 se/hsapiens/hsa_art_se_100/hsa_art_se_100.ke.tsv | sed 1d \
| sort -k 1b,1 | join -1 1 -2 1 hsa_art_se_100.tsv - | tr ' ' '\t' > hsa_art_se_100_kalt.tsv

#sim vs kal genes
cut -f1,4 se/hsapiens/hsa_art_se_100/hsa_art_se_100.ke.tsv | sed 1d | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 tx2gene.tsv ) | awk '{print $3"\t"$2}' \
| sort | agg | sort -k 1b,1 \
| join -1 1 -2 1 hsa_art_se_100_g.tsv - > hsa_art_se_100_kalg.tsv


##PAIRED END
#50 bp PE reads
join -1 1 -2 1 hsa_art_se_50_g.tsv <(sort -k 1b,1 pe/hsapiens/hsa_art_pe_50_1/hsa_art_pe_50_1.se.tsv) \
> hsa_art_pe_50_star.tsv

#sim vs kal transcripts
cut -f1,4 pe/hsapiens/hsa_art_pe_50_1/hsa_art_pe_50_1.ke.tsv | sed 1d \
| sort -k 1b,1 | join -1 1 -2 1 hsa_art_pe_50_tx.tsv - | tr ' ' '\t' > hsa_art_pe_50_kalt.tsv

#sim vs kal genes
cut -f1,4 pe/hsapiens/hsa_art_pe_50_1/hsa_art_pe_50_1.ke.tsv | sed 1d | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 tx2gene.tsv ) | awk '{print $3"\t"$2}' \
| sort | agg | sort -k 1b,1 \
| join -1 1 -2 1 hsa_art_pe_50_g.tsv - > hsa_art_pe_50_kalg.tsv

#100 bp SE reads
#sim vs star
join -1 1 -2 1 hsa_art_pe_100_g.tsv <(sort -k 1b,1 pe/hsapiens/hsa_art_pe_100_1/hsa_art_pe_100_1.se.tsv) \
> hsa_art_pe_100_star.tsv

#sim vs kal transcripts
cut -f1,4 pe/hsapiens/hsa_art_pe_100_1/hsa_art_pe_100_1.ke.tsv | sed 1d \
| sort -k 1b,1 | join -1 1 -2 1 hsa_art_pe_100_tx.tsv - | tr ' ' '\t' > hsa_art_pe_100_kalt.tsv

#sim vs kal genes
cut -f1,4 pe/hsapiens/hsa_art_pe_100_1/hsa_art_pe_100_1.ke.tsv | sed 1d | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 tx2gene.tsv ) | awk '{print $3"\t"$2}' \
| sort | agg | sort -k 1b,1 \
| join -1 1 -2 1 hsa_art_pe_100_g.tsv - > hsa_art_pe_100_kalg.tsv

