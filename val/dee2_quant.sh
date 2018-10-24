#!/bin/bash

agg(){
awk  '{array[$1]+=$2} END { for (i in array) {print i, array[i]}}'
}
export -f agg

main(){
set -x
SPECIES=$1
ORG=$(echo $SPECIES | cut -c-3)
################################
echo Run some simulations
################################

if [ $ORG == "hsa" ] ; then
  REF=Homo_sapiens.GRCh38.cdna.all.fa
  wget -N ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
  gunzip -kf $REF.gz
elif [ $ORG == "ath" ] ; then
  REF=Arabidopsis_thaliana.TAIR10.cdna.all.fa
  wget -N ftp://ftp.ensemblgenomes.org/pub/release-35/plants/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
  gunzip -kf $REF.gz
fi

#check ART version art_bin_MountRainier is present
MD5=363ff45aa04df1ec28c7be46b9c437dd
ART=./art_illumina

if [ ! -r $ART ] ; then
  echo ART simulator not present. Get version ART-MountRainier-2016-06-05
  echo https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
  exit
else
  ARTMD5=$(md5sum $ART | cut -d ' ' -f1)
  if [ $MD5 != $ARTMD5 ] ; then
    echo ART simulator not present. Get version ART-MountRainier-2016-06-05
    echo https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
    exit
  fi
fi

#check docker is installed
if [ -z $(which docker) ] ; then
  echo docker not installed. Quitting.
  exit
fi

#Run the simulation
PFX=${ORG}_art
SEED=1540165885
for LEN in 50 100 ; do
  $ART --rndSeed $SEED -sam -i $REF -l $LEN -f 2 -o ${PFX}_se_${LEN} &
done

for LEN in 50 100 ; do
  $ART --rndSeed $SEED -sam -i $REF -l $LEN -p -m 300 -s 30 -f 2 -o ${PFX}_pe_${LEN}_ &
done

wait
pigz -f *fq
rm *sam *aln

################################
echo establish simulated counts
################################
grep '>' $REF | sed 's/gene:/\n/' \
| cut -d ' ' -f1 | paste - - | sed 's/>//' | cut -d '.' -f-2 > tx2gene.tsv

#tx counts
zcat ${ORG}_art_se_50.fq.gz | sed -n '1~4p' | cut -d '-' -f1 \
| sed 's/@//' | sort | uniq -c | awk '{print $2"\t"$1}' \
| sort -k 1b,1 > ${ORG}_art_se_50_tx.tsv

zcat ${ORG}_art_se_100.fq.gz | sed -n '1~4p' | cut -d '-' -f1 \
| sed 's/@//' | sort | uniq -c | awk '{print $2"\t"$1}' \
| sort -k 1b,1 > ${ORG}_art_se_100_tx.tsv

zcat ${ORG}_art_pe_50_1.fq.gz | sed -n '1~4p' | cut -d '-' -f1 \
| sed 's/@//' | sort | uniq -c | awk '{print $2"\t"$1}' \
| sort -k 1b,1 > ${ORG}_art_pe_50_tx.tsv

zcat ${ORG}_art_pe_100_1.fq.gz | sed -n '1~4p' | cut -d '-' -f1 \
| sed 's/@//' | sort | uniq -c | awk '{print $2"\t"$1}' \
| sort -k 1b,1 > ${ORG}_art_pe_100_tx.tsv

#gene counts
join -1 1 -2 1 <(sort -k 1b,1 tx2gene.tsv)  ${ORG}_art_se_50_tx.tsv | sort -k 2b,2 \
| cut -d ' ' -f2- | sort -k 1b,1 | agg | sort -k 1b,1 > ${ORG}_art_se_50_g.tsv

join -1 1 -2 1 <(sort -k 1b,1 tx2gene.tsv)  ${ORG}_art_se_100_tx.tsv | sort -k 2b,2 \
| cut -d ' ' -f2- | sort -k 1b,1 | agg | sort -k 1b,1 > ${ORG}_art_se_100_g.tsv

join -1 1 -2 1 <(sort -k 1b,1 tx2gene.tsv)  ${ORG}_art_pe_50_tx.tsv | sort -k 2b,2 \
| cut -d ' ' -f2- | sort -k 1b,1 | agg | sort -k 1b,1 > ${ORG}_art_pe_50_g.tsv

join -1 1 -2 1 <(sort -k 1b,1 tx2gene.tsv)  ${ORG}_art_pe_100_tx.tsv | sort -k 2b,2 \
| cut -d ' ' -f2- | sort -k 1b,1 | agg | sort -k 1b,1 > ${ORG}_art_pe_100_g.tsv


################################
echo Run mapping
################################

if [ $(docker images | grep -c mziemann/tallyup_${ORG}) -gt 0 ] ; then
  docker run -v $(pwd):/dee2/mnt mziemann/tallyup_${ORG} ${SPECIES} \
  -f ${ORG}_art_se_100.fq.gz,${ORG}_art_se_50.fq.gz
  rm -rf se
  docker cp $(docker ps -alq):/dee2/data/ se

  docker run -v $(pwd):/dee2/mnt mziemann/tallyup_${ORG} ${SPECIES} \
  -f ${ORG}_art_pe_100_1.fq.gz,${ORG}_art_pe_50_1.fq.gz  ${ORG}_art_pe_100_2.fq.gz,${ORG}_art_pe_50_2.fq.gz
  rm -rf pe
  docker cp $(docker ps -alq):/dee2/data/ pe
else
  docker run -v $(pwd):/dee2/mnt mziemann/tallyup ${SPECIES} \
  -f ${ORG}_art_se_100.fq.gz,${ORG}_art_se_50.fq.gz
  rm -rf se
  docker cp $(docker ps -alq):/dee2/data/ se

  docker run -v $(pwd):/dee2/mnt mziemann/tallyup ${SPECIES} \
  -f ${ORG}_art_pe_100_1.fq.gz,${ORG}_art_pe_50_1.fq.gz  ${ORG}_art_pe_100_2.fq.gz,${ORG}_art_pe_50_2.fq.gz
  rm -rf pe
  docker cp $(docker ps -alq):/dee2/data/ pe
fi


################################
echo join simulated and processed data
################################
#SINGLE END
#50 bp SE reads
#sim vs star
join -1 1 -2 1 ${ORG}_art_se_50_g.tsv <(sort -k 1b,1 se/${SPECIES}/${ORG}_art_se_50/${ORG}_art_se_50.se.tsv) \
> ${ORG}_art_se_50_star.tsv

#sim vs kal transcripts
cut -f1,4 se/${SPECIES}/${ORG}_art_se_50/${ORG}_art_se_50.ke.tsv | sed 1d \
| sort -k 1b,1 | join -1 1 -2 1 ${ORG}_art_se_50_tx.tsv - | tr ' ' '\t' > ${ORG}_art_se_50_kalt.tsv

#sim vs kal genes
cut -f1,4 se/${SPECIES}/${ORG}_art_se_50/${ORG}_art_se_50.ke.tsv | sed 1d | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 tx2gene.tsv ) | awk '{print $3"\t"$2}' \
| sort | agg | sort -k 1b,1 \
| join -1 1 -2 1 ${ORG}_art_se_50_g.tsv - > ${ORG}_art_se_50_kalg.tsv

#100 bp SE reads
#sim vs star
join -1 1 -2 1 ${ORG}_art_se_100_g.tsv <(sort -k 1b,1 se/${SPECIES}/${ORG}_art_se_100/${ORG}_art_se_100.se.tsv) \
> ${ORG}_art_se_100_star.tsv

#sim vs kal transcripts
cut -f1,4 se/${SPECIES}/${ORG}_art_se_100/${ORG}_art_se_100.ke.tsv | sed 1d \
| sort -k 1b,1 | join -1 1 -2 1 ${ORG}_art_se_100_tx.tsv - | tr ' ' '\t' > ${ORG}_art_se_100_kalt.tsv

#sim vs kal genes
cut -f1,4 se/${SPECIES}/${ORG}_art_se_100/${ORG}_art_se_100.ke.tsv | sed 1d | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 tx2gene.tsv ) | awk '{print $3"\t"$2}' \
| sort | agg | sort -k 1b,1 \
| join -1 1 -2 1 ${ORG}_art_se_100_g.tsv - > ${ORG}_art_se_100_kalg.tsv


##PAIRED END
#50 bp PE reads
join -1 1 -2 1 ${ORG}_art_se_50_g.tsv <(sort -k 1b,1 pe/${SPECIES}/${ORG}_art_pe_50_1/${ORG}_art_pe_50_1.se.tsv) \
> ${ORG}_art_pe_50_star.tsv

#sim vs kal transcripts
cut -f1,4 pe/${SPECIES}/${ORG}_art_pe_50_1/${ORG}_art_pe_50_1.ke.tsv | sed 1d \
| sort -k 1b,1 | join -1 1 -2 1 ${ORG}_art_pe_50_tx.tsv - | tr ' ' '\t' > ${ORG}_art_pe_50_kalt.tsv

#sim vs kal genes
cut -f1,4 pe/${SPECIES}/${ORG}_art_pe_50_1/${ORG}_art_pe_50_1.ke.tsv | sed 1d | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 tx2gene.tsv ) | awk '{print $3"\t"$2}' \
| sort | agg | sort -k 1b,1 \
| join -1 1 -2 1 ${ORG}_art_pe_50_g.tsv - > ${ORG}_art_pe_50_kalg.tsv

#100 bp SE reads
#sim vs star
join -1 1 -2 1 ${ORG}_art_pe_100_g.tsv <(sort -k 1b,1 pe/${SPECIES}/${ORG}_art_pe_100_1/${ORG}_art_pe_100_1.se.tsv) \
> ${ORG}_art_pe_100_star.tsv

#sim vs kal transcripts
cut -f1,4 pe/${SPECIES}/${ORG}_art_pe_100_1/${ORG}_art_pe_100_1.ke.tsv | sed 1d \
| sort -k 1b,1 | join -1 1 -2 1 ${ORG}_art_pe_100_tx.tsv - | tr ' ' '\t' > ${ORG}_art_pe_100_kalt.tsv

#sim vs kal genes
cut -f1,4 pe/${SPECIES}/${ORG}_art_pe_100_1/${ORG}_art_pe_100_1.ke.tsv | sed 1d | sort -k 1b,1 \
| join -1 1 -2 1 - <(sort -k 1b,1 tx2gene.tsv ) | awk '{print $3"\t"$2}' \
| sort | agg | sort -k 1b,1 \
| join -1 1 -2 1 ${ORG}_art_pe_100_g.tsv - > ${ORG}_art_pe_100_kalg.tsv

}
export -f main

main athaliana
main hsapiens
