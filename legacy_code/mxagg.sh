#!/bin/bash
#set -x
#aggregate individual counts files into a mega matrix
ORG=$1
#ORG=rnorvegicus
VERS=v1
MX=${ORG}_gene.mx
EMX=`echo $MX | sed 's/gene.mx/exon.mx/'`
GNAMES=$2
MXG=`echo $MX | sed 's/.mx$/_n.mx/'`
EMXG=`echo $EMX | sed 's/.mx$/_n.mx/'`
HD=header.txt
IND=individual_files

#####################################################################
#Assemble the gene-wise matrix
#####################################################################
#Make the header
Z=`find .. | grep gene.cnt.bz2 | head -1`
pbzip2 -dc $Z | sed 1d | cut -f1 \
| gawk -f tmx.awk | tr ' ' '\t' > $HD

#Get read counts
#for Z in `find .. | grep gene.cnt.bz2 ` ; do
#pbzip2 -dc $Z | sed 1d | cut -f7 | gawk -f tmx.awk | tr ' ' '\t'
#done | sed "s/fastq.subjunc.bam/${VERS}cs/" | sed "s/bam/${VERS}/" > tmp

#make the genewise unnamed flat matrix
find .. | grep _gene.cnt.bz2 \
| sort \
| parallel -j -4 " pbzip2 -dc {} | sed 1d \
| awk '{if(NR>1) print | \"sort -k 1b,1 \"; else print}' \
| cut -f7 | gawk -f tmx.awk | tr ' ' '\t' " \
| sed "s/fastq.subjunc.bam/${VERS}cs/" | sed "s/bam/${VERS}/" \
| cat $HD - > $MX
#(read -r; printf "%s\n" "$REPLY"; sort)
#awk 'NR == 1; NR > 1 {print $0 | \"sort -n\"}' \

#prep the gene names
echo "Geneid GeneName" > tmp2
gawk -f tmx.awk $HD | sed 1d | sort -k 1b,1 | join -1 1 -2 1 $GNAMES - \
| cat tmp2 - | tr ' ' '\t' | gawk -f tmx.awk | tr ' ' '\t' > tmp
mv tmp $HD
sed 1d $MX | cat $HD - > $MXG

mkdir $IND
find .. | grep _gene.cnt.bz2 \
| sort \
| parallel -j -4 " pbzip2 -dc {} | sed 1d \
| awk '{if(NR>1) print | \"sort -k 1b,1 \"; else print}' \
| cut -f7 | sed 's/fastq.subjunc.bam/v1cs/' | sed 's/.bam/v1/' > ${IND}/{/.} "

gawk -f tmx.awk $HD | sed 's/ $//' | tr ' ' '\t' | cut -f-2 > ${IND}/rownames.txt

#scp ${IND}/ mziemann
#${ORG}












exit
#####################################################################
#Assemble the exon-wise matrix
#####################################################################
#Make the header
Z=`find .. | grep exon.cnt.bz2 | head -1 `
pbzip2 -dc $Z | sed 1d | cut -f-4 | gawk -f tmx.awk | tr ' ' '\t' > $HD

#Get exon counts and make the exonwise unnamed flat matrix
find .. | grep exon.cnt.bz2 \
| sort \
| parallel -j -4 "pbzip2 -dc {} | sed 1d \
| awk '{if(NR>1) print | \"sort -k 1b,1 \"; else print}' \
| cut -f7 | gawk -f tmx.awk | tr ' ' '\t' " \
| sed "s/fastq.subjunc.bam/${VERS}cs/" | sed "s/bam/${VERS}/" \
| cat $HD - > $EMX

#prep the gene names
echo "Geneid GeneName Chr Start End" > tmp2
gawk -f tmx.awk $HD | sed 1d | sort -k 1b,1 | join -1 1 -2 1 $GNAMES - \
| cat tmp2 - | tr ' ' '\t' | gawk -f tmx.awk | tr ' ' '\t' > tmp
mv tmp $HD
sed '1,4d' $EMX | cat $HD - > $EMXG
rm tmp2 $HD


exit

gawk '$0 !~ /#/{arr[$1]=arr[$1] " " $2}END{for(i in arr)print i,arr[i]}' *f | tr -s ' ' \
| tr ' ' '\t' > $EXONMX

grep ^Geneid $EXONMX > $HD

grep -v ^Geneid $EXONMX | cat $HD - | sed 's/.bam//g' > tmp

mv tmp $EXONMX

rm *.cnt *.cnt.f

#attach genename
sed 's/:/\t/' $EXONMX | sort -k 1b,1 | join -1 1 -2 1 $GNAMES - | cat $HD - \
| sed 's/^Geneid/Geneid Genename Coordinate/' | tr ' ' '\t' > $EXONMXG

rm $HD

#make a transposed copy
for M in $EXONMX $EXONMXG
do
NAME=`echo $M | sed 's/.mx$/_t.mx/'`
gawk -f tmx.awk $M | tr ' ' '\t' > $NAME &
done

wait

exit

#tmx.R code
x<-read.table("mx",header=TRUE,row.names=1)
t<-t(x)
write.table(t,file="mxt",sep="\t",quote=F,col.names=T,row.names=T)

#tmx awk code
$ cat tmx.awk
# syntax: GAWK -f MATRIX_TRANSPOSITION.AWK filename
{   if (NF > nf) {
      nf = NF
    }
    for (i=1; i<=nf; i++) {
      row[i] = row[i] $i " "
    }
}
END {
    for (i=1; i<=nf; i++) {
      printf("%s\n",row[i])
    }
    exit(0)
}
