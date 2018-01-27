#!/bin/bash
set -x

ORG=$1
CODE_DIR=$(pwd)
PROJ_DIR=$(echo $CODE_DIR | rev | cut -d '/' -f2- | rev )
DATA_DIR=$PROJ_DIR/data/$ORG/
MX_DIR=$PROJ_DIR/mx/
SE_MX=$MX_DIR/${ORG}_se.mx

#use list output by R to generate matrix

SE_LIST=$DATA_DIR/${ORG}_se_list.txt
SE_GENES=$DATA_DIR/${ORG}_se_genes.txt

#divide the list into sets of 1000
cd $DATA_DIR
rm *split
split -l 500 --additional-suffix=split $SE_LIST

paste1(){
  SPLIT=$1
  cut -d '/' -f3 $SPLIT | tr '\n' '\t' | sed 's#\t$#\n#'  > $SPLIT.tsv
  paste  $(cat $SPLIT) \
  | awk '{for(i=2;i<=NF;i=i+2){printf "%s\t", $i}{printf "%s", RS}}' \
  | sed 1d >> $SPLIT.tsv
}
export -f paste1
parallel paste1 ::: *split

#Single core
#for SPLIT in *split ; do
#  cut -d '/' -f3 $SPLIT | tr '\n' '\t' | sed 's#\t$#\n#'  > $SPLIT.tsv
#  paste  $(cat $SPLIT) \
#  | awk '{for(i=2;i<=NF;i=i+2){printf "%s\t", $i}{printf "%s", RS}}' \
#  | sed 1d >> $SPLIT.tsv
#done

echo GeneID | sed 's/$/\t/' > $SE_GENES
cut -f1 $(head -1 $SE_LIST ) | sed 1d >> $SE_GENES

#bring it all together
ls *split*.tsv

paste $SE_GENES *split*.tsv > $SE_MX

#tidy up temporary files
rm $SE_GENES *split *split.tsv
