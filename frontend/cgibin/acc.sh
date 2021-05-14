#!/bin/bash
#set -x #-v
#/usr/lib/cgi-bin
#QUERY_STRING='ORG=ecoli&sub=Submit'
cleanit(){
tr '<>&*?/' ' '
}
export -f cleanit

echo "Content-type: text/html"
echo ""
echo '<html>'
echo '<head>'
echo '<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">'
echo '<title>Tempemon</title>'
echo '</head>'
echo '<body>'
# check whether ORG is known and quit if not in the list
ORG=$(echo "$QUERY_STRING" | cleanit | tr ' ' '\n' | grep -m1 ^ORG=  | cut -d '=' -f2)
ORGLIST='athaliana
celegans
dmelanogaster
drerio
ecoli
hsapiens
mmusculus
rnorvegicus
scerevisiae'

ORG_OK=$(echo $ORGLIST | grep -cw $ORG )

if [ "$ORG_OK" != "1" ] ; then
  echo "Unknown species selected. Quitting."
  echo '<br>'
  echo "Help available from the developer mark.ziemann@gmail.com"
  echo '<br>'
  echo '</body>'
  echo '</html>'
  exit 1
fi

DATA=/usr/lib/cgi-bin/acc_data
TODO_NEW=/home/ubuntu/Public/${ORG}.queue.txt
ALLOC=${DATA}/${ORG}.alloc.txt
SHORTLIST=${DATA}/${ORG}.shortlist.txt

if [ ! -r $SHORTLIST ] ; then
  awk '{OFS="\t"}{print $0,$0}' $TODO_NEW | cut -c4-  | sort -k1 -gr | cut -f2 > tmp
  mv tmp $TODO_NEW
  head -1000 $TODO_NEW > $SHORTLIST
  ACCESSION=$(head -1 $SHORTLIST)
  echo $ACCESSION > $ALLOC
else
  PREV=$(tail -1 $ALLOC)
  ACCESSION=$(grep -w -A1 $PREV $SHORTLIST | tail -1)
  echo $ACCESSION >> $ALLOC
fi

DIFF=$(( $(wc -l < $SHORTLIST )  - $(wc -l < $ALLOC ) ))

if [ $DIFF -le 50 ] ; then
  grep -w -A1000 $ACCESSION $TODO_NEW | sed 1d > $SHORTLIST
  > $ALLOC
  if [  $(wc -l < $SHORTLIST) -le 1 ] ; then rm $SHORTLIST ; fi
fi

echo '<br>'
echo "ACCESSION=$ACCESSION"
echo '<br>'
echo "Thank you for your contribution!"
echo '<br>'
echo '</body>'
echo '</html>'
