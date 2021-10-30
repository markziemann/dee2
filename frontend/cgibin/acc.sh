#!/bin/bash
#set -x #-v
#/usr/lib/cgi-bin
#QUERY_STRING='ORG=drerio&sub=Submit'
cleanit(){
tr '<>&*?/' ' '
}
export -f cleanit

echo "Content-type: text/html"
echo ""
echo '<html>'
echo '<head>'
echo '<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">'
echo '<title>DEE2 job queue</title>'
echo '</head>'
echo '<body>'
# check whether ORG is known and quit if not in the list
ORG=$(echo "$QUERY_STRING" | cleanit | tr ' ' '\n' | grep -m1 ^ORG=  | cut -d '=' -f2)
#ORG="ecoli"
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
Q=$DATA/${ORG}.queue.txt
NEWQ=/home/ubuntu/Public/${ORG}.queue.txt

if [ -r $NEWQ ] ; then
  awk '{OFS="\t"}{print $0,$0}' $NEWQ | cut -c4-  | sort -k1 -gr | cut -f2 > tmp
  mv tmp $Q
  rm $NEWQ
fi

ACCESSION=$(head -1 $Q)
sed -i '1d' $Q

echo '<br>'
echo "ACCESSION=$ACCESSION"
echo '<br>'
echo "Thank you for your contribution!"
echo '<br>'
echo '</body>'
echo '</html>'
