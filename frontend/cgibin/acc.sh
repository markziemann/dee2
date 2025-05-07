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

SPEC=$(echo "$QUERY_STRING"  | egrep -c '(\!|\@|\#|\$|\%|\^|\*|\:|\;|\}|\{|\[|\]|\/|\\|\@)')
if [ "$SPEC" -gt 0 ] ; then
  echo "$QUERY_STRING"
  echo "<br>"
  echo "Avoid special characters"
  exit
fi
QUERY_STRING=$(echo "$QUERY_STRING" | tr -d ':;{}()[]\/<>@#$%^*' )

STRING_LEN=$(echo "$QUERY_STRING" | wc -c)
if [ "$STRING_LEN" -gt 50 ] ; then
  echo "$QUERY_STRING"
  echo "<br>"
  echo "Invalid input. Please keep queries to less than 50 characters."
  exit
fi

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
osativa
rnorvegicus
scerevisiae
zmays'

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
  awk '{OFS="\t"}{print $0,$0}' $NEWQ | cut -c4-  | sort -k1 -gr | cut -f2 > $DATA/tmp
  mv $DATA/tmp $Q
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
