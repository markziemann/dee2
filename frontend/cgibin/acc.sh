#!/bin/bash
#set -x #-v
#/usr/lib/cgi-bin
#QUERY_STRING='ORG=athaliana&sub=Submit'
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
TODO=${DATA}/${ORG}.queue.txt
TODO_NEW=/home/ubuntu/Public/${ORG}.queue.txt
ALLOCATED=${DATA}/${ORG}.allocated.txt
WORKLIST=${DATA}/${ORG}.worklist.txt

if [ ! -r $WORKLIST ] ; then
  cut -f1 $TODO_NEW | shuf | sed 's/$/\t0/' > $WORKLIST
fi

#ACCESSION=$(trap - PIPE ; shuf $WORKLIST | sort -k2g | awk 'NR==1 {print $1}' )
ACCESSION=$(sort -k2n -k1R $WORKLIST | awk 'NR==1 {print $1}' )
COUNT=$(grep -w $ACCESSION $WORKLIST | cut -f2)
INCREMENT=$((COUNT+1))
echo $COUNT $INCREMENT

TIME=$(date +%s)
EDIT_TIME=$(stat --format "%Y" $TODO_NEW)
DIFF=$((TIME-EDIT_TIME))
sed -i "s/${ACCESSION}\t${COUNT}$/${ACCESSION}\t${INCREMENT}/" $WORKLIST

comm -23 <(cut -f1 $TODO_NEW | sort ) <(cut -f1 $WORKLIST | sort) | sed 's/$/\t0/' >> $WORKLIST
RAND=$RANDOM
grep -wFf $TODO_NEW $WORKLIST > /tmp.$RAND
mv /tmp.$RAND $WORKLIST



echo '<br>'
echo "ACCESSION=$ACCESSION"
echo '<br>'
echo "Thank you for your contribution!"
echo '<br>'
echo '</body>'
echo '</html>'

#detail joblist length for each organism
TIME=$(date +%s)
ACC_HTML=/var/www/html/acc.html
EDIT_TIME=$(stat --format "%Y" $ACC_HTML)
DIFF=$((TIME-EDIT_TIME))
if [ "$DIFF" -gt "600" ] ; then
  touch $ACC_HTML
  ORGLIST="athaliana celegans dmelanogaster drerio ecoli hsapiens mmusculus rnorvegicus scerevisiae"
  for ORG in $ORGLIST ; do
    TODO=${DATA}/${ORG}.queue.txt
    NUMJOBS=$(wc -l < $TODO)
    sed -i "/${ORG}/s/[0-9]//g;/${ORG}/s/()/(${NUMJOBS})/" $ACC_HTML
  done
fi

