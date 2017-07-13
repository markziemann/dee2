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
TODO_NEW=/home/pi/Public/${ORG}.queue.txt
cp -u $TODO_NEW $TODO
ALLOCATED=${DATA}/${ORG}.allocated.txt
ACCESSION=$(tail -n +100 $TODO | head -1)
TIME=$(date +%s)
echo $ACCESSION $TIME >> $ALLOCATED
sed -i "/${ACCESSION}/d" $TODO

echo '<br>'
echo "ACCESSION=$ACCESSION"
echo '<br>'
echo "Thank you for your contribution!"
echo '<br>'
echo '</body>'
echo '</html>'

TIME=$(date +%s)
sed "s/$/\t${TIME}/" $ALLOCATED | awk '($3-$2)>86400 {print $1}' >> $TODO
sort -u -o $TODO $TODO

exit 1

