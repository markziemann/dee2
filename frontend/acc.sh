#!/bin/bash
set -x #-v
#/usr/lib/cgi-bin
QUERY_STRING='ORG=athaliana&sub=Submit'

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
if [ -r $TODO_NEW ] ; then
  mv $TODO_NEW $TODO
fi
#cp -u $TODO_NEW $TODO
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
sed "s/$/\t${TIME}/" $ALLOCATED | awk '($3-$2)<86400' > $ALLOCATED.tmp
mv $ALLOCATED.tmp $ALLOCATED

#detail joblist length for each organism
TIME=$(date +%s)
ACC_HTML=/var/www/html/acc.html
EDIT_TIME=$(stat --format "%Y" $ACC_HTML)
DIFF=$((TIME-EDIT_TIME))
if [ "$DIFF" -gt "86400" ] ; then
#if [ "$DIFF" -gt "0" ] ; then
  touch $ACC_HTML
  ORGLIST="athaliana celegans dmelanogaster drerio ecoli hsapiens mmusculus rnorvegicus scerevisiae"
  for ORG in $ORGLIST ; do
    TODO=${DATA}/${ORG}.queue.txt
    NUMJOBS=$(wc -l < $TODO)
    sed -i "/${ORG}/s/[0-9]//g;/${ORG}/s/()/(${NUMJOBS})/" $ACC_HTML
  done
fi

#the below needs to be integrated to automate incorporation of volunteer data
if [ ! -r started ] ; then
  touch started
  SFTP_INCOMING=/sftp/guestuser/incoming
validate_zip(){
ZIP=$1
#zipinfo $ZIP
SRR=$(echo $ZIP | cut -d '.' -f1)
SRR=$(basename $SRR)
# check names of the contents by checking the md5sum of the names after stripping
# SRR number and sorting
INVALID=0
MD5SUM=$(unzip -ql $ZIP | grep RR | cut -d '/' -f2 | sed "s#${SRR}.##" \
| sort | md5sum | awk '{print $1}')

# Check the md5sum of the current volunteer pipeline and check it against incoming
# datasets
REFERENCE_PIPELINE_MD5SUM=$(scp mdz@Z620:/home/mdz/bfx/dee2/code/volunteer_pipeline.sh  /dev/stdout \
| md5sum | awk '{print $1}')

PIPELINE_MD5SUM=$(unzip -p $ZIP $SRR/volunteer_pipeline.sh | md5sum | awk '{print $1}')

if [ "$REFERENCE_PIPELINE_MD5SUM" != "PIPELINE_MD5SUM" ] ; then INVALID=$((INVALID+1)) ; fi

REPORTED_MD5SUM=$(unzip -p $ZIP $SRR/volunteer_pipeline.sh | md5sum | awk '{print $1}')
if [ "$MD5SUM" != "579c3e8fe9c178e57199c135a5e29de4" ] ; then INVALID=$((INVALID+1)) ; fi

}
export -f validate_zip


#stick this script on the acc.sh script
if [ "$(ls -A ${SFTP_INCOMING})" ]; then

  for FILE in ${SFTP_INCOMING}/* ; do

    if [ "$(echo $FILE | grep -c .zip$)" -ne "1" ] ; then

      echo "now rm $FILE"
      sudo rm $FILE

    else

      TIME=$(date +%s)
      FILETIME=$(stat --format "%Y" $FILE)

      if [ $((TIME-FILETIME)) -gt "60" ] ; then

        echo process $FILE $TIME $FILETIME
        validate_zip $FILE
        if [ "$INVALID" -eq "0" ] ; then
          unzip $FILE
          BASE=$(echo $FILE | sed 's/.zip$//')
          ORG=$(echo $BASE | tr '.' ' ' | awk '{print $NF}')
          SRR=$(echo $BASE | sed "s/.${ORG}//")
          scp -r $SRR mdz@Z620:~/bfx/dee/data/$ORG
          sudo rm -rf $SRR
          sudo rm $FILE
        else
          sudo rm $FILE
        fi
      else
        echo "wait for $FILE $TIME $FILETIME "
      fi
    fi
  done
fi

  rm started
fi
