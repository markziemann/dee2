#!/bin/bash
set -x
echo "Content-type: text/html"
echo ''
echo '<!DOCTYPE html>
<html xml:lang="en" lang="en">
<link rel='shortcut icon' type='image/x-icon' href='favicon.ico' />
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Digital Expression Explorer 2
</title>
<style type="text/css">
 body{
  -webkit-text-size-adjust:none;
  margin: 0;
  max-width:720px;
  line-height:1.6;
  font:normal normal normal 100% verdana,sans-serif;
  font-size:22px;
  color:#3a3a3a;
  padding:10px 10px
 }
 h1,h2,h3{
  line-height:1.2;
  color:#406fef ;
 }
 aside,p,ul{
  color:#3a3a3a;
 }
 table,th,tr,td{
   line-height:1.2;
   font-family:sans-serif;
   font-size:20px;
   background-color:transparent;
   padding: 2px;
  }

input[type=checkbox] {
    zoom: 1.5;
}


.tooltip {
    position: relative;
    display: inline-block;
    border-bottom: 1px dotted black;
}

.tooltip .tooltiptext {
    visibility: hidden;
    width: 380px;
    background-color: black;
    color: #fff;
    text-align: center;
    border-radius: 6px;
    padding: 5px 0;
    position: absolute;
    z-index: 1;
    top: -5px;
    left: 110%;
}

.tooltip .tooltiptext::after {
    content: "";
    position: absolute;
    top: 2%;
    right: 100%;
    margin-top: -5px;
    border-width: 5px;
    border-style: solid;
    border-color: transparent black transparent transparent;
}
.tooltip:hover .tooltiptext {
    visibility: visible;
}

</style>
</head>
<body>
'

DIR=/var/www/dee2.io/srpqueue

#SPEC=$(echo $QUERY_STRING  | egrep -c '(:|;|}|\{|\[|\]|\/|\\|\@|\<|\>)')
SPEC=$(echo $QUERY_STRING  | egrep -c '(:|;|}|\{|\[|\]|\/|\\|\@)')

if [ $SPEC -gt 0 ] ; then
  echo $QUERY_STRING
  echo "<br>"
  echo "<br>"
  echo "Avoid special characters"
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;"></FORM>'
  exit
fi
QUERY_STRING=$(echo $QUERY_STRING | tr -d ':;{}()[]\/<>' )

ORG=$(echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2)
ACC=$(echo $QUERY_STRING | cut -d '&' -f2 | cut -d '=' -f2)
EMAIL=$(echo $QUERY_STRING | cut -d '&' -f3 | cut -d '=' -f2 | sed 's/%40/@/')

#echo $QUERY_STRING
echo Your request: "$ORG $ACC $EMAIL"
echo "<br>"
echo "<br>"

#Error handling if no input provided
if [ -z "$ACC" ] ; then
  echo 'No search term provided.'
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

# check previously analysed data
PRESENT=$(find /dee2_data/huge/hsapiens/ | grep -c ${ACC}_)
if [ $PRESENT -eq 1 ] ; then
  FILEPATH=$(find /dee2_data/huge/hsapiens/ | grep ${ACC}_)
  DATALINK=$(echo $FILEPATH | sed 's#/dee2_data/#https://dee2.io/#')
  echo "Data set is already available at the URL $DATALINK"
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;"></FORM>'
else
  echo "Data set not found on our server. Let's see if it is in the queue."
  echo "<br>"
  echo "<br>"
fi

# check local queue for metadata
SRPQUEUE=/dee2_data/srpqueue/${ORG}_srpqueue.txt
QHIT=$(grep -wc $ACC $SRPQUEUE)
if [ $QHIT -eq 1 ] ; then
  echo $ACC was listed in the queue.
  echo "<br>"
else
  echo $ACC was not found in the queue.
  echo "<br>"
fi

echo "Fetching and checking $ACC metadata with pysradb."
echo "<br>"
echo "<br>"

# check SRA for metadata
TMPFILE=$(mktemp)
export PATH="$PATH:/var/www/.local/bin"
pysradb metadata $ACC > $TMPFILE
NLINES=$(wc -l < $TMPFILE)
echo The metadata was found and has $NLINES lines.
echo "<br>"
echo "<br>"

if [ $NLINES -lt 2 ] ; then
  echo "Error: SRA project not found with pysradb! Check your accession numbers"
  echo "<br>"
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;"></FORM>'
  exit
fi

RUNS=$(cut -f22 $TMPFILE | sed 1d | paste -s -d ',')
echo RUNS:$RUNS
echo "<br>"
echo "<br>"

ORG2=$(cut -f7 $TMPFILE | head -2 | tail -1)
L1=$(echo $ORG2 | cut -c1 | tr  '[:upper:]' '[:lower:]')
W2=$(echo $ORG2 | cut -d ' ' -f2)
ORG2=$L1$W2

if [ "$ORG2" == "$ORG" ] ; then
  echo "Organism species checked okay."
  echo "<br>"
  echo "<br>"
else
  echo "Error: Organism submitted doesn't match that on SRA."
  echo "<br>"
  echo "<br>"
  echo ORG2:$ORG2
  echo "<br>"
  echo "<br>"
  echo ORG:$ORG
  echo "<br>"
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;"></FORM>'
  exit
fi

mv $TMPFILE /home/ubuntu/dee2_data/newrequests
TMPFILE=$(basename $TMPFILE)
TMPFILE=/home/ubuntu/dee2_data/newrequests/$TMPFILE
echo $EMAIL >> $TMPFILE
chmod 664 $TMPFILE

echo "We've sent an email to your nominated address.
Please follow the information in there to initiate the data processing."
echo "<br>"
echo "<br>"

#echo "BODY:test2" | mailx -s "SUBJECT:test2" $EMAIL && echo "Email sent!"

sudo -u www-data sh -c 'echo "BODY:test3" | mailx -s "SUBJECT:test3" mark.ziemann@gmail.com' \
  && echo "Email sent!" \
  || echo "Email failed!"

echo '</body>
</html>'
