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

#QUERY_STRING="accession=SRP514223&confirmcode=xt26"

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

ACC=$(echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2)
COMBO=$(echo $QUERY_STRING | cut -d '&' -f2 | cut -d '=' -f2)

#echo $QUERY_STRING
echo Your request: "Accession:$ACC Combination:$COMBO"
echo "<br>"
echo "<br>"

#Error handling if no input provided
if [ -z "$ACC" ] ; then
  echo 'Error: An SRA project accession (SRP) is required.'
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Go Back" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

#Error handling if no input provided
if [ -z "$COMBO" ] ; then
  echo 'Error: A confirmation combination is required.'
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Go Back" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

# check that the SRP has been requested
METADATAFILE=/usr/lib/cgi-bin/newrequests/$ACC
if [ ! -r "$METADATAFILE" ] ; then
  echo "Error: We don't seem to have any record of a request for that SRA project accession (SRP). Would you like to request it? If so follow the <a href=newrequest2.html>link</a>."
  echo "<br>"
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Go Back" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

# check that the combos match
FILECOMBO=$(grep "COMBO=" "$METADATAFILE" | cut -d '=' -f2)
if [ "$FILECOMBO" != "$COMBO" ] ; then
  echo "Error: It seems the combinations don't match what we sent by email. Check the spelling and try again."
  echo "<br>"
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Go Back" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

#check the number of requests in the last 370 days
EMAIL_ADDRESS=$(grep "EMAIL=" "$METADATAFILE" | cut -d '=' -f2)
USERFILE="/usr/lib/cgi-bin/newrequests/$EMAIL_ADDRESS"
CURRENT_DATE=$(date +%s | awk '{print $1/86400}' | cut -d '.' -f1)
LAST_YEAR_DATE=$(echo "$CURRENT_DATE" | awk '{print $1 - 370}')
echo "$CURRENT_DATE" "$ACC" >> "$USERFILE"
USER_REQUEST_COUNT_IN_LAST_YEAR=$(awk -v d=$LAST_YEAR_DATE '$1>=d {print $1}' "$USERFILE" | wc -l )
echo "Number of requests in past year:$USER_REQUEST_COUNT_IN_LAST_YEAR"
echo "<br>"
echo "<br>"

if [ "$USER_REQUEST_COUNT_IN_LAST_YEAR" -le 10 ] ; then
  echo "Starting the analysis!"
  echo "<br>"
  echo "<br>"
else
  echo "Sorry, can't start the analysis, user has requested >10 datasets alredy this year."
  echo "<br>"
  echo "<br>"
fi

echo '</body>
</html>'

echo "CONFIRMED=TRUE" >> $METADATAFILE
