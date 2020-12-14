#!/bin/bash
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

DIR=/var/www/html/metadata/

#QUERY_STRING="org=ecoli&accessionsearch=SRP253578&emailaddress=mark.ziemann%40gmail.com"
ORG=$(echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2)
ACC=$(echo $QUERY_STRING | cut -d '&' -f2 | cut -d '=' -f2)
EMAIL=$(echo $QUERY_STRING | cut -d '&' -f3 | cut -d '=' -f2 | sed 's/%40/@/')

#echo $QUERY_STRING
echo Your request: "$ORG $ACC $EMAIL"
echo "<br>"
echo "<br>"

QUEUE=$DIR/${ORG}_srpqueue.txt

#Error handling if no input provided
if [ -z "$ACC" ] ; then
  echo 'No search term provided.'
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

#Accession number workflow
if [ -n "$ACC" ] ; then
  ACC=$(echo $ACC | tr '[:lower:]' '[:upper:]'  )
  CNT=$(egrep -wc "$ACC" $QUEUE)

  if [ $CNT -eq 0 ] ; then
    echo "SRA Project accession $ACC wasn't found in the queue. Please check that it exists and has not already been included in DEE2."
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;"></FORM>'
    exit
  fi

  if [ $CNT -gt 0 ]; then
    echo SRA project accession $ACC was found in the normal queue. It will be moved to the express queue!
    echo We will contact you by email when data processing is completed.
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
    echo "<br><br>"
    echo $ORG $ACC $EMAIL >> /var/www/html/request.txt
    exit
  fi

fi

rm -f $TMP
echo '</body>
</html>'


