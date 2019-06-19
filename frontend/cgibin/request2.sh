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
  max-width:1000px;
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

CNT=`echo "$QUERY_STRING" | grep -c '&'`
if [ $CNT -eq "0" ] ; then
  echo "Content-type: text/html"
  echo ""
  echo No datasets selected
  echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;"></FORM>'
  exit
fi


#ID=$(< /dev/urandom tr -dc A-Za-z0-9 | head -c${1:-32};echo;)
ID=${RANDOM}_${RANDOM}
USRDIR=/mnt/tmp/$ID/
mkdir -p $USRDIR
LOGDIR=$USRDIR/logs
mkdir -p $LOGDIR
ORG=`echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2 | tr 'A-Z' 'a-z'`
DATA_DIR=/dee2_data/data/${ORG}/
ACCESSIONS=/mnt/dee2_data/metadata/${ORG}_metadata.tsv.cut
METADATA=/mnt/dee2_data/metadata/${ORG}_metadata.tsv

cd $DATA_DIR

#Format query string to be compatible with egrep multi-query search
QS=`echo $QUERY_STRING | cut -d '&' -f2- | sed 's/DataSetList=on&//' \
| sed 's/x=/\ /g' | sed 's/ //' | tr -d '&' | sed 's/|//' \
| sed 's/ /_gene.cnt /g' | sed 's/$/_gene.cnt/'`

##remove SRRs with files absent
QS2=`echo $QS | tr ' ' '\n' | cut -d '_' -f1 | tr '\n' ' '`

##get the 3 digit species prefix
PFX=$(echo $ORG | cut -c-3)

#################################################
# Gene info - names and length
#################################################
cp ${PFX}_gene_info.tsv $USRDIR/GeneInfo.tsv

#################################################
# STAR Gene counts
#################################################
GENECOUNTS=$(echo $QS2 | tr ' ' '\n' | awk '{print $1"/"$1"_gene.cnt"}' | tr '\n' ' ' | sed 's/$/\n/')
#add gene symbol
cut -f-2 $USRDIR/GeneInfo.tsv | tr '\t' '_' \
| paste - $GENECOUNTS > $USRDIR/GeneCountMatrix.tsv

cp $LOGS $LOGDIR
cd $USRDIR

TOKEN=$(curl -b cookies.txt -c cookies.txt 'http://degust.erc.monash.edu/upload' | sort | grep 'name="authenticity_token"' | cut -d '"' -f8)
#echo "TOKEN $TOKEN"
#echo "<br>"
curl -b cookies.txt -v -L 'http://degust.erc.monash.edu/upload' -F "authenticity_token=$TOKEN" -F "filename=@/$PWD/GeneCountMatrix.tsv" 2> tmp \
| grep -v 'front-loader.gif'

URL=$(grep -m1 'http://degust.erc.monash.edu/degust/compare.html?code=' tmp | cut -d ' ' -f3)
echo "Sending your GeneCountMatrix.tsv to Degust."
echo "Here is your Degust <a href=$URL target=\"_blank\"> link</a>. It will open a new tab."
echo "<br><br>"
echo '<FORM><INPUT Type="button" VALUE="Go back to DEE2" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'


