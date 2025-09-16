#!/bin/bash
#set -x
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

#SPEC=$(echo $QUERY_STRING  | egrep -c '(:|;|}|\{|\[|\]|\/|\\|\@)')
SPEC=$(echo "$QUERY_STRING"  | egrep -c '(\!|\@|\#|\$|\^|\*|\:|\;|\}|\{|\[|\]|\/|\\|\@)')
if [ "$SPEC" -gt 0 ] ; then
  echo "$QUERY_STRING"
  echo "<br>"
  echo "Avoid special characters"
  exit
fi
QUERY_STRING=$(echo "$QUERY_STRING" | tr -d ':;{}()[]\/<>' )

STRING_LEN=$(echo "$QUERY_STRING" | wc -c)
if [ "$STRING_LEN" -gt 200 ] ; then
  echo "$QUERY_STRING"
  echo "<br>"
  echo "Please keep queries to less than 200 characters."
  exit
fi

#QUERY_STRING="org=bdistachyon&accessionsearch=&keywordsearch=abiotic"
#QUERY_STRING="org=scerevisiae&accessionsearch=&keywordsearch=metaboli"

ORG=`echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2`
ACC=`echo $QUERY_STRING | cut -d '&' -f2 | cut -d '=' -f2`
KEY=`echo $QUERY_STRING | cut -d '&' -f3 | cut -d '=' -f2`

#echo $QUERY_STRING
#echo "<br>"
#echo "$ORG $ACC $KEY"
#echo "<br>"

MD=$DIR/${ORG}_srp.tsv
MDCUT=$DIR/${ORG}_srp.tsv.cut

#save some awk functions for tabulation
# note that the shell functions are used when less than 500 results as the tooltip can be used
tblx(){
echo "<table border="1"><tr><th>SRA Project Accession</th><th> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\">QC summary </a> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\" target=\"_blank\"> <img src=\"/images/question.png\" alt=\"alttext\" title=\"Learn more about the quality metrics\" style=\"width:30px;height:30px;\"> </a> </th><th>Project Title</th><th>Project Description</th><th>Study Type</th><th>GEO Series</th></tr>"
while read line ; do
  C1=$(echo "$line" | cut -f1 | tr -d '"')
  ZIPURL=https://dee2.io/$(find /dee2_data/huge/$ORG | cut -d '/' -f3- | grep ${C1}_ )
  C3=$(echo "$line" | cut -f3 | tr -d '"' )
  C4=$(echo "$line" | cut -f4 | tr -d '"' )
  C6=$(echo "$line" | cut -f6 | tr -d '"' )
  C9=$(echo "$line" | cut -f9 | tr -d '"' )
  C10=$(echo "$line" | cut -f10 | tr -d '"' )
  echo "<tr><td>" $C1 "<br><br> <a href=$ZIPURL target=_blank >" DEE2 data bundle link "</a> <br><br> <a href=https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=$C1 target=_blank >"SRA link" </a> </td><td>" Put $C1 QC here "</td><td>" $C3 "</td><td>" $C4 "</td><td>" $C6"</td><td>" $C10 "</td></tr>"
done
}
export -f tblx

tbl1x(){
echo "<table border=1><tr> <th> SRA run accession </th><th> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\">QC summary </a> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\" target=\"_blank\"> <img src=\"/images/question.png\" alt=\"alttext\" title=\"Learn more about the quality metrics\" style=\"width:30px;height:30px;\"> </a> </th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>Sample name / GEO sample accession</th><th>GEO series accession</th><th>Experiment name</th></tr>"
while read line ; do
  C1=$(echo "$line" | cut -f1)
  C2=$(echo "$line" | cut -f2)
  C3=$(echo "$line" | cut -f3)
  C4=$(echo "$line" | cut -f4)
  C5=$(echo "$line" | cut -f5)
  C6=$(echo "$line" | cut -f6)
  C7=$(echo "$line" | cut -f7)
  C8=$(echo "$line" | cut -f8)
  C9=$(echo "$line" | cut -f9)
  C10=$(echo "$line" | cut -f10)
  C11=$(echo "$line" | cut -f11)

  echo "<tr><td> <a href=http://www.ncbi.nlm.nih.gov/sra/?term="$C1" target=\"_blank\"  >"$C1"</a>  </td><td> <a href=/data/"$ORG"/"$C1"/"$C1".qc target=_blank > <div class=tooltip>"$C2"<span class=tooltiptext > $(cat /dee2_data/data/"$ORG"/"$C1"/"$C1".qc) </span> </div> </a>          </td><td>" $C3 "</td><td>" $C4 "</td><td>" $C5 "</td><td>" $C6 "</td><td>" $C7 "</td><td>" $C8 "</td></tr>"
done
}
export -f tbl1x

tbl2(){
awk -v o=$ORG ' {OFS="\t";FS="\t"} BEGIN { print "<table border="1"><tr><th> <input type=\"checkbox\" name=\"DataSetList\" onClick=\"toggle(this)\" />Select all</th><th> SRA run accession </th><th> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\">QC summary </a> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\" target=\"_blank\"> <img src=\"/images/question.png\" alt=\"alttext\" title=\"Learn more about the quality metrics\" style=\"width:30px;height:30px;\"> </a> </th><th>Keyword context</th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>Sample name / GEO sample accession</th><th>GEO series accession</th><th>Experiment name</th></tr>" }
        { print "<tr><td> <input type='checkbox' name='x' value="$1">  </td><td>  <a href=http://www.ncbi.nlm.nih.gov/sra/"$1" target=\"_blank\" >"$1"</a> </td><td> <a href=/data/"o"/"$1"/"$1".qc  target=\"_blank\" > "$3" </a> </td><td>..." $2 "...</td><td>" $4 "</td><td>" $5 "</td><td>" $6 "</td><td>" $7 "</td><td>" $8  "</td></tr>" }
     END   { print "</table>" }'
}
export -f tbl2

tbl2x(){
echo "<table border=1> <tr><th> <input type=checkbox name=DataSetList onClick=\"toggle(this)\" />Select all</th><th> SRA run accession </th><th> Keyword context </th><th> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\">QC summary </a> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\" target=\"_blank\"> <img src=\"/images/question.png\" alt=\"alttext\" title=\"Learn more about the quality metrics\" style=\"width:30px;height:30px;\"> </a> </th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>Sample name / GEO sample accession</th><th>GEO series accession</th><th>Experiment.title</th></tr>"
while read line ; do
  C1=$(echo "$line" | cut -f1)
  C2=$(echo "$line" | cut -f2)
  C3=$(echo "$line" | cut -f3)
  C4=$(echo "$line" | cut -f4)
  C5=$(echo "$line" | cut -f5)
  C6=$(echo "$line" | cut -f6)
  C7=$(echo "$line" | cut -f7)
  C8=$(echo "$line" | cut -f8)
  C9=$(echo "$line" | cut -f9)
  echo "<tr><td> <input type=checkbox name=x value="$C1"> </td><td> <a href=http://www.ncbi.nlm.nih.gov/sra/"$C1" target=_blank >"$C1"</a> </td><td> ..."$C2"...</td><td> <a href=/data/"$ORG"/"$C1"/"$C1".qc target=_blank > <div class=tooltip>"$C3"<span class=tooltiptext > $(cat /dee2_data/data/"$ORG"/"$C1"/"$C1".qc) </span> </div> </a> </td><td>"$C4"</td><td>"$C5"</td><td>"$C6"</td><td>"$C7"</td><td>"$C8"</td><td>"$C9"</td></tr>"
done
}
export -f tbl2x

tbl3(){
awk -F'\t' -v o=$ORG '{OFS=FS} BEGIN { print "<table border="1"><tr><th> SRA run accession </th><th>Keyword context </th><th> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\">QC summary </a> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\" target=\"_blank\"> <img src=\"/images/question.png\" alt=\"alttext\" title=\"Learn more about the quality metrics\" style=\"width:30px;height:30px;\"> </a> </th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>Sample name / GEO series accession</th><th>GEO sample accession </th><th>Experiment.title </th></tr>" }
	{ print "<tr><td>  <a href=http://www.ncbi.nlm.nih.gov/sra/"$1" target=\"_blank\" >"$1"  </a>  </td><td> ..." $2 "... </td><td> <a href=/data/"o"/"$1"/"$1".qc  target=\"_blank\" > "$3" </a></td><td>" $4 "</td><td>" $5 "</td><td>" $6 "</td><td>" $7 "</td><td>" $8 "</td><td>" $9   "</td></tr>" }
     END   { print "</table>" }'
}
export -f tbl3

tbl3x(){
echo "<table border=1><tr><th> SRA run accession </th><th> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\">QC summary </a> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\" target=\"_blank\"> <img src=\"/images/question.png\" alt=\"alttext\" title=\"Learn more about the quality metrics\" style=\"width:30px;height:30px;\"> </a> </th><th>Keyword context</th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>Sample name / GEO sample accession</th><th>GEO series accession  </th><th>Experiment name</th></tr>"
while read line ; do
  C1=$(echo "$line" | cut -f1)
  C2=$(echo "$line" | cut -f2)
  C3=$(echo "$line" | cut -f3)
  C4=$(echo "$line" | cut -f4)
  C5=$(echo "$line" | cut -f5)
  C6=$(echo "$line" | cut -f6)
  C7=$(echo "$line" | cut -f7)
  C8=$(echo "$line" | cut -f8)
  C9=$(echo "$line" | cut -f9)
  echo "<tr><td>  <a href=http://www.ncbi.nlm.nih.gov/sra/"$C1" target=_blank >"$C1"</a> \
  </td><td> <a href=/data/"o"/"$C1"/"$C1".qc target=_blank> <div class=tooltip>"$C3"<span class=tooltiptext > $(cat /dee2_data/data/"$ORG"/"$C1"/"$C1".qc) </span></div></a> \
  </td><td>..."$C2"...</td><td>"$C4"</td><td>"$C5"</td><td>"$C6"</td><td>"$C7"</td><td>"$C8"</td><td>"$C9"</td></tr>"
done
}
export -f tbl3x

#Error handling if no input provided
if [ -z "$ACC" -a -z "$KEY" ] ; then
  echo 'No search terms provided.'
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

#Error handling if both keyword and accession provided
if [ -n "$ACC" -a -n "$KEY" ] ; then
  echo 'Please enter an accession number OR keyword, not both.'
  echo "<br>"
+  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

#Accession number workflow
if [ -n "$ACC" -a -z "$KEY" ] ; then
  Q=$(echo $ACC | sed 's/\%2C/\|/g' | sed 's/^/\(/' | sed 's/$/\)/' | tr '+' ' ')
  CNT=$(cut -f-9 $MD | awk '!arr[$1]++' | egrep -iwc "$Q")

  echo "<script type=\"text/javascript\"> function toggle(source) { checkboxes = document.getElementsByName('x'); for(var i=0, n=checkboxes.length;i<n;i++) { checkboxes[i].checked = source.checked; } } </script>"
  echo '<form action="request.sh" method="get">'
  echo '<input type="hidden" name="org" value="ORG">' | sed "s/ORG/${ORG}/"

  if [ $CNT -eq 0 ] ; then
    echo No results found
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    exit
  fi

  if [ $CNT -gt 3000 ]; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter accession number search, or consider a '<a href="/bulk">bulk data download</a>.'
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
    exit
  fi

  if [ $CNT -gt 500 ] ; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter accession number search, or a '<a href="/bulk">bulk data download</a>.'
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
    #display all results
    cut -f-8 $MD | egrep -iw "$Q" | awk '!arr[$1]++' | sort -k1 | tbl1
    echo '</table>'
    exit
  fi

  echo ${CNT} datasets found.
  cut -f-10 $MD | egrep -iw "$Q" | awk '!arr[$1]++' | sort -k1 | tblx
  echo '</table>'

  echo '<input type="submit" value="Get Counts" class="tfbutton" style="font-size : 22px;" >'
  echo '<br>'
  echo 'Send data to the Degust analysis tool: '
  echo '<input type="submit" formaction="request2.sh" method="get" value="Degust (STAR gene counts)" class="tfbutton" style="font-size : 22px;" >'
  echo '<input type="submit" formaction="request3.sh" method="get" value="Degust (Kallisto tx counts)" class="tfbutton" style="font-size : 22px;" >'
  echo '<input type="submit" formaction="request4.sh" method="get" value="Degust (Kallisto tx2gene counts)" class="tfbutton" style="font-size : 22px;" >'
  echo '<br>'
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  echo Please hit the submit button just once. Retrieval time is about 1 dataset per second.
  exit
fi

#keyword workflow
if [ -n "$KEY" -a -z "$ACC" ] ; then
  Q=$(echo $KEY | sed 's/\%2C/\|/g' | sed 's/^/\(/' | sed 's/$/\)/' | tr '+' ' ')
  CNT=$(cut -f-9 $MD | awk '!arr[$1]++' | egrep -ic "$Q")

  echo "<script type=\"text/javascript\"> function toggle(source) { checkboxes = document.getElementsByName('x'); for(var i=0, n=checkboxes.length;i<n;i++) { checkboxes[i].checked = source.checked; } } </script>"
  echo '<form action="request.sh" method="get">'
  echo '<input type="hidden" name="org" value="ORG">' | sed "s/ORG/${ORG}/"

  if [ $CNT -eq 0 ] ; then
    echo No results found
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    exit
  fi

  if [ $CNT -gt 3000 ]; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter accession number search, or consider a '<a href="/bulk">bulk data download</a>.'
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
    exit
  fi

  if [ $CNT -gt 500 ] ; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter accession number search, or a '<a href="/bulk">bulk data download</a>.'
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
    #display all results
    cut -f-8 $MD | egrep -i "$Q" | awk '!arr[$1]++' | sort -k1 | tbl1
    echo '</table>'
    exit
  fi

  echo ${CNT} datasets found.
  cut -f-10 $MD | egrep -i "$Q" | awk '!arr[$1]++' | sort -k1 | tblx
  echo '</table>'

  echo '<input type="submit" value="Get Counts" class="tfbutton" style="font-size : 22px;" >'
  echo '<br>'
  echo 'Send data to the Degust analysis tool: '
  echo '<input type="submit" formaction="request2.sh" method="get" value="Degust (STAR gene counts)" class="tfbutton" style="font-size : 22px;" >'
  echo '<input type="submit" formaction="request3.sh" method="get" value="Degust (Kallisto tx counts)" class="tfbutton" style="font-size : 22px;" >'
  echo '<input type="submit" formaction="request4.sh" method="get" value="Degust (Kallisto tx2gene counts)" class="tfbutton" style="font-size : 22px;" >'
  echo '<br>'
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  echo Please hit the submit button just once. Retrieval time is about 1 dataset per second.
  exit
fi
rm -f $TMP




echo '</body>
</html>'


