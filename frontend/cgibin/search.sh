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

SPEC=$(echo $QUERY_STRING  | egrep -c '(:|;|}|\{|\[|\]|\/|\\|\@|\<|\>)')
if [ $SPEC -gt 0 ] ; then
  echo $QUERY_STRING
  echo "<br>"
  echo "Avoid special characters"
  exit
fi
QUERY_STRING=$(echo $QUERY_STRING | tr -d ':;{}()[]\/<>' )

#QUERY_STRING="org=ecoli&accessionsearch=&keywordsearch=chaperone"
#QUERY_STRING="org=scerevisiae&accessionsearch=&keywordsearch=metaboli"

ORG=`echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2`
ACC=`echo $QUERY_STRING | cut -d '&' -f2 | cut -d '=' -f2`
KEY=`echo $QUERY_STRING | cut -d '&' -f3 | cut -d '=' -f2`

#echo $QUERY_STRING
#echo "<br>"
#echo "$ORG $ACC $KEY"
#echo "<br>"

MD=$DIR/${ORG}_metadata.tsv
MDCUT=$DIR/${ORG}_metadata.tsv.cut

#save some awk functions for tabulation
# note that the shell functions are used when less than 500 results as the tooltip can be used
tblx(){
echo "<table border="1"><tr><th> <input type=\"checkbox\" name=\"DataSetList\" onClick=\"toggle(this)\" />Select all</th><th> SRA run accession</th><th> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\">QC summary </a> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\" target=\"_blank\"> <img src=\"/images/question.png\" alt=\"alttext\" title=\"Learn more about the quality metrics\" style=\"width:30px;height:30px;\"> </a> </th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>Sample Name / GEO sample accession</th><th>GEO series accession</th><th>Experiment name</th></tr>"
while read line ; do
  C1=$(echo "$line" | cut -f1)
  C2=$(echo "$line" | cut -f2)
  C3=$(echo "$line" | cut -f3)
  C4=$(echo "$line" | cut -f4)
  C5=$(echo "$line" | cut -f5)
  C6=$(echo "$line" | cut -f6)
  C7=$(echo "$line" | cut -f7)
  C8=$(echo "$line" | cut -f8)
  echo "<tr><td> <input type=checkbox name=x value="$C1">  </td><td> <a href=http://www.ncbi.nlm.nih.gov/sra/?term="$C1" target=_blank >"$C1" </a> </td><td> <a href=/data/"$ORG"/"$C1"/"$C1".qc target=_blank > <div class=tooltip>"$C2"<span class=tooltiptext > $(cat /dee2_data/data/"$ORG"/"$C1"/"$C1".qc) </span> </div> </a> </td><td>" $C3 "</td><td>" $C4 "</td><td>" $C5 "</td><td>" $C6 "</td><td>" $C7 "</td><td>" $C8 "</td></tr>"
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
  echo "<tr><td> <a href=http://www.ncbi.nlm.nih.gov/sra/?term="$C1" target=\"_blank\"  >"$C1"</a>  </td><td> <a href=/data/"$ORG"/"$C1"/"$C1".qc target=_blank > <div class=tooltip>"$C2"<span class=tooltiptext > $(cat /dee2_data/data/"$ORG"/"$C1"/"$C1".qc) </span> </div> </a>  </td><td>" $C3 "</td><td>" $C4 "</td><td>" $C5 "</td><td>" $C6 "</td><td>" $C7 "</td><td>" $C8 "</td></tr>"
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
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

#Accession number workflow
if [ -n "$ACC" -a -z "$KEY" ] ; then
  Q=`echo $ACC | sed 's/\%2C/\|/g' | sed 's/^/\(/' | sed 's/$/\)/'`
  CNT=$(cut -f-9 $MD | awk '!arr[$1]++' | egrep -wc "$Q")

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
    cut -f-8 $MD | egrep -w "$Q" | awk '!arr[$1]++' | sort -k1 | tbl1
    echo '</table>'
    exit
  fi

  echo ${CNT} datasets found. Use the checkboxes to select ones of interest.
  cut -f-10 $MD | egrep -w "$Q" | awk '!arr[$1]++' | sort -k1 | tblx
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
  echo "<script type=\"text/javascript\"> function toggle(source) { checkboxes = document.getElementsByName('x'); for(var i=0, n=checkboxes.length;i<n;i++) { checkboxes[i].checked = source.checked; } } </script>"
  echo '<form action="request.sh" method="get">'
  echo '<input type="hidden" name="org" value="ORG">' | sed "s/ORG/${ORG}/"

  Q=`echo $KEY | tr '+' ' '`
  STR=`< /dev/urandom tr -dc A-Za-z0-9 | head -c${1:-32};echo;`
  TMP=/tmp/TMP_${STR}.txt
  TMP_TSV=/tmp/TMP_${STR}.tsv

  Q=$(echo $Q | tr '[:upper:]' '[:lower:]' | sed 's/ and / AND /')
  QQ=$(echo $Q | sed 's/ and /\n/gI' | tail -1)
  AWK_CMD=$(echo $Q | tr -s ' ' | sed 's#^#/#' | sed 's#AND#/ \&\& /#g' | sed 's#$#/#' )
  CMD="awk ' tolower(\$0) ~ "${AWK_CMD}" '"
  CNT=$(eval $CMD $MD | awk '!arr[$1]++'| wc -l)

  if [ $CNT -eq 0 ]; then
    echo No results found
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size:22px;" ></FORM>'
    exit
  fi

  if [ $CNT -gt 5000 ]; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter keyword or accession number search, or a '<a href="/bulk">bulk download</a>.' "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size:22px;" ></FORM>'
    exit
  fi

  AWK_CMD=$(echo $Q | tr -s ' ' | sed 's#^#/#' | sed 's#AND#/ \&\& /#g' | sed 's#$#/#' )
  CMD="awk ' tolower(\$0) ~ "${AWK_CMD}" '"
  eval $CMD $MD > $TMP

  if [ $CNT -gt 500 ]; then
    HL=$(echo "<a href=\"${TMP_TSV}\"> here</a>")
    echo "SRA_run_accession Keyword_context QC_summary SRA_experiment_accession SRA_sample_accession SRA_project_accession GEO_series_accession GEO_sample_accession Experiment.title" | tr ' ' '\t' >$TMP_TSV
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter keyword or accession number search, or a '<a href="/bulk">bulk download</a>.' \
    Click $HL to download the search results as a tab separated file \(TSV\). "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    #display all results
    sed "s#${QQ}#x@x#I" $TMP |  egrep -io ".{0,30}x@x.{0,30}" | tr '\t' ' ' | sed "s/x@x/${QQ}/g" \
    | tr '\t' ' ' | paste - $TMP | cut -f-9 \
    | awk -F'\t' 'BEGIN{OFS=FS} {print $2,$1,$3,$4,$5,$6,$7,$8,$9}' \
    | sort -k1 | awk '!arr[$1]++'| tee -a $TMP_TSV | tbl3
    exit
  fi

  echo ${CNT} datasets found. Use the checkboxes to select ones of interest.
  sed "s#${QQ}#x@x#I" $TMP |  egrep -io ".{0,30}x@x.{0,30}" | tr '\t' ' ' | sed "s/x@x/${QQ}/g" \
  | tr '\t' ' ' | paste - $TMP | cut -f-9 \
  | awk -F'\t' 'BEGIN{OFS=FS} ;{print $2,$1,$3,$4,$5,$6,$7,$8,$9}' \
  | sort -k1 | awk '!arr[$1]++' | tbl2x
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
  rm -f $TMP

fi
rm -f $TMP
echo '</body>
</html>'


