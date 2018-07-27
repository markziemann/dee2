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
</style>
</head>
<body>
'




#save some awk functions for tabulation
tbl(){
awk 'BEGIN { print "<table border="1"><tr><th> <input type=\"checkbox\" name=\"DataSetList\" onClick=\"toggle(this)\" />Select all</th><th> SRA run accession</th><th> QC summary </th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>SRA submission accession</th><th>GEO series accession</th><th>GEO sample accession</th></tr>" }
	{ print "<tr><td> <input type='checkbox' name='x' value="$1">  </td><td>  <a href=http://www.ncbi.nlm.nih.gov/sra/?term="$1" target=\"_blank\"   >"$1"  </a>  </td><td>  <a href=/qc/"$1".qc > "$2" </a> </td><td>" $3 "</td><td>" $4 "</td><td>" $5 "</td><td>" $6 "</td><td>" $7 "</td><td>" $8 "</td></tr>" }
     END   { print "</table>" }'
}
export -f tbl

tbl1(){
awk 'BEGIN { print "<table border="1"><tr> <th> SRA run accession </th><th> QC summary </th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>SRA submission accession</th><th>GEO series accession</th><th>GEO sample accession</th></tr>" }
        { print "<tr><td> <a href=http://www.ncbi.nlm.nih.gov/sra/?term="$1" target=\"_blank\"  >"$1"</a>  </td><td>" $2 "</td><td>" $3 "</td><td>" $4 "</td><td>" $5 "</td><td>" $6 "</td><td>" $7 "</td><td>" $8 "</td></tr>" }
     END   { print "</table>" }'
}
export -f tbl1

tbl2(){
awk ' {OFS="\t";FS="\t"} BEGIN { print "<table border="1"><tr><th> <input type=\"checkbox\" name=\"DataSetList\" onClick=\"toggle(this)\" />Select all</th><th> SRA run accession </th><th> QC summary </th><th>Keyword context</th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>SRA submission accession</th><th>GEO series accession</th><th>GEO sample accession</th></tr>" }
        { print "<tr><td> <input type='checkbox' name='x' value="$1">  </td><td>  <a href=http://www.ncbi.nlm.nih.gov/sra/"$1" target=\"_blank\" >"$1"</a> </td><td> <a href=/qc/"$1".qc > "$3" </a> </td><td>..." $2 "...</td><td>" $4 "</td><td>" $5 "</td><td>" $6 "</td><td>" $7 "</td><td>" $8 "</td><td>"$9 "</td></tr>" }
     END   { print "</table>" }'
}
export -f tbl2

tbl3(){
awk ' {OFS="\t";FS="\t"} BEGIN { print "<table border="1"><tr><th> SRA run accession </th><th> QC summary </th><th>Keyword context</th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>SRA submission accession</th><th>GEO series accession</th><th>GEO sample accession</th></tr>" }
	{ print "<tr><td>  <a href=http://www.ncbi.nlm.nih.gov/sra/"$1" target=\"_blank\" >"$1"  </a>  </td><td> <a href=/qc/"$1".qc > "$3" </a>  </td><td>..." $2 "...</td><td>" $4 "</td><td>" $5 "</td><td>" $6 "</td><td>" $7 "</td><td>" $8 "</td><td>" $9 "</td></tr>" }
     END   { print "</table>" }'
}
export -f tbl3


#DIR=/var/www/metadata
DIR=/var/www/html/metadata/

#QUERY_STRING="org=ecoli&accessionsearch=&keywordsearch=chaperone"
#QUERY_STRING="org=ecoli&accessionsearch=GSE33671&keywordsearch="

ORG=`echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2`
ACC=`echo $QUERY_STRING | cut -d '&' -f2 | cut -d '=' -f2`
KEY=`echo $QUERY_STRING | cut -d '&' -f3 | cut -d '=' -f2`

#echo $QUERY_STRING
#echo "<br>"
#echo "$ORG $ACC $KEY"
#echo "<br>"

MD=$DIR/${ORG}_metadata.tsv
MDCUT=$DIR/${ORG}_metadata.tsv.cut

if [ -z "$ACC" -a -z "$KEY" ] ; then
  echo 'No search terms provided.'
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

if [ -n "$ACC" -a -n "$KEY" ] ; then
  echo 'Please enter an accession number OR keyword, not both.'
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

if [ -n "$ACC" -a -z "$KEY" ] ; then
  echo "<script type=\"text/javascript\"> function toggle(source) { checkboxes = document.getElementsByName('x'); for(var i=0, n=checkboxes.length;i<n;i++) { checkboxes[i].checked = source.checked; } } </script>"
  echo '<form action="request.sh" method="get">'
  echo '<input type="hidden" name="org" value="ORG">' | sed "s/ORG/${ORG}/"

  Q=`echo $ACC | sed 's/\%2C/\|/g' | sed 's/^/\(/' | sed 's/$/\)/'`

  if [ $CNT -eq 0 ] ; then
    echo No results found
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    exit
  fi

  if [ $CNT -gt 3000 ]; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter accession number search, or consider a bulk data download.
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
    exit
  fi

  if [ $CNT -gt 500 ] ; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter accession number search, or consider a bulk data download.
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
    #display all results
    cut -f-8 $MD | egrep -w "$Q" | sort -k1 | tbl1
    exit
  fi

  cut -f-8 $MD | egrep -w "$Q" | sort -k1 | tbl
  echo '<input type="submit" value="Get Counts" class="tfbutton" style="font-size : 22px;" >'
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  echo Please hit the submit button just once. Retrieval time is about 1 dataset per second.
  exit
fi

if [ -n "$KEY" -a -z "$ACC" ] ; then
  echo "<script type=\"text/javascript\"> function toggle(source) { checkboxes = document.getElementsByName('x'); for(var i=0, n=checkboxes.length;i<n;i++) { checkboxes[i].checked = source.checked; } } </script>"
  echo '<form action="request.sh" method="get">'
  echo '<input type="hidden" name="org" value="ORG">' | sed "s/ORG/${ORG}/"


  Q=`echo $KEY | tr '+' ' '`
  #echo $Q
  STR=`< /dev/urandom tr -dc A-Za-z0-9 | head -c${1:-32};echo;`
  TMP=/tmp/TMP_${STR}.txt
  CNT=`egrep -ic "$Q" $MD `

  if [ $CNT -eq 0 ]; then
    echo No results found
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size:22px;" ></FORM>'
    exit
  fi

  if [ $CNT -gt 3000 ]; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter keyword or accession number search, or consider a bulk data download.
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size:22px;" ></FORM>'
    exit
  fi

  egrep -i "$Q" $MD | tr '\t' '\n' | egrep -i "$Q" $MD > $TMP

  if [ $CNT -gt 500 ]; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter keyword or accession number search, or consider a bulk data download.
    echo "<br>"
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    #display all results
    sed "s/${Q}/x@x/I" $TMP | egrep -o ".{0,30}x@x.{0,30}" | sed "s/x@x/${Q}/" | tr '\t ' '_' \
    | paste - $TMP | cut -f-9 \
    | tr -d ' ' | awk '{FS="\t";OFS="\t"} {print $2,$1,$3,$4,$5,$6,$7,$8,$9}' \
    | sort -k1 | tbl3
    exit
  fi

  echo Use the checkboxes to select data sets of interest.
  sed "s/${Q}/x@x/I" $TMP | egrep -o ".{0,30}x@x.{0,30}" | sed "s/x@x/${Q}/" | tr '\t ' '_' \
  | paste - $TMP | cut -f-9 \
  | tr -d ' ' | awk '{FS="\t";OFS="\t"} {print $2,$1,$3,$4,$5,$6,$7,$8,$9}' \
  | sort -k1 | tbl2
  echo '<input type="submit" value="Get Counts" class="tfbutton" style="font-size:22px;" >'
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"  style="font-size : 22px;"  ></FORM>'
  echo Please hit the submit button just once. Retrieval time is about 1 dataset per second.
  rm -f $TMP
fi
rm -f $TMP
echo '</body>
</html>'

