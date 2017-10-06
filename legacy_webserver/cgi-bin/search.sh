#!/bin/bash
echo "Content-type: text/html"
echo ''
#echo 'Welcome to The RNA Expression Compendium'

DIR=/var/www/metadata

#QUERY_STRING="org=ecoli&accessionsearch=&keywordsearch=pathogenicity+island"

ORG=`echo $QUERY_STRING | cut -d '&' -f1 | cut -d '=' -f2`
ACC=`echo $QUERY_STRING | cut -d '&' -f2 | cut -d '=' -f2`
KEY=`echo $QUERY_STRING | cut -d '&' -f3 | cut -d '=' -f2`

#echo "$ORG $ACC $KEY"

MD=/var/www/metadata/${ORG}_metadata.txt
MDCUT=/var/www/metadata/${ORG}_metadata.txt.cut

if [ -z "$ACC" -a -z "$KEY" ] ; then
  echo 'No search terms provided.'
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
  exit
fi

if [ -n "$ACC" -a -n "$KEY" ] ; then
  echo 'Please enter an accession number OR keyword, not both.'
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
  exit
fi

if [ -n "$ACC" -a -z "$KEY" ] ; then
  echo "<script type=\"text/javascript\"> function toggle(source) { checkboxes = document.getElementsByName('x'); for(var i=0, n=checkboxes.length;i<n;i++) { checkboxes[i].checked = source.checked; } } </script>"
  echo '<form action="extract.sh" method="get">'
  echo '<input type="hidden" name="org" value="ORG">' | sed "s/ORG/${ORG}/"

  Q=`echo $ACC | sed 's/\%2C/\|/g' | sed 's/^/\(/' | sed 's/$/\)/'`

  if [ -r "$MDCUT" ] ; then
    CNT=`egrep -cw "$Q" $MDCUT`
  else
    CNT=`cut -f-8 $MD | tee $MDCUT | egrep -cw "$Q"`
  fi

  if [ $CNT -eq 0 ] ; then
    echo No results found
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    exit
  fi

  if [ $CNT -gt 3000 ]; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter accession number search, or consider a bulk data download.
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    #don't display all results
    #cut -f-7 $MD | egrep -w "$Q" | sort -k1 | ./tbl1.awk
    exit
  fi

  if [ $CNT -gt 500 ] ; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter accession number search, or consider a bulk data download.
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    #display all results
    cut -f-8 $MD | egrep -w "$Q" | sort -k1 | ./tbl1.awk
    exit
  fi

  cut -f-8 $MD | egrep -w "$Q" | sort -k1 | ./tbl.awk
  echo '<input type="submit" value="Get Counts" class="tfbutton">'
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
  exit
fi

if [ -n "$KEY" -a -z "$ACC" ] ; then
  echo "<script type=\"text/javascript\"> function toggle(source) { checkboxes = document.getElementsByName('x'); for(var i=0, n=checkboxes.length;i<n;i++) { checkboxes[i].checked = source.checked; } } </script>"
  echo '<form action="extract.sh" method="get">'
  echo '<input type="hidden" name="org" value="ORG">' | sed "s/ORG/${ORG}/"

  Q=`echo $KEY | tr '+' ' '`
  #echo $Q
  STR=`< /dev/urandom tr -dc A-Za-z0-9 | head -c${1:-32};echo;`
  TMP=/tmp/TMP_${STR}.txt
  CNT=`egrep -ic "$Q" $MD `

  if [ $CNT -eq 0 ]; then
    echo No results found
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    exit
  fi

  if [ $CNT -gt 3000 ]; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter keyword or accession number search, or consider a bulk data download.
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    #display all results
    #sed "s/${Q}/x@x/I" $TMP | egrep -o ".{0,30}x@x.{0,30}" | sed "s/x@x/${Q}/" | tr '\t ' '_' \
    #| paste - $TMP | cut -f-8 \
    #| tr -d ' ' | awk '{FS="\t";OFS="\t"} {print $2,$1,$3,$4,$5,$6,$7,$8}' \
    #| sort -k1 | ./tbl3.awk
    exit
  fi

  egrep -i "$Q" $MD | tr '\t' '\n' | egrep -i "$Q" $MD > $TMP

  if [ $CNT -gt 500 ]; then
    echo Too many results found \(${CNT}\). The webserver is limited to 500 datasets per search. \
    Try a stricter keyword or accession number search, or consider a bulk data download.
    echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
    #display all results
    sed "s/${Q}/x@x/I" $TMP | egrep -o ".{0,30}x@x.{0,30}" | sed "s/x@x/${Q}/" | tr '\t ' '_' \
    | paste - $TMP | cut -f-9 \
    | tr -d ' ' | awk '{FS="\t";OFS="\t"} {print $2,$1,$3,$4,$5,$6,$7,$8,$9}' \
    | sort -k1 | ./tbl3.awk
    exit
  fi

  sed "s/${Q}/x@x/I" $TMP | egrep -o ".{0,30}x@x.{0,30}" | sed "s/x@x/${Q}/" | tr '\t ' '_' \
  | paste - $TMP | cut -f-9 \
  | tr -d ' ' | awk '{FS="\t";OFS="\t"} {print $2,$1,$3,$4,$5,$6,$7,$8,$9}' \
  | sort -k1 | ./tbl2.awk
  echo '<input type="submit" value="Get Counts" class="tfbutton">'
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;"></FORM>'
  rm -f $TMP
fi
rm -f $TMP

