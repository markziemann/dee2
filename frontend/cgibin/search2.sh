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

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            color: #e0e0e0;
            background: #190136;
            padding: 0;
            overflow-x: hidden;
            transition: background 0.3s, color 0.3s;
        }

        body.light-mode {
            background: #f5f5f5;
            color: #333;
        }

        .container {
            max-width: 900px;
            margin: 0 auto;
            background: #0f0f0f;
            position: relative;
            transition: background 0.3s;
        }

        body.light-mode .container {
            background: #ffffff;
        }

        .navbar {
            background: rgba(10, 10, 10, 0.95);
            backdrop-filter: blur(10px);
            border-bottom: 1px solid rgba(0, 255, 194, 0.2);
            padding: 20px 40px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            position: sticky;
            top: 0;
            z-index: 100;
            transition: background 0.3s, border-color 0.3s;
        }

        body.light-mode .navbar {
            background: rgba(255, 255, 255, 0.95);
            border-bottom: 1px solid rgba(168, 85, 247, 0.2);
        }

        .theme-toggle {
            background: linear-gradient(135deg, #00ffc2 0%, #a855f7 100%);
            border: none;
            color: #000;
            font-weight: 700;
            padding: 10px 20px;
            border-radius: 50px;
            cursor: pointer;
            font-size: 0.85em;
            letter-spacing: 1px;
            text-transform: uppercase;
            transition: transform 0.3s, box-shadow 0.3s;
        }

        .theme-toggle:hover {
            transform: scale(1.05);
            box-shadow: 0 5px 20px rgba(0, 255, 194, 0.4);
        }

        .nav-logo {
            font-size: 1.5em;
            font-weight: 700;
            background: linear-gradient(135deg, #00ffc2 0%, #a855f7 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .nav-links {
            display: flex;
            gap: 30px;
            list-style: none;
        }

        .nav-links a {
            color: #888;
            text-decoration: none;
            font-weight: 600;
            text-transform: uppercase;
            font-size: 0.85em;
            letter-spacing: 1px;
            transition: color 0.3s;
            position: relative;
        }

        body.light-mode .nav-links a {
            color: #666;
        }

        .nav-links a::after {
            content: '';
            position: absolute;
            bottom: -5px;
            left: 0;
            width: 0;
            height: 2px;
            background: linear-gradient(90deg, #00ffc2, #a855f7);
            transition: width 0.3s;
        }

        .nav-links a:hover {
            color: #00ffc2;
        }

        .nav-links a:hover::after {
            width: 100%;
        }

        .header {
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            color: white;
            padding: 80px 40px 60px;
            text-align: center;
            position: relative;
            overflow: hidden;
            transition: background 0.3s;
        }

        body.light-mode .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        }

        .header::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background:
                radial-gradient(circle at 20% 50%, rgba(0, 255, 194, 0.1) 0%, transparent 50%),
                radial-gradient(circle at 80% 80%, rgba(138, 43, 226, 0.1) 0%, transparent 50%);
            animation: pulse 8s ease-in-out infinite;
        }

        @keyframes pulse {
            0%, 100% { opacity: 0.5; }
            50% { opacity: 1; }
        }

        .header-content {
            position: relative;
            z-index: 1;
        }

        .logo {
            font-size: 3.5em;
            margin-bottom: 10px;
            filter: drop-shadow(0 0 20px rgba(0, 255, 194, 0.5));
        }

        .header h1 {
            font-size: 3em;
            margin-bottom: 15px;
            font-weight: 800;
            background: linear-gradient(135deg, #00ffc2 0%, #a855f7 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .header .subtitle {
            font-size: 1.2em;
            opacity: 0.8;
            font-weight: 300;
            letter-spacing: 2px;
            text-transform: uppercase;
        }

        .header .issue {
            margin-top: 20px;
            font-size: 0.9em;
            opacity: 0.6;
            letter-spacing: 1px;
        }

        .content {
            padding: 60px 40px;
            background: #0f0f0f;
            transition: background 0.3s;
        }

        body.light-mode .content {
            background: #ffffff;
        }

        .section {
            margin-bottom: 60px;
        }

        .section-title {
            font-size: 2em;
            background: linear-gradient(135deg, #00ffc2 0%, #a855f7 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 30px;
            font-weight: 700;
            position: relative;
            padding-left: 20px;
        }

        .section-title::before {
            content: '';
            position: absolute;
            left: 0;
            top: 50%;
            transform: translateY(-50%);
            width: 4px;
            height: 80%;
            background: linear-gradient(180deg, #00ffc2 0%, #a855f7 100%);
            border-radius: 2px;
        }

        .featured-article {
            background: linear-gradient(135deg, rgba(26, 26, 46, 0.6) 0%, rgba(15, 52, 96, 0.6) 100%);
            padding: 40px;
            border-radius: 20px;
            margin-bottom: 20px;
            border: 1px solid rgba(0, 255, 194, 0.2);
            backdrop-filter: blur(10px);
            position: relative;
            overflow: hidden;
            transition: transform 0.3s, box-shadow 0.3s, background 0.3s, border-color 0.3s;
        }

        body.light-mode .featured-article {
            background: linear-gradient(135deg, rgba(102, 126, 234, 0.1) 0%, rgba(168, 85, 247, 0.1) 100%);
            border: 1px solid rgba(168, 85, 247, 0.3);
        }

        .featured-article::before {
            content: '';
            position: absolute;
            top: -50%;
            right: -50%;
            width: 200%;
            height: 200%;
            background: radial-gradient(circle, rgba(0, 255, 194, 0.1) 0%, transparent 70%);
            animation: rotate 20s linear infinite;
        }

        @keyframes rotate {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }

        .featured-article:hover {
            transform: translateY(-5px);
            box-shadow: 0 20px 40px rgba(0, 255, 194, 0.2);
        }

        .featured-content {
            position: relative;
            z-index: 1;
        }

        .featured-article h3 {
            color: #00ffc2;
            font-size: 1.8em;
            margin-bottom: 15px;
            font-weight: 700;
            transition: color 0.3s;
        }

        body.light-mode .featured-article h3 {
            color: #667eea;
        }

        .featured-article .meta {
            font-size: 0.9em;
            color: #888;
            margin-bottom: 20px;
            display: flex;
            align-items: center;
            gap: 15px;
        }

        .featured-article p {
            color: #b0b0b0;
            margin-bottom: 25px;
            font-size: 1.05em;
            line-height: 1.8;
            transition: color 0.3s;
        }

        body.light-mode .featured-article p {
            color: #555;
        }

        .read-more {
            display: inline-block;
            padding: 14px 30px;
            background: linear-gradient(135deg, #00ffc2 0%, #a855f7 100%);
            color: #000;
            text-decoration: none;
            border-radius: 50px;
            font-weight: 700;
            transition: transform 0.3s, box-shadow 0.3s;
            letter-spacing: 1px;
            text-transform: uppercase;
            font-size: 0.85em;
        }

        .read-more:hover {
            transform: scale(1.05);
            box-shadow: 0 10px 30px rgba(0, 255, 194, 0.4);
        }

        .article-list {
            list-style: none;
        }

        .article-list li {
            padding: 25px;
            border-bottom: 1px solid rgba(255, 255, 255, 0.1);
            background: rgba(26, 26, 46, 0.3);
            margin-bottom: 15px;
            border-radius: 12px;
            border-left: 3px solid transparent;
            transition: border-color 0.3s, background 0.3s, transform 0.3s;
        }

        body.light-mode .article-list li {
            background: rgba(102, 126, 234, 0.05);
            border-bottom: 1px solid rgba(0, 0, 0, 0.05);
        }

        .article-list li:hover {
            border-left-color: #a855f7;
            background: rgba(26, 26, 46, 0.5);
            transform: translateX(5px);
        }

        body.light-mode .article-list li:hover {
            background: rgba(102, 126, 234, 0.1);
        }

        .article-list li:last-child {
            border-bottom: none;
        }

        .article-list h4 {
            color: #00ffc2;
            margin-bottom: 8px;
            font-size: 1.3em;
            font-weight: 600;
            transition: color 0.3s;
        }

        body.light-mode .article-list h4 {
            color: #667eea;
        }

        .article-list .description {
            color: #888;
            font-size: 0.95em;
            transition: color 0.3s;
        }

        body.light-mode .article-list .description {
            color: #666;
        }

        .tools-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 25px;
        }

        .tool-card {
            background: linear-gradient(135deg, rgba(26, 26, 46, 0.5) 0%, rgba(15, 52, 96, 0.3) 100%);
            padding: 30px;
            border-radius: 16px;
            border: 1px solid rgba(168, 85, 247, 0.2);
            transition: transform 0.3s, box-shadow 0.3s, border-color 0.3s, background 0.3s;
            position: relative;
            overflow: hidden;
        }

        body.light-mode .tool-card {
            background: linear-gradient(135deg, rgba(102, 126, 234, 0.08) 0%, rgba(168, 85, 247, 0.08) 100%);
            border: 1px solid rgba(168, 85, 247, 0.2);
        }

        .tool-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 2px;
            background: linear-gradient(90deg, transparent, #a855f7, transparent);
            transform: translateX(-100%);
            transition: transform 0.5s;
        }

        .tool-card:hover::before {
            transform: translateX(100%);
        }

        .tool-card:hover {
            transform: translateY(-8px);
            box-shadow: 0 15px 35px rgba(168, 85, 247, 0.3);
            border-color: #a855f7;
        }

        .tool-card h4 {
            color: #a855f7;
            margin-bottom: 12px;
            font-size: 1.2em;
            font-weight: 700;
            transition: color 0.3s;
        }

        body.light-mode .tool-card h4 {
            color: #764ba2;
        }

        .tool-card p {
            color: #888;
            font-size: 0.95em;
            line-height: 1.6;
            transition: color 0.3s;
        }

        body.light-mode .tool-card p {
            color: #666;
        }

        .footer {
            background: linear-gradient(135deg, #1a1a2e 0%, #0a0a0a 100%);
            color: white;
            padding: 50px 40px;
            text-align: center;
            border-top: 1px solid rgba(0, 255, 194, 0.2);
            transition: background 0.3s, border-color 0.3s;
        }

        body.light-mode .footer {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            border-top: 1px solid rgba(168, 85, 247, 0.2);
        }

        .footer p {
            margin-bottom: 12px;
            color: #888;
        }

        a {
            color: #00ffc2;
            text-decoration: none;
            transition: color 0.3s;
        }

        .footer a {
            color: #00ffc2;
            text-decoration: none;
            transition: color 0.3s;
        }

        .footer a:hover {
            color: #a855f7;
        }

        .social-links {
            margin-top: 30px;
            display: flex;
            justify-content: center;
            gap: 20px;
        }

        .social-links a {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 50px;
            height: 50px;
            border-radius: 50%;
            background: rgba(0, 255, 194, 0.1);
            color: #00ffc2;
            text-decoration: none;
            font-weight: 700;
            transition: all 0.3s;
            border: 1px solid rgba(0, 255, 194, 0.2);
        }

        .social-links a:hover {
            background: #00ffc2;
            color: #000;
            transform: scale(1.1);
            box-shadow: 0 5px 20px rgba(0, 255, 194, 0.4);
        }

        @media (max-width: 600px) {
            .header h1 {
                font-size: 2em;
            }

            .content {
                padding: 40px 20px;
            }

            .tools-grid {
                grid-template-columns: 1fr;
            }

            .header {
                padding: 60px 20px 40px;
            }

            .navbar {
                padding: 15px 20px;
                flex-direction: column;
                gap: 15px;
            }

            .nav-links {
                gap: 20px;
            }
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
MDCUT=$DIR/${ORG}_metadata.tsv.cut

#save some awk functions for tabulation
# note that the shell functions are used when less than 500 results as the tooltip can be used
tblx(){
echo "<table border=\"1\"><tr><th style=\"width:10%\" >SRA Project Accession</th><th> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\">QC summary </a> <a href=\"https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md\" target=\"_blank\"> <img src=\"/images/question.png\" alt=\"alttext\" title=\"Learn more about the quality metrics\" style=\"width:30px;height:30px;\"> </a> </th><th style=\"width:10%\" >Project Title</th><th style=\"width:500%\" >Project Description </th><th style=\"width:10%\" >Study Type</th><th style=\"width:10%\" >GEO Series</th></tr>"
while read line ; do
  C1=$(echo "$line" | cut -f1 | tr -d '"')
  ZIPURL=https://dee2.io/$(find /dee2_data/huge/$ORG | cut -d '/' -f3- | grep ${C1}_ )
  QCINFO=$(grep -w $C1 $MDCUT | cut -f2 | cut -d '(' -f1 | sort | uniq -c | awk '{print $2":", $1}')
  C3=$(echo "$line" | cut -f3 | tr -d '"' )
  C4=$(echo "$line" | cut -f4 | tr -d '"' )
  C6=$(echo "$line" | cut -f6 | tr -d '"' )
  C9=$(echo "$line" | cut -f9 | tr -d '"' )
  C10=$(echo "$line" | cut -f10 | tr -d '"' )
  echo "<tr><td>" $C1 "<br><br> <a href=$ZIPURL target=_blank >" DEE2 data bundle link "</a> <br><br> <a href=https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=$C1 target=_blank >"SRA link" </a> </td><td>" $QCINFO "</td><td>" $C3 "</td><td>" $C4 "</td><td>" $C6"</td><td>" $C10 "</td></tr>"
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
  echo "<p>Accesion Query: $ACC</p>"
  Q=$(echo $ACC | sed 's/\%2C/\|/g' | sed 's/^/\(/' | sed 's/$/\)/' | tr '+' ' ')
  CNT=$(cut -f-10 $MD | awk '!arr[$1]++' | egrep -wic "$ACC")

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

  echo "<h4>${CNT} datasets found.</h4>"
  cut -f-10 $MD | egrep -iw "$Q" | awk '!arr[$1]++' | sort -k1 | tblx
  echo '</table>'

  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
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

  echo "<p>Keyword Query: $KEY</p>"
  echo "<h4>${CNT} datasets found.</h4>"
  cut -f-10 $MD | egrep -i "$Q" | awk '!arr[$1]++' | sort -k1 | tblx
  echo '</table>'

  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi
rm -f $TMP




echo '</body>
</html>'


