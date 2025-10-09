#!/bin/bash
echo "Content-type: text/html"
echo ''
echo '<!DOCTYPE html>
<html xml:lang="en" lang="en">
<link rel='shortcut icon' type='image/x-icon' href='favicon.ico' />
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DEE2
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
    <div class="container">
        <div class="content">
            <div class="section">

'

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
echo "<h2 class='section-title'>Your request: $ORG $ACC $EMAIL</h2>"
echo "<br>"
echo "<br>"

#Error handling if no input provided
if [ -z "$ACC" ] ; then
  echo 'Error: No accession number provided.'
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Search again" onClick="history.go(-1);return true;" style="font-size : 22px;" ></FORM>'
  exit
fi

#Error handlingif no email address provided
if [ -z "$EMAIL" ] ; then
  echo 'Error: No email address provided.'
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
  exit
fi

# check the previously requested data
PRESENT=$(find /dee2_data/requests | grep -c ${ACC}.zip )
if [ $PRESENT -eq 1 ] ; then
  FILEPATH=$(find /dee2_data/requests | grep ${ACC}.zip )
  DATALINK=$(echo $FILEPATH | sed 's#/dee2_data/#https://dee2.io/#')
  echo "Data set is already available at the URL $DATALINK"
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;"></FORM>'
  exit
fi

# check local queue for metadata
SRPQUEUE=/dee2_data/srpqueue/${ORG}_srpqueue.txt
QHIT=$(grep -wc $ACC $SRPQUEUE)
if [ $QHIT -eq 1 ] ; then
  echo "$ACC was listed in the queue."
  echo "<br>"
fi

echo "Fetching and checking $ACC metadata with pysradb.<br>"

# check SRA for metadata
TMPFILE=$(mktemp)
export PATH="$PATH:/var/www/.local/bin"
pysradb metadata $ACC > $TMPFILE
NLINES=$(wc -l < $TMPFILE)
echo "The metadata was found and has $NLINES lines."
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
echo "Here are the runs that will beprocessed: $RUNS" | sed 's/,/, /g'
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
  echo pysradb: $ORG2
  echo "<br>"
  echo "<br>"
  echo Species you submitted: $ORG
  echo "<br>"
  echo "<br>"
  echo '<FORM><INPUT Type="button" VALUE="Go back" onClick="history.go(-1);return true;"></FORM>'
  exit
fi

mv $TMPFILE /usr/lib/cgi-bin/newrequests
TMPFILE=$(basename $TMPFILE)
TMPFILE=/usr/lib/cgi-bin/newrequests/$TMPFILE
echo "EMAIL=$EMAIL" >> $TMPFILE
chmod 664 $TMPFILE

echo "We've sent an email to your nominated address.
Please follow the information in there to initiate the data processing.
If you didn't receive it, then check your spam folder and add www-data@dee2.io and ubuntu@dee2.io to your
contacts."
echo "<br>"
echo "<br>"

COMBO=$(for i in {1..10} ; do echo {a..z} {0..9} | tr ' ' '\n' | shuf | head -1 ; done | paste -s -d '')

echo "COMBO=$COMBO" >> $TMPFILE

EMAIL_BODY=$(cat <<'EOF'
Hello DEE2 user,

A request was made to process dataset ACC from this email address.

Please confirm it by visiting the link https://dee2.io/cgi-bin/confirm.sh?accession=ACC&confirmcode=COMBO

Data processing won't begin unless it is validated.

Best wishes,

Mark and the DEE2 Team

EOF
)

echo "$EMAIL_BODY" \
| sed "s/ACC/$ACC/" \
| sed "s/COMBO/$COMBO/" \
| mailx -s "Confirm $ACC data processing for DEE2" "$EMAIL" \
&& echo "email success" \
|| echo "email failed"

echo '
</div></div></div>
</body>
</html>'

cp $TMPFILE /usr/lib/cgi-bin/newrequests/$ACC
