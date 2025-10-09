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
echo "<h2 class='section-title'>Your request: Accession: $ACC Combination: $COMBO </h2>"
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

FREE_LIM=5
PREMIUM_LIM=200
SUPER_LIM=5000

PREMIUM_MEMBER=$(grep -c $EMAIL_ADDRESS /usr/lib/cgi-bin/newrequests/premium)
SUPER_MEMBER=$(grep -c $EMAIL_ADDRESS /usr/lib/cgi-bin/newrequests/super)

if [ $SUPER_MEMBER -eq 1 ] ; then
  echo "Welcome back SUPER member. Thanks for your support!"
  if [ "$USER_REQUEST_COUNT_IN_LAST_YEAR" -le "$SUPER_LIM" ] ; then
  echo "Starting the analysis for you now. You will receive a notification email when the data is ready."
    echo "<br>"
    echo "<br>"
    CONFIRMED="TRUE"
  else
    echo "Sorry, can't start the analysis. Premium membership has a maximum of $SUPER_LIM datasets per year."
    echo "Get in touch with the DEE2 team to make other arrangements"
    echo "<br>"
    echo "<br>"
    CONFIRMED="FALSE"
  fi
fi

if [ $PREMIUM_MEMBER -eq 1 ] ; then
  echo "Welcome back premium member. Thanks for your support!"
  if [ "$USER_REQUEST_COUNT_IN_LAST_YEAR" -le "$PREMIUM_LIM" ] ; then
  echo "Starting the analysis for you now. You will receive a notification email when the data is ready."
    echo "<br>"
    echo "<br>"
    CONFIRMED="TRUE"
  else
    echo "Sorry, can't start the analysis. Premium membership has a maximum of $PREMIUM_LIM datasets per year."
    echo "Get in touch with the DEE2 team to arrange Super membership, giving $SUPER_LIM datasets per year."
    echo "<br>"
    echo "<br>"
    CONFIRMED="FALSE"
  fi
fi

if [ $PREMIUM_MEMBER -eq 0 ] ;then
  if [ $SUPER_MEMBER -eq 0 ] ; then
    if [ "$USER_REQUEST_COUNT_IN_LAST_YEAR" -le $FREE_LIM ] ; then
      echo "Starting the analysis! You will receive a notification email when the data is ready."
      echo "<br>"
      echo "<br>"
      CONFIRMED="TRUE"
    else
      echo "Sorry, can't start the analysis. Free level has a maximum of $FREE_LIM datasets per year."      echo "Get in touch with the DEE2 team to arrange premium membership."
      echo "<br>"
      echo "<br>"
      CONFIRMED="FALSE"
    fi
  fi
fi

if [ $CONFIRMED == "TRUE" ] ; then
  mv $METADATAFILE $METADATAFILE.confirmed
else
  rm $METADATAFILE
fi

echo '</body>
</html>'

