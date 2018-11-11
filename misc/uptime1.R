# The idea is to catch most common type of failures to DEE2
# This script shuld be run every 10 minuts via a cron job

# 1 The webserver is responding

# 2 The search is working

# 3 The download doesnt work properly

# All 9 species can be tested, one species every 10 minutes via a cron job
# Use "touch" and ifelse to determine which species should be tested next

#Dependancy in ubuntu 16
#sudo apt install libcurl4-openssl-dev
#sudo apt-get install mailutils
#sudo apt-get install html2text
#install.packges("RCurl")

library(RCurl)
library(sendmailR)

###########################################
# Part 1 check that URL is up
###########################################
ALERTSENT=0
if ( ! url.exists("dee2.io") ) { 
 message("URL does not exist") 
 DATE=date()
 from <- "uptime_bot@dee2.io"
 to <- "mark.ziemann@gmail.com"
 subject <- "DEE2 server is down"
 body <- paste("The URL dee2.io is not available as of",date(),". Next test in 10 minutes") 
 sendmail(from=from,to=to,subject=subject,msg=body)
 ALERTSENT=1
} else {
 message ("URL OK")
}

###########################################
# Part 2 check that the search is working
# Select a few searches that should give a minimum number of hits for each species
###########################################

#athaliana,petiole,21
CMD="curl \'http://dee2.io/cgi-bin/search.sh?org=athaliana&accessionsearch=&keywordsearch=petiole\' | tee petiole.html | html2text | awk \'NR==2 {print $1}\' "
CNT_ATH<-system(CMD,intern=TRUE)
#21 or greater
if ( CNT_ATH < 21 ) {
 message("ERROR too few datasets were identified. there may e something wrong with the webserver")
 DATE=date()
 from <- "uptime_bot@dee2.io"
 to <- "mark.ziemann@gmail.com"
 subject <- "DEE2 search error"
 body <- paste("The DEE2 search function is not working properly as of ",date(),". Less than 21 results for Arabidopsis keyword search for petioles. Next test in 10 minutes")
 if ( ALERTSENT!=1 ) { 
  sendmail(from=from,to=to,subject=subject,msg=body)
  ALERTSENT=1
 }
}

###########################################
# Part 3 check that the download is working
###########################################
source("https://raw.githubusercontent.com/markziemann/dee2/master/getDEE2.R")
SRRvec<-system("html2text petiole.html | cut -d '|' -f3 | grep RR | tr -d ' ' " ,intern=TRUE)
x<-getDEE2("athaliana",SRRvec)
CNT_ATH_DL=dim(x$GeneCounts)[2]
if ( CNT_ATH != CNT_ATH_DL ) {
  message("ERROR too few datasets were identified. there may be something wrong with the webserver")
 DATE=date()
 from <- "uptime_bot@dee2.io"
 to <- "mark.ziemann@gmail.com"
 subject <- "DEE2 search error"
 body <- paste("The DEE2 download function is not working properly as of ",date(),". The obtained datasets did not match the search results. Next test in 10 minutes")
 if ( ALERTSENT!=1 ) {
  sendmail(from=from,to=to,subject=subject,msg=body)
  ALERTSENT=1
 }
}
