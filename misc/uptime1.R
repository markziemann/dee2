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

#source("https://raw.githubusercontent.com/markziemann/dee2/master/getDEE2.R")
source("getDEE2.R")


#################################################
# This is an example of try
#################################################
#message("Should work")
#res<-try(getDEE2("drerio", "ERR1759699") )
#length(res)

#message("Wont work")
#res<-try(getDEE2("drerio","SRR7252513") )
#length(res)

###########################################
# Part 1 check that URL is up
###########################################
message ("Starting check 1: checking that URL exists")
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
 message ("Check 1 result: URL OK")
}

###########################################
# Part 2 check that the search is working
# Select a few searches that should give a minimum number of hits for each species
###########################################
message ("Starting check 2: checking that search is working properly")
#athaliana,petiole,21
ATH<-t(list("athaliana","petiole",21))
CEL<-t(list("celegans","sexual",52))
DME<-t(list("dmelanogaster","pathogenic",145))
DRE<-t(list("drerio","cardiac",256))
ECO<-t(list("ecoli","pathogen",334))
HSA<-t(list("hsapiens","mettl3",128))
MMU<-t(list("mmusculus","saha",71))
RNO<-t(list("rnorvegicus","epilepsy",22))
SCE<-t(list("scerevisiae","pathogen",238))

QUERIES<-rbind(ATH,CEL,DME,DRE,ECO,HSA,MMU,RNO,SCE)

for ( i in 1:nrow(QUERIES)) {
#for ( i in 4) {
 ORG=QUERIES[i,][1]
 KEYWORD=QUERIES[i,][2]
 NUM=QUERIES[i,][3]

 CMD=paste("curl \'http://dee2.io/cgi-bin/search.sh?org=",ORG,"&accessionsearch=&keywordsearch=",KEYWORD,"\' | tee search_res.html | html2text | awk \'NR==2 {print $1}\' ",sep="")
 CNT<-system(CMD,intern=TRUE)

 CNT>=NUM
 #21 or greater
 if ( NUM > CNT ) {
  message("ERROR too few datasets were identified. there may be something wrong with the webserver")
  DATE=date()
  from <- "uptime_bot@dee2.io"
  to <- "mark.ziemann@gmail.com"
  subject <- "DEE2 search error"
  body <- paste("The DEE2 search function is not working properly as of ",date(),". Less than 21 results for ",ORG,"keyword search for ",KEYWORD,". Next test in 10 minutes")
  if ( ALERTSENT!=1 ) { 
   sendmail(from=from,to=to,subject=subject,msg=body)
   ALERTSENT=1
  }
 } else {
  message("Check 2 result: Search working OK")
 }

###########################################
# Part 3 check that the download is working
###########################################
 message ("Starting check 3: checking that download is working properly")
 SRRvec<-system("html2text search_res.html | cut -d '|' -f3 | grep RR | tr -d ' ' " ,intern=TRUE)
 x<-try(getDEE2(ORG,SRRvec,quiet=TRUE))
 if ( length(x)>3 ) {
 
  CNT_DL=dim(x$GeneCounts)[2]
  if ( CNT == CNT_DL ) {
   message("Check 3 result: Download working OK")
  } else {
   message("ERROR too few datasets were identified. there may be something wrong with the webserver")
   DATE=date()
   from <- "uptime_bot@dee2.io"
   to <- "mark.ziemann@gmail.com"
   subject <- "DEE2 search error"
   body <- paste("The DEE2 download function is not working properly as of ",date(),". The obtained datasets for ",ORG,"keyword search for ",KEYWORD,"did not match the search results. Next test in 10 minutes")
   if ( ALERTSENT!=1 ) {
    sendmail(from=from,to=to,subject=subject,msg=body)
    ALERTSENT=1
   }
  }
 }
}
