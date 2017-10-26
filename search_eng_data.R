library(SRAdb)
library("elasticsearchr")
#library(dplyr)
#library(readr)
#library(jsonlite)
#library(curl)

setwd("../sradb")
sqlfile <- 'SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),sqlfile)
study<-dbGetQuery(sra_con,"select * from study limit 100000000")
experiment<-dbGetQuery(sra_con,"select * from experiment limit 100000000")
run<-dbGetQuery(sra_con,"select * from run limit 100000000")
sample<-dbGetQuery(sra_con,"select * from sample limit 100000000")

#a<-merge(head(experiment),study,by="study_accession")
a<-merge(experiment,study,by="study_accession")
a<-merge(a,sample,by="sample_accession")
a<-merge(a,run,by="experiment_accession")

a<-a[,grep("\\.y$",names(a),invert=T)]
a<-a[,grep("\\.1$",names(a),invert=T)] 
names(a)=gsub("\\.x","",names(a))

b<-head(a[which(a$scientific_name=="Saccharomyces cerevisiae"),])

#match completed SRR numbers from the 

elastic("http://localhost:9200", "bdata", "data") %index% b


q()
#This is how to query
#curl -XPOST "http://localhost:9200/_search" -d'{"query": {"query_string": {"query": "DRX000006" } }}'
#http://joelabrahamsson.com/elasticsearch-101/
#https://cran.r-project.org/web/packages/elasticsearchr/vignettes/quick_start.html
