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

for (ORG in c("Arabidopsis thaliana","Caenorhabditis elegans","Dario rerio","Drosophila melanogaster","Escherichia coli","Homo sapiens","Mus musculus","Rattus norvegicus","Saccharomyces cerevisiae")) {
#for (ORG in c("Arabidopsis thaliana","Caenorhabditis elegans")) {
  s<-sample[which(sample$scientific_name==ORG),]
  a<-merge(experiment,study,by="study_accession")
  a<-merge(a,s,by="sample_accession")
  a<-merge(a,run,by="experiment_accession")
  a<-a[,grep("\\.y$",names(a),invert=T)]
  a<-a[,grep("\\.1$",names(a),invert=T)] 
  names(a)=gsub("\\.x","",names(a))
  org=gsub(' ','.',ORG)
  org=tolower(org)
  elastic("http://localhost:9200", org, "data") %index% a
}

q()
#This is how to query
#curl -XPOST "http://localhost:9200/_search" -d'{"query": {"query_string": {"query": "DRX000006" } }}'
#http://joelabrahamsson.com/elasticsearch-101/
#https://cran.r-project.org/web/packages/elasticsearchr/vignettes/quick_start.html


#for (ORG in c("Arabidopsis thaliana","Caenorhabditis elegans","Dario rerio","Drosophila melanogaster","Escherichia coli","Homo sapiens","Mus musculus","Rattus norvegicus","Saccharomyces cerevisiae") {


