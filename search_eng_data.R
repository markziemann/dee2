library(SRAdb)
library(dplyr)
library(readr)
library(jsonlite)

setwd("sradb")

sqlfile <- 'SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),sqlfile)


study<-dbGetQuery(sra_con,"select * from study limit 100000000")
experiment<-dbGetQuery(sra_con,"select * from experiment limit 100000000")
run<-dbGetQuery(sra_con,"select * from run limit 100000000")
sample<-dbGetQuery(sra_con,"select * from sample limit 100000000")


a<-merge(head(experiment),study,by="study_accession")
a<-merge(a,sample,by="sample_accession")
a<-merge(a,run,by="experiment_accession")


a %>% 
    toJSON() %>%
    write_lines("test.json")

