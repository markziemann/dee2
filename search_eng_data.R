library(SRAdb)
library(dplyr)
library(readr)
library(jsonlite)

setwd("../sradb")

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



q()
http://joelabrahamsson.com/elasticsearch-101/
#need to write json to file then use curl to index
sed 's/{/\n\{\n/g;s/\}/\n\}\n/g;s/","/",\n"/g' test.json  \
| tr -s ' ' | tr -d '\n' | cut -c2- | rev | cut -c2- | rev | sed 's/},{/}\n{/g' | tr '.' ' ' \
| while read line ; do
  echo $line > tmp.json
  curl -XPUT "http://localhost:9200/dat/dat/1" -d"$(cat tmp.json)"
done

