library(SRAdb)
sqlfile <- '../SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),sqlfile)
write.table(dbGetQuery(sra_con,"select * from study limit 100000000"),file="study.txt",sep="\t")
write.table(dbGetQuery(sra_con,"select * from experiment limit 100000000"),file="experiment.txt",sep="\t")
write.table(dbGetQuery(sra_con,"select * from run limit 100000000"),file="run.txt",sep="\t")
write.table(dbGetQuery(sra_con,"select * from sample limit 100000000"),file="sample.txt",sep="\t")
sracompl<-read.table("queue/SRR_complete.txt")
write.table(sraConvert(sracompl$V1, sra_con = sra_con),file="SRXaccessions.txt",sep="\t",quote=FALSE)

