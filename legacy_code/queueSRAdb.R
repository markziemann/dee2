library(SRAdb)
sqlfile <- '/data/projects/mziemann/geo2mx_project/v1/metadata/dlSRADB/SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),sqlfile)
write.table(dbGetQuery(sra_con,"select * from study limit 10000000"),file="study.txt",sep="\t")
write.table(dbGetQuery(sra_con,"select * from experiment limit 10000000"),file="experiment.txt",sep="\t")
write.table(dbGetQuery(sra_con,"select * from run limit 10000000"),file="run.txt",sep="\t")
write.table(dbGetQuery(sra_con,"select * from sample limit 10000000"),file="sample.txt",sep="\t")
srx<-read.table("SRX.txt")
write.table(sraConvert(srx$V1, sra_con = sra_con),file="SRXaccessions.txt",sep="\t",quote=FALSE)
