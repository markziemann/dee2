#!/usr/bin/env Rscript
setwd("/scratch/mziemann/dee2/code/") 

list.of.packages <- c("SRAdb")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) {
  #install.packages(new.packages)
  source("https://bioconductor.org/biocLite.R") ; biocLite("SRAdb")
}

library(SRAdb)

for (org in c("athaliana", "celegans", "dmelanogaster", "drerio", "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae") ) {

  #create a list of species full names
  species_list <- c("'Arabidopsis thaliana'", "'Caenorhabditis elegans'", "'Drosophila melanogaster'", "'Danio rerio'", "'Escherichia coli'", "'Homo sapiens'", "'Mus musculus'", "'Rattus norvegicus'", "'Saccharomyces cerevisiae'")
  #now annotate the short names 
  names(species_list)<- c("athaliana", "celegans", "dmelanogaster", "drerio", "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae")

  species_name<-species_list[[org]]

  ###Set some directories
  CODEWD=getwd()
  DATAWD=paste(normalizePath("../data/"),org,sep="/")
  SRADBWD=normalizePath("../sradb/")
  MXDIR=normalizePath("../mx/")
  QUEUEWD=normalizePath("../queue/")

  setwd(SRADBWD)

  sqlfile <- 'SRAmetadb.sqlite'
  #sqlfile<-paste(SRADBWD,"/SRAmetadb.sqlite",sep="")
  if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()

  #creation time
  ct<-as.numeric(file.info('SRAmetadb.sqlite')$ctime)
  #current time
  now<-as.numeric(Sys.time())
  #calculate difference between
  dif<-now-ct
  #define 4week in seconds
  wk<-60*60*24*7*4*12
  #Test if ctime older than a week, redownload if neccessary
  if (dif>wk) {print("sqlfile too old - getting a new one now..") ; unlink(sqlfile) ; sqlfile <<- getSRAdbFile() } else {print("sqlfile still OK!")}

  #connect to DB
  sra_con <- dbConnect(SQLite(),sqlfile)
  sra_tables <- dbListTables(sra_con)

  #make a list of RNA-seq samples
  rs<-dbGetQuery(sra_con,"select study_accession, sample_accession, experiment_accession, library_source from experiment where library_strategy like 'RNA-seq' ")

  #make a list of samples from a particular species
  sqlStatement <- paste("select sample_accession , scientific_name FROM sample WHERE scientific_name LIKE ",species_name)
  species<-dbGetQuery(sra_con, sqlStatement)

  #filter accession list by
  accessions<-merge(rs,species,by="sample_accession")
  accessions<-sraConvert(accessions$experiment_accession,sra_con = sra_con)
  runs<-accessions$run

  ###### Here need to determine which datasets already processed
  setwd(DATAWD)
  finished_files<-list.files(path = ".", pattern = "finished" , full.names = FALSE, recursive = TRUE, no.. = FALSE)
  runs_done<-basename(as.character(strsplit(finished_files,".finished")))
  runs_todo<-setdiff(runs, runs_done)
  setwd(CODEWD)
  queue_name=paste(QUEUEWD,"/",org,".queue.txt",sep="")
  write.table(runs_todo,queue_name,quote=F,row.names=F,col.names=F)
  SCP_COMMAND=paste("scp -i ~/.ssh/cloud/id_rsa ",queue_name ," ubuntu@118.138.241.34:~/Public")
  system(SCP_COMMAND)

  setwd(DATAWD)

  #define read.tsv
  read.tsv<-function(file, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", ...){
  read.table(file = file, header = header, sep = sep, quote = quote, dec = dec, fill = fill, row.names=1, comment.char = comment.char, ...)}

  print("generate a list of se tsv files to import and merge")

  filterempty<-function(file){CNT=length(readLines(file));if(CNT>0){print(file)}}

  se_list<-list.files(path = ".", pattern = "se.tsv" , full.names = TRUE , recursive = TRUE, no.. = FALSE)
  se_list<-as.character(mclapply(se_list,filterempty))
  se<-do.call("cbind", lapply(se_list, read.tsv))
  print("se list finished OK")

  ke_list<-list.files(path = ".", pattern = "ke.tsv" , full.names = TRUE , recursive = TRUE, no.. = FALSE)
  ke<-do.call("cbind", lapply(ke_list, read.tsv))
  COLS=paste(grep("_est_counts", names(ke), value = TRUE))
  ke_counts<-ke[,COLS]
  COLS=NULL
  COLS=paste(grep("_tpm", names(ke), value = TRUE))
  ke_tpm<-ke[,COLS]
  print("ke list finished OK")

  print("writing SE mx")
  OUT=paste(MXDIR,"/",org,"_se.tsv",sep="")
  write.table(se,file=OUT,sep="\t",quote=F,row.names=T)

  print("writing KE counts mx")
  OUT=paste(MXDIR,"/",org,"_ke_counts.tsv",sep="")
  write.table(ke_counts,file=OUT,sep="\t",quote=F,row.names=T)

  print("writing KE tpm mx")
  OUT=paste(MXDIR,"/",org,"_ke_tpm.tsv",sep="")
  write.table(ke_tpm,file=OUT,sep="\t",quote=F,row.names=T)

  setwd(CODEWD)

}
