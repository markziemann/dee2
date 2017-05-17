#!/usr/bin/env Rscript
#read the args get species name in shorthand dmelanogaster
args = commandArgs(trailingOnly=TRUE)
org<-as.character(strsplit(args, " ")[1])
org
if (org=="NULL") {org="dmelanogaster"}
org

list.of.packages <- c("SRAdb")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
#install.packages(new.packages)
source("https://bioconductor.org/biocLite.R") ; biocLite("SRAdb")
}

#if (!require("pacman")) install.packages("pacman")
#pacman::p_load(SRAdb)
#pacman::p_load(package1, package2, package_n), package2, package_n)
library(SRAdb)

#create a list of species full names
species_list <- c("'Arabidopsis thaliana'", "'Caenorhabditis elegans'", "'Drosophila melanogaster'", "'Danio rerio'", "'Escherichia coli'", "'Homo sapiens'", "'Mus musculus'", "'Rattus norvegicus'", "'Saccharomyces cerevisiae'")
#now annotate the short names 
names(species_list)<- c("athaliana", "celegans", "dmelanogaster", "drerio", "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae")
species_name<-species_list[[org]]

###Set some directories
CODEWD=getwd()
PIPELINE=normalizePath("pipeline.sh")
DATAWD=paste(normalizePath("../data/"),org,sep="/")
SRADBWD=normalizePath("../sradb/")
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
wk<-60*60*24*7*4
#Test if ctime older than a week, redownload if neccessary
if (dif>wk) {print("sqlfile too old - getting a new one now..") ; unlink(sqlfile) ; sqlfile <<- getSRAdbFile() } else {print("sqlfile still OK!")}

#connect to DB
sra_con <- dbConnect(SQLite(),sqlfile)
sra_tables <- dbListTables(sra_con)

#make a list of RNA-seq samples
rs<-dbGetQuery(sra_con,"select study_accession, sample_accession, experiment_accession, library_source from experiment where library_strategy like 'RNA-seq' ")
#rs<-dbGetQuery(sra_con,"select study_accession, sample_accession, experiment_accession, library_source from experiment where library_strategy like 'miRNA-seq' ")

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

##check system df disk usage in while loop
##start pipelines

#specify ascp command
#download aspera license here: https://github.com/PRIDE-Toolsuite/pride-inspector/blob/master/aspera/etc/asperaweb_id_dsa.openssh
ascpCMD <- 'ascp -QT -l 300m -i ~/.ascp/aspera-license'

#to get 1 file
#ascpSRA ( runs_todo[1] , sra_con, ascpCMD, fileType = 'sra', destDir = getwd() )

dee.pipeline<-function(SRR){
#for (SRR in runs_todo){
  setwd(DATAWD)
  print(paste("run",SRR))
  #build SRR path for later pipeline
  SRR_FILENAME=paste(SRR,".sra",sep="")
  SRR_PATH=paste(DATAWD,SRR_FILENAME,sep="/")

  #create directory if not present
  if (!dir.exists(file.path(DATAWD, SRR))) {
    dir.create(file.path(DATAWD, SRR))
  }

  #enter srr directory
  setwd(file.path(DATAWD, SRR))

  #download SRR file if not present and disk space is available
  if(!file.exists(SRR_PATH)){
    print(SRR_PATH)
    DISK_AVAIL=system('df . | awk \'NR=="2"{print $4}\'',intern=TRUE)
    DISK_LIM=100000000
    if (DISK_AVAIL>DISK_LIM) {
      #best to download SRA into the SRR subdirectory
      ascpSRA ( SRR , sra_con, ascpCMD, fileType = 'sra', destDir = getwd() )
      #kick off pipeline
    }
  }
  COMMAND_ARGS=paste(SRR_PATH,org)
  COMMAND=paste(PIPELINE,COMMAND_ARGS)
  system(COMMAND)
}
#lapply(runs_todo,dee.pipeline)
#mclapply(runs_todo,dee.pipeline,mc.cores=3)
mclapply(runs_todo[1001:1105],dee.pipeline,mc.cores=1)









q()
#for (SRR in runs_todo[1000:2000]){
for (SRR in runs_todo){
  setwd(DATAWD)
  print(paste("run",SRR))
  #build SRR path for later pipeline
  SRR_FILENAME=paste(SRR,".sra",sep="")
  SRR_PATH=paste(DATAWD,SRR_FILENAME,sep="/")

  #create directory if not present
  if (!dir.exists(file.path(DATAWD, SRR))) {
    dir.create(file.path(DATAWD, SRR))
  }

  #enter srr directory
  setwd(file.path(DATAWD, SRR))

  #download SRR file if not present and disk space is available
  if(!file.exists(SRR_PATH)){
    print(SRR_PATH)
    DISK_AVAIL=system('df . | awk \'NR=="2"{print $4}\'',intern=TRUE)
    DISK_LIM=100000000
    if (DISK_AVAIL>DISK_LIM) {
      #best to download SRA into the SRR subdirectory
      ascpSRA ( SRR , sra_con, ascpCMD, fileType = 'sra', destDir = getwd() )
      #kick off pipeline
    }
  }
  COMMAND_ARGS=paste(SRR_PATH,org)
  COMMAND=paste(PIPELINE,COMMAND_ARGS)
  system(COMMAND)
}





q()
#in_acc<-dm[,5]
#ascpSRA(in_acc, sra_con, ascpCMD, fileType = 'sra', destDir = getwd())

#this works
#ascp -QT -l 300m -i  ~/.ascp/aspera-license  anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR536/SRR5366387/SRR5366387.sra .
-----BEGIN DSA PRIVATE KEY-----
MIIBuwIBAAKBgQDkKQHD6m4yIxgjsey6Pny46acZXERsJHy54p/BqXIyYkVOAkEp
KgvT3qTTNmykWWw4ovOP1+Di1c/2FpYcllcTphkWcS8lA7j012mUEecXavXjPPG0
i3t5vtB8xLy33kQ3e9v9/Lwh0xcRfua0d5UfFwopBIAXvJAr3B6raps8+QIVALws
yeqsx3EolCaCVXJf+61ceJppAoGAPoPtEP4yzHG2XtcxCfXab4u9zE6wPz4ePJt0
UTn3fUvnQmJT7i0KVCRr3g2H2OZMWF12y0jUq8QBuZ2so3CHee7W1VmAdbN7Fxc+
cyV9nE6zURqAaPyt2bE+rgM1pP6LQUYxgD3xKdv1ZG+kDIDEf6U3onjcKbmA6ckx
T6GavoACgYEAobapDv5p2foH+cG5K07sIFD9r0RD7uKJnlqjYAXzFc8U76wXKgu6
WXup2ac0Co+RnZp7Hsa9G+E+iJ6poI9pOR08XTdPly4yDULNST4PwlfrbSFT9FVh
zkWfpOvAUc8fkQAhZqv/PE6VhFQ8w03Z8GpqXx7b3NvBR+EfIx368KoCFEyfl0vH
Ta7g6mGwIMXrdTQQ8fZs
-----END DSA PRIVATE KEY-----





#dbGetQuery(sra_con, "SELECT library_strategy where 'Library Strategy' like RNA-Seq")
