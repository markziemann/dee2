#!/usr/bin/env Rscript
setwd("/scratch/mziemann/dee2/code/") 

library(SRAdb)
library(parallel)
library(data.table)

CORES=ceiling(detectCores()/2)
for (org in c("ecoli" , "scerevisiae" , "celegans", "athaliana",  "rnorvegicus" , "dmelanogaster", "drerio", "hsapiens", "mmusculus" ) ) {
#org="ecoli"
  #create a list of species full names
  species_list<-c("'Arabidopsis thaliana'","'Caenorhabditis elegans'","'Drosophila melanogaster'","'Danio rerio'",
  "'Escherichia coli'","'Homo sapiens'", "'Mus musculus'", "'Rattus norvegicus'", "'Saccharomyces cerevisiae'")
  #now annotate the short names 
  names(species_list)<- c("athaliana", "celegans", "dmelanogaster", "drerio",
  "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae")
  species_name<-species_list[[org]]

  print(org)

  ###Set some directories
  CODEWD=getwd()
  DATAWD=paste(normalizePath("../data/"),org,sep="/")
  SRADBWD=normalizePath("../sradb/")
  MXDIR=normalizePath("../mx/")
  QUEUEWD=normalizePath("../queue/")

########################
# Get info from sradb
########################

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
  #define 1 week in seconds
  wk<-60*60*24*7*1

  #Test if ctime older than a week, redownload if neccessary
  if (dif>wk) {
    print("sqlfile too old - getting a new one now..")
    unlink(sqlfile)
    sqlfile <<- getSRAdbFile()
  } else {
    print("sqlfile still OK!")
  }

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

########################
# Now determine which datasets have already been processed and completed
########################
  setwd(CODEWD)
  finished_files<-list.files(path = DATAWD, pattern = "finished" , full.names = FALSE, recursive = TRUE, no.. = FALSE)

  if ( length(finished_files) > 0 ) { 
   system(paste("./dee_pipeline.sh",org))
   runs_done<-unique( read.table(paste(DATAWD,"/",org,"_val_list.txt",sep=""),stringsAsFactors=F)[,1] )
   print(paste(length(runs_done),"new runs completed"))
   runs_todo<-base::setdiff(runs, runs_done)
   print(paste(length(runs_todo),"requeued runs"))

   #Update queue on webserver
   queue_name=paste(QUEUEWD,"/",org,".queue.txt",sep="")
   write.table(runs_todo,queue_name,quote=F,row.names=F,col.names=F)
   SCP_COMMAND=paste("scp -i ~/.ssh/cloud/id_rsa ",queue_name ," ubuntu@118.138.240.228:~/Public")
   system(SCP_COMMAND)

   #Update metadata on webserver
   accessions_done<-accessions[which(accessions$run %in% runs_done),]
   write.table(accessions_done,file=paste(SRADBWD,"/",org,"_accessions.tsv",sep=""),quote=F,row.names=F)

   #now the metadata   experiment submission     study    sample        run
   x<-paste("('",paste(accessions_done$sample,collapse="','"),"')",sep="")
   sql_statement<-paste("select * FROM sample WHERE sample_accession IN ",x)
   sample_metadata<-dbGetQuery(sra_con, sql_statement)
 
   x<-paste("('",paste(accessions_done$run,collapse="','"),"')",sep="")
   sql_statement<-paste("select * FROM run WHERE run_accession IN ",x)
   run_metadata<-dbGetQuery(sra_con, sql_statement)

   x<-paste("('",paste(accessions_done$experiment,collapse="','"),"')",sep="")
   sql_statement<-paste("select * FROM experiment WHERE experiment_accession IN ",x)
   experiment_metadata<-dbGetQuery(sra_con, sql_statement)

   x<-paste("('",paste(accessions_done$submission,collapse="','"),"')",sep="")
   sql_statement<-paste("select * FROM submission WHERE submission_accession IN ",x)
   submission_metadata<-dbGetQuery(sra_con, sql_statement)

   x<-paste("('",paste(accessions_done$study,collapse="','"),"')",sep="")
   sql_statement<-paste("select * FROM study WHERE study_accession IN ",x)
   study_metadata<-dbGetQuery(sra_con, sql_statement)

   #now merging this data together to make a big metadata table
   x<-merge(run_metadata,experiment_metadata,by='experiment_accession')
   x<-x[,-grep(".y$",colnames(x))]
   colnames(x)<-gsub(".x$","",colnames(x))

   x<-merge(x,sample_metadata,by='sample_accession')
   x<-x[,-grep(".y$",colnames(x))]
   colnames(x)<-gsub(".x$","",colnames(x))

   x<-merge(x,submission_metadata,by='submission_accession')
   x<-x[,-grep(".y$",colnames(x))] 
   colnames(x)<-gsub(".x$","",colnames(x))

   x<-merge(x,study_metadata,by='study_accession')
   x<-x[,-grep(".y$",colnames(x))]
   colnames(x)<-gsub(".x$","",colnames(x))

   #collect QC info - this is temporary and logic will be incorporated in future
   x$QC_summary="PASS"

   #Need to rearrange columns
   GSE<-function(i) {res=grepl("GEO:",i) ; if (res == FALSE) { j="NA" } else { j<-gsub("GEO: ","",i) ; j} }
   x$GSE_accession<-as.vector(sapply(x[,66],GSE))
   GSM<-function(i) {res=grepl("GEO Accession:",i) ; if (res == FALSE) { j="NA" } else { j<-gsub("GEO Accession: ","",i) ; j} }
   x$GSM_accession<-as.vector(sapply(x[,53],GSM))

   #extract out the important accessions in order
   x2<-as.data.frame(cbind(x$run_accession,x$QC_summary,x$experiment_accession,x$sample_accession,
   x$study_accession,x$submission_accession, x$GSE_accession, x$GSM_accession))

   colnames(x2)<-c("SRR_accession","QC_summary","SRX_accession","SRS_accession",
   "SRP_accession","SRA_accession","GSE_accession","GSM_accession")

   #now remove the moved cols from x
   x<-x[, !(colnames(x) %in% c("run_accession", "QC_summary","experiment_accession","sample_accession",
   "study_accession","submission_accession","GSE_accession","GSM_accession"))]

   x<-cbind(x2,x)
  
   #save metadata as mysql
   connection <- dbConnect(SQLite())
   dbWriteTable(connection, value = x, name = org, append = TRUE )

   write.table(x,file=paste(SRADBWD,"/",org,"_metadata.tsv",sep=""),quote=F,sep="\t",row.names=F)

   #upload 
   SCP_COMMAND=paste("scp -i ~/.ssh/cloud/id_rsa ", paste(SRADBWD,"/",org,"_metadata.tsv",sep="") ," ubuntu@118.138.240.228:/mnt/dee2_data/metadata")
   system(SCP_COMMAND)

  }

  setwd("/scratch/mziemann/dee2/code/")

  rowcnt2<-function( file) { z<-system(paste("wc -l < ",file) , intern=TRUE) ; z}

  png("dee_datasets.png",width=580,height=580)

  FILES1<-list.files(pattern="*queue.txt$",path="/scratch/mziemann/dee2/queue/",full.names=T)
  x<-as.data.frame(sapply(FILES1,rowcnt2),stringsAsFactors=FALSE)
  rownames(x)=c("A. thaliana","C. elegans","D. melanogaster","D. rerio","E. coli","H. sapiens","M. musculus","R. norvegicus","S. cerevisiae")
  colnames(x)="queued"

  FILES2<-list.files(pattern="*accessions.tsv$",path="/scratch/mziemann/dee2/sradb/",full.names=T)
  y<-as.data.frame(sapply(FILES2,rowcnt2),stringsAsFactors=FALSE)
  rownames(y)=c("A. thaliana","C. elegans","D. melanogaster","D. rerio","E. coli","H. sapiens","M. musculus","R. norvegicus","S. cerevisiae")
  colnames(y)="completed"

  z<-merge(x,y,by=0)
  rownames(z)=z$Row.names
  z$Row.names=NULL

  DATE=strsplit(as.character(file.info(FILES2[6])[,4])," ",fixed=T)[[1]][1]
  HEADER=paste("Updated",DATE)
  #colnames(z)=HEADER
  z<-z[order(rownames(z),decreasing=T ), ,drop=F]
  par(las=2) ; par(mai=c(1,2.5,1,0.5))
  MAX=max(as.numeric(z[,1]))+80000

  bb<-barplot( rbind( as.numeric(z$queued) , as.numeric(z$completed) ) ,
   names.arg=rownames(z) ,xlim=c(0,MAX),beside=T, main=HEADER, col=c("darkblue","red") ,
   horiz=T , las=1, cex.axis=1.3, cex.names=1.4,cex.main=1.4 )

  legend("topright", colnames(z), fill=c("darkblue","red") , cex=1.2)

  text( cbind(as.numeric(z[,1])+20000 ,as.numeric(z[,2])+20000 )  ,t(bb),labels=c(z[,1],z[,2]) ,cex=1.2)
  dev.off()
  system("scp -i ~/.ssh/cloud/cloud2.key dee_datasets.png ubuntu@118.138.240.228:/mnt/dee2_data/mx")

}
