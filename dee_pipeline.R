#!/usr/bin/env Rscript
setwd("/scratch/mziemann/dee2/code/") 

library(SRAdb)
library(parallel)
library(data.table)

CORES=ceiling(detectCores()/2)
for (org in c("ecoli" , "scerevisiae" ,  "rnorvegicus" , "athaliana", "celegans", "dmelanogaster", "drerio", "hsapiens", "mmusculus" ) ) {
#org="ecoli"
  #create a list of species full names
  species_list<-c("'Arabidopsis thaliana'","'Caenorhabditis elegans'","'Drosophila melanogaster'","'Danio rerio'",
  "'Escherichia coli'","'Homo sapiens'", "'Mus musculus'", "'Rattus norvegicus'", "'Saccharomyces cerevisiae'")
  #now annotate the short names 
  names(species_list)<- c("athaliana", "celegans", "dmelanogaster", "drerio",
  "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae")
  species_name<-species_list[[org]]

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
  setwd(DATAWD)
  finished_files<-list.files(path = ".", pattern = "finished" , full.names = FALSE, recursive = TRUE, no.. = FALSE)

  if ( length(finished_files) > 0 ) { 

   runs_done<-basename(as.character(strsplit(finished_files,".finished")))
   runs_todo<-setdiff(runs, runs_done)
   setwd(CODEWD)
   queue_name=paste(QUEUEWD,"/",org,".queue.txt",sep="")

   setwd(DATAWD)

   #define read.tsv
   read.tsv<-function(file, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", ...){
   read.table(file = file, header = header, sep = sep, quote = quote, dec = dec, fill = fill, row.names=1, comment.char = comment.char, ...)}

   print("generate a list of se tsv files to import and merge")
   se_list<-list.files(path = ".", pattern = "\\.se.tsv$" , full.names = TRUE , recursive = TRUE, no.. = FALSE)

   #rowcnt<-function(file){nrow(fread(file,fill=T)) }
   rowcnt<-function(file){nrow(read.table(file,fill=T)) }

   #add the failed datasets to the todo list
   LEN=length(rownames(subset(file.info(se_list),size==0)))
   if (LEN!=0){
    p<-rownames(subset(file.info(se_list),size==0))
    p<-unlist(strsplit(rownames(p),'/'))
    p<-p[seq(1,length(p),2)]
    runs_todo<-unique(union(runs_todo,p))
   }

   se_list<-rownames(subset(file.info(se_list),size!=0))

  #Need to ensure matrix will be square
   x<-as.matrix(table(factor(as.numeric(mclapply(se_list,rowcnt, mc.cores = CORES )))))
   expected_len=as.numeric(rownames(tail(x,n=1)))
   y<-t(as.data.frame(mclapply(se_list,rowcnt, mc.cores = CORES  )))
   rownames(y)=as.character(se_list)
 
   LEN=length(rownames(subset(y,y!=expected_len)))
   if (LEN!=0){
    p<-rownames(subset(y,y!=expected_len))
    p<-unlist(strsplit(p,'/'))
    p<-p[seq(1,length(p),2)]
    runs_todo<-unique(union(runs_todo,p))
   }

   se_list<-rownames(subset(y,y=!expected_len))
   se_name=paste(org,"_se_list.txt",sep="")
   write.table(se_list,file=se_name,quote=F,row.names=F)
 
############################
#  Matrix generation is simply too memory hungry and slow not doing for now
#############################

   print("se list finished OK")

   #now focus on kallisto data
  
   ke_list<-list.files(path = ".", pattern = "\\.ke.tsv$" , full.names = TRUE , recursive = TRUE, no.. = FALSE)

   LEN=length(rownames(subset(file.info(se_list),size==0)))
   if (LEN!=0){
    p<-rownames(subset(file.info(ke_list),size==0))
    p<-unlist(strsplit(rownames(subset(p,p!=expected_len)),'/'))
    p<-p[seq(1,length(p),2)]
    runs_todo<-unique(union(runs_todo,p))
   }

   ke_list<-rownames(subset(file.info(ke_list),size!=0))
   x<-as.matrix(table(factor(as.numeric(mclapply(ke_list,rowcnt , mc.cores = CORES )))))
   expected_len=as.numeric(rownames(tail(x,n=1)))

   y<-t(as.data.frame(mclapply(ke_list,rowcnt, mc.cores = CORES  )))
   rownames(y)=as.character(ke_list)

   LEN=length(rownames(subset(y,y!=expected_len)))
   if (LEN!=0){
    p<-rownames(subset(y,y!=expected_len))
    p<-unlist(strsplit(p,'/'))
    p<-p[seq(1,length(p),2)]
    runs_todo<-unique(union(runs_todo,p))
   }

   y<-t(as.data.frame(mclapply(ke_list,rowcnt , mc.cores = CORES )))

   rownames(y)=as.character(ke_list)
   ke_list<-rownames(subset(y,y==expected_len))
   ke_name=paste(org,"_ke_list.txt",sep="")
   write.table(ke_list,file=ke_name,quote=F,row.names=F)

############################
#  Matrix generation is simply too memory hungry and slow not doing for now
############################
   print("ke list finished OK")

   #QC data
   #create the qc matrix
   system("for I in `find . | grep qc$` ; do grep ':' $I > tmp ; mv tmp $I ; done")
   print("generate a list of qc files to import and merge")
   qc_list<-list.files(path = ".", pattern = "\\.qc$" , full.names = TRUE , recursive = TRUE, no.. = FALSE)

  # rowcnt<-function(file){nrow(fread(file,fill=T)) }

   #add the failed datasets to the todo list
   LEN=length(rownames(subset(file.info(qc_list),size==0)))
   if (LEN!=0){
    p<-rownames(subset(file.info(qc_list),size==0))
    p<-strsplit(p,'/')
    p<-sapply(p, "[[", 2)
    runs_todo<-unique(union(runs_todo,p))
   }

   qc_list<-rownames(subset(file.info(qc_list),size!=0))

   #Need to ensure matrix will be square
   x<-as.matrix(table(factor(as.numeric(mclapply(qc_list,rowcnt, mc.cores = CORES  )))))
   expected_len=as.numeric(rownames(tail(x,n=1)))
   y<-t(as.data.frame(mclapply(qc_list,rowcnt, mc.cores = CORES  )))
   rownames(y)=as.character(qc_list)

   LEN=length(rownames(subset(y,y!=expected_len)))
   if (LEN!=0){
    p<-rownames(subset(y,y!=expected_len))
    p<-unlist(strsplit(p,'/'))
    p<-p[seq(1,length(p),2)]
    runs_todo<-unique(union(runs_todo,p))
   }

   qc_list<-rownames(subset(y,y==expected_len))
   qc_name=paste(org,"_qc_list.txt",sep="")
   write.table(qc_list,file=qc_name,quote=F,row.names=F)

   #redefine read.colon.tsv
   read.colon.tsv<-function(file, header = FALSE, sep = ":", quote = "\"", dec = ".", fill = TRUE, comment.char = "", ...){
   read.table(file = file, sep = sep, quote = quote, dec = dec, fill = fill, row.names=1, comment.char = comment.char, ...)}

   #bringing in some QC data
   qc<-do.call("cbind", mclapply(qc_list, read.colon.tsv,mc.cores = CORES ))
   colnames(qc)<-sapply(strsplit(sub("./","",qc_list), "/"), head, 1)
   print("qc list finished OK")

   #need to do some logic to pass fail data based on QC info and fill in the field
   #for col in numcols do blah blah
   #if then if then
   #qc$COL[29]=PASS
   #qc$COL[29]=FAIL:FASTQ_LENGTH

   #Moved this part to the end due to a small number of files posessing the wrong number of lines 
   write.table(runs_todo,queue_name,quote=F,row.names=F,col.names=F)
   SCP_COMMAND=paste("scp -i ~/.ssh/cloud/id_rsa ",queue_name ," ubuntu@118.138.240.228:~/Public")
   system(SCP_COMMAND)

   setwd(CODEWD)
 
   #PASTE_MX_CMD=paste("./pastemx.sh",org)
   #system(PASTE_MX_CMD)

   #new section to obatin and update webserver metadata for completed runs
   #start with accession information
   runs_done<-setdiff(runs_done,runs_todo)

  if (length(runs_done)>0) {

   #TODO: runs done .finished files renamed to .validated
   #TODO: runs done become write only
   se=NULL
   ke=NULL
   qc=NULL

   se_func<-function(file){
    y<-read.table(file)
    y$run=strsplit(file,"/")[[1]][7]
    y$gene=rownames(y)
    rownames(y)=NULL
    colnames(y)=c("count","run","gene")
    y2<-y[,c(2,3,1)]
    return(y2)
  }

  message("rbinding new se data")
  se.new<-mclapply(paste("/scratch/mziemann/dee2/data/",org,"/",runs_done,"/",runs_done,".se.tsv",sep=""),se_func, mc.cores = CORES)
  se.new<-rbindlist(se.new)

  se_name=paste("/scratch/mziemann/dee2/mx/",org,"_se.tsv",sep="")
  if(file.exists(se_name)) {
   se<-read.table(se_name,row.names=F)
   se<-rbind(se,se.new)
  } else {
   se<-se.new
  }
  write.table(se,file=se_name,quote=F,row.names=F)
  system(paste("pbzip2 -fk " ,se_name))

#  message("appending new se data to db")
#  con <- dbConnect(RSQLite::SQLite(), dbname = paste(org,"_se",sep="") )
#  dbWriteTable(con, paste(org,"_se",sep="") , se, append = TRUE)

  message("rbinding new ke data")
  ke_func<-function(file){
    y<-read.table(file,header=T)
    y$run=strsplit(file,"/")[[1]][7]
    y2<-y[,c(6,1,4)]
    colnames(y2)=c("run","gene","count")
    return(y2)
  }

  ke.new<-mclapply(paste("/scratch/mziemann/dee2/data/",org,"/",runs_done,"/",runs_done,".ke.tsv",sep=""),ke_func, mc.cores = CORES)
  ke.new<-rbindlist(ke.new)

  ke_name=paste("/scratch/mziemann/dee2/mx/",org,"_ke.tsv",sep="")
  if(file.exists(ke_name)) {
   ke<-read.table(ke_name,row.names=F)
   ke<-rbind(ke,ke.new)
  } else {
   ke<-ke.new
  }
  write.table(ke,file=ke_name,quote=F,row.names=F)
  system(paste("pbzip2 -fk " ,ke_name))

#  message("appending new ke data to db")
#  con <- dbConnect(RSQLite::SQLite(), dbname = paste(org,"_ke",sep="") )
#  dbWriteTable(con, paste(org,"_ke",sep="") , ke, append = TRUE)

  message("rbinding new qc data")
  qc_func<-function(file){
    y<-read.colon.tsv(file)
    y$run=strsplit(file,"/")[[1]][7]
    y$metric=rownames(y)
    rownames(y)=NULL
    y2<-y[,c(2,1,3)]
    colnames(y2)=c("run","metric","count")
    return(y2)
  }

  qc.new<-mclapply(paste("/scratch/mziemann/dee2/data/",org,"/",runs_done,"/",runs_done,".qc",sep=""),qc_func, mc.cores = CORES)
  qc.new<-rbindlist(qc.new)

  qc_name=paste("/scratch/mziemann/dee2/mx/",org,"_qc.tsv",sep="")
  if(file.exists(qc_name)) {
   qc<-read.table(qc_name,row.names=F)
   qc<-rbind(qc,qc.new)
  } else {
   qc<-qc.new
  }
  write.table(qc,file=qc_name,quote=F,row.names=F)
  system(paste("pbzip2 -fk " ,qc_name))

#  message("appending new qc data to db")
#  con <- dbConnect(RSQLite::SQLite(), dbname = paste(org,"_qc",sep="") )
#  dbWriteTable(con, paste(org,"_qc",sep="") , qc, append = TRUE)

  for (run in runs_done) {
    FIN=paste("/scratch/mziemann/dee2/data/",org,"/",run,"/",run,".finished",sep="")
    VAL=paste("/scratch/mziemann/dee2/data/",org,"/",run,"/",run,".validated",sep="")
    DIR=paste("/scratch/mziemann/dee2/data/",org,"/",run,sep="")

    file.rename(FIN,VAL)
    system(paste("chmod 0544 ",DIR) )
   }

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
  }

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

  DATE=strsplit(as.character(file.info(FILES2[1])[,4])," ",fixed=T)[[1]][1]
  HEADER=paste("Updated",DATE)
  #colnames(z)=HEADER
  z<-z[order(rownames(z),decreasing=T ), ,drop=F]
  par(las=2) ; par(mai=c(1,2.5,1,0.5))
  MAX=max(as.numeric(z[,1]))+30000

  bb<-barplot( rbind( as.numeric(z$queued) , as.numeric(z$completed) ) ,
   names.arg=rownames(z) ,xlim=c(0,MAX),beside=T, main=HEADER, col=c("darkblue","red") ,
   horiz=T , las=1, cex.axis=1.3, cex.names=1.4,cex.main=1.4 )

  legend("topright", colnames(z), fill=c("darkblue","red") , cex=1.2)

  text( cbind(as.numeric(z[,1])+20000 ,as.numeric(z[,2])+20000 )  ,t(bb),labels=c(z[,1],z[,2]) ,cex=1.2)
  dev.off()
  system("scp -i ~/.ssh/cloud/cloud2.key dee_datasets.png ubuntu@118.138.240.228:/mnt/dee2_data/mx")

}
