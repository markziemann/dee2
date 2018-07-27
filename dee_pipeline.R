#!/usr/bin/env Rscript
setwd("/scratch/mziemann/dee2/code/") 

#library(SRAdb)
library(parallel)
library(data.table)
library(SRAdbV2)

CORES=ceiling(detectCores()/2)
for (org in c("ecoli" , "scerevisiae" , "celegans", "athaliana",  "rnorvegicus" , "celegans", "dmelanogaster", "drerio", "hsapiens", "mmusculus" ) ) {

  #create a list of NCBI taxa full names
  species_list<-c("3702","6239","7227","7955","562","9606", "10090", "10116", "4932")
 
 #now annotate the short names 
  names(species_list)<- c("athaliana", "celegans", "dmelanogaster", "drerio",
  "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae")

  taxa_name<-species_list[[org]]

  print(org)

  ###Set some directories
  CODEWD=getwd()
  DATAWD=paste(normalizePath("../data/"),org,sep="/")
  SRADBWD=normalizePath("../sradb/")
  MXDIR=normalizePath("../mx/")
  QUEUEWD=normalizePath("../queue/")

########################
# Get info from sradb vers 2
########################

  setwd(SRADBWD)

  oidx = Omicidx$new()
  query=paste( paste0('sample_taxon_id:', taxa_name), 'AND experiment_library_strategy : "rna-seq"')
  z = oidx$search(q=query,entity='full',size=100L)
  s = z$scroll()
  res = s$collate(limit = Inf)

  accessions<-as.data.frame(cbind(res$experiment_accession,res$study_accession,res$sample_accession,res$run_accession))
  colnames(accessions)=c("experiment","study","sample","run")
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

   #collect QC info - this is temporary and logic will be incorporated in future
   QC_summary="PASS"

   #Need to rearrange columns
   GSE<-function(i) {res=grepl("GSE",i) ; if (res == FALSE) {j="NA"} else {j=i } ; j }
   GSE_accession<-as.vector(sapply(res$study_GEO,GSE))

   GSM<-function(i) {res=grepl("GSM",i) ; if (res == FALSE) {j="NA"} else { j=i} ; j }
   GSM_accession<-as.vector(sapply(res$sample_GEO,GSM))

   #extract out the important accessions in order
   x2<-as.data.frame(cbind(res$run_accession,QC_summary,res$experiment_accession,res$sample_accession,
   res$study_accession, GSE_accession, GSM_accession))

   colnames(x2)<-c("SRR_accession","QC_summary","SRX_accession","SRS_accession",
   "SRP_accession","GSE_accession","GSM_accession")

   #write out the accession number info and upload to webserver
   write.table(x2,file=paste(SRADBWD,"/",org,"_metadata.tsv.cut",sep=""),quote=F,sep="\t",row.names=F)
   SCP_COMMAND=paste("scp -i ~/.ssh/cloud/id_rsa ", paste(SRADBWD,"/",org,"_metadata.tsv.cut",sep="") ," ubuntu@118.138.240.228:/mnt/dee2_data/metadata")
   system(SCP_COMMAND)

   #now attach the additional metadata and upload
   x<-res[, !(colnames(res) %in% c("run_accession", "QC_summary","experiment_accession","sample_accession","study_accession","submission_accession","GSE_accession","GSM_accession"))]
   x<-as.data.frame(cbind(x2,x))
 
   x<-x[which(x$SRR_accession %in% runs_done),]
   x <- apply(x,2,as.character)
   x<-gsub("\r?\n|\r", " ", x)
   write.table(x,file=paste(SRADBWD,"/",org,"_metadata.tsv",sep=""),quote=F,sep="\t",row.names=F)
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

  DATE=strsplit(as.character(file.info(FILES2[1])[,4])," ",fixed=T)[[1]][1]
  HEADER=paste("Updated",DATE)
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
