#!/usr/bin/env Rscript

setwd("/mnt/md0/dee2/code")

#library(SRAdb)
library(parallel)
library(data.table)
library(SRAdbV2)
library(R.utils)

IPADD="118.138.234.131"
#simple rowcount function
rowcnt2<-function( file) { z<-system(paste("wc -l < ",file) , intern=TRUE) ; z}

CORES=ceiling(detectCores()/2)

qc_analysis<-function(org,srr) {
 QCFILE=paste("../data/",org,"/",srr,"/",srr,".qc",sep="")
 QCLFILE=paste("../data/",org,"/",srr,"/",srr,".qcl",sep="")
 q<-read.table(QCFILE,stringsAsFactors=F)
 q<-as.data.frame(t(   data.frame(strsplit(as.character(q$V1),split=":"),stringsAsFactors=F) ))
 rownames(q)=1:nrow(q)
 q<-q[1:28,]
 #get number of genes to determine required number of reads
 NumReadsQcPass=as.numeric(as.character(q[10,2]))
 NumReadsQcPass_PerGene=NumReadsQcPass/numgenes
 QcPassRate=as.numeric(as.character(gsub("%","",q[11,2])))
 STAR_UniqMapRate=as.numeric(as.character(gsub("%","",q[24,2])))
 STAR_AssignRate=as.numeric(as.character(gsub("%","",q[25,2])))
 Kallisto_MapRate=as.numeric(as.character(gsub("%","",q[28,2])))
 Kallisto_MappedReads=as.numeric(as.character(q[27,2]))
 Kallisto_MappedReads_PerGene=Kallisto_MappedReads/numgenes
 STAR_AssignedReads=as.numeric(as.character(q[23,2]))
 STAR_AssignedReads_PerGene=STAR_AssignedReads/numgenes
 WARN=FAIL=NULL
  # QC checks
  if ( is.na(NumReadsQcPass) == TRUE ) {
    FAIL=paste(FAIL,"1",sep=",")
  } else {
    if ( NumReadsQcPass_PerGene < 500 ) {
      if ( NumReadsQcPass_PerGene < 50 ) {
        FAIL=paste(FAIL,"1",sep=",")
      } else {
        WARN=paste(WARN,"1",sep=",")
      }
    }
  }

  if ( is.na(QcPassRate) == TRUE ) {
    FAIL=paste(FAIL,"2",sep=",")
    } else {
    if ( QcPassRate < 80 ) {
      if ( QcPassRate < 60 ) {
        FAIL=paste(FAIL,"2",sep=",")
      } else {
        WARN=paste(WARN,"2",sep=",")
      }
    }
  }

  if ( is.na(STAR_UniqMapRate) == TRUE ) {
    FAIL=paste(FAIL,"3",sep=",")
  } else {
    if ( STAR_UniqMapRate < 70 ) {
      if (STAR_UniqMapRate < 50 ) {
        FAIL=paste(FAIL,"3",sep=",")
      } else {
         WARN=paste(WARN,"3",sep=",")
      }
    }
  }

  if ( is.na(STAR_AssignRate) == TRUE ) {
    FAIL=paste(FAIL,"4",sep=",")
  } else {
    if ( STAR_AssignRate < 60 ) {
      if (STAR_AssignRate < 40 ) {
        FAIL=paste(FAIL,"4",sep=",")
      } else {
         WARN=paste(WARN,"4",sep=",")
      }
    }
  }

  if ( is.na(STAR_AssignedReads) == TRUE ) {
    FAIL=paste(FAIL,"5",sep=",")
  } else {
    if ( STAR_AssignedReads_PerGene < 500 ) {
      if ( STAR_AssignedReads_PerGene < 50 ) {
        FAIL=paste(FAIL,"5",sep=",")
      } else {
        WARN=paste(WARN,"5",sep=",")
      }
    }
  }

  if ( is.na(Kallisto_MapRate) == TRUE ) {
    FAIL=paste(FAIL,"6",sep=",")
  } else {
    if ( Kallisto_MapRate < 60 ) {
      if (Kallisto_MapRate < 40 ) {
        FAIL=paste(FAIL,"6",sep=",")
      } else {
         WARN=paste(WARN,"6",sep=",")
      }
    }
  }

  if ( is.na(Kallisto_MappedReads) == TRUE ) {
    FAIL=paste(FAIL,"7",sep=",")
  } else {
    if ( Kallisto_MappedReads_PerGene < 500 ) {
      if ( Kallisto_MappedReads_PerGene < 50 ) {
        FAIL=paste(FAIL,"7",sep=",")
      } else {
        WARN=paste(WARN,"7",sep=",")
      }
    }
  }

  #incorporate a warning about low correlation data
  mycor<-cors[which(rownames(cors) %in% srr),]
  if ( mycor < 0.5 ) {
    WARN=paste(WARN,"8",sep=",")
  }
  mycor=c("DatasetCorrel",round(mycor,digits=2))

  #attach mycor to the qc file
  q<-rbind(q,t(mycor))


  if (!is.null(FAIL)) {
    FAIL=sub(',','',FAIL)
    QCRES=paste("FAIL(",FAIL,")",sep="")
  } else if (!is.null(WARN)) {
    WARN=sub(',','',WARN)
    QCRES=paste("WARN(",WARN,")",sep="")
  } else {
    QCRES="PASS"
  }

  qcres=c("QC_SUMMARY",QCRES)
  q<-rbind(q,t(qcres))
  write.table(q,file=QCFILE,quote=F,sep=':',col.names=F,row.names=F)
  write.table(q$V2,file=QCLFILE,quote=F,sep=':',col.names=F,row.names=F)
  QCRES
}


getmean<-function(org) {
# Generate an "average" sample that can be used for correlation QC analysis
# Make the solution scalable from 100 datasets to 1M.

# load libraries
library(data.table)
library(reshape2)


meanfile=paste("../mx/",org,"_means.tsv",sep="")

DOIT=0
if( !file.exists(meanfile) ) {
  DOIT=1
} else {
  MODTIME=as.numeric(difftime(Sys.time() ,file.mtime(meanfile),units="s"))
  if (MODTIME>60*60*24*30) {
    DOIT=1
  }
}

if ( DOIT==1) {
  # read in metadata
  mdat=paste("../sradb/",org,"_metadata.tsv.cut",sep="")
  m<-read.table(mdat,header=T,sep="\t",quote="",fill=F)

  # make a list of samples to use
  p<-m[grep("PASS",m$QC_summary),1]

  # make blocks of datasets to analyse
  num_blocks=ceiling(length(p)/1000)
  block_size=floor(length(p)/num_blocks)

  # grow a data frame with colums
  df=NULL
  TSV=paste("../mx/",org,"_se.tsv.bz2",sep="")
  for ( i in 1:num_blocks) {
    s<-sample(p,block_size)
    p<-setdiff(p,s)
    d<-fread(TSV)[which(fread(TSV)$V1 %in% s),]
    d<-as.matrix(acast(d, V2~V1, value.var="V3"))
    d<-d/colSums(d)*1000000
    d<-rowMeans(d)
    df<-cbind(df,d)
  }

  #get the means and format as dataframe
  df<-as.data.frame(rowMeans(df))
  colnames(df)="mean"
  write.table(df,file=meanfile,quote=F,sep="\t")

} else {
  df<-read.table(meanfile,sep="\t",header=T,row.names=1)
}
df
}

#check the folder contents and validate
check_contents<-function(d,gre,tre) {
SRR=sapply(strsplit(d,"/"),"[[",7)
DELETE=0
SE=paste(d,"/",SRR,".se.tsv",sep="")
G=paste(d,"/",SRR,"_gene.cnt",sep="")
KE=paste(d,"/",SRR,".ke.tsv",sep="")
TX=paste(d,"/",SRR,"_tx.cnt",sep="")
QC=paste(d,"/",SRR,".qc",sep="")
LOG=paste(d,"/",SRR,".log",sep="")
FIN=paste(d,"/",SRR,".finished",sep="")
VAL=paste(d,"/",SRR,".validated",sep="")

if ( !file.exists( SE ) ) {
  DELETE=1
} else {
  se<-fread( SE ,header=F)
  gro<-dim( se )[1]
  write(SRR,file=G)
  write(se$V2,file=G,append=T,ncolumns=1)    
  gzip(SE)
  if ( gro!=gre ) {
    DELETE=1
  }
}

if ( !file.exists( KE ) ) {
  DELETE=1
} else {
  ke<-fread( KE ,header=T)
  tro<-dim( ke )[1]
  write(SRR,file=TX)
  write.table(ke[,4],file=TX,append=T,row.names = F,col.names=F)
  gzip(KE)
  if ( tro!=tre ) {
    DELETE=1
  }
}

if ( !file.exists( QC ) ) {
  DELETE=1
} else {
  qc<-fread( QC ,header=F)
  qcro<-dim( fread( QC ,header=F) )[1]
  if ( qcro<29 ) {
    DELETE=1
  }
}

if ( !file.exists( LOG ) ) {
  DELETE=1
} else {
  LOGOK=length(grep("completed mapping pipeline successfully",readLines( LOG ) ) )
  if ( LOGOK!=1 ) {
    DELETE=1
  }
}

if ( DELETE==1 ) {
  unlink(d,recursive=TRUE)
} else {
  file.rename(FIN, VAL)
}
}


#start the analysis
for (org in c("scerevisiae") ) {
#for (org in c("ecoli", "scerevisiae" , "athaliana",  "rnorvegicus" , "celegans", "dmelanogaster", "drerio", "hsapiens", "mmusculus" ) ) {
  #create a list of NCBI taxa full names
  species_list<-c("3702","6239","7227","7955","562","9606", "10090", "10116", "4932")
 
  #now annotate the short names 
  names(species_list)<- c("athaliana", "celegans", "dmelanogaster", "drerio",
  "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae")

  taxa_name<-species_list[[org]]

  print(org)

  #number of protein coding genes from ensembl
  numgenes<-c(27655,20362,13918,25903,4140,20338,22598,22250,6692)
  names(numgenes)<- c("athaliana", "celegans", "dmelanogaster", "drerio",
  "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae")
  numgenes=numgenes[[org]]

  ###Set some directories
  CODEWD=getwd()
  DATAWD=paste(normalizePath("../data/"),org,sep="/")
  SRADBWD=normalizePath("../sradb/")
  MXDIR=normalizePath("../mx/")
  QUEUEWD=normalizePath("../queue/")

  ########################
  # Get metadata mod date
  ########################

  RDA=paste(SRADBWD,"/",org,".RData",sep="")
  GETNEW=0

  if(!file.exists(RDA)){ 
    GETNEW=1
  } else {
    MODTIME=as.numeric(difftime(Sys.time() ,file.mtime(RDA),units="s"))
    if (MODTIME>60*60*24*30) {
      GETNEW=1
    }
  }

  if ( GETNEW==1 ) { 
    # Get info from sradb vers 2
    message("part A")
    if ( exists("s") ) { s$reset() }

    oidx=z=s=res=accessions=runs=NULL
    message("part B")
    oidx = Omicidx$new()
    message("part C")
    query=paste( paste0('sample.taxon_id:', taxa_name), 'AND experiment.library_strategy : "rna-seq"')
    message("part D")
    z = oidx$search(q=query,entity='full',size=100L)
    message("part E")
    s = z$scroll()
    message("part F")
    res = s$collate(limit = Inf)
    message("part G")
    accessions<-as.data.frame(cbind(res$experiment.accession,res$study.accession,res$sample.accession,res$accession))
    colnames(accessions)=c("experiment","study","sample","run")
    runs<-accessions$run
    s$reset()
    save.image(file = RDA)
  } else {
    message("using existing metadata")
    load(RDA)
  }

  ########################
  # Now determine which datasets have already been processed and completed
  ########################

  folders<-list.files(DATAWD,full.names=T,pattern="RR")
  fin_new<-paste(DATAWD,sapply(strsplit(list.files(path = DATAWD, pattern = "finished" , full.names = T, recursive = T),"/"),"[[",7),sep="/")
  fin_new<-fin_new[grep("RR",fin_new)]
  val_old<-paste(DATAWD,sapply(strsplit(list.files(path = DATAWD, pattern = "validated" , full.names = T, recursive = T),"/"),"[[",7),sep="/")

  # delete folders that are not expected as they are not in known runs
  expected_folders<-paste(DATAWD,runs,sep="/")
  unexpected_folders<-setdiff(folders,expected_folders)
  write(unexpected_folders,file=paste(DATAWD,"/unexpected_folders.txt",sep=""))

  # delete folders without finished or validated files
  allocated<-union(fin_new,val_old)
  unalloc<-setdiff(folders,allocated)
  unlink(unalloc,recursive=TRUE)

  #Expected rows of se and ke files
  gre<-dim( read.table(paste(DATAWD,"/rownames_gene.txt",sep=""),header=F) )[1]
  tre<-dim( read.table(paste(DATAWD,"/rownames_tx.txt",sep=""),header=T) )[1]

  # run the check of new datasets
  lapply(fin_new,check_contents,gre,tre)

  val<-paste(DATAWD,sapply(strsplit(list.files(path = DATAWD, pattern = "validated" , full.names = T, recursive = T),"/"),"[[",7),sep="/")
  runs_done<-sapply(strsplit(val,"/"),"[[",7)

  av<-getmean(org)

  se_files<-paste( DATAWD , "/" , runs_done , "/" , runs_done , ".se.tsv.gz",sep="")

  corav<-function(x) {
    xx<-read.table(x)
    cor(merge(xx,av,by=0)[2:3],method="p")[1,2]
  }

  cors<-data.frame(as.numeric(  mclapply( se_files , corav , mc.cores=CORES) ) )
  rownames(cors)<-sapply(strsplit(se_files,"/"),"[[",7)

  runs_done<-sapply(strsplit(val,"/"),"[[",7)
  print(paste(length(runs_done),"runs completed"))
  runs_todo<-base::setdiff(runs, runs_done)
  print(paste(length(runs_todo),"requeued runs"))

  #Update queue on webserver
  queue_name=paste(QUEUEWD,"/",org,".queue.txt",sep="")
  write.table(runs_todo,queue_name,quote=F,row.names=F,col.names=F)
  SCP_COMMAND=paste("scp -i ~/.ssh/monash/cloud2.key ",queue_name ," ubuntu@118.138.234.131:~/Public")
  system(SCP_COMMAND)

  #Update metadata on webserver
  accessions_done<-accessions[which(accessions$run %in% runs_done),]
  write.table(accessions_done,file=paste(SRADBWD,"/",org,"_accessions.tsv",sep=""),quote=F,row.names=F)

  save.image(file = paste(org,".RData",sep=""))

  #collect QC info - this is temporary and logic will be incorporated in future
  QC_summary="BLANK"

  #Need to rearrange columns
  GSE<-function(i) {res=grepl("GSE",i) ; if (res == FALSE) {j="NA"} else {j=i } ; j }
  GSE_accession<-as.vector(sapply(res$study.GEO,GSE))

  GSM<-function(i) {res=grepl("GSM",i) ; if (res == FALSE) {j="NA"} else { j=i} ; j }
  GSM_accession<-as.vector(sapply(res$sample.GEO,GSM))

  #extract out the important accessions in order
  ##x2<-as.data.frame(cbind(res$run.accession,QC_summary,res$experiment.accession,res$sample.accession,res$study.accession, GSE_accession, GSM_accession))
  x2<-as.data.frame(cbind(res$accession,QC_summary,res$experiment.accession,res$sample.accession,res$study.accession, GSE_accession, GSM_accession,res$experiment.title))

  colnames(x2)<-c("SRR_accession","QC_summary","SRX_accession","SRS_accession",
  "SRP_accession","GSE_accession","GSM_accession","experiment_title")

  #write out the accession number info and upload to webserver
  write.table(x2,file=paste(SRADBWD,"/",org,"_metadata.complete.tsv.cut",sep=""),quote=F,sep="\t",row.names=F)
  x2<-x2[which(x2$SRR_accession %in% runs_done),]

  x2$QC_summary<-unlist(lapply(x2$SRR_accession , qc_analysis, org=org))

  #here we rsync files to server in chunks of 1000
  rsync<-function(d,org) {
    while ( length(d)>0 ) {
      CHUNKSIZE=1000
      if ( length(d)>CHUNKSIZE ) {
        chunk<-paste(d[1:CHUNKSIZE],collapse=" ")
        d<-setdiff(d,d[1:1000])
      } else {
        chunk<-paste(d[1:length(d)],collapse=" ")
        d<-setdiff(d,d[1:1000])
      }
      CMD=paste('rsync -avzh -e \"ssh -i  ~/.ssh/monash/cloud2.key \" ', chunk ,' ubuntu@118.138.234.131:/dee2_data/data/',org,sep="")
      system(CMD)
    }
  } 
  d<-val
  rsync(d,org)

  write.table(x2,file=paste(SRADBWD,"/",org,"_metadata.tsv.cut",sep=""),quote=F,sep="\t",row.names=F)

  SCP_COMMAND=paste("scp -i ~/.ssh/monash/cloud2.key", paste(SRADBWD,"/",org,"_metadata.tsv.cut",sep="") ," ubuntu@118.138.234.131:/mnt/dee2_data/metadata")
  system(SCP_COMMAND)

  save.image(file = paste(org,".RData",sep=""))

  #now attach the additional metadata and upload
  x<-res[, !(colnames(res) %in% c("QC_summary","experiment.accession","sample.accession","study.accession","submission.accession","GSE_accession","GSM_accession"))]
  x<-merge(x2,x,by.x="SRR_accession",by.y="accession")
 
  x<-x[which(x$SRR_accession %in% runs_done),]
  x <- apply(x,2,as.character)
  x<-gsub("\r?\n|\r", " ", x)
  write.table(x,file=paste(SRADBWD,"/",org,"_metadata.tsv",sep=""),quote=F,sep="\t",row.names=F)
  SCP_COMMAND=paste("scp -i ~/.ssh/monash/cloud2.key ", paste(SRADBWD,"/",org,"_metadata.tsv",sep="") ," ubuntu@118.138.234.131:/mnt/dee2_data/metadata")
  system(SCP_COMMAND)

  save.image(file = paste(org,".RData",sep=""))

  png("dee_datasets.png",width=580,height=580)

  FILES1<-list.files(pattern="*queue.txt$",path="/mnt/md0/dee2/queue/",full.names=T)
  x<-as.data.frame(sapply(FILES1,rowcnt2),stringsAsFactors=FALSE)
  rownames(x)=c("A. thaliana","C. elegans","D. melanogaster","D. rerio","E. coli","H. sapiens","M. musculus","R. norvegicus","S. cerevisiae")
  colnames(x)="queued"

  FILES2<-list.files(pattern="*accessions.tsv$",path="/mnt/md0/dee2/sradb/",full.names=T)
  y<-as.data.frame(sapply(FILES2,rowcnt2),stringsAsFactors=FALSE)
  rownames(y)=c("A. thaliana","C. elegans","D. melanogaster","D. rerio","E. coli","H. sapiens","M. musculus","R. norvegicus","S. cerevisiae")
  colnames(y)="completed"

  z<-merge(x,y,by=0)
  rownames(z)=z$Row.names
  z$Row.names=NULL

  DATE=strsplit(as.character(file.info(FILES2[1])[,6])," ",fixed=T)[[1]][1]
  HEADER=paste("Updated",DATE)
  z<-z[order(rownames(z),decreasing=T ), ,drop=F]
  par(las=2) ; par(mai=c(1,2.5,1,0.5))
  MAX=max(as.numeric(z[,1]))+100000

  bb<-barplot( rbind( as.numeric(z$queued) , as.numeric(z$completed) ) ,
   names.arg=rownames(z) ,xlim=c(0,MAX), beside=T, main=HEADER, col=c("darkblue","red") ,
   horiz=T , las=1, cex.axis=1.3, cex.names=1.4, cex.main=1.4 )

  legend("topright", colnames(z), fill=c("darkblue","red") , cex=1.2)

  text( cbind(as.numeric(z[,1])+50000 ,as.numeric(z[,2])+50000 )  ,t(bb),labels=c(z[,1],z[,2]) ,cex=1.2)
  dev.off()
  system("scp -i ~/.ssh/monash/cloud2.key dee_datasets.png ubuntu@118.138.234.131:/mnt/dee2_data/mx")

}
