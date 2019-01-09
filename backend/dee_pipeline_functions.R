CORES=ceiling(detectCores())

rowcnt2<-function( file) { z<-system(paste("wc -l < ",file) , intern=TRUE) ; z}

qc_analysis<-function(srr,org) {
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
  if ( is.na(mycor) ) {
    WARN=paste(WARN,"8",sep=",")
  } else if ( mycor < 0.5 ) {
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


getmean_old_old<-function(org) {
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


 

getmean_old<-function(org) {
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
  if (MODTIME>60*60*24*30*6) {
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

  dt<-fread(TSV)

  setkey(dt, V1)

  for ( i in 1:num_blocks) {
    s<-sample(p,block_size)
    p<-setdiff(p,s)
    d<-dt[which(dt$V1 %in% s)]
    d<-as.matrix(acast(d, V2~V1, value.var="V3"))
    d<-d/colSums(d)*1000000
    d<-rowMeans(d)
    df<-cbind(df,d)
  }


  #get the means and format as dataframe
  df<-as.data.frame(rowMeans(df))
  colnames(df)="mean"
  write.table(df,file=meanfile,quote=F,sep="\t")
  rm(dt)
  gc()

} else {
  df<-read.table(meanfile,sep="\t",header=T,row.names=1)
}
df
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
  if (MODTIME>60*60*24*30*6) {
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
  num_blocks=ceiling(length(p)/100)
  block_size=floor(length(p)/num_blocks)

  chunks=NULL
  for (n in 1:num_blocks) {
    chunk<-sample(p,block_size)
    p<-setdiff(p,chunk)
    chunks<-c(list(chunk),chunks)
  }

  chunkmean<-function(chunk) {
  df=NULL
  for ( SRR in chunk ) {
    sefile=paste("../data/",org,"/",SRR,"/",SRR,".se.tsv.gz",sep="")
    if ( file.exists(sefile) ) {
      se<-fread(sefile,header=F)
      x<-as.matrix(se$V2/sum(se$V2)*1000000)
      colnames(x)=SRR
      df = cbind(df,x)
    }
  }
  rowMeans(df)
  }

  chunkmeans<-mclapply(chunks,chunkmean,mc.cores=CORES)
  means<-as.data.frame(matrix(unlist(chunkmeans), nrow=length(unlist(chunkmeans[1]))))
  means<-as.data.frame(rowMeans(means))
  colnames(means)="means"
  gene_names_file=paste("../data/",org,"/rownames_gene.txt",sep="")
  genes<-readLines(gene_names_file)
  genes<-genes[2:length(genes)]
  rownames(means)=genes
  write.table(means,file=meanfile,quote=F,sep="\t")
} else {
  means<-read.table(meanfile,sep="\t",header=T,row.names=1)
}
means
}


# check the folder contents and validate
check_contents<-function(d,gre,tre) {
SRR=sapply(strsplit(d,"/"),"[[",7)
DELETE=0
SE=paste(d,"/",SRR,".se.tsv",sep="")
SEZ=paste(d,"/",SRR,".se.tsv.gz",sep="")
G=paste(d,"/",SRR,"_gene.cnt",sep="")
KE=paste(d,"/",SRR,".ke.tsv",sep="")
KEZ=paste(d,"/",SRR,".ke.tsv.gz",sep="")
TX=paste(d,"/",SRR,"_tx.cnt",sep="")
QC=paste(d,"/",SRR,".qc",sep="")
LOG=paste(d,"/",SRR,".log",sep="")
FIN=paste(d,"/",SRR,".finished",sep="")
VAL=paste(d,"/",SRR,".validated",sep="")

# compress if necessary
if ( file.exists( SE ) ) {
  gzip(SE,overwrite=T)
}
if ( file.exists( KE ) ) {
  gzip(KE,overwrite=TRUE)
}

if ( !file.exists( SEZ ) ) {
  DELETE=1
} else {
  se<-fread( SEZ ,header=F)
  gro<-dim( se )[1]
  write(SRR,file=G)
  write(se$V2,file=G,append=T,ncolumns=1)    
  if ( gro!=gre ) {
    DELETE=1
  }
}

if ( !file.exists( KEZ ) ) {
  DELETE=1
} else {
  ke<-fread( KEZ ,header=T)
  tro<-dim( ke )[1]
  write(SRR,file=TX)
  write.table(ke[,4],file=TX,append=T,row.names = F,col.names=F)
  if ( tro!=tre ) {
    DELETE=1
  }
}

if ( !file.exists( QC ) ) {
  DELETE=1
} else {
  qc<-fread( QC ,header=F,sep=":")
  qcro=nrow(qc)
  if ( qcro<29 ) {
    DELETE=1
  }
}

if ( !file.exists( LOG ) ) {
  DELETE=1
} else {
  LOGOK=length(grep("completed mapping pipeline successfully",  tail(readLines( LOG ) )    ) )
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

