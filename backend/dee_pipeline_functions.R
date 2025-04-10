CORES=ceiling(8)

rowcnt2<-function( file) { z<-system(paste("cat ",file, "| cut -f1 | sed 1d | sort -u | wc -l" ) , intern=TRUE) ; z}

qc_analysis<-function(srr,org) {
 QCFILE=paste("../data/",org,"/",srr,"/",srr,".qc",sep="")
 QCLFILE=paste("../data/",org,"/",srr,"/",srr,".qcl",sep="")

 q <- read.csv(QCFILE,sep=":",stringsAsFactors=FALSE,header=FALSE)

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
} else if ( file.info(SEZ)[1]<1000 ) {
  DELETE=1
} else {
  se<-NULL
  tryCatch({ se<-fread(SEZ,header=F) },
    error = function(e) { DELETE=1 },
    finally = {
      gro<-dim( se )[1]
      write(SRR,file=G)
      write(se$V2,file=G,append=T,ncolumns=1)    
      if ( is.null( se )) {
        DELETE=1
      } else if ( gro!=gre ) {
        DELETE=1
      }
    }
  )
}

if ( !file.exists( KEZ ) ) {
  DELETE=1
} else if ( file.info(KEZ)[1]<1000 ) {
  DELETE=1
} else {
  ke<-NULL                 
  tryCatch({ ke<-fread(KEZ,header=T) },
    error = function(e) { DELETE=1 },
    finally = {
      tro<-dim( ke )[1]
      write(SRR,file=TX)
      write.table(ke[,4],file=TX,append=T,row.names = F,col.names=F)
      if ( is.null( ke )) {
        DELETE=1
      } else if ( tro!=tre ) {
        DELETE=1
      }
    }
  )
}

if ( !file.exists( QC ) ) {
  DELETE=1
} else if (file.info(QC)[1]<100 ) {
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

