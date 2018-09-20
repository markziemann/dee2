getDee2Metadata<-function(species,outfile=NULL){
  metadataURL=paste("http://dee2.io/metadata/",species,"_metadata.tsv.cut",sep="")
  if(is.null(outfile)){
    metadataname=tempfile()
  } else {
    metadataname=outfile
    if(!grepl(".tsv$",metadataname)){metadataname=paste0(metadataname,".tsv")}
  }
  download.file(metadataURL, destfile=metadataname)
  mdat<-read.table(metadataname,header=T)

  if(is.null(outfile)){unlink(metadataname)}
  return(mdat)
}


queryDee2<-function(species, SRRvec) {
  present<-SRRvec[which(SRRvec %in% dee2metadata.celegans$SRR_accession)]
  absent<-SRRvec[-which(SRRvec %in% dee2metadata.celegans$SRR_accession)]
  dat <- list("present" = present, "absent" = absent)
  return(dat)
}


loadGeneCounts<-function(zipname){
  CM="GeneCountMatrix.tsv"
  TF=tempfile()
  unzip(zipname, files = CM, exdir = tempdir() )
  mxname<-paste0(tempdir(),"/",CM)
  file.rename(mxname,TF)
  dat <- read.table(TF,row.names=1,header=T)
  unlink(TF)
  return(dat)
}


loadTxCounts<-function(zipname){
  CM="TxCountMatrix.tsv"
  TF=tempfile()
  unzip(zipname, files = CM, exdir = tempdir() )
  mxname<-paste0(tempdir(),"/",CM)
  file.rename(mxname,TF)
  dat <- read.table(TF,row.names=1,header=T)
  unlink(TF)
  return(dat)
}

loadQcMx<-function(zipname){
  CM="QC_Matrix.tsv"
  TF=tempfile()
  unzip(zipname, files = CM, exdir = tempdir() )
  mxname<-paste0(tempdir(),"/",CM)
  file.rename(mxname,TF)
  dat <- read.table(TF,row.names=1,header=T,fill=T)
  unlink(TF)
  return(dat)
}

getDEE2<-function(species, SRRvec, outfile=NULL ,
  baseURL="http://dee2.io/cgi-bin/request.sh?"){
  SRRvec<-gsub(" ","",SRRvec)
  llist<-paste0("&x=",paste(SRRvec,collapse = "&x="))
  murl <- paste0(baseURL,"org=",species, llist)
  if(is.null(outfile)){
        zipname=tempfile()
  } else {
        zipname=outfile
        if(!grepl(".zip$",zipname)){zipname=paste0(zipname,".zip")}
  }
  download.file(murl, destfile=zipname)

  GeneCounts<-loadGeneCounts(zipname)
  TxCounts<-loadTxCounts(zipname)
  QcMx<-loadQcMx(zipname)
  dat <- list("GeneCounts" = GeneCounts, "TxCounts" = TxCounts, "QcMx" = QcMx)

  if(is.null(outfile)){unlink(zipname)}
  return(dat)
}

#mytable<-getDEE("Ecoli",c("SRR922260","SRR922261"))
#data is returned as a list of three dataframes
