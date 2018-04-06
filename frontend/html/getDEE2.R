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
  TMP="/tmp/tmp.html"
  download.file(baseURL,destfile=TMP)
  a<-readChar("/tmp/tmp.html", file.info(TMP)$size)
  a<-unlist(strsplit(a," "))
  a<-a[grep("http",a)][2]
  a<-strsplit(a,"/")[[1]][3]
  ip<-gsub('\"',"",a)
  baseURL=gsub("dee2.io",ip,baseURL)

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
