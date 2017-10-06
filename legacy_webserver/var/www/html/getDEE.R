loadDEE<-function(zipname){
  CM="CountMatrix.xls"
  TF=tempfile()
  unzip(zipname, files = CM, exdir = tempdir() )
  mxname<-paste0(tempdir(),"/",CM)
  file.rename(mxname,TF)
  dat <- read.table(TF,row.names=1,header=T)
  unlink(TF)
  return(dat)
}


getDEE<-function(species, SRRvec, outfile=NULL ,
	baseURL="http://dee.bakeridi.edu.au/cgi-bin/extract.sh?"){
  llist<-paste0("&x=",paste(SRRvec,collapse = "&x="))
  murl <- paste0(baseURL,"org=",species, llist)
  if(is.null(outfile)){
	zipname=tempfile()
  } else {
	zipname=outfile
	if(!grepl(".zip$",zipname)){zipname=paste0(zipname,".zip")}
  }
  download.file(murl, destfile=zipname)
  
  dat<-loadDEE(zipname)
  if(is.null(outfile)){unlink(zipname)}
  return(dat)  
}
#mytable<-getDEE("Ecoli",c("SRR922260","SRR922261"))
