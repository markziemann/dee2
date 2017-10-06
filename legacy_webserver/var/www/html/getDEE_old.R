loadDEE<-function(tarname){
  CM="CountMatrix.xls"
  TF=tempfile()
  untar(tarname, files = CM, compressed = "gzip", exdir = tempdir() )
  mxname<-paste0(tempdir(),"/",CM)
  file.rename(mxname,TF)
  dat <- read.table(TF,row.names=1,header=T)
  unlink(TF)
  return(dat)
}


getDEE<-function(species, SRRvec, outfile=NULL ,
	baseURL="http://172.16.115.17/cgi-bin/extract.sh?"){
  llist<-paste0("&x=",paste(SRRvec,collapse = "&x="))
  murl <- paste0(baseURL,"org=",species, llist)
  if(is.null(outfile)){
	tarname=tempfile()
  } else {
	tarname=outfile
	if(!grepl(".tar.gz$",tarname)){tarname=paste0(tarname,".tar.gz")}
  }
  download.file(murl, destfile=tarname)
  
  dat<-loadDEE(tarname)
  if(is.null(outfile)){unlink(tarname)}
  return(dat)  
}
#mytable<-getDEE("Ecoli",c("SRR922260","SRR922261"))
