#!/usr/bin/env Rscript

# if an error occurs on here: mclapply(fin_new,check_contents,gre,tre,mc.cores=CORES)
# then need to run this rm -rf `find ../data/org | grep tmp$ | head -1000`

setwd("/mnt/md0/dee2/code")

library("parallel")
library("data.table")
library("R.utils")
library("reutils")
library("XML")
library("rjson")

IPADD="118.138.235.221"

CORES=5

#start the analysis
for ( org in c( "ecoli" )) {

#for ( org in c(   "athaliana", "celegans", "dmelanogaster", "drerio",
#"rnorvegicus", "scerevisiae" , "mmusculus", "ecoli", "hsapiens" )) {

#args = commandArgs(trailingOnly=TRUE)
#org=args[1]

species_list <- c("3702","6239","7227","7955","562","9606", "10090", "10116", "4932")
 
#now annotate the short names 
names(species_list) <- c("athaliana", "celegans", "dmelanogaster", "drerio",
"ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae")

taxa_name <- species_list[[org]]

species_names <- c("'Arabidopsis thaliana'","'Caenorhabditis elegans'","'Drosophila melanogaster'","'Danio rerio'",
  "'Escherichia coli'","'Homo sapiens'", "'Mus musculus'", "'Rattus norvegicus'", "'Saccharomyces cerevisiae'")
 #now annotate the short names 
names(species_names)<- c("athaliana", "celegans", "dmelanogaster", "drerio",
  "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae")
species_name <- species_names[[org]]

print(org)

#number of protein coding genes from ensembl
numgenes <- c(27655,20362,13918,25903,4140,20338,22598,22250,6692)
names(numgenes) <- c("athaliana", "celegans", "dmelanogaster", "drerio",
"ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae")
numgenes = numgenes[[org]]

###Set some directories
CODEWD = getwd()
DATAWD = paste(normalizePath("../data/"),org,sep="/")
SRADBWD = normalizePath("../sradb/")
MXDIR = normalizePath("../mx/")
QUEUEWD = normalizePath("../queue/")

########################
# Get metadata 
########################

CSV=paste(SRADBWD,"/",org,".csv",sep="")

if (!file.exists(CSV)) {
  file.rename(CSV,paste(CSV,".old",sep=""))
}

SRA_DL_CMDS <- c(
' wget -O ../sradb/athaliana.csv \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="transcriptomic"[Source] AND cluster_public[prop] AND "Arabidopsis thaliana"[Organism]\' ' ,
' wget -O ../sradb/celegans.csv \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="transcriptomic"[Source] AND cluster_public[prop] AND "Caenorhabditis elegans"[Organism]\' ' ,
' wget -O ../sradb/dmelanogaster.csv \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="transcriptomic"[Source] AND cluster_public[prop] AND "Drosophila melanogaster"[Organism]\' ' ,
' wget -O ../sradb/drerio.csv \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="transcriptomic"[Source] AND cluster_public[prop] AND "Danio rerio"[Organism]\' ' ,
' wget -O ../sradb/ecoli.csv \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="transcriptomic"[Source] AND cluster_public[prop] AND "Escherichia coli"[Organism]\' ' ,
' wget -O ../sradb/hsapiens.csv \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="transcriptomic"[Source] AND cluster_public[prop] AND "Homo sapiens"[Organism]\' ' ,
' wget -O ../sradb/mmusculus.csv \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="transcriptomic"[Source] AND cluster_public[prop] AND "Mus musculus"[Organism]\' ' ,
' wget -O ../sradb/rnorvegicus.csv \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="transcriptomic"[Source] AND cluster_public[prop] AND "Rattus norvegicus"[Organism]\' ' ,
' wget -O ../sradb/scerevisiae.csv \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="transcriptomic"[Source] AND cluster_public[prop] AND "Saccharomyces cerevisiae"[Organism]\' '
)

MY_SRA_CMD <- SRA_DL_CMDS[grep(gsub("'","",species_name),SRA_DL_CMDS)]
MY_SRA_CMD

system(MY_SRA_CMD)

if (!file.exists(CSV)) {
    stop("Error: the metadata CSV file does not exist")
}

res <-read.csv(CSV,stringsAsFactors=FALSE)

# remove tabs
res <- apply(res,2,function(x)gsub('\t', ' ',x))

res <-as.data.frame(res,stringsAsFactors=FALSE)

res <-res[order(res$Run),]

accessions <- as.data.frame(cbind(res$Experiment,res$SRAStudy,res$Sample,res$Run))
nrow(accessions)
colnames(accessions) <- c("experiment","study","sample","run")

runs <- res$Run

source("dee_pipeline_functions.R")
########################
# Now determine which datasets have already been processed and completed
########################

unlink(list.files(DATAWD,pattern="_STARtmp",recursive=T))

folders<-list.files(DATAWD,full.names=T,pattern="RR")
fin_new<-paste(DATAWD,sapply(strsplit(list.files(path = DATAWD, pattern = "finished" , 
  full.names = T, recursive = T),"/"),"[[",7),sep="/")
fin_new<-fin_new[grep("RR",fin_new)]
val_old<-paste(DATAWD,sapply(strsplit(list.files(path = DATAWD, pattern = "validated" , 
  full.names = T, recursive = T),"/"),"[[",7),sep="/")

# delete folders that are not expected as they are not in known runs
expected_folders<-paste(DATAWD,runs,sep="/")
unexpected_folders<-setdiff(folders,expected_folders)
write(unexpected_folders,file=paste(DATAWD,"/unexpected_folders.txt",sep=""))

# delete folders without finished or validated files
allocated<-union(fin_new,val_old)
unalloc<-setdiff(folders,allocated)
unlink(unalloc,recursive=TRUE)

#Expected rows of se and ke files
gre<-dim( read.table(paste(DATAWD,"/rownames_gene.txt",sep=""),header=TRUE, comment.char="") )[1]
tre<-dim( read.table(paste(DATAWD,"/rownames_tx.txt",sep=""),header=TRUE, comment.char="") )[1]

message("run the check of new datasets")
source("dee_pipeline_functions.R")

mclapply(fin_new,check_contents,gre,tre,mc.cores=5)
message("done")

val<-paste(DATAWD,sapply(strsplit(list.files(path = DATAWD, pattern = "validated" , 
  full.names = T, recursive = T),"/"),"[[",7),sep="/")
runs_done<-sapply(strsplit(val,"/"),"[[",7)

se_files<-paste( DATAWD , "/" , runs_done , "/" , runs_done , ".se.tsv.gz",sep="")

print(paste(length(runs_done),"runs completed"))
runs_todo<-base::setdiff(runs, runs_done)
print(paste(length(runs_todo),"requeued runs"))

#Update queue on webserver if older than
queue_name=paste(QUEUEWD,"/",org,".queue.txt",sep="")
DIFFTIME=difftime ( ( Sys.Date()-90 ) , file.mtime(queue_name,units="s") )[1]

if ( DIFFTIME > 0 ) {
  write.table(runs_todo,queue_name,quote=F,row.names=F,col.names=F)
  SCP_COMMAND=paste("scp -i ~/.ssh/dee2 ",queue_name ," ubuntu@118.138.235.221:~/Public")
  system(SCP_COMMAND)
}

#Update metadata on webserver
accessions_done<-accessions[which(accessions$run %in% runs_done),]
write.table(accessions_done,file=paste(SRADBWD,"/",org,"_accessions.tsv",sep=""),quote=F,row.names=F)

save.image(file = paste(org,".RData",sep=""))

#collect QC info - this is temporary and logic will be incorporated in future
QC_summary="BLANK"

# here extract any samples that have GEO sample IDs ie GSM
# res$GEO_Accession

# here we use reutils use https://www.rdocumentation.org/packages/reutils/versions/0.2.2
GEO_QUERY_TERMS <- c(
'"Arabidopsis thaliana"[porgn:__txid3702] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Caenorhabditis elegans"[porgn:__txid6239] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Drosophila melanogaster"[porgn:__txid7227] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Danio rerio"[porgn:__txid7955] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Escherichia coli"[porgn:__txid562] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Homo sapiens"[porgn:__txid9606] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Mus musculus"[porgn:__txid10090] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Rattus norvegicus"[porgn:__txid10116] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Saccharomyces cerevisiae"[porgn:__txid4932] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]'
)

MY_GEO_QUERY_TERM <- GEO_QUERY_TERMS[grep(gsub("'","",species_name),GEO_QUERY_TERMS)]


# so now we are going with JSON format because XML was failing due to 
# angled brackets in the GEO data
ESEARCH <- esearch(term = MY_GEO_QUERY_TERM , db = "gds", rettype = "uilist", 
  retmode = "json", retstart = 0, retmax = 500000000, usehistory = TRUE, 
  webenv = NULL, querykey = NULL, sort = NULL, field = NULL,
  datetype = NULL, reldate = NULL, mindate = NULL, maxdate = NULL)

j <- fromJSON(ESEARCH$content)
COUNT <- j$esearchresult$count
myrange <- seq(0,COUNT,500)

gsel <- lapply(myrange, function(i) {
  ESUMMARY <- esummary(ESEARCH,retstart=i,retmax=500,retmode="json")
  myjson <- fromJSON(ESUMMARY$content)
  myjson <- myjson$result
  myjson[1] = NULL
  geodf <- do.call(rbind, myjson)
  geodf <- geodf[,c("gse","accession")]
  return(geodf)
})
gsel <- do.call(rbind, gsel)
GSE <- unlist(gsel[,1])
GSM <- unlist(gsel[,2])
gse <- data.frame(GSE,GSM)
colnames(gse) <- c("GEO_series","GEO_sample")
gse$GEO_series[which(gse$GEO_series=="")] <- "NA"
gse$GEO_series <- paste("GSE",sapply(strsplit(gse$GEO_series,";"),"[[",1),sep="")

resx <- merge(res,gse,by.x="SampleName",by.y="GEO_sample",all.x=TRUE)

# here is a good opportunity to check that the join has worked
res<-resx
res<-res[order(res$Run),]

#extract out the important accessions in order
x2<-as.data.frame(cbind(res$Run, QC_summary, res$Experiment, res$Sample, 
  res$SRAStudy,res$SampleName, res$GEO_series, res$LibraryName), stringsAsFactors=FALSE)

colnames(x2)<-c("SRR_accession","QC_summary","SRX_accession","SRS_accession",
  "SRP_accession","Sample_name","GEO_series","Library_name")

# NA values replaced with blank
x2[is.na(x2)] <- ""

#write out the accession number info and upload to webserver
write.table(x2,file=paste(SRADBWD,"/",org,"_metadata.complete.tsv.cut",sep=""),
  quote=F,sep="\t",row.names=F)
x2<-x2[which(x2$SRR_accession %in% runs_done),]

# write the srpqueue to enable new request from users
srpqueue <- accessions[which(! accessions$run %in% x2$SRR_accession),2]
srpqueue <- unique(srpqueue)
srpqueue <- srpqueue[order(srpqueue)]
srpqueuename = paste(SRADBWD,"/",org,"_srpqueue.txt",sep="")
writeLines(srpqueue,con=srpqueuename)
SCP_COMMAND=paste("scp -i ~/.ssh/dee2", srpqueuename ,
    " ubuntu@118.138.235.221:/dee2_data/srpqueue")
system(SCP_COMMAND)


source("dee_pipeline_functions.R")
library("parallel")
# this line is giving some errors with mouse data - now executes in batches
myaccessions <- x2$SRR_accession
qc_res=NULL
chunksize=1000
while (length(myaccessions) > chunksize ) {
    qc_res_part<-unlist(mclapply(myaccessions[1:chunksize], qc_analysis, org=org,
        mc.cores=4))
    qc_res <- c(qc_res,qc_res_part)
    myaccessions <- myaccessions[(chunksize+1):length(myaccessions)]
}

if (length(myaccessions) > 0 ) {
    qc_res_part<-unlist(mclapply(myaccessions, qc_analysis, org=org,mc.cores=4))
    qc_res <- c(qc_res,qc_res_part)
}
x2$QC_summary <- qc_res

# write metadata
write.table(x2,file=paste(SRADBWD,"/",org,"_metadata.tsv.cut",sep=""),
  quote=F,sep="\t",row.names=F)

#aggregate se ke and qc data
CMD=paste("./dee_pipeline.sh",org)
system(CMD)


#delete *e.tsv.gz after each chunk of 10000
#not needed if --exclude \"*.gz\" is used
CMD2=paste('ssh -i ~/.ssh/dee2 ubuntu@118.138.235.221 "find /dee2_data/data/',org,' | grep e.tsv.gz | parallel -j1 rm {}"',sep="")

#here we rsync files to server in chunks of 10000
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
    CMD=paste('rsync -azh -e \"ssh -i  ~/.ssh/dee2 \" --exclude \"*.gz\" ',
      chunk ,' ubuntu@118.138.235.221:/dee2_data/data/',org,sep="")
    system(CMD)
    system(CMD2)
  }
}
d<-val
rsync(d,org)

#upload metadata
SCP_COMMAND=paste("scp -i ~/.ssh/dee2", 
  paste(SRADBWD,"/",org,"_metadata.tsv.cut",sep="") ,
    " ubuntu@118.138.235.221:/dee2_data/metadata")
system(SCP_COMMAND)

save.image(file = paste(org,".RData",sep=""))

#now attach the additional metadata and upload
x<-res[, !(colnames(res) %in% 
  c("QC_summary","Experiment","sample_acc","SRA.Study","GEO_series","GEO_Accession"))]
x<-merge(x2,x,by.x="SRR_accession",by.y="Run")

x<-x[which(x$SRR_accession %in% runs_done),]
x <- apply(x,2,as.character)
x<-gsub("\r?\n|\r", " ", x)
write.table(x,file=paste(SRADBWD,"/",org,"_metadata.tsv",sep=""),quote=F,sep="\t",row.names=F)
SCP_COMMAND=paste("scp -i ~/.ssh/dee2 ", 
  paste(SRADBWD,"/",org,"_metadata.tsv",sep="") ,
  " ubuntu@118.138.235.221:/dee2_data/metadata")
system(SCP_COMMAND)

save.image(file = paste(org,".RData",sep=""))


}

png("dee_datasets.png",width=600,height=600)
options(bitmapType="cairo")
FILES1<-list.files(pattern="*queue.txt$",path="/mnt/md0/dee2/queue/",full.names=T)
x<-as.data.frame(sapply(FILES1,rowcnt2),stringsAsFactors=FALSE)
rownames(x)=c("A. thaliana","C. elegans","D. melanogaster","D. rerio",
  "E. coli","H. sapiens","M. musculus","R. norvegicus","S. cerevisiae")
colnames(x)="queued"

FILES2<-list.files(pattern="*accessions.tsv$",path="/mnt/md0/dee2/sradb/",full.names=T)
y<-as.data.frame(sapply(FILES2,rowcnt2),stringsAsFactors=FALSE)
rownames(y)=c("A. thaliana","C. elegans","D. melanogaster","D. rerio",
  "E. coli","H. sapiens","M. musculus","R. norvegicus","S. cerevisiae")
colnames(y)="completed"

z<-merge(x,y,by=0)
rownames(z)=z$Row.names
z$Row.names=NULL

DATE=date()
HEADER=paste("Updated",DATE)
z<-z[order(rownames(z),decreasing=T ), ,drop=F]
par(las=2) ; par(mai=c(1,2.5,1,0.5))
MAX=1000000

bb<-barplot( rbind( as.numeric(z$queued) , as.numeric(z$completed) ) ,
  names.arg=rownames(z) ,xlim=c(0,MAX), beside=T, main=HEADER, col=c("darkblue","red") ,
  horiz=T , las=1, cex.axis=1.2, cex.names=1.4, cex.main=1.4 ,
  xlab="number of SRA runs")

legend("topright", colnames(z), fill=c("darkblue","red") , cex=1.2)

text( cbind(as.numeric(z[,1])+70000 ,as.numeric(z[,2])+70000 )  ,t(bb),labels=c(z[,1],z[,2]) ,cex=1.2)
dev.off()
system("scp -i ~/.ssh/dee2 dee_datasets.png ubuntu@118.138.235.221:/home/ubuntu/dee2/frontend/html/images")
