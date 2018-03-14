library(ggplot2)
rowcnt<-function(file){nrow(read.table(file,fill=T)) }
rowcnt2<-function( file) { z<-system(paste("wc -l < ",file) , intern=TRUE) ; z}


png("/home/mziemann/dee_datasets.png",width=480,height=480)
FILES<-list.files(pattern="*a.tsv$",path="/scratch/mziemann/dee2/sradb/",full.names=T)
z<-as.data.frame(sapply(FILES,rowcnt2),stringsAsFactors=FALSE)
rownames(z)=c("A. thaliana","C. elegans","D. melanogaster","D. rerio","E. coli","H. sapiens","M. musculus","R. norvegicus","S. cerevisiae")
DATE=strsplit(as.character(file.info(FILES[1])[,4])," ",fixed=T)[[1]][1]
HEADER=paste("Hosted datasets",DATE)
colnames(z)=HEADER
z<-z[order(rownames(z),decreasing=T ), ,drop=F]
par(las=2) ; par(mai=c(1,2.5,1,0.5))
MAX=max(as.numeric(z[,1]))+20000
options(scipen=10000)
bb<-barplot(as.numeric(z[,1]),names.arg=rownames(z),xlim=c(0,MAX),beside=T,main=HEADER,horiz=T , las=1, cex.axis=1.5, cex.names=1.5,cex.main=1.4)
#barplot(z,beside=T,names.arg=rownames(z),ylim=c(0,250000),main=HEADER)
text(as.numeric(z[,1])+10000,bb,labels=z[,1],cex=1.4)	
dev.off()



x<-read.table("dee_progress.txt")
colnames(x)=c("species","date","num_datasets")
attach(x)


png("/home/mziemann/dee_progress.png",width=650,height=480)
ggplot(x, aes(x = date, y = num_datasets, group = species, colour = species)) + geom_line() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

library(parallel)
l<-list.files(path = "/scratch/mziemann/dee2/queue", pattern="*txt",full.names =T)
y<-t(as.data.frame(mclapply(l,rowcnt, mc.cores = detectCores()  )))
ll<-list.files(path = "/scratch/mziemann/dee2/queue", pattern="*txt",full.names =F)
rownames(y)=sub(".queue.txt","",as.character(ll))
png(filename = "/home/mziemann/dee_queue.png",width=600,height=480)
par(las=2) ; par(mai=c(2,1,1,1))
bb<-barplot(y,beside=T,names.arg=rownames(y),ylim=c(0,250000),main="Number of datasets in queue")
text(bb,as.numeric(y)+5000,labels=as.numeric(y),cex=0.9)
dev.off()

