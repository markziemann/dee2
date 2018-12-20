library(reshape2)

for (QCZ in list.files(pattern="_qc.tsv.bz2") ) {

ORG=sapply(strsplit(QCZ,"_"),"[[",1)
QC=gsub(".bz2","",QCZ)
CMD=paste("pbzip2 -dc",QCZ,"| awk 'NF==3' > ", QC)
system(CMD)
tmp<-read.table(QC)
x<-as.matrix(acast(tmp, V1~V2, value.var="V3"))
#ROWNAMES=rownames(x)
#heax<-as.data.frame(x)
print(ORG)
PDF=paste(ORG,"_qc.pdf",sep="")
pdf(PDF)
par(mfrow=c(3,3))
par(mar=c(5,6,3,3))

# Seq format
barplot(summary(x[,"SequenceFormat"])[1:2],main="Sequence format",ylab="no. SRA runs",cex.main=0.8)
text(  ((1:2)-0.5)*1.3, ( summary(x[,"SequenceFormat"])[1:2] *1.0 ) +400, labels=summary(x[,"SequenceFormat"])[1:2], xpd=TRUE)

# Quality encoding
z<-summary(x[,"QualityEncoding"])[1:2]
names(z)=c("1.9","1.5")
barplot(z,main="Quality encoding",ylab="no. SRA runs",cex.main=0.8)
text(  ((1:2)-0.5)*1.3, ( summary(x[,"QualityEncoding"])[1:2] *1.0 )+400 , labels=summary(x[,"QualityEncoding"])[1:2], xpd=TRUE,cex.lab=0.5)

# Read length
hist(as.numeric(as.character(x[,"Read1MedianLength"])),main="Read length",xlab="bp",cex.main=0.8)

# Read QC pass rate
hist(as.numeric(gsub("%","",as.character(x[,"QcPassRate"]))),main="Read QC Pass Rate",xlab="Percent",cex.main=0.8)

# Number of QC passed reads
hist(log10(as.numeric(as.character(x[,"NumReadsQcPass"]))),main="Number of QC passed reads per run",xlab="log10(reads)",cex.main=0.8)

#Star map rate
z<-as.numeric(as.character(x[,"STAR_UniqMappedReads"])) / as.numeric(as.character(x[,"NumReadsQcPass"])) *100
z<-z[which(z<100)]
hist(z ,main="STAR unique mapping percent",xlab="Percent",cex.main=0.8)

#Star assignment rate
z<-as.numeric(as.character(x[,"STAR_AssignedReads"])) / as.numeric(as.character(x[,"NumReadsQcPass"])) *100
z<-z[which(z<100)]
hist(z ,main="STAR reads assigned percent",xlab="Percent",cex.main=0.8)

#STAR_Strandedness
z<-summary(x[,"STAR_Strandedness"])[1:3]
print(names(z))
names(z)=c("pos","none","neg")
barplot(z ,main="STAR Strandedness",ylab="no. SRA runs",cex.main=0.8)
text(  ((1:3)-0.5)*1.3, ( z *1.0 ) +400, labels=z, xpd=TRUE)

#Kallisto map rate
z<-as.numeric(as.character(x[,"Kallisto_MappedReads"])) / as.numeric(as.character(x[,"NumReadsQcPass"])) *100
z<-z[which(z<100)]
hist(z ,main="Kallisto mapped reads percent",xlab="Percent",cex.main=0.8)

dev.off()

}
