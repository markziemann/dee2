for ( ORG in c("ath","cel","dme","dre","eco","hsa","mmu","rno","sce")) {

PDFNAME=paste("sim_",ORG,".pdf",sep="")
pdf(PDFNAME,width=12.6,height=16.8)
par(mfrow=c(4,3))

##50bp SE
#star
TSV=paste(ORG,"_art_se_50_star.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 50 bp SE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd", "Rho=", RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")

#Kallisto transcripts
TSV=paste(ORG,"_art_se_50_kalt.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 50 bp SE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd", "Rho=", RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")

#Kallisto genes
TSV=paste(ORG,"_art_se_50_kalg.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 50 bp SE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd", "Rho=",RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")

#100bp SE
#star
TSV=paste(ORG,"_art_se_100_star.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 100 bp SE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd","Rho=",RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")

#Kallisto transcripts
TSV=paste(ORG,"_art_se_100_kalt.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 100 bp SE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd","Rho=",RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")

#Kallisto genes
TSV=paste(ORG,"_art_se_100_kalg.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 100 bp SE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd","Rho=",RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")


#50bp PE
#star
TSV=paste(ORG,"_art_pe_50_star.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 50 bp PE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd","Rho=",RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")

#Kallisto transcripts
TSV=paste(ORG,"_art_pe_50_kalt.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 50 bp PE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd","Rho=",RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")

#Kallisto genes
TSV=paste(ORG,"_art_pe_50_kalg.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 50 bp PE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd","Rho=",RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")


#100bp PE
#star
TSV=paste(ORG,"_art_pe_100_star.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 100 bp PE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd","Rho=",RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")

#Kallisto transcripts
TSV=paste(ORG,"_art_pe_100_kalt.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 100 bp PE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd","Rho=",RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")

#Kallisto genes
TSV=paste(ORG,"_art_pe_100_kalg.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 100 bp PE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd","Rho=",RHO)
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")

dev.off()

}
