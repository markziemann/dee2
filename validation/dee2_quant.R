for ( ORG in c("ath","cel","dme","dre","eco","hsa","mmu","rno","sce")) {

PDFNAME=paste("sim_",ORG,".pdf",sep="")
pdf(PDFNAME,width=12.6,height=16.8)
par(mfrow=c(4,3))
par(cex.lab=1.4)
par(mar=c(5.1,5.1,5.1,2.1))

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
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")


#Kallisto transcripts
TSV=paste(ORG,"_art_se_50_kalt.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 50 bp SE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")


#Kallisto genes
TSV=paste(ORG,"_art_se_50_kalg.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 50 bp SE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

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
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

#Kallisto transcripts
TSV=paste(ORG,"_art_se_100_kalt.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 100 bp SE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

#Kallisto genes
TSV=paste(ORG,"_art_se_100_kalg.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 100 bp SE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2) ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")


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
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

#Kallisto transcripts
TSV=paste(ORG,"_art_pe_50_kalt.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 50 bp PE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

#Kallisto genes
TSV=paste(ORG,"_art_pe_50_kalg.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 50 bp PE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")


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
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

#Kallisto transcripts
TSV=paste(ORG,"_art_pe_100_kalt.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 100 bp PE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

#Kallisto genes
TSV=paste(ORG,"_art_pe_100_kalg.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 100 bp PE reads"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

dev.off()

}



#now for the figure
ORG="hsa"
PDFNAME=paste("fig_sim_",ORG,".pdf",sep="")
pdf(PDFNAME,width=12.6,height=16.8)
par(mfrow=c(4,3))
par(cex.lab=1.4)
par(mar=c(5.1,5.1,5.1,2.1))

#100bp SE
#star
TSV=paste(ORG,"_art_se_100_star.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

#Kallisto transcripts
TSV=paste(ORG,"_art_se_100_kalt.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2) ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

#Kallisto genes
TSV=paste(ORG,"_art_se_100_kalg.tsv",sep="")
x<-read.table(TSV,row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM"
RHO=round(cor(x,method="s")[1,2] , digits=3)
SUBHEADER=paste(dim(x)[1] ,"contigs", signif(dim(th)[1]/dim(x)[1]*100,digits=2)  ,"%"  , "overest'd", signif(dim(tl)[1]/dim(x)[1]*100,digits=2), "%", "underest'd", "Rho=", RHO)
plot(log2(x),col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),main=HEADER,pch=19,cex=0.5,xlab="log2(ground truth)",ylab="log2(DEE2-pipeline-inferred expression)",cex.axis=1.4,cex.main=1.6)
mtext(SUBHEADER,cex=0.8)
abline(1,1,lty = 5,lwd=3,col="darkgray")
abline(-1,1,lty = 5,lwd=3,col="darkgray")

dev.off()
