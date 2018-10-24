pdf("sim_hsapiens.pdf",width=12.6,height=16.8)
par(mfrow=c(4,3))

##50bp SE
#star
x<-read.table("hsa_art_se_50_star.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 50 bp SE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,12,LABEL)

#Kallisto transcripts
x<-read.table("hsa_art_se_50_kalt.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 50 bp SE reads"
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,9,LABEL)

#Kallisto genes
x<-read.table("hsa_art_se_50_kalg.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 50 bp SE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,9,LABEL)

#100bp SE
#star
x<-read.table("hsa_art_se_100_star.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 100 bp SE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,12,LABEL)

#Kallisto transcripts
x<-read.table("hsa_art_se_100_kalt.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 100 bp SE reads"
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,9,LABEL)

#Kallisto genes
x<-read.table("hsa_art_se_100_kalg.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 100 bp SE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,9,LABEL)


#50bp PE
#star
x<-read.table("hsa_art_pe_50_star.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 50 bp PE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,8,LABEL)

#Kallisto transcripts
x<-read.table("hsa_art_pe_50_kalt.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 50 bp PE reads"
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,5,LABEL)

#Kallisto genes
x<-read.table("hsa_art_pe_50_kalg.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 50 bp PE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,9,LABEL)


#100bp PE
#star
x<-read.table("hsa_art_pe_100_star.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 100 bp PE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,8,LABEL)

#Kallisto transcripts
x<-read.table("hsa_art_pe_100_kalt.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 100 bp PE reads"
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,5,LABEL)

#Kallisto genes
x<-read.table("hsa_art_pe_100_kalg.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 100 bp PE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(0,9,LABEL)

dev.off()

pdf("sim_athaliana.pdf",width=12.6,height=16.8)
par(mfrow=c(4,3))

##50bp SE
#star
x<-read.table("ath_art_se_50_star.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 50 bp SE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,8,LABEL)

#Kallisto transcripts
x<-read.table("ath_art_se_50_kalt.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 50 bp SE reads"
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,5,LABEL)

#Kallisto genes
x<-read.table("ath_art_se_50_kalg.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 50 bp SE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,9,LABEL)

#100bp SE
#star
x<-read.table("ath_art_se_100_star.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 100 bp SE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,8,LABEL)

#Kallisto transcripts
x<-read.table("ath_art_se_100_kalt.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 100 bp SE reads"
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,5,LABEL)

#Kallisto genes
x<-read.table("ath_art_se_100_kalg.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 100 bp SE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,9,LABEL)


#50bp PE
#star
x<-read.table("ath_art_pe_50_star.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 50 bp PE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,8,LABEL)

#Kallisto transcripts
x<-read.table("ath_art_pe_50_kalt.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 50 bp PE reads"
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,5,LABEL)

#Kallisto genes
x<-read.table("ath_art_pe_50_kalg.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 50 bp PE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,9,LABEL)


#100bp PE
#star
x<-read.table("ath_art_pe_100_star.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","STAR")
x<-x/colSums(x)*1000000
th<-x[which(x$STAR>(x$Simulated*2)),]
tl<-x[which(x$STAR<(x$Simulated/2)),]
HEADER="STAR gene RPM - 100 bp PE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,8,LABEL)

#Kallisto transcripts
x<-read.table("ath_art_pe_100_kalt.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto transcript RPM - 100 bp PE reads"
SUBHEADER=paste( dim(x)[1] ,"transcripts", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,5,LABEL)

#Kallisto genes
x<-read.table("ath_art_pe_100_kalg.tsv",row.names=1,header=F)
colnames(x)=c("Simulated","Kallisto")
x<-x/colSums(x)*1000000
th<-x[which(x$Kallisto>(x$Simulated*2)),]
tl<-x[which(x$Kallisto<(x$Simulated/2)),]
HEADER="Kallisto gene RPM - 100 bp PE reads"
SUBHEADER=paste(dim(x)[1] ,"genes", dim(th)[1], "overest'd", dim(tl)[1], "underest'd")
plot(log2(x),main=HEADER,pch=19,cex=0.5)
mtext(SUBHEADER,cex=0.8)
points(log2(th),main="Ensembl gene RPM",pch=19,cex=0.5,col="red")
points(log2(tl),main="Ensembl gene RPM",pch=19,cex=0.5,col="blue")
LABEL=paste("Rho=",round(cor(x,method="s")[1,2] , digits=3))
text(1,9,LABEL)

dev.off()
