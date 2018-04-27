library("plyr")
library("statmod")
library("edgeR")
library("locfit")
library("reshape2")
library("parallel")

source("https://raw.githubusercontent.com/markziemann/dee2/master/frontend/html/getDEE2.R")

#ATHALIANA
x<-getDEE2("athaliana",c("SRR1044945","SRR1044946","SRR1044947","SRR1044948")) 
x$GeneCounts<-x$GeneCounts[which(rowMeans(x$GeneCounts)>10),]
group<-c(1,1,2,2)
y <- DGEList(counts=x$GeneCounts, group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,robust=TRUE,prior.df=1)
fit <- glmFit(y)
lrt <- glmLRT(fit)
dge<-as.data.frame(topTags(lrt,n=1000000))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,y$counts,by='row.names')
dee_res<-dge[order(dge$PValue),]

system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1281nnn/GSM1281703/suppl/GSM1281703_Col_1.txt.gz | gunzip > SRR1044945.tsv")
SRR1044945<-read.table("SRR1044945.tsv",header=T)
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1281nnn/GSM1281704/suppl/GSM1281704_Col_2.txt.gz | gunzip > SRR1044946.tsv")
SRR1044946<-read.table("SRR1044946.tsv",header=T)
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1281nnn/GSM1281705/suppl/GSM1281705_HBI1-Ox_1.txt.gz | gunzip > SRR1044947.tsv")
SRR1044947<-read.table("SRR1044947.tsv",header=T)
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1281nnn/GSM1281706/suppl/GSM1281706_HBI1-Ox_2.txt.gz | gunzip > SRR1044948.tsv")
SRR1044948<-read.table("SRR1044948.tsv",header=T)

b<-cbind(SRR1044945[,2],SRR1044946[,2],SRR1044947[,2],SRR1044948[,2])
rownames(b)=gsub("-TAIR-G","",SRR1044945[,1])
colnames(b)=c("SRR1044945","SRR1044946","SRR1044947","SRR1044948")
b<-b[which(rowMeans(b)>10),]
group<-c(1,1,2,2)
y <- DGEList(counts=b, group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,robust=TRUE,prior.df=1)
fit <- glmFit(y)
lrt <- glmLRT(fit)
dge<-as.data.frame(topTags(lrt,n=1000000))
dge2$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,y$counts,by='row.names')
geo_res<-dge[order(dge$PValue),]

colnames(x$GeneCounts)=gsub("$","_dee",colnames(x$GeneCounts))
head(merge(x$GeneCounts,b,by=0))
 cor(d[,2:9],method="s")[1,5]
 cor(d[,2:9],method="s")[2,6]
 cor(d[,2:9],method="s")[3,7]
 cor(d[,2:9],method="s")[4,8]

