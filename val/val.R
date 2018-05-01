library("plyr")
library("statmod")
library("edgeR")
library("locfit")
library("reshape2")
library("parallel")

source("https://raw.githubusercontent.com/markziemann/dee2/master/frontend/html/getDEE2.R")

#############
# A. thaliana GSE53078 ctrl=c(“SRR1044945”,”SRR1044946”), trt=c(“SRR1044947”,”SRR1044948”)
#############
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
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,y$counts,by='row.names')
geo_res<-dge[order(dge$PValue),]

dee_geo_res<-merge(dee_res,geo_res,by="Row.names")

#contrast wise correlation
dee_geo_res$dee_metric=dee_geo_res$logFC.x/-log10(dee_geo_res$PValue.x)
dee_geo_res$geo_metric=dee_geo_res$logFC.y/-log10(dee_geo_res$PValue.y)
ath_cor=cor(dee_geo_res$dee_metric,dee_geo_res$geo_metric,method="s")

#sample wise correlation
colnames(x$GeneCounts)=gsub("$","_dee",colnames(x$GeneCounts))
colnames(x$GeneCounts)=gsub("$","_geo",colnames(b))
d<-merge(x$GeneCounts,b,by=0)
SRR1044945=cor(d[,2:9],method="s")[1,5]
SRR1044946=cor(d[,2:9],method="s")[2,6]
SRR1044947=cor(d[,2:9],method="s")[3,7]
SRR1044948=cor(d[,2:9],method="s")[4,8]

#save result
ath_res=c(ath_cor,SRR1044945,SRR1044946,SRR1044947,SRR1044948)

#########
#C. elegans GSE46344  ctrl=c(“SRR834594”,”SRR834595”,”SRR834596”), trt=c(“SRR834600”,”SRR834601”,”SRR834602”)
#########
x<-getDEE2("celegans",c("SRR834594","SRR834595","SRR834596","SRR834600","SRR834601","SRR834602"))
x$GeneCounts<-x$GeneCounts[which(rowMeans(x$GeneCounts)>10),]
group<-c(1,1,1,2,2,2)
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

system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46344/suppl/GSE46344_ce_11samples_DOG_counts_and_RPKM.txt.gz | gunzip | cut -f1,4-6,10-12 > GSE46344.tsv")
GSE46344<-read.table("GSE46344.tsv",header=T,row.names=1)
colnames(GSE46344)=c("SRR834594","SRR834595","SRR834596","SRR834600","SRR834601","SRR834602")

system("curl ftp://ftp.ensembl.org/pub/release-90/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.90.gtf.gz | gunzip | grep -w gene | cut -d '\"' -f2,4 | tr '\"' '\t' > celegans_genenames.tsv")
celegans_genenames<-read.table("celegans_genenames.tsv")

b<-merge(celegans_genenames,GSE46344,by.x="V2",by.y=0)
rownames(b)=b$V1
b$V1=NULL
b$V2=NULL
b<-b[which(rowMeans(b)>10),]
group<-c(1,1,1,2,2,2)
y <- DGEList(counts=b, group=group)
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
geo_res<-dge[order(dge$PValue),]

dee_geo_res<-merge(dee_res,geo_res,by="Row.names")

#contrast wise correlation
dee_geo_res$dee_metric=dee_geo_res$logFC.x/-log10(dee_geo_res$PValue.x)
dee_geo_res$geo_metric=dee_geo_res$logFC.y/-log10(dee_geo_res$PValue.y)
cel_cor=cor(dee_geo_res$dee_metric,dee_geo_res$geo_metric,method="s")

#sample wise correlation
colnames(x$GeneCounts)=gsub("$","_dee",colnames(x$GeneCounts))
colnames(b)=gsub("$","_geo",colnames(b))
d<-merge(x$GeneCounts,b,by=0)
SRR834594=cor(d[,2:12])[1,7]
SRR834595=cor(d[,2:12])[2,8]
SRR834596=cor(d[,2:12])[3,9]
SRR834600=cor(d[,2:12])[4,10]
SRR834601=cor(d[,2:12])[5,11]
SRR834602=cor(d[,2:12])[6,12]

cel_res=c(cel_cor,SRR834594,SSRR834595,SRR834596,SRR834600,SRR834601,SRR834602)

