library("plyr")
library("statmod")
library("edgeR")
library("locfit")
library("reshape2")
library("parallel")

source("https://raw.githubusercontent.com/markziemann/dee2/master/getDEE2.R")


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
dee_geo_res$dee_metric=sign(dee_geo_res$logFC.x)/-log10(dee_geo_res$PValue.x)
dee_geo_res$geo_metric=sign(dee_geo_res$logFC.y)/-log10(dee_geo_res$PValue.y)
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
dee_geo_res$dee_metric=sign(dee_geo_res$logFC.x)/-log10(dee_geo_res$PValue.x)
dee_geo_res$geo_metric=sign(dee_geo_res$logFC.y)/-log10(dee_geo_res$PValue.y)
cel_cor=cor(dee_geo_res$dee_metric,dee_geo_res$geo_metric,method="s")

#sample wise correlation
colnames(x$GeneCounts)=gsub("$","_dee",colnames(x$GeneCounts))
colnames(b)=gsub("$","_geo",colnames(b))
d<-merge(x$GeneCounts,b,by=0)
SRR834594=cor(d[,2:13],method="s")[1,7]
SRR834595=cor(d[,2:13],method="s")[2,8]
SRR834596=cor(d[,2:13],method="s")[3,9]
SRR834600=cor(d[,2:13],method="s")[4,10]
SRR834601=cor(d[,2:13],method="s")[5,11]
SRR834602=cor(d[,2:13],method="s")[6,12]

cel_res=c(cel_cor,SRR834594,SRR834595,SRR834596,SRR834600,SRR834601,SRR834602)

################
#D. melanogaster GSE43180 ctrl=c(“SRR641382”,”SRR641383”), trt=(“SRR641384”,”SRR641385”)
################

x<-getDEE2("dmelanogaster",c("SRR641382","SRR641383","SRR641384","SRR641385"))
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

system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1057nnn/GSM1057982/suppl/GSM1057982_1906.counts.txt.gz | gunzip > SRR641382.tsv")
SRR641382<-read.table("SRR641382.tsv",header=F)
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1057nnn/GSM1057983/suppl/GSM1057983_1907.counts.txt.gz | gunzip > SRR641383.tsv")
SRR641383<-read.table("SRR641383.tsv",header=F)
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1057nnn/GSM1057984/suppl/GSM1057984_1904.counts.txt.gz | gunzip > SRR641384.tsv")
SRR641384<-read.table("SRR641384.tsv",header=F)
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1057nnn/GSM1057985/suppl/GSM1057985_1905.counts.txt.gz | gunzip > SRR641385.tsv")
SRR641385<-read.table("SRR641385.tsv",header=F)

b<-cbind(SRR641382[,2],SRR641383[,2],SRR641384[,2],SRR641384[,2])
rownames(b)=SRR641382[,1]
colnames(b)=c("SRR641382","SRR641383","SRR641384","SRR641385")
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
dee_geo_res$dee_metric=sign(dee_geo_res$logFC.x)/-log10(dee_geo_res$PValue.x)
dee_geo_res$geo_metric=sign(dee_geo_res$logFC.y)/-log10(dee_geo_res$PValue.y)
dme_cor=cor(dee_geo_res$dee_metric,dee_geo_res$geo_metric,method="s")

#sample wise correlation
colnames(x$GeneCounts)=gsub("$","_dee",colnames(x$GeneCounts))
colnames(b)=gsub("$","_geo",colnames(b))
d<-merge(x$GeneCounts,b,by=0)
SRR641382=cor(d[,2:9],method="s")[1,5]
SRR641383=cor(d[,2:9],method="s")[2,6]
SRR641384=cor(d[,2:9],method="s")[3,7]
SRR641385=cor(d[,2:9],method="s")[4,8]

dme_res=c(dme_cor,SRR641382,SRR641383,SRR641384,SRR641385)

##########
# D. rerio GSE59683 ctrl=c(“SRR1523211”,”SRR1523212”,”SRR1523213”), trt=c(“SRR1523214”,”SRR1523215”,”SRR1523216”)
##########
#HERENOW
x<-getDEE2("drerio",c("SRR1523211","SRR1523212","SRR1523213","SRR1523214","SRR1523215","SRR1523216"))
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

system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1442nnn/GSM1442819/suppl/GSM1442819_1.txt.gz | gunzip | cut -f3- > SRR1523211.tsv")
SRR1523211<-read.table("SRR1523211.tsv",header=T)
dre_genes<-SRR1523211[-nrow(SRR1523211),1]
SRR1523211<-as.numeric(SRR1523211[-nrow(SRR1523211),2])
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1442nnn/GSM1442820/suppl/GSM1442820_2.txt.gz | gunzip | cut -f3- > SRR1523212.tsv")
SRR1523212<-read.table("SRR1523212.tsv",header=T)
SRR1523212<-as.numeric(SRR1523212[-nrow(SRR1523212),2])
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1442nnn/GSM1442821/suppl/GSM1442821_3.txt.gz | gunzip | cut -f3- > SRR1523213.tsv")
SRR1523213<-read.table("SRR1523213.tsv",header=T)
SRR1523213<-as.numeric(SRR1523213[-nrow(SRR1523213),2])
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1442nnn/GSM1442822/suppl/GSM1442822_4.txt.gz | gunzip | cut -f3- > SRR1523214.tsv")
SRR1523214<-read.table("SRR1523214.tsv",header=T)
SRR1523214<-as.numeric(SRR1523214[-nrow(SRR1523214),2])
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1442nnn/GSM1442823/suppl/GSM1442823_5.txt.gz | gunzip | cut -f3- > SRR1523215.tsv")
SRR1523215<-read.table("SRR1523215.tsv",header=T)
SRR1523215<-as.numeric(SRR1523215[-nrow(SRR1523215),2])
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1442nnn/GSM1442824/suppl/GSM1442824_6.txt.gz | gunzip | cut -f3- > SRR1523216.tsv")
SRR1523216<-read.table("SRR1523216.tsv",header=T)
SRR1523216<-as.numeric(SRR1523216[-nrow(SRR1523216),2])

b<-cbind(SRR1523211,SRR1523212,SRR1523213,SRR1523214,SRR1523215,SRR1523216)
rownames(b)=dre_genes
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
dee_geo_res$dee_metric=sign(dee_geo_res$logFC.x)/-log10(dee_geo_res$PValue.x)
dee_geo_res$geo_metric=sign(dee_geo_res$logFC.y)/-log10(dee_geo_res$PValue.y)
dre_cor=cor(dee_geo_res$dee_metric,dee_geo_res$geo_metric,method="s")

#sample wise correlation
colnames(x$GeneCounts)=gsub("$","_dee",colnames(x$GeneCounts))
colnames(b)=gsub("$","_geo",colnames(b))
d<-merge(x$GeneCounts,b,by=0)
SRR1523211=cor(d[,2:13],method="s")[1,7]
SRR1523212=cor(d[,2:13],method="s")[2,8]
SRR1523213=cor(d[,2:13],method="s")[3,9]
SRR1523214=cor(d[,2:13],method="s")[4,10]
SRR1523215=cor(d[,2:13],method="s")[5,11]
SRR1523216=cor(d[,2:13],method="s")[6,12]

dre_res=c(dre_cor,SRR1523211,SRR1523212,SRR1523213,SRR1523214,SRR1523215,SRR1523216)

###########
# E. coli GSE80251
###########
# DEE2
x<-getDEE2("ecoli",c("SRR3379590","SRR3379591","SRR3379592","SRR3379593","SRR3379594","SRR3379595"))
#x<-getDEE2("ecoli",c("SRR933983","SRR933984","SRR933985","SRR933989","SRR933990","SRR933991"))

#Add gene names
system("wget -N ftp://ftp.ensemblgenomes.org/pub/release-36/bacteria//gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.36.gtf.gz  ")
system("gunzip -f Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.36.gtf.gz  ")
system("grep -w gene Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.36.gtf | cut -d '\"' -f2,4 | tr '\"' '\t' | sort > ecoli_genes.tsv")
ecoli_genes<-read.table("ecoli_genes.tsv",row.names=1)
x$GeneCounts<-merge(ecoli_genes,x$GeneCounts,by=0)
x$GeneCounts$Row.names=NULL

#Aggregate as there are a few redundant names
xx<-aggregate(x=x$GeneCounts[,2:ncol(x$GeneCounts)], by = list(unique.values = x$GeneCounts$V2 ) , FUN=sum)
rownames(xx)=xx$unique.values
xx$unique.values=NULL

#run edgeR
xx<-xx[which(rowMeans(xx)>10),]
group<-c(1,1,1,2,2,2)
y <- DGEList(counts=xx, group=group)
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
rownames(dee_res)=dee_res$Row.names
dee_res$Row.names=NULL

#GEO
system("curl 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE80251&format=file&file=GSE80251%5Fprocessed%5FRNA%5Fexpression%5Fmnfyap%2Exlsx' > GSE80251.xlsx ")
system("ssconvert -S --export-type Gnumeric_stf:stf_assistant -O 'separator=\"''\t''\"' GSE80251.xlsx GSE80251.xlsx.txt")
GSE80251<-read.table("GSE80251.xlsx.txt.1",header=T,stringsAsFactors = FALSE)
GSE80251<-GSE80251[,-c(1:3,5)]
GSE80251<-aggregate(x=GSE80251[,2:ncol(GSE80251)], by = list(unique.values = GSE80251$gene ) , FUN=sum)
rownames(GSE80251)=GSE80251$unique.values
GSE80251$unique.values=NULL
b<-GSE80251[,1:6]
#The authors dof the study report the basepair coverage, not the number of reads. Here I'm scaling the matrix down by the observed median readlength 
b<-round(b/108)
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
rownames(dge)=dge$Row.names
dge$Row.names=NULL
geo_res<-dge[order(dge$PValue),]

dee_geo_res<-merge(dee_res,geo_res,by=0)
rownames(dee_geo_res)=dee_geo_res$Row.names
dee_geo_res$Row.names=NULL

#contrast wise correlation
dee_geo_res$dee_metric=sign(dee_geo_res$logFC.x)/-log10(dee_geo_res$PValue.x)
dee_geo_res$geo_metric=sign(dee_geo_res$logFC.y)/-log10(dee_geo_res$PValue.y)
eco_cor=cor(dee_geo_res$dee_metric,dee_geo_res$geo_metric,method="s")

#sample wise correlation
colnames(xx)=gsub("$","_dee",colnames(xx))
colnames(b)=gsub("$","_geo",colnames(b))
d<-merge(xx,b,by=0)
SRR933983=cor(d[,2:13],method="s")[1,7]
SRR933984=cor(d[,2:13],method="s")[2,8]
SRR933985=cor(d[,2:13],method="s")[3,9]
SRR933989=cor(d[,2:13],method="s")[4,10]
SRR933990=cor(d[,2:13],method="s")[5,11]
SRR933991=cor(d[,2:13],method="s")[6,12]

eco_res=c(eco_cor,SRR933983,SRR933984,SRR933985,SRR933989,SRR933990,SRR933991)

###########
# H. sapiens
###########
# DEE2
x<-getDEE2("hsapiens",c("SRR1692137","SRR1692138","SRR1692139","SRR1692140","SRR1692141","SRR1692142"))

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

#GEO
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63776/suppl/GSE63776_PANC1_counts.txt.gz | gunzip | sed 's/\\r/\\n/g' | egrep -v '(2-Mar|1-Mar)'> GSE63776.tsv")
GSE63776<-read.table("GSE63776.tsv",header=T)
colnames(GSE63776)=c("SRR1692137","SRR1692138","SRR1692139","SRR1692140","SRR1692141","SRR1692142")

system("curl ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz | gunzip | grep -w gene | cut -d '\"' -f2,6 | tr '\"' '\t' > hsapiens_genenames.tsv")
hsapiens_genenames<-read.table("hsapiens_genenames.tsv")
b<-merge(hsapiens_genenames,GSE63776,by.x="V2",by.y=0)
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
dee_geo_res$dee_metric=sign(dee_geo_res$logFC.x)/-log10(dee_geo_res$PValue.x)
dee_geo_res$geo_metric=sign(dee_geo_res$logFC.y)/-log10(dee_geo_res$PValue.y)
hsa_cor=cor(dee_geo_res$dee_metric,dee_geo_res$geo_metric,method="s")

#sample wise correlation
colnames(x$GeneCounts)=gsub("$","_dee",colnames(x$GeneCounts))
colnames(b)=gsub("$","_geo",colnames(b))
d<-merge(x$GeneCounts,b,by=0)
SRR1692137=cor(d[,2:13],method="s")[1,7]
SRR1692138=cor(d[,2:13],method="s")[2,8]
SRR1692139=cor(d[,2:13],method="s")[3,9]
SRR1692140=cor(d[,2:13],method="s")[4,10]
SRR1692141=cor(d[,2:13],method="s")[5,11]
SRR1692142=cor(d[,2:13],method="s")[6,12]

hsa_res=c(hsa_cor,SRR1692137,SRR1692138,SRR1692139,SRR1692140,SRR1692141,SRR1692142)

###########
# M. musculus
###########

#DEE2
x<-getDEE2("mmusculus",c("SRR1533761","SRR1533762","SRR1533763","SRR1533764","SRR1533765","SRR1533766"))

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

#GEO
system("curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59970/suppl/GSE59970_mRNA_P1_P14_f.txt.gz | gunzip > GSE59970.tsv")
GSE59970<-read.table("GSE59970.tsv",header=T)
colnames(GSE59970)=c("SRR1533761","SRR1533763","SRR1533765","SRR1533762","SRR1533764","SRR1533766")
rownames(GSE59970)<-(gsub("_.*","",rownames(GSE59970)))
b<-GSE59970
b<-b[which(rowMeans(b)>10),]
group<-c(1,2,1,2,1,2)
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
dee_geo_res$dee_metric=sign(dee_geo_res$logFC.x)/-log10(dee_geo_res$PValue.x)
dee_geo_res$geo_metric=sign(dee_geo_res$logFC.y)/-log10(dee_geo_res$PValue.y)
mmu_cor=cor(dee_geo_res$dee_metric,dee_geo_res$geo_metric,method="s")

#sample wise correlation
colnames(x$GeneCounts)=gsub("$","_dee",colnames(x$GeneCounts))
colnames(b)=gsub("$","_geo",colnames(b))
d<-merge(x$GeneCounts,b,by=0)
SRR1533761=cor(d[,2:13],method="s")[1,7]
SRR1533762=cor(d[,2:13],method="s")[2,9]
SRR1533764=cor(d[,2:13],method="s")[3,11]
SRR1533766=cor(d[,2:13],method="s")[4,8]
SRR1533763=cor(d[,2:13],method="s")[5,10]
SRR1533765=cor(d[,2:13],method="s")[6,12]

mmu_res=c(mmu_cor,SRR1533761,SRR1533762,SRR1533763,SRR1533764,SRR1533765,SRR1533766)

###########
# R. norvigicus
###########
#R. norvegicus GSE65715 ctrl=c(“SRR1793792”,”SRR1793793”,”SRR1793794”), trt=c(“SRR1793795”,”SRR1793796”,”SRR1793797”)

###########
# S. cerevisiae
###########
#S. cerevisiae GSE19685 ctrl=c(“SRR039177”,”SRR039178”), trt=c(“SRR039179”,”SRR039179”)

#DEE2
x<-getDEE2("scerevisiae",c("SRR039177","SRR039178","SRR039179","SRR039180"))

x$GeneCounts<-x$GeneCounts[which(rowMeans(x$GeneCounts)>10),]
group<-c(1,2,1,2)
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

#GEO
system("curl \"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE19685&format=file\" > GSE19685.tar")
system("tar xf GSE19685.tar")
system("gunzip GSM49151[2-5].txt.gz")
SRR039177<-read.table("GSM491512.txt",header=F,row.names=1)
SRR039178<-read.table("GSM491513.txt",header=F,row.names=1)
SRR039179<-read.table("GSM491514.txt",header=F,row.names=1)
SRR039180<-read.table("GSM491515.txt",header=F,row.names=1)

cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
GSE19685<-as.data.frame(cbind.fill(SRR039177,SRR039178,SRR039179,SRR039180))

colnames(GSE19685)=c("SRR039177","SRR039178","SRR039179","SRR039180")
b<-GSE19685
b<-b[which(rowMeans(b)>10),]
group<-c(1,2,1,2)
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
dee_geo_res$dee_metric=sign(dee_geo_res$logFC.x)/-log10(dee_geo_res$PValue.x)
dee_geo_res$geo_metric=sign(dee_geo_res$logFC.y)/-log10(dee_geo_res$PValue.y)
sce_cor=cor(dee_geo_res$dee_metric,dee_geo_res$geo_metric,method="s")

#there is something very wrone with the geo entry for SRR039179 and SRR039180 which doesn't match the raw data.

#sample wise correlation
colnames(x$GeneCounts)=gsub("$","_dee",colnames(x$GeneCounts))
colnames(b)=gsub("$","_geo",colnames(b))
d<-merge(x$GeneCounts,b,by=0)

SRR039177=cor(d[,2:9],method="s")[1,5]
SRR039178=cor(d[,2:9],method="s")[2,6]
SRR039179=cor(d[,2:9],method="s")[3,7]
SRR039180=cor(d[,2:9],method="s")[4,8]

sce_res=c(sce_cor,SRR039177,SRR039178,SRR039179,SRR039180)

