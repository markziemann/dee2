library(reshape2)
library(gplots)
#install_github("GuangchuangYu/bitr")
library(clusterProfiler)
library("org.Sc.sgd.db")
#library("viridis")
library(RColorBrewer)

tmp<-read.table("scerevisiae_se.tsv")
x<-as.matrix(acast(tmp, V2~V1, value.var="V3"))

pass<-read.table("scerevisiae_pass.txt")
warn<-read.table("scerevisiae_warn.txt")

xp<-x[,which(colnames(x) %in% pass$V1)]
yp<-scale(xp)
cyp<-cor(yp,method="s")
tcyp<-cor(t(yp),method="s")
zp<-cmdscale(dist(t(yp)))
zzp<-cmdscale(dist(yp))

pdf("metaanalysis.pdf") 
#calculate the row mean and 
avcor<-cor(rowMeans(yp),yp,method="p")
hist(avcor,main="Pearson correlation to dataset average: Passed data")
mtext(paste("mean=",signif(mean(avcor),3)),cex=1.2)
avcor<-cor(rowMeans(yp),yp,method="s")
hist(avcor,main="Spearman correlation to dataset average: Passed data")
mtext(paste("mean=",signif(mean(avcor),3)),cex=1.2)


#noabel
plot(zp, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS by sample: Passed data", cex=0.6, pch=19)
#labeled
plot(zp, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS by sample: Passed data", cex=0.6, pch=19)
text( zp , labels=rownames(zp) ,cex=0.8) 
#nolabel
plot(zzp, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS by gene: Passed data", cex=0.6, pch=19) 
#label
plot(zzp, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS by gene: Passed data", cex=0.6, pch=19)
text(zzp , labels=rownames(zzp) ,cex=0.8)
colfunc <- colorRampPalette(c("blue", "white", "red"))

# Dataset level correlation
heatmap.2(  cyp, col=colfunc(25),scale="none", trace="none",margins = c(6,6), cexRow=.4, main="Sample correlation: Passed data")

c<-as.dist(1-cor(yp, method="spearman"))
hr <- hclust(c, method="complete")
mycl <- cutree(hr, h=max(hr$height/1.6))
#clusterCols <- rainbow(length(unique(mycl)))
clusterCols <- brewer.pal(length(unique(mycl)),"Paired")
myClusterSideBar <- clusterCols[mycl]
colfunc <- colorRampPalette(c("blue", "white", "red"))
write.table(mycl,file="DatasetClusters1.txt",quote=F,sep="\t")

heatmap.2(cyp, main="Dataset Clustering 1", Rowv=as.dendrogram(hr), Colv=as.dendrogram(hr),
 dendrogram="row", scale="none", col = colfunc(25), trace="none",
 RowSideColors= myClusterSideBar)

#Gene level correlation analysis
heatmap.2( tcyp, col=colfunc(25), scale="none",trace="none",margins = c(6,6), cexRow=.4, main="Gene correlation: Passed data")

c<-as.dist(1-cor(t(yp), method="spearman"))
hr <- hclust(c , method="complete")
mycl <- cutree(hr, h=max(hr$height/1.2))
clusterCols <- brewer.pal(length(unique(mycl)),"Paired")
myClusterSideBar <- clusterCols[mycl]
colfunc <- colorRampPalette(c("blue", "white", "red"))
write.table(mycl,file="GeneClusters1.txt",quote=F,sep="\t")

heatmap.2(tcyp, main="Gene Clustering 1", Rowv=as.dendrogram(hr), Colv=as.dendrogram(hr), 
 dendrogram="row", scale="none", col = colfunc(25), trace="none", 
 RowSideColors= myClusterSideBar)

warnpass<-rbind(pass,warn)
xw<-x[,which(colnames(x) %in% warnpass$V1)]
yw<-scale(xw)
cyw<-cor(yw,method="s")
tcyw<-cor(t(yw),method="s")
zw<-cmdscale(dist(t(yw)))
zzw<-cmdscale(dist(yw))

avcor<-cor(rowMeans(yw),yw,method="p")
hist(avcor,main="Pearson correlation to dataset average: Warn data")
mtext(paste("mean=",signif(mean(avcor),3)),cex=1.2)
avcor<-cor(rowMeans(yw),yw,method="s")
hist(avcor,main="Spearman correlation to dataset average: Warn data")
mtext(paste("mean=",signif(mean(avcor),3)),cex=1.2)
#nolabel
plot(zw, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS by sample: incl warn data", cex=0.6, pch=19)
#lebelled
plot(zw, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS by sample: incl warn data", cex=0.6, pch=19)
text(zw , labels=rownames(zw) ,cex=0.8)
#nolabel
plot(zzw, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS by gene: incl warn data", cex=0.6, pch=19)
#label
plot(zzw, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS by gene: incl warn data", cex=0.6, pch=19)
text(zzw , labels=rownames(zzw) ,cex=0.8)
#dataset correlation
heatmap.2(  cyw, col=colfunc(25),scale="none", trace="none",margins = c(6,6), cexRow=.4, main="Sample correlation: Pass and warn data")

#dataset clustering
c<-as.dist(1-cor(yw, method="spearman"))
hr <- hclust(c, method="complete")
mycl <- cutree(hr, h=max(hr$height/1.6))
clusterCols <- brewer.pal(length(unique(mycl)),"Paired")
myClusterSideBar <- clusterCols[mycl]
colfunc <- colorRampPalette(c("blue", "white", "red"))
write.table(mycl,file="DatasetClusters2.txt",quote=F,sep="\t")

heatmap.2(cyw, main="Dataset Clustering 2", Rowv=as.dendrogram(hr), Colv=as.dendrogram(hr),
 dendrogram="row", scale="none", col = colfunc(25), trace="none",
 RowSideColors= myClusterSideBar)

#gene correlation
heatmap.2( tcyw, col=colfunc(25), scale="none",trace="none",margins = c(6,6), cexRow=.4, main="Gene correlation: Pass and warn data")

# gene clustering
c<-as.dist(1-cor(t(yw), method="spearman"))
hr <- hclust(c , method="complete")
mycl <- cutree(hr, h=max(hr$height/1.2))
clusterCols <- brewer.pal(length(unique(mycl)),"Paired")
myClusterSideBar <- clusterCols[mycl]
colfunc <- colorRampPalette(c("blue", "white", "red"))
write.table(mycl,file="GeneClusters2.txt",quote=F,sep="\t")

heatmap.2(tcyw, main="Gene Clustering 2", Rowv=as.dendrogram(hr), Colv=as.dendrogram(hr),
 dendrogram="row", scale="none", col = colfunc(25), trace="none",
 RowSideColors= myClusterSideBar)
dev.off()



g<-read.table("GeneClusters2.txt")

bg<-unique(rownames(g))

bge<-bitr(as.character(bg), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Sc.sgd.db")

bge<-as.character(bge$ENTREZID)

res=NULL
for (n in unique(g$x) ) {
  c<-rownames(subset(g,x==n))
  ce<-bitr(as.character(c), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Sc.sgd.db")
  ce<-as.character(ce$ENTREZID)
  r<-enrichGO(ce,OrgDb="org.Sc.sgd.db",pvalueCutoff = 0.5,pAdjustMethod = "fdr",universe=bge)
  #r<-enrichPathway(gene=ce,organism = "yeast",pvalueCutoff = 0.5,pAdjustMethod = "BH", universe=bge)
  r<-r@result
  r<-r[which(r$p.adjust<0.1),]
  r<-r[1:3,]
  r$cluster<-n
  res<-rbind(res,r)
}

res$sign<- -log10(res$p.adjust)
cols<-as.data.frame(brewer.pal(10, "Paired"))
#cols<-as.data.frame(magma(10))
#cols<-as.data.frame(brewer.pal(10, "Spectral"))
#cols<-as.data.frame(rainbow(10))
cols$num<-1:10
colnames(cols)=c("colour","num")

res<-merge(res,cols,by.x="cluster",by.y="num")

res$Description<-strtrim(res$Description,80)


pdf("GO_chart.pdf")
par(mai=c(1,4,1,0.5))
barplot(rev(res$sign),main="Enriched molecular function", horiz=TRUE,names.arg=rev(res$Description),las=1,cex.axis=1,cex.names=0.6,xlab="-log10(FDR)",col=as.character(rev(res$colour)))
dev.off()

