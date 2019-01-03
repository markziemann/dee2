x<-read.table("smry.tsv",sep="\t",header=T)

cols<-as.data.frame(rainbow(10))
cols$num<-1:10
colnames(cols)=c("colour","num")

x<-merge(x,cols,by.x="Cluster",by.y="num")


pdf("GO_chart.pdf")
par(mai=c(1,4,1,0.5))
barplot(rev(x$Significance),main="Enriched molecular function", horiz=TRUE,names.arg=rev(x$TERM),las=1,cex.axis=1,cex.names=0.6,xlab="-log10(FDR)",col=as.character(rev(x$colour)))
dev.off()
