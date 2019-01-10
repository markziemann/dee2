s1<-read.table("s1.txt")
s2<-read.table("s2.txt")

s1<-subset(s1,V2>29)
s2<-subset(s2,V2>29)

s1mean=mean(s1$V2)
s2mean=mean(s2$V2)

s<-rbind(s1,s2)
smean=mean(s$V2)
n=nrow(s)

pdf("time.pdf")

layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))

par(mar=c(0, 4.1, 1.1, 2.1))

boxplot(s$V2 , horizontal=TRUE , ylim=c(0,5000), xaxt="n" , col="gray" , frame=F)

par(mar=c(4, 4.1, 1.1, 2.1))

hist(s$V2,breaks=seq(1,50000,by=200),xlim=c(0,5000),xlab="seconds",main=paste("elapsed time per run;","mean=",signif(smean,4),"s;","n=",n) )

abline(v=smean,col="red",lwd=2, lty=2)

dev.off()
