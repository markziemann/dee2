# UPDATE WEBPAGE

library(xtable)

DATE = date()

rowcnt2 <- function( file) { z<-system(paste("cat ",file, "| cut -f1 | sed 1d | sort -u | wc -l" ) , intern=TRUE) ; z}

# Collect run data stats
FILES1 <- list.files(pattern="*metadata.complete.tsv.cut$",path="/mnt/hdd1/dee2/sradb/",full.names=T)
x <- as.data.frame(sapply(FILES1,rowcnt2),stringsAsFactors=FALSE)

orgnames <- c("A. thaliana", "B. distachyon", "C. elegans",
  "D. melanogaster", "D. rerio", "E. coli", "G. max",
  "H. sapiens", "H. vulgare", "M. musculus", "O. sativa",
  "P. trichocarpa", "R. norvegicus", "S. bicolor",
  "S. cerevisiae", "S. lycopersicum", "S. tuberosum",
  "T. aestivum", "V. vinifera", "Z. mays")

rownames(x) <- orgnames

colnames(x)="total"

FILES2 <- list.files(pattern="*_metadata.tsv$",path="/mnt/hdd1/dee2/sradb/",full.names=T)
y <- as.data.frame(sapply(FILES2,rowcnt2),stringsAsFactors=FALSE)

rownames(y) <- orgnames

colnames(y) = "completed"

z <- merge(x,y,by=0)
rownames(z) = z$Row.names
z$Row.names = NULL

zz <- as.data.frame(apply(z,2,as.integer))
zz$queued <- zz$total - zz$completed
rownames(zz) <- rownames(z)
z <- zz
z$pc_complete <- signif(z$completed / z$total * 100 ,3)

# Collect bundle data stats

bundles <- list.files("../sradb/big_proj/",pattern="*zip$",recursive=TRUE)
bundle_tbl <- table(sapply(strsplit(bundles,"/"),"[[",1))
z$n_bundles <- bundle_tbl

plants <- c("Z. mays", "V. vinifera", "T. aestivum","S. tuberosum","S. lycopersicum", "S. bicolor",
  "P. trichocarpa","O. sativa","H. vulgare","G. max","B. distachyon","A. thaliana")

animals <- c("R. norvegicus","M. musculus", "H. sapiens", "D. rerio", "D. melanogaster", "C. elegans")

microbes <- c("S. cerevisiae","E. coli")

plantdf <- z[rownames(z) %in% plants,]
animaldf <- z[rownames(z) %in% animals,]
microbedf <- z[rownames(z) %in% microbes,]

animaldf$group="Animal"
plantdf$group="Plant"
microbedf$group="Microbe"
df <- rbind(animaldf,plantdf,microbedf)

df <- df[,c(6,1:5)]

colnames(df) <- c("Group","Runs total","Runs completed","Runs queued","% Runs completed","No. bundles completed")

## Chart of 
png("dee_datasets.png",width=600,height=520)
options(bitmapType="cairo")
layout.matrix <- matrix(c(1, 3, 2, 2), nrow = 2, ncol = 2)

layout(mat=layout.matrix, heights=c(2,1),widths=c(2,2))

par(las=2) ; par(mai=c(0.5,2,0.5,0.2))
MAX=1800000
z2 <- animaldf
z2 <- z2[seq(nrow(z2),1),]
bb <- barplot( t(z2[,c("completed","queued")]), beside=TRUE , xlim=c(0,MAX), main="Animals", col=c("#FFDD00","#0057B7") ,
  horiz=T , las=1, cex.axis=1.2, cex.names=1.3, cex.main=1.4, xlab="number of SRA runs")
text( cbind( as.numeric(z2[,2])+220000 , as.numeric(z2[,3])+220000 ) * 1.1 , t(bb),labels=c(z2[,2],z2[,3]) ,cex=1.15)
legend("topright", c("Runs queued", "Runs completed")  , fill=c("#0057B7","#FFDD00") , cex=0.9)

par(las=2) ; par(mai=c(0.5,2,0.5,0.2))
MAX=180000
z2 <- plantdf
z2 <- z2[seq(nrow(z2),1),]
bb <- barplot( t(z2[,c("completed","queued")]), beside=TRUE , xlim=c(0,MAX), main="Plants", col=c("#FFDD00","#0057B7") ,
  horiz=T , las=1, cex.axis=1.2, cex.names=1.3, cex.main=1.4, xlab="number of SRA runs")
text( cbind( as.numeric(z2[,2])+18000 , as.numeric(z2[,3])+18000 ) * 1.1 , t(bb),labels=c(z2[,2],z2[,3]) ,cex=1.15)

par(las=2) ; par(mai=c(0.5,2,0.5,0.2))
MAX=70000
z2 <- microbedf
z2 <- z2[seq(nrow(z2),1),]
bb <- barplot( t(z2[,c("completed","queued")]), beside=TRUE , xlim=c(0,MAX), main="Microbes", col=c("#FFDD00","#0057B7") ,
  horiz=T , las=1, cex.axis=1.2, cex.names=1.3, cex.main=1.4, xlab="number of SRA runs")
text( cbind( as.numeric(z2[,2])+6000 , as.numeric(z2[,3])+6000 ) * 1.1 , t(bb),labels=c(z2[,2],z2[,3]) ,cex=1.15)

dev.off()

system("scp -i ~/.ssh/dee2_2026 dee_datasets.png ubuntu@dee2.io:/home/ubuntu/dee2/frontend/html/images")

# export the table
print(xtable(df), type = "html",file = 'table1.html')
system("scp -i ~/.ssh/dee2_2026 table1.html ubuntu@dee2.io:/home/ubuntu/dee2/frontend/html/table1.html")
