
# UPDATE WEBPAGE

FILES1 <- list.files(pattern="*metadata.complete.tsv.cut$",path="/mnt/hdd1/dee2/sradb/",full.names=T)
x <- as.data.frame(sapply(FILES1,rowcnt2),stringsAsFactors=FALSE)

rownames(x)=c("A. thaliana", "B. distachyon", "C. elegans",
  "D. melanogaster", "D. rerio", "E. coli", "G. max",
  "H. sapiens", "H. vulgare", "M. musculus", "O. sativa",
  "P. trichocarpa", "R. norvegicus", "S. bicolor",
  "S. cerevisiae", "S. lycopersicum", "S. tuberosum",
  "T. aestivum", "V. vinifera", "Z. mays")

colnames(x)="queued"

FILES2 <- list.files(pattern="*_metadata.tsv$",path="/mnt/hdd1/dee2/sradb/",full.names=T)
y <- as.data.frame(sapply(FILES2,rowcnt2),stringsAsFactors=FALSE)

rownames(y)=c("A. thaliana", "B. distachyon", "C. elegans",
  "D. melanogaster","D. rerio", "E. coli", "G. max",
  "H. sapiens", "H. vulgare", "M. musculus", "O. sativa",
  "P. trichocarpa", "R. norvegicus", "S. bicolor",
  "S. cerevisiae", "S. lycopersicum", "S. tuberosum",
  "T. aestivum", "V. vinifera", "Z. mays")

colnames(y) = "completed"

z <- merge(x,y,by=0)
rownames(z) = z$Row.names
z$Row.names = NULL

zz <- as.data.frame(apply(z,2,as.numeric))
zz$queued <- zz$queued - zz$completed
rownames(zz) <- rownames(z)
z <- zz

DATE = date()
HEADER = paste("Updated",DATE)
z <- z[order(rownames(z),decreasing=T ), ,drop=F]

plants <- c("Z. mays", "V. vinifera", "T. aestivum","S. tuberosum","S. lycopersicum", "S. bicolor",
  "P. trichocarpa","O. sativa","H. vulgare","G. max","B. distachyon","A. thaliana")

animals <- c("R. norvegicus","M. musculus", "H. sapiens", "D. rerio", "D. melanogaster", "C. elegans")

microbes <- c("S. cerevisiae","E. coli")

plantdf <- z[rownames(z) %in% plants,]
animaldf <- z[rownames(z) %in% animals,]
microbedf <- z[rownames(z) %in% microbes,]

png("dee_datasets.png",width=600,height=600)
options(bitmapType="cairo")
layout.matrix <- matrix(c(1, 3, 2, 2), nrow = 2, ncol = 2)

layout(mat=layout.matrix, heights=c(2,1),widths=c(2,2))

par(las=2) ; par(mai=c(0.5,2,0.5,0.2))
MAX=1000000
z2 <- animaldf
bb <- barplot( t(z2), beside=TRUE , xlim=c(0,MAX), main="Animals", col=c("#FFDD00","#0057B7") ,
  horiz=T , las=1, cex.axis=1.2, cex.names=1.3, cex.main=1.4, xlab="number of SRA runs")
legend("topright", rev(colnames(z2)), fill=c("#0057B7","#FFDD00") , cex=1.1)
text( cbind(as.numeric(z2[,1])+150000 ,as.numeric(z2[,2])+150000 ),
  t(bb),labels=c(z2[,1],z2[,2]) ,cex=1.15)

par(las=2) ; par(mai=c(0.5,2,0.5,0.2))
MAX=150000
z2 <- plantdf
bb <- barplot( t(z2), beside=TRUE , xlim=c(0,MAX), main="Plants", col=c("#FFDD00","#0057B7") ,
  horiz=T , las=1, cex.axis=1.2, cex.names=1.3, cex.main=1.4, xlab="number of SRA runs")
mtext(HEADER, las=1, adj=1, cex=0.8)
text( cbind(as.numeric(z2[,1])+20000 ,as.numeric(z2[,2])+20000 ),
  t(bb),labels=c(z2[,1],z2[,2]) ,cex=1.15)

par(las=2) ; par(mai=c(0.5,2,0.5,0.2))
MAX=55000
z2 <- microbedf
bb <- barplot( t(z2), beside=TRUE , xlim=c(0,MAX), main="Microbes", col=c("#FFDD00","#0057B7"),
horiz=T , las=1, cex.axis=1.2, cex.names=1.3, cex.main=1.4, xlab="number of SRA runs")
text( cbind(as.numeric(z2[,1])+8000 ,as.numeric(z2[,2])+8000 ),
  t(bb),labels=c(z2[,1],z2[,2]) ,cex=1.15)

dev.off()

system("scp -i ~/.ssh/dee2_2026 dee_datasets.png ubuntu@dee2.io:/home/ubuntu/dee2/frontend/html/images")

animaldf$group="Animal"
plantdf$group="Plant"
microbedf$group="Microbe"
df <- rbind(animaldf,plantdf,microbedf)
df$pc_complete <- signif(df$completed /  (df$completed +  df$queued) *100 ,3)

df <- df[,c(3,1,2,4)]
colnames(df) <- c("Group","Queued","Completed","% Complete")

library(gridExtra)
library(gtable)
library(grid)

# Sample data
data <- data.frame(
  System = c("Neural Link", "Cyber Deck", "Bio Monitor", "Optic Array"),
  Status = c("ACTIVE", "STANDBY", "ACTIVE", "OFFLINE"),
  Capacity = c("98%", "76%", "100%", "0%"),
  Uptime = c("247h", "180h", "312h", "0h")
)

# Cyberpunk color scheme
bg_dark <- "#0a0e27"
bg_header <- "#1a1f3a"
text_cyan <- "#00ffff"
text_magenta <- "#ff00ff"
text_green <- "#39ff14"
text_red <- "#ff0055"
border_neon <- "#00ffff"

# Custom theme function
cyberpunk_theme <- ttheme_minimal(
  core = list(
    fg_params = list(
      col = text_cyan,
      fontface = "bold",
      fontfamily = "mono",
      cex = 1.1
    ),
    bg_params = list(
      fill = c(rep(c(bg_dark, "#0f1435"), length.out = nrow(data))),
      col = border_neon,
      lwd = 2
    )
  ),
  colhead = list(
    fg_params = list(
      col = text_magenta,
      fontface = "bold",
      fontfamily = "mono",
      cex = 1.2
    ),
    bg_params = list(
      fill = bg_header,
      col = border_neon,
      lwd = 3
    )
  ),
  rowhead = list(
    fg_params = list(col = text_green),
    bg_params = list(
      fill = bg_header,
      col = border_neon
    )
  )
)

# Draw the table
png("dee_datasets2.png",width=480,height=550,bg=NA)
grid.newpage()
grid.rect(gp = gpar(fill = bg_dark, col = NA))
g <- tableGrob(df, theme = cyberpunk_theme)

# Add neon glow effect (outer border)
g <- gtable_add_grob(
  g,
  grobs = rectGrob(gp = gpar(col = text_magenta, fill = NA, lwd = 4)),
  t = 1, l = 1, b = nrow(g), r = ncol(g),
  z = 0
)
grid.draw(g)

# Add title with glow effect
grid.text(
  "⚡ DEE2 Data by SRA Run ⚡",
  x = 0.5, y = 0.95,
  gp = gpar(
    col = text_cyan,
    fontsize = 20,
    fontface = "bold",
    fontfamily = "mono"
  )
)

dev.off()
system("scp -i ~/.ssh/dee2_2026 dee_datasets2.png ubuntu@dee2.io:/home/ubuntu/dee2/frontend/html/images")

## UPDATE BUNDLES

bundles <- list.files("../sradb/big_proj/",pattern="*zip$",recursive=TRUE)

bundle_tbl <- table(sapply(strsplit(bundles,"/"),"[[",1))

names(bundle_tbl) <- c("A. thaliana", "B. distachyon", "C. elegans",
  "D. melanogaster","D. rerio", "E. coli", "G. max",
  "H. sapiens", "H. vulgare", "M. musculus", "O. sativa",
  "P. trichocarpa", "R. norvegicus", "S. bicolor",
  "S. cerevisiae", "S. lycopersicum", "S. tuberosum",
  "T. aestivum", "V. vinifera", "Z. mays")

animal_bundles <- bundle_tbl[names(bundle_tbl) %in% animals]
plant_bundles <- bundle_tbl[names(bundle_tbl) %in% plants]
microbe_bundles <- bundle_tbl[names(bundle_tbl) %in% microbes]

animal_bundles <- data.frame(animal_bundles,row.names=names(animal_bundles))
plant_bundles <- data.frame(plant_bundles,row.names=names(plant_bundles))
microbe_bundles <- data.frame(microbe_bundles,row.names=names(microbe_bundles))

animal_bundles$Var1 <- "Animal"
plant_bundles$Var1 <- "Plant"
microbe_bundles$Var1 <- "Microbe"

bundle_tbl2 <- rbind(animal_bundles,plant_bundles,microbe_bundles)
colnames(bundle_tbl2) <- c("Group","No. bundles")

df2 <- bundle_tbl2

# Draw the table
png("dee_bundles.png",width=480,height=550,bg=NA)
grid.newpage()
grid.rect(gp = gpar(fill = bg_dark, col = NA))
g <- tableGrob(df2, theme = cyberpunk_theme)

# Add neon glow effect (outer border)
g <- gtable_add_grob(
  g,
  grobs = rectGrob(gp = gpar(col = text_magenta, fill = NA, lwd = 4)),
  t = 1, l = 1, b = nrow(g), r = ncol(g),
  z = 0
)
grid.draw(g)

# Add title with glow effect
grid.text(
  "⚡ DEE2 Data Bundles ⚡",
  x = 0.5, y = 0.95,
  gp = gpar(
    col = text_cyan,
    fontsize = 20,
    fontface = "bold",
    fontfamily = "mono"
  )
)

dev.off()
system("scp -i ~/.ssh/dee2_2026 dee_bundles.png ubuntu@dee2.io:/home/ubuntu/dee2/frontend/html/images")
