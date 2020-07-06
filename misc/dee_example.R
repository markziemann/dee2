library("getDEE2")

mdat <- getDEE2Metadata("hsapiens")
head(mdat)
dim(mdat)

mdat[which(mdat$SRP_accession == "SRP096177"),]
mdat1 <- mdat[which(mdat$SRP_accession == "SRP096177"),]
SRRvec <- as.vector(mdat1$SRR_accession)
SRRvec

## Fetching DEE2 data using SRA run accession numbers
x <- getDEE2("hsapiens", SRRvec, metadata=mdat, legacy=TRUE)
names(x)
head(x$GeneCounts)
colSums(x$GeneCounts)

# aggregate transcript counts to genes
x <- Tx2Gene(x)
head(x$Tx2Gene)
colSums(x$Tx2Gene)
