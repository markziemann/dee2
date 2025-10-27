library("parallel")
library("data.table")
library("R.utils")
library("reutils")
library("XML")
library("rjson")
library("rhdf5")
library("reutils")

IPADD="118.138.235.221"

CORES=5

#start the analysis

#for ( org in c( "athaliana", "osativa", "zmays",
# "celegans", "dmelanogaster", "drerio",
#"rnorvegicus", "scerevisiae" , "mmusculus", "ecoli", "hsapiens" )) {

args = commandArgs(trailingOnly=TRUE)
org=args[1]

species_list <- c("3702", "4530", "4577", "6239", "7227", "7955",
  "562", "9606", "10090", "10116", "4932",
  "15368", "3847", "1753", "3694", "4558",
  "1753", "4113", "4565", "1753")

#now annotate the short names 
names(species_list) <- c("athaliana", "osativa", "zmays",
  "celegans", "dmelanogaster", "drerio",
  "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae",
  "bdistachyon", "gmax", "hvulgare", "ptrichocarpa", "sbicolor",
  "slycopersicum", "stuberosum", "taestivum", "vvinifera")

taxa_name <- species_list[[org]]

species_names <- c("'Arabidopsis thaliana'", "'Oryza sativa'",
  "'Zea mays'", "'Caenorhabditis elegans'",
  "'Drosophila melanogaster'", "'Danio rerio'",
  "'Escherichia coli'", "'Homo sapiens'", "'Mus musculus'",
  "'Rattus norvegicus'", "'Saccharomyces cerevisiae'",
  "'Brachypodium distachyon'", "'Glycine max'",
  "'Hordeum vulgare'", "'Populus trichocarpa'",
  "'Sorghum bicolor'","'Solanum lycopersicum'",
  "'Solanum tuberosum'", "'Triticum aestivum'",
  "'Vitis vinifera'")

#now annotate the short names 
names(species_names)<- c("athaliana", "osativa", "zmays",
  "celegans", "dmelanogaster", "drerio",
  "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae",
  "bdistachyon", "gmax", "hvulgare", "ptrichocarpa", "sbicolor",
  "slycopersicum", "stuberosum", "taestivum", "vvinifera")

species_name <- species_names[[org]]

print(org)

#number of protein coding genes from ensembl
numgenes <- c(27655, 35806, 39756, 20362, 13918, 25903, 4140, 20338, 22598, 22250, 6692,
  34310, 55897, 35826, 34699, 34118, 34658, 39021, 107891, 37480)

names(numgenes) <- c("athaliana", "osativa", "zmays",
  "celegans", "dmelanogaster", "drerio",
  "ecoli", "hsapiens", "mmusculus", "rnorvegicus", "scerevisiae",
  "bdistachyon", "gmax", "hvulgare", "ptrichocarpa", "sbicolor",
  "slycopersicum", "stuberosum", "taestivum", "vvinifera")

numgenes = numgenes[[org]]

###Set some directories
CODEWD = getwd()
DATAWD = paste(normalizePath("../data/"),org,sep="/")
SRADBWD = normalizePath("../sradb/")
MXDIR = normalizePath("../mx/")
QUEUEWD = normalizePath("../queue/")

### Set years range
YEARS=2000:2025

#collect QC info - this is temporary and logic will be incorporated in future
QC_summary="BLANK"

# here we use reutils use https://www.rdocumentation.org/packages/reutils/versions/0.2.2
GEO_QUERY_TERMS <- c(
'"Arabidopsis thaliana"[porgn:__txid3702] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Oryza sativa"[porgn:__txid4530] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Zea mays"[porgn:__txid4577] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Caenorhabditis elegans"[porgn:__txid6239] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Drosophila melanogaster"[porgn:__txid7227] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Danio rerio"[porgn:__txid7955] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Escherichia coli"[porgn:__txid562] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Homo sapiens"[porgn:__txid9606] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Mus musculus"[porgn:__txid10090] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Rattus norvegicus"[porgn:__txid10116] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Saccharomyces cerevisiae"[porgn:__txid4932] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Brachypodium distachyon"[porgn:__txid15368] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Glycine max"[porgn:__txid3847] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Hordeum vulgare"[porgn:__txid4513] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Populus trichocarpa"[porgn:__txid3694] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Sorghum bicolor"[porgn:__txid4558] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Solanum lycopersicum"[porgn:__txid4081] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Solanum tuberosum"[porgn:__txid4113] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Triticum aestivum"[porgn:__txid4565] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]',
'"Vitis vinifera"[porgn:__txid29760] AND "gsm"[Filter])) AND "high throughput sequencing"[Platform Technology Type]')

MY_GEO_QUERY_TERM <- GEO_QUERY_TERMS[grep(gsub("'","",species_name),GEO_QUERY_TERMS)]

# helper funciton
# as the geo resource is unreliable I need to use trycatch to prevent errors
# from terminating the script
fetch_esummary <- function(i) {
  out <- tryCatch(
    {
      message(i)
      esummary(ESEARCH,retstart=i,retmax=500,retmode="json")
    },error=function(cond) {
      message(paste("There was an error processing request",i))
      NULL
    }
  )
}

## Data for earlier years is already archived

YEAR=2025
#YEAR <- sapply(strsplit(as.character(Sys.Date()),"-"),"[[",1)
#MINDATE=paste(YEAR,"/01/01",sep="")
MINDATE=paste("2024","/01/01",sep="")
MAXDATE=paste(YEAR,"/12/31",sep="")
GEOFILENAME=paste(SRADBWD,"/",org,"_geo_",YEAR,".csv",sep="")
Sys.sleep(1)

# so now we are going with JSON format because XML was failing due to 
# angled brackets in the GEO data
ESEARCH <- esearch(term = MY_GEO_QUERY_TERM , db = "gds", rettype = "uilist",
  retmode = "json", retstart = 0, retmax = 500000000, usehistory = TRUE,
  webenv = NULL, querykey = NULL, sort = NULL, field = NULL,
  datetype = "pdat", reldate = NULL, mindate = MINDATE, maxdate = MAXDATE)

j <- fromJSON(ESEARCH$content)
COUNT <- as.numeric(j$esearchresult$count)
message(paste("number of GEO results in year",YEAR,":",COUNT))

if (COUNT>1) {
  myrange <- seq(0,COUNT,500)
  gsel <- lapply(myrange, function(i) {
    message(i)
    Sys.sleep(1)
    ESUMMARY <- NULL
    attempt <- 1
    while ( is.null(ESUMMARY) ) {
      message(paste("attempt",attempt))
      attempt <- attempt + 1
      ESUMMARY <- fetch_esummary(i)
    }
    myjson <- fromJSON(ESUMMARY$content)
    myjson <- myjson$result
    myjson[1] = NULL
    geodf <- do.call(rbind, myjson)
    geodf <- geodf[,c("gse","accession")]
    return(geodf)
  })
}

gsel <- do.call(rbind, gsel)
GSE <- unlist(gsel[,1])
GSM <- unlist(gsel[,2])
gse <- data.frame(GSE,GSM)
colnames(gse) <- c("GEO_series","GEO_sample")
gse$GEO_series[which(gse$GEO_series=="")] <- "NA"
gse$GEO_series <- paste("GSE",sapply(strsplit(gse$GEO_series,";"),"[[",1),sep="")
write.table(gse,GEOFILENAME,sep=",")

# now join old and new
PATTERN=paste(org,"_geo_20",sep="")
geofiles <- list.files("../sradb/", pattern=PATTERN,full.names=TRUE)
gf <- lapply(geofiles,read.csv)
gf <- do.call(rbind,gf)
gse <- unique(gf)

# write a backup of GEO data
GEONEW=paste(SRADBWD,"/",org,"_all.csv",sep="")
write.table(gse,GEONEW,sep=",")
