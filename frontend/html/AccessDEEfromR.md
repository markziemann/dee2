# Incorporating dee2 data into your RNA-seq workflow

This tutorial provides a walkthrough for how to work with dee2 expression data,
starting with data searches, obtaining the data from dee2.io and then performing
a differential analysis with DESeq or edgeR.

To use search and obtain dee2 data in R, you will need to source the 'getDEE2' script as follows:

`> source("https://raw.githubusercontent.com/markziemann/dee2/master/frontend/html/getDEE2.R")`

## Searching for datasets of interest starting with accession numbers

The first step is to download the list of accession numbers of available datasets with the 'getDee2Metadata' function, specifying a species name. The options for species currently are:

* athaliana
* celegans
* dmelanogaster
* drerio
* ecoli
* hsapiens
* mmusculus
* rnorvegicus
* scerevisiae

If the species name is incorrect, an error will be thrown.

```
> mdat<-getDee2Metadata("celegans")
trying URL 'http://dee2.io/metadata/celegans_metadata.tsv.cut'
Content type 'text/tab-separated-values' length 549206 bytes (536 KB)
==================================================
downloaded 536 KB

> head(metadata)
  SRR_accession QC_summary SRX_accession SRS_accession SRP_accession
1    SRR1557807       PASS     SRX686598     SRS689414     SRP045778
2    SRR2537196       PASS    SRX1295772    SRS1094603     SRP064324
3    SRR1658975       PASS     SRX765271     SRS749590     SRP050116
4     SRR087429       PASS     SRX035773     SRS150525     SRP004901
5    SRR7443620       PASS    SRX4314169    SRS3473340     SRP151477
6    SRR3535786       PASS    SRX1770003    SRS1442527     SRP075276
  GSE_accession GSM_accession
1      GSE60755    GSM1487404
2      GSE73589    GSM1898659
3      GSE63528    GSM1551798
4      GSE26165     GSM642428
5          <NA>          <NA>
6      GSE81522    GSM2155077
```

If we have a GEO series accession number in mind already (eg: GSE33569) then we can see if the datasets are present.

```> mdat[which(mdat$GSE_accession %in% "GSE33569"),]
     SRR_accession QC_summary SRX_accession SRS_accession SRP_accession
2127     SRR363796       PASS     SRX105188     SRS270025     SRP009256
4608     SRR363797       PASS     SRX105189     SRS270026     SRP009256
5976     SRR363798       PASS     SRX105190     SRS270027     SRP009256
8072     SRR363799       PASS     SRX105191     SRS270028     SRP009256
     GSE_accession GSM_accession
2127      GSE33569     GSM829554
4608      GSE33569     GSM829555
5976      GSE33569     GSM829556
8072      GSE33569     GSM829557
```

Dee2 data is centred around SRA run accessions numbers, these SRR_accessions can be obtained like this:

```
> mdat1<-mdat[which(mdat$GSE_accession %in% "GSE33569"),]
> SRRlist<-as.vector(mdat1$SRR_accession)
> SRRlist
[1] "SRR363796" "SRR363797" "SRR363798" "SRR363799"
```

## Fetching expression data from dee2.io

The general syntax for obtaining dee2 data is this:

`> getDEE2(species,SRRlist,outfile="NULL")`

First, the function downloads a copy of the metadata, then runs a query to make sure that the requested datasets are present. It then downloads the requested expression data as a zip archive that contains the following:
1. a gene-wise count expression matrix, 
2. a transcript-wise count expression matrix,
3. a matrix of quality metrics, and
4. a folder of run logs detailing the processing of the data including base quality scores, alignment rates, etc.

If 'outfile' is defined, then files will be downloaded to the current working directory. If it is not defined, then the files are downloaded to the temporary directory of R and deleted immediately after use.

The SRR numbers need to exactly match those in SRA.

Here is an example of using the SRRlist as defined above. 
```
> x<-getDEE2("celegans",SRRlist)
trying URL 'http://dee2.io/metadata/celegans_metadata.tsv.cut'
Content type 'text/tab-separated-values' length 549206 bytes (536 KB)
==================================================
downloaded 536 KB

trying URL 'http://dee2.io/cgi-bin/request.sh?org=celegans&x=SRR363796&x=SRR363797&x=SRR363798&x=SRR363799'
downloaded 479 KB

> x$
x$GeneCounts  x$TxCounts    x$QcMx        x$absent      
> head(x$GeneCounts)
               SRR363796 SRR363797 SRR363798 SRR363799
WBGene00197333         0         0         0         0
WBGene00198386         0         0         0         0
WBGene00015153         4        16         6         4
WBGene00002061        44       100       217        77
WBGene00255704         0         5         1         1
WBGene00235314         0         0         0         1
> head(x$TxCounts)
           SRR363796 SRR363797 SRR363798 SRR363799
Y110A7A.10        11        23        48        45
F27C8.1            0         0         0         0
F07C3.7            0         2         0        21
F52H2.2a           0         0         7         0
F52H2.2b           0         4         0         5
T13A10.10a         0         0         0         0
> head(x$QcMx)
                                   SE               SE.1               SE.2
SequenceFormat     Sanger/Illumina1.9 Sanger/Illumina1.9 Sanger/Illumina1.9
QualityEncoding                    36                 36                 36
Read1MinimumLength                 36                 36                 36
Read1MedianLength                  36                 36                 36
Read1MaxLength                   NULL               NULL               NULL
Read2MinimumLength               NULL               NULL               NULL
                                 SE.3
SequenceFormat     Sanger/Illumina1.9
QualityEncoding                    36
Read1MinimumLength                 36
Read1MedianLength                  36
Read1MaxLength                   NULL
Read2MinimumLength               NULL
```

You can directly specify the SRR accessions in the command line, but be sure to type them correctly. In case SRR accessions are not present in the database, there will be a warning message.

```
> x<-getDEE2("celegans",c("SRR363798","SRR363799","SRR3581689","SRR3581692"))
trying URL 'http://dee2.io/metadata/celegans_metadata.tsv.cut'
Content type 'text/tab-separated-values' length 549206 bytes (536 KB)
==================================================
downloaded 536 KB

trying URL 'http://dee2.io/cgi-bin/request.sh?org=celegans&x=SRR363798&x=SRR363799'
downloaded 369 KB

Warning, datasets not found: 'SRR3581689,SRR3581692'
```

In this case the accessions SRR3581689 and SRR3581692 are A. thaliana accessions and therefore not present in the C. elegans accession list. The list of absent accessions is provided in case you need these for your records.
```
> x$absent
[1] "SRR3581689" "SRR3581692"
```
## Keyword searching metadata

* Use the SRAdb package to make queries of the SRA metadata
* Intersect SRAdb hits with dee2 accession numbers. 




## Stand-alone functions for downloading and loading dee2 data

In case you want to download the data once, then reuse it many times, using the standalone scripts can be faster. Data is saved in zip format as follows:

```
x<-getDEE2("celegans",SRRlist,outfile="DEE_count_data.zip")
```

Loading data from previously downloaded zip files is done as follows
```
> myGeneCounts<-loadGeneCounts("DEE_count_data.zip")
> myTxCounts<-loadTxCounts("DEE_count_data.zip")
> myQcMx<-loadQcMx("DEE_count_data.zip")
> head(myGeneCounts)
               SRR363796 SRR363797 SRR363798 SRR363799
WBGene00197333         0         0         0         0
WBGene00198386         0         0         0         0
WBGene00015153         4        16         6         4
WBGene00002061        44       100       217        77
WBGene00255704         0         5         1         1
WBGene00235314         0         0         0         1
```





## Aggregating runs data

In case you need to aggregate runs that are technical replicates, read this
thread (http://stackoverflow.com/questions/26046776/sum-two-columns-in-r). Here
is an example:

```
> head(mytable)
             SRR922260v1 SRR922261v1
ECDH10B_0001         470         452
ECDH10B_0002        6046        8138
ECDH10B_0003        1758        2918
ECDH10B_0004        3082        4188
ECDH10B_0005          95         268
ECDH10B_0006         169         734
> mytable$SRR922260v1_SRR922261v1<-mytable$SRR922260v1 + mytable$SRR922261v1
> head(mytable)
             SRR922260v1 SRR922261v1 SRR922260v1_SRR922261v1
ECDH10B_0001         470         452                     922
ECDH10B_0002        6046        8138                   14184
ECDH10B_0003        1758        2918                    4676
ECDH10B_0004        3082        4188                    7270
ECDH10B_0005          95         268                     363
ECDH10B_0006         169         734                     903
```
## Running a differential analysis
* edgeR
* DESeq

## Report bugs, issues and suggestions at the github page
https://github.com/markziemann/dee2
