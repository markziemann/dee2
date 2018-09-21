# Incorporating dee2 data into your RNA-seq workflow

Copyright, Mark Ziemann, 2018.

## Table of Contents
1. [Getting started](#getting-started)
2. [Searching for datasets of interest starting with accession numbers](#searching-for-datasets-of-interest-starting-with-accession-numbers)
3. [Fetching dee2 data using SRA run accession numbers](#fetching-dee2-data-using-sra-run-accession-numbers)
4. [Keyword searching metadata](#keyword-searching-metadata)
5. [Stand-alone functions for downloading and loading dee2 data](#stand-alone-functions-for-downloading-and-loading-dee2-data)
6. [Aggregating runs data](#aggregating-runs-data)
7. [Running a differential analysis](#running-a-differential-analysis)
8. [Report bugs, issues and suggestions](#report-bugs-issues-and-suggestions)

## Getting started

This tutorial provides a walkthrough for how to work with dee2 expression data,
starting with data searches, obtaining the data from dee2.io and then performing
a differential analysis with DESeq or edgeR.

To use search and obtain dee2 data in R, you will need to source the 'getDEE2' script as follows:

`> source("https://raw.githubusercontent.com/markziemann/dee2/master/getDEE2.R")`

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

```
> mdat[which(mdat$GSE_accession %in% "GSE33569"),]
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

## Fetching dee2 data using SRA run accession numbers

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

> names(x)
[1] "GeneCounts" "TxCounts"   "QcMx"       "absent"   
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

Here I'll demonstrate two ways to query SRAdbv2 ([Docs](https://seandavi.github.io/SRAdbV2/articles/SRAdbv2.html)).

You will need to install, then load the library.

```
#installation
install.packages('BiocManager')
BiocManager::install('seandavi/SRAdbV2')
#load library
library(SRAdbV2)
```

The first method is to download all the transcriptome metadata for a particular species (taxid:6239 is for C. elegans).

```
#new query
> oidx = Omicidx$new()
> query=paste(
>  paste0('sample_taxon_id:', 6239),
>  'AND experiment_library_source:transcriptomic')
> z = oidx$search(q=query,entity='full',size=100L)
> s = z$scroll()
> res = s$collate(limit = 10000L)
 working: ( - ) [================================================] 100% eta:  0s
> head(res)
# A tibble: 6 x 85
  experiment_Insdc experiment_LastMet… experiment_LastUpd… experiment_Publish…
  <lgl>            <dttm>              <dttm>              <dttm>             
1 TRUE             2014-08-26 14:00:12 2016-02-18 21:16:28 2016-02-18 21:16:28
2 TRUE             2015-09-29 20:32:16 2016-06-01 19:55:21 2016-06-01 19:55:21
3 TRUE             2014-11-21 13:38:14 2016-01-07 20:57:05 2016-01-07 20:57:05
4 TRUE             2013-09-23 16:45:14 2013-09-23 16:45:14 2011-02-09 16:35:05
5 TRUE             2018-06-27 19:10:23 2018-06-27 19:21:11 2018-06-27 19:21:11
6 TRUE             2017-09-04 18:59:05 2018-05-11 20:51:14 2018-05-11 20:51:14
# ... with 81 more variables: experiment_Received <dttm>,
#   experiment_Status <chr>, experiment_accession <chr>,
#   experiment_alias <chr>, experiment_attributes <list>,
#   experiment_center_name <chr>, experiment_identifiers <list>,
#   experiment_instrument_model <chr>,
#   experiment_library_construction_protocol <chr>,
#   experiment_library_layout <chr>, experiment_library_selection <chr>,
#   experiment_library_source <chr>, experiment_library_strategy <chr>,
#   experiment_platform <chr>, experiment_title <chr>, experiment_xrefs <list>,
#   run_FileDate <dbl>, run_FileMd5 <chr>, run_FileSize <int>, run_Insdc <lgl>,
#   run_LastMetaUpdate <dttm>, run_LastUpdate <dttm>, run_Published <dttm>,
#   run_Received <dttm>, run_Status <chr>, run_accession <chr>,
#   run_alias <chr>, run_bases <dbl>, run_center_name <chr>,
#   run_identifiers <list>, run_nreads <int>, run_reads <list>,
#   run_spot_length <int>, run_spots <int>, sample_BioSample <chr>,
#   sample_GEO <chr>, sample_Insdc <lgl>, sample_LastMetaUpdate <dttm>,
#   sample_LastUpdate <dttm>, sample_Published <dttm>, sample_Received <dttm>,
#   sample_Status <chr>, sample_accession <chr>, sample_alias <chr>,
#   sample_attributes <list>, sample_identifiers <list>, sample_organism <chr>,
#   sample_taxon_id <int>, sample_title <chr>, sample_xrefs <list>,
#   study_BioProject <chr>, study_GEO <chr>, study_Insdc <lgl>,
#   study_LastMetaUpdate <dttm>, study_LastUpdate <dttm>,
#   study_Published <dttm>, study_Received <dttm>, study_Status <chr>,
#   study_abstract <chr>, study_accession <chr>, study_alias <chr>,
#   study_attributes <list>, study_center_name <chr>, study_identifiers <list>,
#   study_title <chr>, study_type <chr>, study_xrefs <list>,
#   experiment_library_name <chr>, run_attributes <list>,
#   sample_center_name <chr>, experiment_design <chr>,
#   sample_description <chr>, experiment_broker_name <chr>,
#   run_broker_name <chr>, run_center <chr>, sample_broker_name <chr>,
#   study_broker_name <chr>, study_description <chr>, run_file_addons <list>,
#   experiment_library_layout_length <dbl>,
#   experiment_library_layout_sdev <chr>
```
Next, you will want to filter these results for those that have dee2 datasets available. Couple of ways to do this but the best might be to attach a logical value for whether there is data available for the runs.

```
#get the C. elegans metadata if not already obtained
> mdat<-getDee2Metadata("celegans")
#next attach logical value to SRAdb metadata
> res$dee2data<-res$run_accession %in% mdat$SRR_accession
#Check that it has been added
> head(cbind(res$run_accession,res$dee2data))
     [,1]         [,2]   
[1,] "SRR1557807" "TRUE" 
[2,] "SRR2537196" "TRUE" 
[3,] "SRR1658975" "TRUE" 
[4,] "SRR087429"  "TRUE" 
[5,] "SRR7443620" "TRUE" 
[6,] "SRR6002322" "FALSE"
```

You can also see how complete the dee2 coverage is:

```
> length(which(res$dee2data==T))
[1] 8721
> length(which(res$dee2data==F))
[1] 1246
```

You can then browse and make targeted searches of this metadata locally with grepl. In this case I will search for studies with "PAR-CLIP" in the abstract.

```
> res[grep("PAR-CLIP",res$study_abstract),]
# A tibble: 21 x 86
   experiment_Insdc experiment_LastMet… experiment_LastUpd… experiment_Publish…
   <lgl>            <dttm>              <dttm>              <dttm>             
 1 TRUE             2014-06-17 23:15:06 2014-06-17 23:15:29 2011-12-09 15:17:33
 2 TRUE             2014-02-25 18:52:07 2014-11-28 16:04:14 2014-11-21 21:37:23
 3 TRUE             2014-03-28 19:34:16 2014-05-23 20:22:05 2014-05-23 11:06:30
 4 TRUE             2013-09-23 16:45:14 2013-09-23 16:45:14 2011-12-09 15:17:15
 5 TRUE             2013-09-23 16:45:14 2013-09-23 16:45:14 2011-12-09 15:17:15
 6 TRUE             2014-03-28 19:34:16 2014-05-23 20:22:05 2014-05-23 11:06:30
 7 TRUE             2014-03-28 19:34:16 2014-05-23 20:22:05 2014-05-23 11:06:30
 8 TRUE             2014-03-28 19:34:16 2015-07-22 21:07:18 2015-07-22 21:07:18
 9 TRUE             2014-02-25 18:52:41 2014-11-28 16:05:03 2014-11-21 21:38:04
10 TRUE             2014-06-17 23:15:06 2014-06-17 23:15:29 2011-12-09 15:17:33
# ... with 11 more rows, and 82 more variables: experiment_Received <dttm>,
#   experiment_Status <chr>, experiment_accession <chr>,
#   experiment_alias <chr>, experiment_attributes <list>,
#   experiment_center_name <chr>, experiment_identifiers <list>,
#   experiment_instrument_model <chr>,
#   experiment_library_construction_protocol <chr>,
#   experiment_library_layout <chr>, experiment_library_selection <chr>,
#   experiment_library_source <chr>, experiment_library_strategy <chr>,
#   experiment_platform <chr>, experiment_title <chr>, experiment_xrefs <list>,
#   run_FileDate <dbl>, run_FileMd5 <chr>, run_FileSize <int>, run_Insdc <lgl>,
#   run_LastMetaUpdate <dttm>, run_LastUpdate <dttm>, run_Published <dttm>,
#   run_Received <dttm>, run_Status <chr>, run_accession <chr>,
#   run_alias <chr>, run_bases <dbl>, run_center_name <chr>,
#   run_identifiers <list>, run_nreads <int>, run_reads <list>,
#   run_spot_length <int>, run_spots <int>, sample_BioSample <chr>,
#   sample_GEO <chr>, sample_Insdc <lgl>, sample_LastMetaUpdate <dttm>,
#   sample_LastUpdate <dttm>, sample_Published <dttm>, sample_Received <dttm>,
#   sample_Status <chr>, sample_accession <chr>, sample_alias <chr>,
#   sample_attributes <list>, sample_identifiers <list>, sample_organism <chr>,
#   sample_taxon_id <int>, sample_title <chr>, sample_xrefs <list>,
#   study_BioProject <chr>, study_GEO <chr>, study_Insdc <lgl>,
#   study_LastMetaUpdate <dttm>, study_LastUpdate <dttm>,
#   study_Published <dttm>, study_Received <dttm>, study_Status <chr>,
#   study_abstract <chr>, study_accession <chr>, study_alias <chr>,
#   study_attributes <list>, study_center_name <chr>, study_identifiers <list>,
#   study_title <chr>, study_type <chr>, study_xrefs <list>,
#   experiment_library_name <chr>, run_attributes <list>,
#   sample_center_name <chr>, experiment_design <chr>,
#   sample_description <chr>, experiment_broker_name <chr>,
#   run_broker_name <chr>, run_center <chr>, sample_broker_name <chr>,
#   study_broker_name <chr>, study_description <chr>, run_file_addons <list>,
#   experiment_library_layout_length <dbl>,
#   experiment_library_layout_sdev <chr>, dee2data <lgl>
```

Now check that the runs have corresponding dee2 datasets:

```
> res2<-res[grep("PAR-CLIP",res$study_abstract),]
> cbind(res2$run_accession,res2$dee2data)
      [,1]         [,2]  
 [1,] "SRR363980"  "TRUE"
 [2,] "SRR1176664" "TRUE"
 [3,] "SRR1207392" "TRUE"
 [4,] "SRR363797"  "TRUE"
 [5,] "SRR363796"  "TRUE"
 [6,] "SRR1207390" "TRUE"
 [7,] "SRR1207391" "TRUE"
 [8,] "SRR1207395" "TRUE"
 [9,] "SRR1176669" "TRUE"
[10,] "SRR363981"  "TRUE"
[11,] "SRR363798"  "TRUE"
[12,] "SRR1176645" "TRUE"
[13,] "SRR1207389" "TRUE"
[14,] "SRR1207394" "TRUE"
[15,] "SRR1176665" "TRUE"
[16,] "SRR1176644" "TRUE"
[17,] "SRR1176646" "TRUE"
[18,] "SRR1176670" "TRUE"
[19,] "SRR1176666" "TRUE"
[20,] "SRR363799"  "TRUE"
[21,] "SRR1207393" "TRUE"
```

The drawback of this approach is that it could be very slow for species like human and mouse with hundreds of thousands of RNA-seq datasets. A more targeted search can be performed with a keyword for one of the fields. In this case, we'll search for datasets with the term "PAR-CLIP" in the study abstract. This is definitely a faster approach.

```
> oidx = Omicidx$new()
> query=paste(
+   paste0('sample_taxon_id:', 6239),
+   'AND experiment_library_source:transcriptomic',
+   'AND study_abstract:PAR-CLIP')
> z = oidx$search(q=query,entity='full',size=100L)
> 
> s = z$scroll()
> 
> res = s$collate(limit = 10000L)
dim(res)
> dim(res)
[1] 38 75
> head(res)
# A tibble: 6 x 75
  experiment_Insdc experiment_LastMet… experiment_LastUpd… experiment_Publish…
  <lgl>            <dttm>              <dttm>              <dttm>             
1 TRUE             2014-02-25 18:47:09 2014-11-28 16:04:13 2014-11-21 21:37:21
2 TRUE             2014-02-25 18:47:09 2014-11-28 16:04:12 2014-11-21 21:37:21
3 TRUE             2014-02-25 18:47:09 2014-11-28 16:04:13 2014-11-21 21:37:21
4 TRUE             2013-09-23 16:45:14 2013-09-23 16:45:14 2011-12-09 15:17:15
5 TRUE             2013-09-23 16:45:14 2013-09-23 16:45:14 2011-12-09 15:17:15
6 TRUE             2013-09-23 16:45:14 2013-09-23 16:45:14 2011-12-09 15:17:15
# ... with 71 more variables: experiment_Received <dttm>,
#   experiment_Status <chr>, experiment_accession <chr>,
#   experiment_alias <chr>, experiment_attributes <list>,
#   experiment_center_name <chr>, experiment_identifiers <list>,
#   experiment_instrument_model <chr>,
#   experiment_library_construction_protocol <chr>,
#   experiment_library_layout <chr>, experiment_library_selection <chr>,
#   experiment_library_source <chr>, experiment_library_strategy <chr>,
#   experiment_platform <chr>, experiment_title <chr>, experiment_xrefs <list>,
#   run_FileDate <dbl>, run_FileMd5 <chr>, run_Insdc <lgl>,
#   run_LastMetaUpdate <dttm>, run_LastUpdate <dttm>, run_Published <dttm>,
#   run_Received <dttm>, run_Status <chr>, run_accession <chr>,
#   run_alias <chr>, run_bases <dbl>, run_center_name <chr>,
#   run_identifiers <list>, run_nreads <int>, run_reads <list>,
#   run_spot_length <int>, run_spots <int>, sample_BioSample <chr>,
#   sample_GEO <chr>, sample_Insdc <lgl>, sample_LastMetaUpdate <dttm>,
#   sample_LastUpdate <dttm>, sample_Published <dttm>, sample_Received <dttm>,
#   sample_Status <chr>, sample_accession <chr>, sample_alias <chr>,
#   sample_attributes <list>, sample_center_name <chr>,
#   sample_identifiers <list>, sample_organism <chr>, sample_taxon_id <int>,
#   sample_title <chr>, sample_xrefs <list>, study_BioProject <chr>,
#   study_GEO <chr>, study_Insdc <lgl>, study_LastMetaUpdate <dttm>,
#   study_LastUpdate <dttm>, study_Published <dttm>, study_Received <dttm>,
#   study_Status <chr>, study_abstract <chr>, study_accession <chr>,
#   study_alias <chr>, study_attributes <list>, study_center_name <chr>,
#   study_identifiers <list>, study_title <chr>, study_type <chr>,
#   study_xrefs <list>, experiment_library_name <chr>, run_FileSize <int>,
#   run_attributes <list>, run_file_addons <list>
```

Lastly, obtain a list of SRR accessions and obtain the dee2 data

```
> SRRvec=res2$run_accession
> x<-getDEE2("celegans",SRRvec)
trying URL 'http://dee2.io/metadata/celegans_metadata.tsv.cut'
Content type 'text/tab-separated-values' length 549206 bytes (536 KB)
==================================================
downloaded 536 KB

trying URL 'http://dee2.io/cgi-bin/request.sh?org=celegans&x=SRR363980&x=SRR1176664&x=SRR1207392&x=SRR363797&x=SRR363796&x=SRR1207390&x=SRR1207391&x=SRR1207395&x=SRR1176669&x=SRR363981&x=SRR363798&x=SRR1176645&x=SRR1207389&x=SRR1207394&x=SRR1176665&x=SRR1176644&x=SRR1176646&x=SRR1176670&x=SRR1176666&x=SRR363799&x=SRR1207393'
downloaded 1.7 MB
```

## Stand-alone functions for downloading and loading dee2 data

In case you want to download the data once, then reuse it many times, using the standalone scripts can be faster. Data is saved in zip format as follows:

```
> x<-getDEE2("celegans",SRRlist,outfile="DEE_count_data.zip")
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

In case you need to aggregate runs that are technical replicates, these can be aggregated easily in R. [This
SO thread](http://stackoverflow.com/questions/26046776/sum-two-columns-in-r) provides some examples. Here
is one way:

```
#Here's a data frame of two datasets
> head(gcounts)
               SRR1176644 SRR1176645
WBGene00197333          0          0
WBGene00198386          0          0
WBGene00015153         67         43
WBGene00002061        131        121
WBGene00255704         11          4
WBGene00235314          0         12
#make a new column sum of the 
> gcounts$sum<-gcounts$SRR1176644+gcounts$SRR1176645
> head(gcounts)
               SRR1176644 SRR1176645 sum
WBGene00197333          0          0   0
WBGene00198386          0          0   0
WBGene00015153         67         43 110
WBGene00002061        131        121 252
WBGene00255704         11          4  15
WBGene00235314          0         12  12

```
## Running a differential analysis

In this example, I'll demonstrate a complete workflow for a GEO series accession (GSE93236) starting with SRAdb metadata query then differential expression analysis. 

```
> library(SRAdbV2)
> oidx = Omicidx$new()
> query=paste(
+ paste0('sample_taxon_id:', 9606),
+ 'AND experiment_library_source:transcriptomic',
+ 'AND study_GEO:GSE93236')

> z = oidx$search(q=query,entity='full',size=100L)
> s = z$scroll()
> res = s$collate(limit = 100L)
> head(res)
# A tibble: 6 x 71
  experiment_Insdc experiment_LastMet… experiment_LastUpd… experiment_Publish…
  <lgl>            <dttm>              <dttm>              <dttm>             
1 TRUE             2017-01-06 15:40:42 2017-01-10 15:19:34 2017-01-10 15:19:34
2 TRUE             2017-01-06 15:40:43 2017-01-10 15:19:34 2017-01-10 15:19:34
3 TRUE             2017-01-06 15:40:42 2017-01-10 15:19:34 2017-01-10 15:19:34
4 TRUE             2017-01-06 15:40:42 2017-01-10 15:19:34 2017-01-10 15:19:34
5 TRUE             2017-01-06 15:40:43 2017-01-10 15:19:34 2017-01-10 15:19:34
6 TRUE             2017-01-06 15:40:42 2017-01-10 15:19:34 2017-01-10 15:19:34
# ... with 67 more variables: experiment_Received <dttm>,
#   experiment_Status <chr>, experiment_accession <chr>,
#   experiment_alias <chr>, experiment_attributes <list>,
#   experiment_identifiers <list>, experiment_instrument_model <chr>,
#   experiment_library_construction_protocol <chr>,
#   experiment_library_layout <chr>, experiment_library_selection <chr>,
#   experiment_library_source <chr>, experiment_library_strategy <chr>,
#   experiment_platform <chr>, experiment_title <chr>, experiment_xrefs <list>,
#   run_FileDate <dbl>, run_FileMd5 <chr>, run_FileSize <int>, run_Insdc <lgl>,
#   run_LastMetaUpdate <dttm>, run_LastUpdate <dttm>, run_Published <dttm>,
#   run_Received <dttm>, run_Status <chr>, run_accession <chr>,
#   run_alias <chr>, run_bases <int>, run_identifiers <list>, run_nreads <int>,
#   run_reads <list>, run_spot_length <int>, run_spots <int>,
#   sample_BioSample <chr>, sample_GEO <chr>, sample_Insdc <lgl>,
#   sample_LastMetaUpdate <dttm>, sample_LastUpdate <dttm>,
#   sample_Published <dttm>, sample_Received <dttm>, sample_Status <chr>,
#   sample_accession <chr>, sample_alias <chr>, sample_attributes <list>,
#   sample_identifiers <list>, sample_numeric_properties <list>,
#   sample_ontology_terms <list>, sample_organism <chr>, sample_taxon_id <int>,
#   sample_title <chr>, sample_type <chr>, sample_xrefs <list>,
#   study_BioProject <chr>, study_GEO <chr>, study_Insdc <lgl>,
#   study_LastMetaUpdate <dttm>, study_LastUpdate <dttm>,
#   study_Published <dttm>, study_Received <dttm>, study_Status <chr>,
#   study_abstract <chr>, study_accession <chr>, study_alias <chr>,
#   study_center_name <chr>, study_identifiers <list>, study_title <chr>,
#   study_type <chr>, study_xrefs <list>

#now extract sample name info into a samplesheet
> sample_sheet<-as.data.frame(cbind(res$run_accession,res$sample_title))
> sample_sheet
          V1          V2
1 SRR5150592 Set7KD_rep1
2 SRR5150597    NTC_rep3
3 SRR5150594 Set7KD_rep3
4 SRR5150593 Set7KD_rep2
5 SRR5150596    NTC_rep2
6 SRR5150595    NTC_rep1
#tidy up the table
> colnames(sample_sheet)=c("run","sample")
> sample_sheet<-sample_sheet[order(sample_sheet$run),]
> sample_sheet
         run      sample
1 SRR5150592 Set7KD_rep1
4 SRR5150593 Set7KD_rep2
3 SRR5150594 Set7KD_rep3
6 SRR5150595    NTC_rep1
5 SRR5150596    NTC_rep2
2 SRR5150597    NTC_rep3

#identify the knockdown samples
> sample_sheet$knockdown<-as.factor(as.numeric(grepl("KD",sample_sheet$sample)))
> sample_sheet
         run      sample knockdown
1 SRR5150592 Set7KD_rep1         1
4 SRR5150593 Set7KD_rep2         1
3 SRR5150594 Set7KD_rep3         1
6 SRR5150595    NTC_rep1         0
5 SRR5150596    NTC_rep2         0
2 SRR5150597    NTC_rep3         0
> sample_sheet$knockdown
[1] 1 1 1 0 0 0
Levels: 0 1
#get expression data
> x<-getDEE2("hsapiens",sample_sheet$run)
trying URL 'http://dee2.io/metadata/hsapiens_metadata.tsv.cut'
Content type 'text/tab-separated-values' length 21089134 bytes (20.1 MB)
==================================================
downloaded 20.1 MB

trying URL 'http://dee2.io/cgi-bin/request.sh?org=hsapiens&x=SRR5150592&x=SRR5150593&x=SRR5150594&x=SRR5150595&x=SRR5150596&x=SRR5150597'
downloaded 3.2 MB

> head(x$GeneCounts)
                SRR5150592 SRR5150593 SRR5150595 SRR5150596 SRR5150597
ENSG00000223972          1          0          0          0          0
ENSG00000227232          0          1          2          1          1
ENSG00000278267          0          0          0          0          0
ENSG00000243485          0          0          0          1          1
ENSG00000284332          0          0          0          0          0
ENSG00000237613          0          0          0          0          0


```
The next step is to run the differential analysis. First with edgeR.

```
> library("edgeR")
#filter counts for genes with 10 or more reads per sample 
> y<-x$GeneCounts
> y<-y[which(rowSums(y)/ncol(y)>=(10)),]
#set knockdown as the treatment
> design<-model.matrix(~sample_sheet$knockdown)
> rownames(design)<-sample_sheet$run
> y<-DGEList(counts=y)
> y<-calcNormFactors(y)
#design<-design[which(rownames(design) %in% colnames(y$counts)),]
> y <- estimateDisp(y, design,robust=TRUE)
> fit <- glmFit(y, design)
> lrt <- glmLRT(fit)
> dge<-as.data.frame(topTags(lrt,n=Inf))
> dge$dispersion<-lrt$dispersion
> dge<-merge(dge,lrt$fitted.values,by='row.names')
> rownames(dge)=dge$Row.names
> dge$Row.names=NULL
> dge<-merge(dge,y$counts,by='row.names')
> dge<-dge[order(dge$PValue),]
> head(dge)
            Row.names     logFC   logCPM       LR        PValue           FDR
12389 ENSG00000168542  2.729477 5.693558 760.7346 1.859234e-167 1.083971e-162
11419 ENSG00000164692  2.259761 8.893806 376.1168  8.717894e-84  2.541353e-79
13332 ENSG00000172531 -1.713319 6.772587 363.1704  5.744332e-81  1.116354e-76
12656 ENSG00000169715 -1.478648 5.482148 225.2870  6.356442e-51  9.264832e-47
3317  ENSG00000106484  1.205877 7.274166 223.5334  1.533471e-50  1.788089e-46
11902 ENSG00000166595 -1.372937 5.681828 208.4336  3.017219e-47  2.931832e-43
        dispersion SRR5150592.x SRR5150593.x SRR5150595.x SRR5150596.x
12389 2.977294e-02    2135.7324    2096.3321     338.4845     331.8332
11419 9.765625e-05   18301.0957   17963.4752    4017.8772    3938.9253
13332 2.977294e-02     931.2977     914.1170    3211.2547    3148.1530
12656 2.977294e-02     434.0033     425.9968    1271.9650    1246.9707
3317  2.977294e-02    4738.4889    4651.0728    2159.7640    2117.3243
11902 2.977294e-02     528.8961     519.1390    1440.5010    1412.1949
      SRR5150597.x SRR5150592.y SRR5150593.y SRR5150595.y SRR5150596.y
12389     337.7491         2312         1923          337          338
11419    4009.1481        20781        15529         4131         3834
13332    3204.2780          902          943         3352         3085
12656    1269.2016          434          426         1328         1195
3317     2155.0718         4936         4457         2172         2129
11902    1437.3714          526          522         1492         1423
      SRR5150597.y
12389          333
11419         4003
13332         3128
12656         1266
3317          2131
11902         1375
```
Next run the DESeq2 analysis 
```
> library("DESeq2")

> sample_sheet$knockdown<-as.factor(as.numeric(grepl("KD",sample_sheet$sample)))
#> sample_sheet<-sample_sheet[which(sample_sheet$run %in% colnames(x1)),]
> x1<-x$GeneCounts[which(rowSums(x$GeneCounts)/ncol(x$GeneCounts)>=(10)),]
> sample_sheet<-sample_sheet[which(sample_sheet$run %in% colnames(x1)),]

> dds <- DESeqDataSetFromMatrix(countData = x1, colData = sample_sheet, design = ~ knockdown)
> res <- DESeq(dds)
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
> z<- results(res)
> vsd <- vst(dds, blind=FALSE)
#stick on the normalised expression values to the table
> zz<-cbind(z,assay(vsd))
#sort by adjusted p-value
> zz<-as.data.frame(zz[order(zz$padj),])
> head(zz)
                 baseMean log2FoldChange      lfcSE      stat        pvalue
ENSG00000168542 1065.1201       2.717235 0.10765778  25.23956 1.474689e-140
ENSG00000164692 9780.1702       2.247032 0.10414641  21.57570 3.038595e-103
ENSG00000172531 2255.1312      -1.726000 0.08077537 -21.36790 2.657985e-101
ENSG00000106484 3186.7105       1.193305 0.07224973  16.51640  2.796178e-61
ENSG00000169715  919.6678      -1.491496 0.09875523 -15.10296  1.548274e-51
ENSG00000166595 1056.7284      -1.385686 0.09346528 -14.82567  9.997384e-50
                         padj SRR5150592 SRR5150593 SRR5150595 SRR5150596
ENSG00000168542 2.173397e-136  11.497073  11.306165   9.641273   9.656953
ENSG00000164692  2.239141e-99  14.402089  14.023006  12.159873  12.090571
ENSG00000172531  1.305779e-97  10.482920  10.545884  11.897193  11.819312
ENSG00000106484  1.030252e-57  12.442404  12.336635  11.375295  11.375280
ENSG00000169715  4.563692e-48   9.853176   9.854604  10.830405  10.741488
ENSG00000166595  2.455691e-46  10.004341  10.014810  10.954333  10.924912
                SRR5150597
ENSG00000168542   9.636455
ENSG00000164692  12.126039
ENSG00000172531  11.817721
ENSG00000106484  11.358812
ENSG00000169715  10.785519
ENSG00000166595  10.872210
```

## Report bugs, issues and suggestions
Submit at the [github page](https://github.com/markziemann/dee2).
