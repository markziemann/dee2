# Incorporating dee2 data into your R-based RNA-seq workflow

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

To use search and obtain dee2 data in R, you can install the DEE2 package as follows:

```
library("devtools")
devtools::install_github("markziemann/dee2/getDEE2")
library("getDEE2")

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

> head(mdat)
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
[1] "GeneCounts" "TxCounts"   "GeneInfo"   "TxInfo"     "QcMx"      
[6] "absent"    
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
> head(x$GeneInfo)
               GeneSymbol mean median longest_isoform merged
WBGene00197333   cTel3X.2  150    150             150    150
WBGene00198386   cTel3X.3  150    150             150    150
WBGene00015153    B0348.5 1051   1178            1178   1178
WBGene00002061      ife-3 1015    949            1107   1107
WBGene00255704   B0348.10  363    363             363    363
WBGene00235314    B0348.9  220    220             220    220
> head(x$TxInfo)
                   GeneID GeneSymbol TxLength
Y110A7A.10 WBGene00000001      aap-1     1787
F27C8.1    WBGene00000002      aat-1     1940
F07C3.7    WBGene00000003      aat-2     1728
F52H2.2a   WBGene00000004      aat-3     1739
F52H2.2b   WBGene00000004      aat-3     1840
T13A10.10a WBGene00000005      aat-4     1734

```

Notice the new objects x$GeneInfo and x$TxInfo. They have information on the gene and transcript lengths
that might be useful for calculating FPKM. Gene symbol information might be useful for downstream analysis. 
The TxInfo dataframe contains the relationships between the transcript and parent genes. There is a function to
aggregate Tx counts to gene level counts called "Tx2Gene", demonstrated below:

```
> x<-Tx2Gene(x)
> names(x)
[1] "Tx2Gene"    "GeneCounts" "TxCounts"   "GeneInfo"   "TxInfo"    
[6] "QcMx"       "absent"    
> head(x$Tx2Gene)
               SRR363796 SRR363797 SRR363798 SRR363799
WBGene00000001        11        23        48        45
WBGene00000002         0         0         0         0
WBGene00000003         0         2         0        21
WBGene00000004         0         4         7         5
WBGene00000005         0         0         0         0
WBGene00000006         2         1         0         2

```

You can directly specify the SRR accessions in the command line, but be sure to type them correctly.
In case SRR accessions are not present in the database, there will be a warning message.

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
                SRR5150592 SRR5150593 SRR5150594 SRR5150595 SRR5150596
ENSG00000223972          1          0          0          0          0
ENSG00000227232          0          1          0          2          1
ENSG00000278267          0          0          1          0          0
ENSG00000243485          0          0          0          0          1
ENSG00000284332          0          0          0          0          0
ENSG00000237613          0          0          0          0          0
                SRR5150597
ENSG00000223972          0
ENSG00000227232          1
ENSG00000278267          0
ENSG00000243485          1
ENSG00000284332          0
ENSG00000237613          0
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
9169 ENSG00000168542  2.719133 5.924633 909.0439 1.061182e-199 1.622017e-195
9720 ENSG00000172531 -1.724929 6.631796 582.9230 8.674970e-129 6.629846e-125
8517 ENSG00000164692  2.280808 9.121734 452.8242 1.751887e-100  8.925865e-97
9333 ENSG00000169715 -1.462740 5.364131 326.6468  5.165698e-73  1.973942e-69
8840 ENSG00000166595 -1.370595 5.567933 313.4334  3.901993e-70  1.192839e-66
2621 ENSG00000106484  1.199930 7.396969 301.3725  1.654854e-67  4.215741e-64
     dispersion SRR5150592.x SRR5150593.x SRR5150594.x SRR5150595.x
9169 0.01268588    2151.7916    2117.7443    2119.1704     338.4094
9720 0.02904578     937.4830     922.6494     923.2707    3210.5625
8517 0.00470866   18844.2254   18546.0577   18558.5468    4016.9643
9333 0.02637610     445.2974     438.2516     438.5467    1271.6812
8840 0.01459492     537.5872     529.0811     529.4374    1440.2005
2621 0.01646871    4788.7056    4712.9351    4716.1088    2159.2906
     SRR5150596.x SRR5150597.x SRR5150592.y SRR5150593.y SRR5150594.y
9169     331.5685     338.0951         2312         1923         2156
9720    3145.6617    3207.5813          902          943          938
8517    3935.7622    4013.2343        20781        15529        19670
9333    1245.9744    1270.5004          434          426          462
8840    1411.0872    1438.8632          526          522          548
2621    2115.6411    2157.2856         4936         4457         4827
     SRR5150595.y SRR5150596.y SRR5150597.y
9169          337          338          333
9720         3352         3085         3128
8517         4131         3834         4003
9333         1328         1195         1266
8840         1492         1423         1375
2621         2172         2129         2131
```
Next run the DESeq2 analysis 
```
> library("DESeq2")

> sample_sheet$knockdown<-as.factor(as.numeric(grepl("KD",sample_sheet$sample)))
> x1<-x$GeneCounts[which(rowSums(x$GeneCounts)/ncol(x$GeneCounts)>=(10)),]
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
ENSG00000168542  1249.2520       2.725565 0.09107000  29.92824 8.447059e-197
ENSG00000164692 11450.2830       2.286642 0.08630846  26.49384 1.141381e-154
ENSG00000172531  2029.4554      -1.718911 0.07115994 -24.15560 6.519887e-129
ENSG00000106484  3462.0694       1.206004 0.06421216  18.78155  1.069084e-78
ENSG00000125148  4424.0310      -1.657291 0.09265727 -17.88624  1.509457e-71
ENSG00000169715   841.1335      -1.456930 0.08785758 -16.58285  9.271849e-62
                         padj SRR5150592 SRR5150593 SRR5150594 SRR5150595
ENSG00000168542 1.266045e-192  11.461774  11.265955  11.394072   9.554022
ENSG00000164692 8.553510e-151  14.389614  14.007636  14.330986  12.135788
ENSG00000172531 3.257336e-125  10.424497  10.488011  10.476727  11.869813
ENSG00000106484  4.005858e-75  12.419366  12.310908  12.407513  11.339600
ENSG00000125148  4.524748e-68  11.216534  11.533624  11.292089  12.849496
ENSG00000169715  2.316108e-58   9.773389   9.773877   9.834018  10.782934
                SRR5150596 SRR5150597
ENSG00000168542   9.570465   9.548311
ENSG00000164692  12.065852  12.100309
ENSG00000172531  11.791029  11.788021
ENSG00000106484  11.339755  11.321661
ENSG00000125148  12.681067  12.913793
ENSG00000169715  10.691897  10.735876
```

## Report bugs, issues and suggestions
Submit at the [github page](https://github.com/markziemann/dee2).
