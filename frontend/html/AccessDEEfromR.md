# Accessing Digital Expression Explorer from R
(March 2018)

To facilitate convenient access to these data from within the R environment, a 
function called 'getDEE2' has been written that interfaces with the web server 
and loads the matrix into R. The function can be 'sourced' as follows:

`> source("https://raw.githubusercontent.com/markziemann/dee2/master/frontend/html/getDEE2.R")`

Function syntax is as follows:

`> getDEE2(organism,SRRlist,outfile="NULL")`

The function will download the data as a zip archive that contains a (1) a gene-wise
count expression matrix, (2) a transcript-wise count expression matrix, (3) a matrix
of quality metrics and (4) a folder of run logs detailing the processing of the 
data including base quality scores, alignment rates, etc. If 'outfile' is 
defined, then files will be downloaded to the current working directory. If it
is not defined, then the files are downloaded to the temporary directory of R
and deleted immediately after use.

The options for organism currently are:

athaliana
celegans
dmelanogaster
drerio
ecoli
hsapiens
mmusculus
rnorvegicus
scerevisiae

The SRR numbers need to exactly match those in SRA.

Here is an example of usage for a small number of E. coli datasets that are 
loaded into a 'mytable' and the zip archive is saved in the current working
directory as "DEE_count_data.zip".

`> x<-getDEE2("celegans",c("SRR029423","SRR029424"))
trying URL 'http://dee2.io/cgi-bin/request.sh?'
Content type 'text/html; charset=utf-8' length 488 bytes
==================================================
downloaded 488 bytes

trying URL 'http://118.138.240.228/cgi-bin/request.sh?org=celegans&x=SRR029423&x=SRR029424'
downloaded 342 KB
`

If all works, then "x" is a list containing 3 data frames (GeneCounts,TxCounts,QcMx).

`> head(x$GeneCounts)
               SRR029423 SRR029424
WBGene00197333         0         0
WBGene00198386         0         0
WBGene00015153         0         2
WBGene00002061        14        30
WBGene00255704         1         0
WBGene00235314         0         1

> head(x$TxCounts)
           SRR029423.ke.tsv SRR029424.ke.tsv
Y110A7A.10               15                7
F27C8.1                   1                0
F07C3.7                   1                1
F52H2.2a                  0                8
F52H2.2b                  0                0
T13A10.10a                0                0

> head(x$QcMx)
                            SRR029423          SRR029424
SequenceFormat                     SE                 SE
QualityEncoding    Sanger/Illumina1.9 Sanger/Illumina1.9
Read1MinimumLength                 36                 36
Read1MedianLength                  36                 36
Read1MaxLength                     36                 36
Read2MinimumLength               NULL               NULL
`

Another example this time starting with a list of SRR accession numbers in a 
plain text file and not saving the zip archive. This approach works well with 
SRAdb integration.

`> SRRlist<-read.table("SRRlist.txt")
> SRRlist
          V1
1 SRR3581689
2 SRR3581692
3 SRR3581858

> x2<-getDEE2("athaliana",SRRlist$V1)
trying URL 'http://dee2.io/cgi-bin/request.sh?'
Content type 'text/html; charset=utf-8' length 488 bytes
==================================================
downloaded 488 bytes

trying URL 'http://118.138.240.228/cgi-bin/request.sh?org=athaliana&x=SRR3581689&x=SRR3581692&x=SRR3581658'
downloaded 593 KB

> head(x2$GeneCounts)
          SRR3581689 SRR3581692
AT1G01010        129         17
AT1G01020        193        206
AT1G03987          0          0
AT1G01030        230        205
AT1G01040       1239        783
AT1G03993          0          0

> head(x2$TxCounts)
            SRR3581689.ke.tsv SRR3581692.ke.tsv
AT1G19090.1            0.0000             0.000
AT1G18320.1           14.3433            41.557
AT5G11100.1           52.0000            37.000
AT4G35335.1          466.0000           607.000
AT1G60930.1          173.0000           103.000
AT1G17000.1            0.0000             0.000

> head(x2$QcMx)
                           SRR3581689         SRR3581692
SequenceFormat                     SE                 SE
QualityEncoding    Sanger/Illumina1.9 Sanger/Illumina1.9
Read1MinimumLength                 50                 50
Read1MedianLength                  50                 50
Read1MaxLength                     50                 50
Read2MinimumLength                 50                 50
`

For convenience, we provide stand alone functions to load each of the data tables 
from previously downloaded zip files.

`> myGeneCounts<-loadGeneCounts("DEE_count_data.zip")
> myTxCounts<-loadTxCounts("DEE_count_data.zip")
> myQcMx<-loadQcMx("DEE_count_data.zip")
> head(mytable3)
             GeneName    SRR922260    SRR922261
ECDH10B_0001     thrL          470          452
ECDH10B_0002     thrA         6046         8138
ECDH10B_0003     thrB         1758         2918
ECDH10B_0004     thrC         3082         4188
ECDH10B_0005     yaaX           95          268
ECDH10B_0006     yaaA          169          734
`

In case you need to aggregate datasets that are technical replicates, read this
thread (http://stackoverflow.com/questions/26046776/sum-two-columns-in-r). Here
is an example:

`> head(mytable)
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
`

## Report bugs, issues and suggestions at the github page
https://github.com/markziemann/dee2
