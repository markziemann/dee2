# DEE2 quality metrics
The DEE2 pipeline gathers a large number of metrics that are critical to determining the validity of 
transcriptome data. Remember that DEE2 is organised by SRA runs. Here, I will describe how these metrics are 
generated and what they mean. Further down, I discuss the logic behind classification of pass, warn and fail
datasets.

## SequenceFormat
Prior to a full analysis, a sample of 4000 reads is obtained from the .sra file to perform some checks. The 
output is checked and the .sra file is classified as single (SE) or paired end (PE).

## QualityEncoding
Fastq files come in a variety of subtle flavours. These vary in the way that quality scores are defined by ASCII
characters. The standard in the field is Sanger/Illumina1.9 but there may be a small number of files with 
different quality encodings.

## Read lengths 
After these initial checks, the entire fastq data are unpacked and FastQC is run. We can see the minimum, median 
and maximum read length of read 1 and read 2 (if paired end). These metrics come under the following headings:
* Read1MinimumLength
* Read1MedianLength
* Read1MaxLength
* Read2MinimumLength
* Read2MedianLength
* Read2MaxLength

## NumReadsTotal
Is the total number of reads or "spots" in the fastq file before filtering.

## NumReadsQcPass
Is the number of reads that pass quality control filtering. This includes 3' quality trimming with Skewer, 
adapter clipping (if an adapter was detected with a frequency of >2.5%), progressive 5' trimming (this is
described in the manuscript) and removal of reads shorted than 18 bp. 

## QcPassRate
The percentage of reads that pass quality control filtering.

## StarMapRateTest
Before performing a full alignment, a sample of 10,000 reads are selected and mapped with STAR. If a dataset is 
paired end, then read 1 and read 2 are mapped independently (PE_Read1_StarMapRateTest and 
PE_Read2_StarMapRateTest). This value ranges from 0 to 100 and is an integer. 

## PE_Read1_Excluded
If the StarMapRateTest number is below 40 for one read, then that read is dropped from the pair. If the value for
PE_Read1_Excluded is TRUE, then read 1 has been excluded. if PE_Read2_Excluded is TRUE then read 2 has been 
excluded.

## MappingFormat
After performing all of the checks, the mapping is performed: SE for single end and PE for paired end. In the
case that SequenceFormat=PE and MappingFormat=SE, one read from the pair has been omitted due to low mapping 
rate.

## STAR_UniqMappedReads
This is the number of uniquely mapped reads as determined by STAR.

## STAR_Strandedness
After mapping, STAR produces three sets of expression counts. These are based on assumptions of the strand 
orientation of the library. The pipeline calculates the strand bias and classifies the strand information if 
there is a bias >5.

## Star read assignment metrics
STAR classifies reads in a number of ways that can be useful to know:
* STAR_UnmappedReads: Number of reads that didn't map to the reference genome
* STAR_MultiMappedReads: Number of reads that mapped to multiple locations in the genome
* STAR_NoFeatureReads: Number of reads that mapped uniquely to a locus without an annotated gene
* STAR_AmbiguousReads: Number of reads that were not assigned because it only partially mapped to an exon
* STAR_AssignedReads: Number of reads that are uniquely mapped to a gene and included in the gene expression 
count data. This information is important because the larger the number of assigned reads, the more accurate the
quantificaiton will be.
* STAR_UniqMapRate: Proportion of reads that were uniquely mapped to the genome. This is a key quality indicator
as a high proportion of unmapped reads could be indicative of a high carryover of rRNA or a sample swap.
* STAR_AssignRate: Proportion of reads that are assigned to genes. This is another important metric as a low value
suggests contamination with gDNA or a high rate of intron mapping. 

## Kallisto_Kmer
This is the k-mer parameter used by Kallisto to map reads to transcripts. Normally this is 31 by default, but 
when the read length is shorter, the kmer used might be smaller.

## Kallisto_MappedReads
This is the number of reads mapped & assigned by Kallisto. This information is important because the larger the 
number of assigned reads, the more accurate the quantificaiton will be. 

## Kallisto_MapRate
This is the proportion of reads that were assigned to transcripts. A low value here could be a sign that there 
are issues with the library preparation.

# QC classifications
Using the QC metrics outlined above, we have classified the datasets as “pass”, “warn” and “fail” according to 
some simple rules as summarised in the table below. Each rule has a numeric code.

| Metric | Fail threshold | Warn threshold | Code |
| ------ | ------ | ------ | ------ | 
| NumReadsQcPass | < 50 reads per gene | < 500 reads per gene | 1 |
| QcPassRate | < 60% | < 80% | 2 |
| STAR_UniqMapRate | < 50% | < 70% | 3 |
| STAR_AssignRate | < 40% | < 60% | 4 |
| STAR_AssignedReads | < 50 reads per gene | < 500 reads per gene | 5 |
| Kallisto_MapRate | < 40% | < 60% | 6 |
| Kallisto_MappedReads | < 50 reads per gene | < 500 reads per gene | 7 |

