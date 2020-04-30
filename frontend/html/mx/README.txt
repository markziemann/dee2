README for DEE2 bulk data

The data are in 'long format' tables with the columns: 'dataset', 'gene',
'count'. Long tables are prefered for loading into databases as compared to
wide matrix format.

STAR counts have the prefix 'se.tsv.bz2'. STAR counts will always be integers.

Kallisto estimated counts have the 'ke.tsv.bz2' prefix and are commonly
floating point numbers.

QC metrics are available in long table format with the columns 'dataset',
'QC metric type', 'QC metric result'. These have the suffix 'qc.tsv.bz2'.

The metadata files corresponding to these bulk data are present at
(http://dee2.io/metadata/). If you are using these data in your studies, please
note the snapshot date in your records.

The checksums.md5 file is available for you to check that the downloads have worked 
properly

