README for DEE2 bulk data

The data are in HDF5 format.
HDF5 is ideal for storing and and accessing large and complex datasets.
Within each of the HDF5 files, there will be three objects.
"Bigmatrix" is the main data object and contains the gene expression data.
In this type of data, the rows represent SRA runs, and the columns represent genes.
The other two objects are the column names and row names.
Changing format to HDF5 provides a few key benefits.
Firstly, the old format was computationally very expensive to convert from long to wide format
for downstream analysis.
HDF5 allows for random access without loading the entire object into memory.
This will be a great benefit for end-users who want to analyse subsets of the bulk data.
Also, HDF5 occupies less disk space.
For kallisto counts with the "ke" suffix, the data are floating point numbers, and the
saving is ~25%, but for STAR counts with suffix "se" the ddadta are integers, so the disk
storage saving can be up to 75%..
This change will allow us to continue adding samples to the database without the need to
purchase additional cloud storage.

STAR counts have the suffix 'se.h5'. STAR counts will always be integers.

Kallisto estimated counts have the 'ke.h5' suffix and are commonly
floating point numbers.

QC metrics are available in the file with the 'qc.h5' suffix.

The metadata files corresponding to these bulk data are present at
(http://dee2.io/metadata/). If you are using these data in your studies, please
download and retain a copy for your records and report the snapshot date in your
publications.


The checksums.md5 file is available for you to check that the downloads have worked 
properly

For guidance working with HDF5 data, I suggest using the <a href=https://bioconductor.org/packages/release/bioc/html/rhdf5.html>rhdf5</a>
R package.
