# use HDF5 format to store bulk data.
# this saves a lot of space. for ecoli se counts 242M down to 85M
# and ke 238M to 111M. QC matrix is slightly bigger.
# it also allows quick random access on our server

library(rhdf5)

org="ecoli"

## Start with the SE TSV data
h5file <- paste(org,"_se.h5",sep="")

if (file.exists(h5file)) { file.remove(h5file) }
h5createFile(h5file)

input_dir <- paste("/mnt/hdd1/dee2/data/",org,sep="")
file_list <- list.files(input_dir, pattern = "se.tsv.gz$", recursive=TRUE, full.names = TRUE)
batch_size <- 1000

num_files <- length(file_list)
rows_per_file <- nrow(read.table(file_list[1],row.names=1,header=TRUE))
num_cols <- 1

h5createDataset(
  file = h5file,
  dataset = "bigmatrix",
  dims = c(num_files, rows_per_file),
  storage.mode = "integer",
  chunk = c(10, rows_per_file),    # 10,000-row chunks
  level = 9                      # compression level (gzip)
)

for (i in seq(1, length(file_list), by = batch_size)) {
  batch_files <- file_list[i:min(i + batch_size - 1, length(file_list))]
  data_list <- lapply(batch_files, function(f) { t(as.matrix(read.table(f,header=TRUE,row.names=1)[,1])) } )
  batch_data <- do.call(rbind, data_list)
  row_start <- (i - 1) + 1
  row_end <- row_start + nrow(batch_data) - 1
  h5write(batch_data, h5file, "bigmatrix", index = list(row_start:row_end, 1:rows_per_file ))
  cat(sprintf("Processed files %d to %d\n", i, i + length(batch_files) - 1))
}

myrownames <- sapply(strsplit(basename(file_list),"\\."),"[[",1) 
h5write(myrownames, h5file,"rownames")

mycolnames <- rownames( read.table(file_list[1]) )
h5write(mycolnames, h5file, "colnames")

H5close()

# NOW FOR THE KE TSV DATASET
h5file <- paste(org,"_ke.h5",sep="")

if (file.exists(h5file)) { file.remove(h5file) }
h5createFile(h5file)

input_dir <- paste("/mnt/hdd1/dee2/data/",org,sep="")
file_list <- list.files(input_dir, pattern = "ke.tsv.gz$", recursive=TRUE, full.names = TRUE)
batch_size <- 1000

num_files <- length(file_list)
rows_per_file <- nrow(read.table(file_list[1],row.names=1,header=TRUE))
num_cols <- 1

h5createDataset(
  file = h5file,
  dataset = "bigmatrix",
  dims = c(num_files, rows_per_file),
  storage.mode = "double",
  chunk = c(10, rows_per_file),    # 10,000-row chunks
  level = 9                      # compression level (gzip)
)

for (i in seq(1, length(file_list), by = batch_size)) {
  batch_files <- file_list[i:min(i + batch_size - 1, length(file_list))]
  data_list <- lapply(batch_files, function(f) { t(as.matrix(read.table(f,header=TRUE,row.names=1)[,3]) ) } )
  batch_data <- do.call(rbind, data_list)
  row_start <- (i - 1) + 1
  row_end <- row_start + nrow(batch_data) - 1
  h5write(batch_data, h5file, "bigmatrix", index = list(row_start:row_end, 1:rows_per_file ))
  cat(sprintf("Processed files %d to %d\n", i, i + length(batch_files) - 1))
}

myrownames <- sapply(strsplit(basename(file_list),"\\."),"[[",1) 
h5write(myrownames, h5file, "rownames")

mycolnames <- rownames( read.table(file_list[1]) )
h5write(mycolnames, h5file, "colnames")

H5close()

## Next is the QC data
h5file <- paste(org,"_qc.h5",sep="")

if (file.exists(h5file)) { file.remove(h5file) }
h5createFile(h5file)

input_dir <- paste("/mnt/hdd1/dee2/data/",org,sep="")
file_list <- list.files(input_dir, pattern = "\\.qc$", recursive=TRUE, full.names = TRUE)
batch_size <- 1000

qc <- read.csv(file_list[1],row.names=1,header=FALSE,sep=":")

num_files <- length(file_list)
rows_per_file <- nrow(read.csv(file_list[1],row.names=1,header=FALSE,sep=":"))
num_cols <- 1

h5createDataset(
  file = h5file,
  dataset = "bigmatrix",
  dims = c(num_files, rows_per_file),
  storage.mode = "character",
  chunk = c(100, rows_per_file),
  level = 9
)

for (i in seq(1, length(file_list), by = batch_size)) {
  batch_files <- file_list[i:min(i + batch_size - 1, length(file_list))]
  data_list <- lapply(batch_files, function(f) { t(as.matrix( read.csv(file_list[1],row.names=1,header=FALSE,sep=":") ) ) } )
  batch_data <- do.call(rbind, data_list)
  row_start <- (i - 1) + 1
  row_end <- row_start + nrow(batch_data) - 1
  h5write(batch_data, h5file, "bigmatrix", index = list(row_start:row_end, 1:rows_per_file ))
  cat(sprintf("Processed files %d to %d\n", i, i + length(batch_files) - 1))
}

myrownames <- sapply(strsplit(basename(file_list),"\\."),"[[",1)
h5write(myrownames, h5file, "rownames")

mycolnames <- rownames( read.table(file_list[1]) )
h5write(mycolnames, h5file, "colnames")

H5close()

#test
h5ls("combined_data.h5")


y <- h5dump("combined_data.h5")

yy <- h5read(file = "ecoli_qc.h5", 
                 name = "/bigmatrix")
