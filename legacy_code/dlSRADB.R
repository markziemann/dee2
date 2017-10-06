library(SRAdb)
sqlfile <- 'SRAmetadb.sqlite'
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()
