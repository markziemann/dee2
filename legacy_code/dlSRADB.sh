#!/bin/bash
cd /data/projects/mziemann/geo2mx_project/v1/metadata/tmp
rm SRAmetadb*
#Rscript dlSRADB.R
wget http://gbnci.abcc.ncifcrf.gov/backup/SRAmetadb.sqlite.gz
gunzip -f SRAmetadb.sqlite.gz
mv SRAmetadb.sqlite ..
cd -
