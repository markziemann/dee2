#!/bin/bash
#use this crontab
#0 1 * * * bash /home/username/dee_progress.sh


cd /scratch/username/dee2/data/
DATE=$(date '+%Y-%m-%d')
for i in * ; do
  CNT=$(ls $i | wc -l)
  echo $i $DATE $CNT
done >> /home/username/dee_progress.txt
cd ~
Rscript dee_progress.R
scp -i ~/.ssh/cloud/cloud2.key dee_progress.png dee_queue.png ubuntu@118.138.240.228:/mnt/dee2_data/mx
