# Under construction
Note that this app is still under development and alpha testing. It may not work and may undergo
major changes in future.

# DEE2
DEE2 is an intelligent RNA-seq pipeline for distributed computing. It uses docker containers to
enhance ease of use, portability and reproducibility. DEE2 currently supports analysis of several 
major species including A. thaliana , C. elegans, D. melanogaster, D. rerio, E. coli, H. sapiens,
M. musculus, R. norvegicus and S. cerevisiae. DEE2 will download data from SRA and process it, 
providing tabulated data that can be used in downstream statistical analysis.

Here are some of the key features:
 * Intelligent adapter detection and clipping
 * Strandedness detection
 * Parallel assignment of reads genes and transcripts with STAR and Kallisto
 * Thorough quality control logs
 * Open source pipeline
 * Distributed approach

## How it works
Using docker allows reproducible analysis across different environments, which allow work to be 
distributed across many machines. The user can provide a species and SRA run accession numbers 
to be processed, and this data will be immediately available to the user when pipline is 
completed. Data is also uploaded to our server by sftp and after sanitation, will be available to
other users via the frontend (still under construction). If users don't provide accession numbers,
then they will receive accessions from the current queue via html request. If the user doesn't
specify a species, then one is selected based on the memory available. 

## Contributions welcome
We welcome contributions to code development as well as machine time. If you have idle compute 
resources and bandwith, consider sheduling DEE2 over the weekend with a cron job.

## Quick start guide
Install docker

`sudo apt install docker.io`

Add user to docker group

`sudo usermod -aG docker $(whoami)`

Log out and log back in or use the below command

`exec su -l $(whoami)`

Now pull image

`docker pull mziemann/tallyup`

Run the pipeline by specifying organism and SRA run accessions

`docker run -it mziemann/tallyup /root/code/volunteer_pipeline.sh ecoli SRR2637695,SRR2637696,SRR2637697,SRR2637698`

return the data directory from the container to the host filesystem

`docker cp $(docker ps -alq):/root/data/ .`


## Donating compute time
If you have a species of interest

`docker run -it mziemann/tallyup /root/code/volunteer_pipeline.sh celegans`

Let the app choose the species and accessions - it will be selected based on available memory

`docker run -it mziemann/tallyup /root/code/volunteer_pipeline.sh`


## Troubleshooting
If there is insufficient space in /var to work, do follow instructions from the
[link](https://stackoverflow.com/a/34731550).
