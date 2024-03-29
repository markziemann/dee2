# DEE2 pipeline 
The DEE2 pipeline currently supports analysis of several major species including A. thaliana , C. elegans, D. melanogaster, D. rerio, E. coli, H. sapiens, M. musculus, R. norvegicus and S. cerevisiae. The DEE2 pipeline downloads data from SRA and processes it, providing tabulated data that can be used in downstream statistical analysis. You can also process your own fastq files.

## Pipeline key features
Here are some of the key features:
 * Intelligent adapter detection and clipping
 * Clipping of non-reference 5' bases (eg UMIs)
 * Strandedness detection
 * Parallel assignment of reads to genes and transcripts with STAR and Kallisto
 * Thorough quality control logs
 * Open source pipeline
 * Distributed approach using containers [link](https://hub.docker.com/r/mziemann/tallyup/)
 * Ability to process own fastq files as well as those from SRA

## How it works
The DEE2 pipeline uses containers (Docker and Singularity) to enhance ease of use, portability and reproducibility.
This means the DEE2 pipeline can be run reproducibly across different environments, making it amenable to distributed computing.
The user can provide a species and SRA run accession numbers to be processed, and this data will be immediately available to the user when completed.
Data is also uploaded to our server by sftp and after sanitation, will be available to other users via the frontend (still under construction).
If users don't provide accession numbers, then they will receive accessions from the current queue via html request.

## Before starting check system requirements

So far this has only been tested on linux but will probably work for any system with sufficient resources capable of running Docker.

Minimum 8GB RAM available. Use the "free" command:

`free -h`

Minimum 32GB of data storage available for Docker. Use the "df" command to check the default Docker data location in /var.

`df -h /var`

## Quick start guide for docker users

If Docker is not intalled, visit [Docker.com](https://www.docker.com/get-started) and follow the instructions.

Now pull image

`docker pull mziemann/tallyup`

Run the pipeline by specifying organism and SRA run accessions

`docker run mziemann/tallyup ecoli SRR2637695,SRR2637696,SRR2637697,SRR2637698`

You can also process your own fastq files in the current working directory. Gzip and bzip2 compressed files are OK. Data remain private and are not transmitted.

`docker run -v $(pwd):/dee2/mnt mziemann/tallyup hsapiens -f sample1_R1.fq.gz,sample2_R1.fq sample1_R2.fq.gz,sample2_R2.fq`

Return the data directory from the container to the host filesystem

`docker cp $(docker ps -alq):/dee2/data/ .`

Please note that SFTP upload is now inoperative, so these processed data cannot be uploaded to DEE2.

## Running on HPC without download ability

The recommended application is to use the downloader.sh script to download a bunch of SRA files in 
advance and keep these in the current working directory.
downloader.sh fetches SRR accession numbers from the queue, so if you want to analyse some different
datasets from SRA, you can download them with SRA toolkit `prefetch`.

Use downloader.sh like this:

`bash downloader.sh hsapiens`

Then you can run the docker image like this:

`docker run -v $(pwd):/dee2/mnt mziemann/tallyup hsapiens -d`

I have included an example script called `run_dee2_hsa1.sh` which I use on my own workstation.

## Donating compute time

Please note that the SFTP upload feature is now inoperative as it is a security weakness.

The downloader.sh script enables the data to be uploaded to dee2, but you will need to request a
unique ssh key first.

## Running into problems?

Raise an issue and we will try to resolve it.
