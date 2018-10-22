# DEE2
The aim of DEE2 is to make all RNA-seq data freely available to everyone. DEE2 consists of three parts:
* Webserver where end-users can search for and obtain data-sets of interest 
* Back-end that collects, filters and organises data provided by contributing worker nodes
* Pipeline that can download and process SRA data as well as users' own fastq files.

DEE2 currently supports analysis of several major species including A. thaliana , C. elegans, D. melanogaster, D. rerio, E. coli, H. sapiens, M. musculus, R. norvegicus and S. cerevisiae. The DEE2 pipeline downloads data from SRA and processes it, providing tabulated data that can be used in downstream statistical analysis.

## How can I access the processed data?
The processed data is available at http://dee2.io but keep in mind that processing is not yet complete. 

Dee2 data can be accessed using our specially developed [R interface](../master/AccessDEEfromR.md)

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

# Want to learn more?
## How it works
The DEE2 pipeline uses containers (Docker and Singularity) to enhance ease of use, portability and reproducibility. This means the DEE2 pipeline can be run reproducibly across different environments, making it amenable to distributed computing. The user can provide a species and SRA run accession numbers to be processed, and this data will be immediately available to the user when completed. Data is also uploaded to our server by sftp and after sanitation, will be available to other users via the frontend (still under construction). If users don't provide accession numbers, then they will receive accessions from the current queue via html request. If the user doesn't specify a species, then one is selected based on the memory available. 

## Before starting check system requirements
So far this has only been tested on linux but will probably work for any system with sufficient resources capable of running Docker.

Minimum 8GB RAM available. Use the "free" command:

`free -h`

Minimum 32GB of data storage available for Docker. Use the "df" command to check the default Docker data location in /var.

`df -h /var`

If you need to change the default Docker data to another location with >32GB free, follow these steps (works for Ubuntu 16.04): 

1) Create a config file for docker data storage

`sudo nano /etc/systemd/system/docker.service.d/docker-storage.conf`

2) Paste in the following to the new file, substituting "/path/to/data/" with the location you'd like to run the image:

```
[Service]
ExecStart=
ExecStart=/usr/bin/dockerd -H fd:// --data-root="/path/to/data/"
```
3) Restart the docker service

`sudo systemctl daemon-reload ; sudo systemctl restart docker`

4) Confirm that the change in location has occurred.

`docker info|grep "Docker Root Dir"`

## Quick start guide for docker users (Singularity users see this guide [link](https://github.com/markziemann/dee2/blob/master/notes/singularity_walkthrough.md))
Install docker

`sudo apt install docker.io`

Add user to docker group

`sudo usermod -aG docker $(whoami)`

Log out and log back in or use the below command

`exec su -l $(whoami)`

Now pull image

`docker pull mziemann/tallyup`

Run the pipeline by specifying organism and SRA run accessions

`docker run mziemann/tallyup ecoli SRR2637695,SRR2637696,SRR2637697,SRR2637698`

You can also process your own fastq files in the current working directory. Gzip and bzip2 compressed files are OK. Data remain private and are not transmitted.

`docker run -v $(pwd):/dee2/mnt mziemann/tallyup hsapiens -f sample1_R1.fq.gz,sample2_R1.fq sample1_R2.fq.gz,sample2_R2.fq`

Return the data directory from the container to the host filesystem

`docker cp $(docker ps -alq):/root/data/ .`

## Donating compute time
If you have a species of interest

`docker run mziemann/tallyup celegans`

Let the app choose the species and accessions - it will be selected based on available memory

`docker run mziemann/tallyup`

If you want to keep it running even if ssh gets terminated, use 'nohup'

`nohup docker run mziemann/tallyup & tailf nohup.out`

If you want to keep one container running contantly, try the 'run-one' utility

`nohup run-one-constantly docker run mziemann/tallyup mmusculus &`

## Contributions welcome
We welcome feedback, bug reports, and contributions to code development. Compute time is also very welcome.

## Troubleshooting
If you're encountering problems, check the system requirements again. If still stuck, contact me on mark.ziemann{at}gmail.com
