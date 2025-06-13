# DEE2 pipeline 

The DEE2 pipeline currently supports analysis of several major species including A. thaliana,
C. elegans, D. melanogaster, D. rerio, E. coli, H. sapiens, M. musculus, R. norvegicus and
S. cerevisiae.
We recently added support for B. distachyon, G. max, H. vulgare, O. sativa, P. trichocarpa,
S. bicolor, S. lycopersicum, S. tuberosum, T. aestivum, V. vinifera and Z. mays.

The DEE2 pipeline downloads data from SRA and processes it, providing tabulated data that
can be used in downstream statistical analysis.
You can also process your own fastq files.

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
The DEE2 pipeline uses containers (Docker and Apptainer) to enhance ease of use, portability
and reproducibility.
This means the DEE2 pipeline can be run reproducibly across different environments, making
it amenable to distributed computing.
The user can provide a species and SRA run accession numbers to be processed, and this data
will be immediately available to the user when completed.
If users don't provide accession numbers, then they will receive accessions from the current
queue via html request.
By default, completed datasets are not uploaded to the dee2.io server, but we welcome
contributions from the community.

## Before starting check system requirements

So far this has only been tested on linux but will probably work for any system with
sufficient resources capable of running Docker.

Minimum 8GB RAM available. Use the "free" command:

`free -h`

Minimum 32GB of data storage available for Docker.
Use the "df" command to check the default Docker data location in /var.

`df -h /var`

## Quick start guide for docker users

If Docker is not installed, visit [Docker.com](https://www.docker.com/get-started) and
follow the instructions.

Then pull the DEE2 pipeline image, which we called "tallyup".
The image tagged "dev" is the currently supported one.

```
docker pull mziemann/tallyup:dev
```

Familiarise yourself with the help page.

```
docker run -v $(pwd):/dee2/mnt mziemann/tallyup:dev --help
```

The most common way of using the pipeline is by specifying organism and SRA run accessions.

```
docker run mziemann/tallyup:dev -s ecoli -a SRR19643580,SRR19643581
```

This approach is okay for a small number of runs, but isn't the most efficient, because
the CPU is idle while downloading data from SRA.

To increase the throughput and efficiency, I recommend runnning another script in parallel
which downloads the SRA files, so the CPUs can be kept busy processing them.
To do this, create a text file called SRR.txt which contains all the SRR identifiers of all
runs you want to process.
Then you can make a script like the one below called `downloader.sh` and execute it

```
#!/bin/bash
ORG=ecoli
for SRR n $(cat SRR.txt) ; do
  prefetch -X 9999999999999 -o ${ORG}_${SRR}.sra $SRR
done
```

After the first SRA file has been successfully downloaded, you can run the mapping pipeline
like this:

```
docker run -v $(pwd):/dee2/mnt mziemann/tallyup:dev -s ecoli -t 8 -d -v
```

The `-d` option tells the pipeline to search for SRA files in the current directory.
It will continue processing the SRA files until they're all completed.

You can also process your own fastq files in the current working directory.
Gzip and bzip2 compressed files are OK.
Data remain private and are not transmitted.

```
docker run -v $(pwd):/dee2/mnt mziemann/tallyup:dev -s ecoli -f1 SRR19643580_pass_1.fastq.gz -f2 SRR19643580_pass_2.fastq.gz

```

Once the job is done, return the data directory from the container to the host filesystem

`docker cp $(docker ps -alq):/dee2/data/ .`

In the above script, I have only demonstrated processing of E. coli data, but it works for
the other supported species too.
For those other species, the genome needs to be downloaded and indexed for each time the
pipeline is launched, which can be inefficient in you're going to be using it regularly.

To keep the indexed genome, I suggest starting with processing just one SRA file.
When it is finished, use `docker ps -a` to list the completed containers.
Find the one that corresponds to the job of interest, where the genome was indexed and copy
the "CONTAINER ID" which will be a hexadecimal code like "532defb08c8a".
Now use the `docker commit` command to save that docker image with a tag that identifies
the species, so you can use it later without having to regenerate the index.


```
docker commit 532defb08c8a mziemann/tallyup:hsa
```

Then when you want to run the pipeline again, you can skip the genome regeneration step.

```
docker run -v $(pwd):/dee2/mnt mziemann/tallyup:hsa -s hsapiens -t 8 -d -v
```

## Running into problems?

Raise an issue and we will try to resolve it.
