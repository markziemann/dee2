# Guide to running DEE2 pipeline using apptainer

Apptainer is useful because it does not require root privileges.
Apptainer is the name for the new and improved container software previously called singularity.

This guide will help you set up the pipeline (called tallyup).
You will need to install apptainer yourself - refer to https://apptainer.org/documentation/

First step is to pull the current dev image from Dockerhub

```
apptainer build --sandbox tallyup/ docker://mziemann/tallyup:dev
```

You will now have a new folder called tallyup that contains the image data.

The normal way of running the pipeline is by downloading data separately.
Here is an example, using a small E. coli dataset.

```
prefetch -X 9999999999999 -o ecoli_SRR1755067.sra SRR1755067
```

Now try out the apptainer based pipeline with this code.

```
apptainer run -w -B ${PWD}:/dee2/mnt/ tallyup -s ecoli -t 10 -d
```

The options used here will make the image writable (-w), and bind the current
working directory to a place where the pipeline will find the .sra files (-B).
The species (-s) is ecoli, if you select a different species then ecoli,
the pipeline will download the genome information and index it which could take
an hour or more, before processing the RNA-seq data.
The pipeline will use 10 parallel CPU threads.
This can be adjusted, but increasing beyond 12-16 does not yield linear speedup.
The -d option simply instructs tallyup to look for *.sra files in /dee2/mnt/
folder.

At the end of the run, the *.sra files will be deleted and the expression data will
be found in zip archives with a name structure like this: SRR1755067.ecoli.zip.

Finally, keep in mind that as we selected the image to be writable, a lot of saved
data will be collected in the folder.
This will cause the folder to increase in size over time, so it is recommended to
flush out the /dee2/data/ folders to prevent it from becoming too oversized.
