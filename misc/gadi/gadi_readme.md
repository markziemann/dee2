# Running DEE2 pipeline on GADI

This is quite a difficult process so I would not recommend it unless you want to process 100k
datasets and have no other facilities.

GADI has a fairly locked down system and to my knowledge it is not possible to run containers
due to user namespace restrictions.
I have tried udocker, singularity and apptainer without luck.
Therefore the general approach is to uncontainerise the pipeline by providing all the
software dependancies in a "sw" folder.
The results of a few runs should be checked/compared to docker results to be confident it
is working OK.

GADI nodes also don't have internet access, so the data needs to be sent to GADI data mover
nodes from the DEE2 workstation computer.

GADI has restrictions on data storage, even with /scratch so we aim to be as lean as
possible, by keeping very few sra files at any time and regularly scrubbing the dee/data
folder.

In the ~/dee2/nocontainer folder I have working folders mmu1, mmu2, etc which continuously
run a job each.
At the moment I'm not sure how to run multiple parallel jobs per folder.
Will work on this later.

In each working folder there is the workflow folder and a folder called tfr which receives
the SRA files.

The DEE2 workstation computer runs the "gadi_dl_tfr.sh" script which sends SRA files to the
tfr directory, and also counts how many files are there.
If there are 20 sra files present, then it submits a job to the PBS queue (via ssh) to
execute the "dee2_pbs.sh" script.

Multiple instances of "gadi_dl_tfr.sh" can be run, each from different folders to provide
SRA files to multiple GADI working folders.

The "dee2_pbs.sh" script just makes sure that the there are some SRA files available before
running the workflow.
It also checks whether there is another job already running, so that it doesn't cause an
avalanche of submitted jobs.

The "volunteer_pipeline.sh" actually does the data processing.
It is a derivation of the main pipeline version which has been adapted to working outside of
a container.


