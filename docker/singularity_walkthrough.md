## Running DEE2 via singularity

Download the image from the dee2 webpage [here (1.6G)](http://dee2.io/images/sing_img.tar.bz2).
Unpack the image with tar:

```
tar xf sing_img.tar.bz2
```

Change into the `img` directory and you will see an img file and tmp folder, these must remain in the same working directory.

Next, it is necessary to alter the configuration of the SRA toolkit so the downloading works properly.
To do this we need to access a shell within the writable container.

```
singularity shell -w -B $(pwd)/tmp:/tmp mziemann_tallyup.img
cd /dee2/sw/sratoolkit.2.8.2-1-ubuntu64/bin/
./vdb-config -i
```

You will be in the configuration tool now. 
Set default import path to `/dee2/ncbi/public` then save and exit.
Then exit from the shell.

```
exit
```

Next run this command which will run an internal test to make sure the pipeline is working on your hardware and establish sftp connection to the dee2 webserver.

```
singularity run -w -B $(pwd)/tmp:/tmp mziemann_tallyup.img
```

If you get no errors, you will be safe now to run the image.
You're now free to run the image.
When analysing another organism that isn't ecoli, be sure to use the `-w` option which will save the index for future use. 
If you're analysing your own fastq files, the `-w` flag is required to retrieve the processed files.
If you are processing pubic RNA-seq for the purpose of sending to dee2.io, the `-w` option is not required.

Here are three examples to try in a normal shell.

```
singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup.img ecoli SRR057751
singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup.img celegans
```

Here are two other examples for long-running instances.

```
nohup run-one-constantly singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup.img hsapiens
nohup bash -c 'while true ; do singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup.img ; sleep 10 ; done'
```

Now try to schedule jobs. Below is an example sbatch file that works for me.
Modify memory requirements and target organism as needed.

```
#!/bin/bash
# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=dee2

# To set a project account for credit charging, 

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1

# Memory usage (MB)
# SBATCH --mem-per-cpu=1000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-10:00:00

#SBATCH --partition=m3a
#SBATCH --reservation=epigenetics

# Set the file for output (stdout)
#SBATCH --output=MyJob-%j.out

# Set the file for error log (stderr)
#SBATCH --error=MyJob-%j.err


module load singularity/2.3.1
cd /path/to/working/image/
>sbatch.log
nohup bash -c 'while true ; do singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup.img >sbatch.log 2>&1 ; done'
```

If this approach doesn't work for you because your HPC doesn't authorise singularity, I have recently had a better experience with [udocker](https://github.com/indigo-dc/udocker).

## Known bugs
Occasionally, wget will hang inside a container inexplicably after downloading a file.
With some investigation, I found that the file ".wget-hsts" in the home directory can cause this behavior inside singularity containers.
Deleting the file will solve the problem.
