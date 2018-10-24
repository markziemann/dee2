# Singularity guide
This guide consists of two parts: (i) creating the image and (ii) running the image on HPC with slurm based scheduling.

## Creating a DEE2 singularity image from a docker image
The strategy is to convert a working docker image and then make some modifications to make it suitable for non-root user on HPC.

Now need a box with root priviledges install singularity:

```
sudo wget -O- http://neuro.debian.net/lists/xenial.us-ca.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
sudo apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9
sudo apt-get update
sudo apt-get install -y singularity-container
```
then check that singularity works
```
singularity --version
```

Great, now change current working directory to somewhere with lots of storage (>100GB).

Now convert docker image to singularity according to this script (link here). Obviously this cant be done on the HPC, so it will need to be done on another PC and then copied over to the HPC. (I'm working on a better way to do this but the solution below works)

```
docker run \
-v /var/run/docker.sock:/var/run/docker.sock \
-v $(pwd):/output \
--privileged -t --rm \
singularityware/docker2singularity \
mziemann/tallyup
```
You'll see that there's a new .img file in the current working directory. Now need to make some modifications as singularity is designed to be run by non-root user as compared to docker.

```
mkdir tmp
sudo singularity shell -w -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img
```

Once inside the container, check disk space folder and shift /dee2 to path with storage, then check that it works, then exit

```
df -h
mv /dee2 /tmp
ln -s /tmp/dee2 /
cd /dee2
chmod -R 777 *
rm test_pass
rm -rf /dee2/data/ecoli/SRR057750
exit
```

The changes made to the container will be saved. Now check that it works from outside the container. 
```
sudo chown mziemann mziemann_tallyup-2018-01-18-xxx.img
singularity run -w -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img
```
Copy the image and tmp folder to the HPC (without root access) and test it as above.

## Running DEE2 via singularity
(Under construction)
Download the current image from [here](https://vm-118-138-241-34.erc.monash.edu.au/images/current_singularity_img.tar.gz).

unpack the image.
```
tar xf current_singularity_img.tar.gz
```
You will see an img file and tmp folder, these must remain in the same working directory.

Run the test in writable mode, ensuring the /tmp directory exists in the present working directory. This will write new ssh keys and confirm the pipeline is working as expected.
```
singularity run -w -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img
```

You're now free to run the image. Here are three examples to try in a normal shell.
```
singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img ecoli SRR057751
singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img celegans
singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img
nohup run-one-constantly singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img
nohup bash -c 'while true ; do singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img ; done'
```

Now try to schedule jobs. Below is an example sbatch file that works for me. Modify memory requirements and target organism as needed.
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
nohup bash -c 'while true ; do singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img >sbatch.log 2>&1 ; done'
```
## Known bugs
Occassionally, wget will hang inside a container inexplicably after downloading a file. With some investigation, I found that the file ".wget-hsts" in the home directory can cause this behavior inside singularity containers. Deleting the file will solve the problem.
