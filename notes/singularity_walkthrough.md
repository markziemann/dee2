# Running DEE2 via singularity
TODO

# Creating a DEE2 singularity image from a docker image
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
You'll see that there's a new .img file in the current working directory

Now time to make some modifications as singularity is designed to be run by non-root user as compared to docker.

```
mkdir tmp
sudo singularity shell -w -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img
```

Once inside the container, check disk space folder and shift /dee2 to path with storage, then check that it works, then exit

```
df -h
mv /dee2 /tmp
ln -s /tmp/dee2 /
cd /dee
chmod -R 777 *
rm test_pass
rm -rf /dee2/data/ecoli/SRR057750
exit
```

The changes made to the container will be saved. Now check that it works from outside the container
```
sudo chown mziemann mziemann_tallyup-2018-01-19-c76488ce3d2d.img
singularity run -w -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-19-c76488ce3d2d.img
```

You're now free to run the image, so long as the /tmp directory exists in the present working directory. Here are three examples.
```
singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img ecoli SRR057751
singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img celegans
singularity run -B $(pwd)/tmp:/tmp mziemann_tallyup-2018-01-18-xxx.img

```
