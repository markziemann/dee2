# Running DEE via singularity
First, install singularity like this:

```
sudo wget -O- http://neuro.debian.net/lists/xenial.us-ca.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
sudo apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9
sudo apt-get update
sudo apt-get install -y singularity-container
```
then check that it work
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

Now check that the image works. Before we do that, consider whether the home directory needs to be bound (this is the default). Many servers have very limited $HOME data storage, so you may need to bind it using the -H option. Similarly the /tmp directory may be very limited and can be bound wuth the -B option

```
singularity run -H /scratch/mziemann/:/home/mziemann/ -B /scratch/mziemann:/tmp/ mziemann_tallyup-2017-12-03-95a2303d5deb.img scerevisiae
```
