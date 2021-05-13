# Essential Docker commands

## Build
docker build -t dee2 .

docker run -it -d dee2 bash

docker ps -a

docker commit container_name mziemann/tallyup

docker push mziemann/tallyup

## Maintenance

alias docker-purge='docker rm $(docker ps -aq) ; docker rmi $(docker images -aq)'

docker stop 4ede2bbff7d9

docker rm $(docker ps -a -q)

docker rmi $(docker images -q)

docker rmi $(docker images -q "dangling=true")

docker system prune

## copy container data to pwd
docker cp $(docker ps -alq):/root/data/ .

## run bash in an entrypoint container 
docker run -it --entrypoint /bin/bash <image>

## how to update test-pass without remaking the genome indexes
#prune system first

docker system prune

#now shell in and then delete the test_pass file. this will force ssh key update. 

docker run -it --entrypoint /bin/bash mziemann/tallyup_dre

#now get the container name

docker ps -a

#commit the container to the image

docker commit flamboyant_albattani mziemann/tallyup_dre

#run this to set the entrypoint and regenerate the ssh keys

docker run -it --entrypoint /dee2/code/volunteer_pipeline.sh mziemann/tallyup ecoli SRR5985593

#get new container name

docker ps -a

#now save the container as the image

docker commit modest_minsky mziemann/tallyup_dre
