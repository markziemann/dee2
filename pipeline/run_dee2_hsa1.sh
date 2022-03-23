#!/bin/bash
while true ; do

  CNT=$(docker ps | grep -cw tallyup )

  if [ $CNT -lt 2 ] ; then

    docker run -v $(pwd):/dee2/mnt mziemann/tallyup hsapiens -d

  fi

  sleep 10

done
