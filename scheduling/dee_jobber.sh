#!/bin/bash
if [ $(ps -aef | grep docker | grep -c tallyup ) -eq 0 ] ; then
  echo "no dockers running, starting tallyup" >> test_jobber.txt
  nohup docker run mziemann/tallyup
else
  echo "docker tallyup job already running - waiting..." >> test_jobber.txt
fi

