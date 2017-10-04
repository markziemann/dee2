#!/bin/bash
set -x #-v
# sanitise data regularly by cron and move validated data to /home/pi/dee2_data/
#add the following tocrontab using "crontab -e"
#*/10 * * * * bash -c "/home/pi/clean.sh"
#put script in the location "/home/pi/clean.sh"

DATA=/home/pi/dee2_data/
#the below needs to be integrated to automate incorporation of volunteer data
cd $DATA
AGE=$(date -r started +%s)
TIME=$(date +%s)
if [ $((TIME-AGE)) -gt 300 ] ; then
  rm started
fi

if [ ! -r started ] ; then
  touch started
  SFTP_INCOMING=/sftp/guestuser/incoming

  if [ -d dee2 ] ; then
    cd dee2 ; git pull ; cd ..
  else
    git clone https://github.com/markziemann/dee2.git
  fi

  REFERENCE_PIPELINE_MD5SUM=$(md5sum dee2/volunteer_pipeline.sh | awk '{print $1}')

  if [ "$(ls -A ${SFTP_INCOMING})" ]; then

    for FILE in ${SFTP_INCOMING}/* ; do

      if [ "$(echo $FILE | grep -c .zip$)" -ne "1" ] ; then
        echo "now rm $FILE"
        sudo rm $FILE
      else
        TIME=$(date +%s)
        FILETIME=$(stat --format "%Y" $FILE)

        if [ $((TIME-FILETIME)) -gt "60" ] ; then
          echo process $FILE $TIME $FILETIME
          BASE=$(basename $(echo $FILE | sed 's/.zip$//') )
          ORG=$(echo $BASE | tr '.' ' ' | awk '{print $NF}')
          SRR=$(echo $BASE | sed "s/.${ORG}//")
          PIPELINE_MD5SUM=$(unzip -p $FILE $SRR/volunteer_pipeline.sh | md5sum | awk '{print $1}')

          INVALID=0
          if [ "$REFERENCE_PIPELINE_MD5SUM" != "$PIPELINE_MD5SUM" ] ; then INVALID=$((INVALID+1)) ; fi

          unzip -t $FILE || INVALID=$((INVALID+1))

          if [ $INVALID -eq "0" ] ; then
            sudo mv $FILE $DATA
          else
            sudo rm $FILE
          fi
        fi
      fi
    done
  fi
  rm started
fi
