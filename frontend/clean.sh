#!/bin/bash
set -x #-v
# sanitise data regularly by cron and move validated data to /home/pi/dee2_data/
#create the file /etc/cron.d/pi and paster in the following
#MAILTO=pi
#*/10 * * * * pi bash -c "/home/pi/clean.sh"
#put this script in the location "/home/pi/clean.sh"

DATA=/home/pi/dee2_data
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

    for FILE in $(find ${SFTP_INCOMING} ) ; do

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

          #INVALID=0 means data is valid and good. 1 or higher is bad
          INVALID=0
          if [ "$REFERENCE_PIPELINE_MD5SUM" != "$PIPELINE_MD5SUM" ] ; then INVALID=$((INVALID+1)) ; fi
          unzip -t $FILE || INVALID=$((INVALID+1))
          if [ $(du -s $FILE | cut -f1) -gt 6000 ] ; then INVALID=$((INVALID+1)) ; fi

          if [ $INVALID -eq "0" ] ; then
            #sudo mv $FILE $DATA
            mkdir $DATA/$ORG
#            unzip -of /sftp/guestuser/incoming/ERR1158067.scerevisiae.zip -d /home/pi/dee2_data/scerevisiae/
            unzip -o $FILE -d $DATA/$ORG && scp -r $DATA/$ORG/$SRR mdz@Z620:~/bfx/dee2/data/$ORG && sudo rm -rf $DATA/$ORG/$SRR $FILE
          else
            sudo rm $FILE
          fi
        fi
      fi
    done
  fi
  rm started
fi

#check time if older than 1 week (604800)
[[ $(( $(date +'%s') - $(stat --format "%Y" /home/pi/dee2_data/mx/*tsv | sort | head -1) )) -gt 120 ]] && \
 scp mdz@Z620:~/bfx/dee2/mx/* /home/pi/dee2_data/mx || \
 echo "not time to refresh"
