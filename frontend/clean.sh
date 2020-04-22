#!/bin/bash
set -x #-v
# sanitise data regularly by cron and move validated data to /home/user/dee2_data/
#create the file /etc/cron.d/user and paster in the following
#MAILTO=user
#*/10 * * * * user bash -c "/home/user/clean.sh"
#put this script in the location "/home/user/clean.sh"

DATA=/home/ubuntu/dee2_data
STARTED_FILE=$DATA/started
SFTP_INCOMING=/sftp/guestuser/incoming
STASH=/mnt/stash
#the below needs to be integrated to automate incorporation of volunteer data
cd $DATA
AGE=$(date -r started +%s)
TIME=$(date +%s)
if [ $((TIME-AGE)) -gt 3600 ] ; then
  rm $STARTED_FILE
fi

if [ ! -r started ] ; then
  touch $STARTED_FILE

  if [ -d dee2 ] ; then
    cd dee2 ; git pull ; cd ..
  else
    git clone https://github.com/markziemann/dee2.git
  fi

  PREV_REFERENCE_PIPELINE_MD5SUM="750d197f551ef3ac98f9e83f8d61ca43"

  REFERENCE_PIPELINE_MD5SUM=$(md5sum dee2/pipeline/volunteer_pipeline.sh | awk '{print $1}')

  if [ "$(ls -A ${SFTP_INCOMING})" ]; then

    for FILE in $(find ${SFTP_INCOMING} ${STASH}) ; do

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
          if [ "$REFERENCE_PIPELINE_MD5SUM" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM" != "$PIPELINE_MD5SUM" ] ; then
              INVALID=$((INVALID+1))
          fi
          unzip -t $FILE || INVALID=$((INVALID+1))
          if [ $(du -s $FILE | cut -f1) -gt 6000 ] ; then INVALID=$((INVALID+1)) ; fi

          if [ $INVALID -eq "0" ] ; then
            #sudo mv $FILE $DATA
            mkdir $DATA/$ORG
#            unzip -o $FILE -d $DATA/$ORG && scp -i ~/.ssh/monash/id_rsa -r -l 8192 $DATA/$ORG/$SRR mziemann@118.138.246.227:/scratch/mziemann/dee2/data/$ORG && sudo rm -rf $DATA/$ORG/$SRR $FILE

            #test the connection
            ssh -i ~/.ssh/monash/id_rsa -p 2211 mdz@localhost "ls" >/dev/null && CONNECT=1 || CONNECT=0

            #test disk space on root is more than 3GB
            DF=$(df / | awk 'NR>1 {print $4}')
            [ $DF -gt 3249932 ] && STORAGE=1 || STORAGE=0

            if [ $CONNECT -eq 1 -a $STORAGE -eq 1 ] ; then
              unzip -o $FILE -d $DATA/$ORG && \
              rsync -Pavz -e "ssh -i ~/.ssh/monash/id_rsa -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -p 2211" $DATA/$ORG/$SRR mdz@localhost:/mnt/md0/dee2/data/$ORG  && \
              sudo rm -rf $DATA/$ORG/$SRR $FILE
              ssh -i ~/.ssh/monash/id_rsa -p 2211 mdz@localhost "chmod -R +w /mnt/md0/dee2/data/$ORG/$SRR"
            else
              mv $FILE $STASH
            fi
          else
            sudo rm $FILE
          fi
        fi
      fi
    done
  fi
  rm $STARTED_FILE
fi

# refresh matrices
# check time if older than 1 week (604800)
#[[ $(( $(date +'%s') - $(stat --format "%Y" /mnt/dee2_data/mx/*.bz2 | sort | head -1) )) -gt 86400 ]] && \
# scp -i ~/.ssh/monash/id_rsa mziemann@118.138.246.227:/scratch/mziemann/dee2/mx/*bz2 /mnt/dee2_data/mx || \
# echo "not time to refresh"
find /var/log/apache2/  -mtime +2 -exec rm {} \;
