#!/bin/bash
set -x #-v
# sanitise data regularly by cron and move validated data to /home/user/dee2_data/
#create the file /etc/cron.d/user and paster in the following
#MAILTO=user
#0 0 * * * user bash -c "/home/user/clean.sh"
#put this script in the location "/home/user/clean.sh"

INCOMING=~/upload/
TMPDATA=~/tmpdata/
STARTED_FILE=$INCOMING/started
STASH=/dee2_data/stash

if [ ! -r started ] ; then
  touch $STARTED_FILE
  if [ -d dee2 ] ; then
    cd dee2 ; git pull ; cd ~
  else
    git clone https://github.com/markziemann/dee2.git
  fi

  PREV_REFERENCE_PIPELINE_MD5SUM1="750d197f551ef3ac98f9e83f8d61ca43"
  PREV_REFERENCE_PIPELINE_MD5SUM2="b658ab07180ba71dce32a255b1c32fa3"
  PREV_REFERENCE_PIPELINE_MD5SUM3="888a4b1e21a693aa53ebc5490bd49053"
  PREV_REFERENCE_PIPELINE_MD5SUM4="0cf85a7690f1a121fe283ebe9467f631"
  PREV_REFERENCE_PIPELINE_MD5SUM5="0b12dedda2f2ab88aa4f098215055ff1"
  PREV_REFERENCE_PIPELINE_MD5SUM6="eba8ca83e14fe4b9484e5d9e0e6ea4e7"
  PREV_REFERENCE_PIPELINE_MD5SUM7="5ef12e0158c2554b8cae2f32c787f942"
  PREV_REFERENCE_PIPELINE_MD5SUM8="7d63a30717e102bf7932ad56a566bfe4"
  PREV_REFERENCE_PIPELINE_MD5SUM9="7610cbde291320df137d9761ccc8090d"
  PREV_REFERENCE_PIPELINE_MD5SUM10="9c7b9b535947214b2b29cacbd12e6ddf"
  PREV_REFERENCE_PIPELINE_MD5SUM11="75a1c95dc37c82ae8aa4078f57fb9360"
  PREV_REFERENCE_PIPELINE_MD5SUM12="73f3388fbedb61bd6b8213d6769ac49c"
  PREV_REFERENCE_PIPELINE_MD5SUM13="fe3bc7b82b416d24cb99c0f85b4381d0"
  PREV_REFERENCE_PIPELINE_MD5SUM14="c55adbf9dac28a62b1fcd832e19a53e7"
  REFERENCE_PIPELINE_MD5SUM=$(md5sum ~/dee2/pipeline/volunteer_pipeline.sh | awk '{print $1}')

  if [ "$(ls -A ${INCOMING})" ]; then

    for FILE in $(find ${INCOMING} ${STASH}) ; do

      # delete files that are not zips
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
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM1" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM2" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM3" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM4" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM5" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM6" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM7" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM8" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM9" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM10" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM11" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM12" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM13" != "$PIPELINE_MD5SUM" ] \
          && [ "$PREV_REFERENCE_PIPELINE_MD5SUM14" != "$PIPELINE_MD5SUM" ] ; then
              INVALID=$((INVALID+1))
          fi
          unzip -t $FILE || INVALID=$((INVALID+1))
          if [ $(du -s $FILE | cut -f1) -gt 6000 ] ; then INVALID=$((INVALID+1)) ; fi

          if [ $INVALID -eq "0" ] ; then
            mkdir -p $TMPDATA/$ORG
            #test the connection
            ssh -p 2210 mdz@localhost "ls" >/dev/null && CONNECT=1 || CONNECT=0
            #test disk space on root is more than 3GB
            DF=$(df / | awk 'NR>1 {print $4}')
            [ $DF -gt 3249932 ] && STORAGE=1 || STORAGE=0

            if [ $CONNECT -eq 1 -a $STORAGE -eq 1 ] ; then
              unzip -o $FILE -d $TMPDATA/$ORG && \
              scp -r -P 2210 $TMPDATA/$ORG/$SRR mdz@localhost:/mnt/md0/dee2/data/$ORG
              sudo rm -rf $TMPDATA/$ORG/$SRR $FILE
              #sleep 1s
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

find /var/log/apache2/  -mtime +2 -exec sudo rm {} \;
#sudo find /tmp -mmin +59 -delete

