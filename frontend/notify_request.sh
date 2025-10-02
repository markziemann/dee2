#!/bin/bash

# crontab line
#*      *     *       *       *       bash /dee2_data/requests/notify_request.sh

cd /dee2_data/requests/

for ZIP in $(find . | grep zip ) ; do
  SRP=$(basename $ZIP .zip)
  echo $ZIP $SRP

  if [ ! -r /usr/lib/cgi-bin/newrequests/$SRP.confirmed ] ; then
    rm $ZIP
    exit
  fi

  if [ ! -r $SRP.notified ] ; then
    touch $SRP.notified
    ACC=$SRP
    EMAIL=$(grep EMAIL= /usr/lib/cgi-bin/newrequests/$SRP.confirmed | cut -d '=' -f2)

    EMAIL_BODY=$(cat <<'EOF'
Hello DEE2 user,

We're happy to notify you that the requested dataset ACC has been processed and is ready
for you to download.
Download link: https://dee2.io/requests/ACC.zip

Thanks for using the DEE2 data service!

Mark and the DEE2 Team

EOF
)

    echo "$EMAIL_BODY" \
    | sed "s/ACC/$ACC/" \
    | mailx -s "DEE2 notification: Dataset $ACC has been completed" "$EMAIL"

  fi
done
