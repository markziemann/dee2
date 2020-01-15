#!/bin/bash
# downloader script for DEE2 pipeline on HPC
#Copyright Mark Ziemann 2015 to 2020 mark.ziemann@gmail.com
set -x

MYPWD=$PWD

#Fix locale issue
export LANGUAGE=C
export LC_ALL=C
export LANG=C
export LC_TYPE=C

NUMVARS=$#

shopt -s expand_aliases

#handling verbosity setting
LASTVAR=$(echo $@ | rev | cut -d ' ' -f1 | rev)
if [ $LASTVAR == "-V" ] ; then
  set -x
  VERBOSE=TRUE
  NUMVARS=$#
  NUMVARS=$((NUMVARS-1))
else
  alias wget='wget -q'
  alias curl='curl -s'
fi

MY_ORG=$1
if [ $NUMVARS -gt 1 ] ; then
  if [ $2 != '-f' ] ; then
    MY_ACCESSIONS=$2
  fi
fi
MEM_FACTOR=2

main(){
#logging all commands
NUMVARS=$#
LASTVAR=$(echo $@ | rev | cut -d ' ' -f1 | rev)
VERBOSE=$(echo $LASTVAR | cut -d '=' -f2)
if [ ! -z $VERBOSE ] ; then
  if [ $VERBOSE == "TRUE" ] ; then
    set -x
  fi
fi

#define bad exit
exit1(){
rm *fastq *.sra *tsv
return 1
}
export -f exit1

#JOB
ORG=$1

if [ $2 != '-f' ] ; then
  SRR_FILE=$2
  SRR=$(basename $SRR_FILE .sra)
  echo $SRR
  wget -O $SRR.html "https://www.ncbi.nlm.nih.gov/sra/${SRR}"
  ORG2=$(echo $ORG | cut -c2-)
  ORG_OK=$(sed 's/class=/\n/g' $SRR.html  | grep 'Organism:' | grep -c $ORG2)
  rm $SRR.html
  if [ $ORG_OK -ne 1 ] ; then
    echo Annotated species name from NCBI SRA does not match user input! Quitting. | tee -a $SRR.log
    exit1 ; exit 1
  else
    echo User input species and SRA metadata match. OK.
  fi
fi

#ENVIRONMENT VARS
DEE_DIR=$MYPWD
cd $DEE_DIR
CODE_DIR=$DEE_DIR/code
PIPELINE=$0
PIPELINE_MD5=$(md5sum $MYPWD/$PIPELINE | cut -d ' ' -f1)
SW_DIR=$DEE_DIR/sw
PATH=$PATH:$SW_DIR
DATA_DIR=$DEE_DIR/data/$ORG
REF_DIR=$DEE_DIR/ref
QC_DIR=$DEE_DIR/qc

#LIMITS
DISKLIM=32000000
DLLIM=1
ALNLIM=2
MEMALNLIM=4
THREADS=$(grep -c ^processor /proc/cpuinfo)
DISK=$(df . | awk 'END{print$4}')
MEM=$(echo $(free | awk '$1 ~ /Mem:/  {print $2-$3}') \
  $(free | awk '$1 ~ /Swap:/  {print $2-$3}') \
  | awk '{print $1+$2}' )

##########################################################################
#Initial disk space check
##########################################################################
if [ $DISK -lt $DISKLIM ] ; then
  echo Error low disk space $DISK available $DISKLIM limit
  exit1 ; return 1
fi

##########################################################################
# Lets test all the input variables
##########################################################################
if [ ! -d "$QC_DIR" ] ; then
  mkdir -p $QC_DIR
fi

##########################################################################
# Lets get started
##########################################################################
if [ ! -d $DATA_DIR ] ; then mkdir -p $DATA_DIR ; fi
cd $DATA_DIR

if [ $2 != '-f' ] ; then
  mkdir $SRR ; cp $PIPELINE $SRR ; cd $SRR
  echo "Starting $PIPELINE $SRR
    current disk space = $DISK
    free memory = $MEM " | tee -a $SRR.log

##########################################################################
# Check number of attempts
##########################################################################
  ATTEMPTS=$SRR.attempts.txt

  if [ -r $SRR.attempts.txt ] ; then
    NUM_ATTEMPTS=$(wc -l < $ATTEMPTS)
    if [ $NUM_ATTEMPTS -gt "2" ] ; then
      echo $SRR has already been tried 3 times, skipping
      exit1 ; return 1
    fi
  fi
  DATE=`date +%Y-%m-%d:%H:%M:%S`
  echo $PIPELINE $PIPELINE_MD5 $DATE >> $ATTEMPTS

##########################################################################
  echo $SRR check if SRA file exists and download if neccessary
##########################################################################
  if [ ! -f $SRR.sra ] ; then
    #build URL
    BASEURL=ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra
    PFX1=$(echo $SRR | cut -c-3)
    PFX2=$(echo $SRR | cut -c-6)
    URL=anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/${PFX1}/${PFX2}/${SRR}/${SRR}.sra
    ID=$DEE_DIR/.ascp/aspera-license

    if [ ! -d $DEE_DIR/.ascp ] ; then
      mkdir -p $DEE_DIR/.ascp
      cat << EOF > $ID
-----BEGIN DSA PRIVATE KEY-----
MIIBuwIBAAKBgQDkKQHD6m4yIxgjsey6Pny46acZXERsJHy54p/BqXIyYkVOAkEp
KgvT3qTTNmykWWw4ovOP1+Di1c/2FpYcllcTphkWcS8lA7j012mUEecXavXjPPG0
i3t5vtB8xLy33kQ3e9v9/Lwh0xcRfua0d5UfFwopBIAXvJAr3B6raps8+QIVALws
yeqsx3EolCaCVXJf+61ceJppAoGAPoPtEP4yzHG2XtcxCfXab4u9zE6wPz4ePJt0
UTn3fUvnQmJT7i0KVCRr3g2H2OZMWF12y0jUq8QBuZ2so3CHee7W1VmAdbN7Fxc+
cyV9nE6zURqAaPyt2bE+rgM1pP6LQUYxgD3xKdv1ZG+kDIDEf6U3onjcKbmA6ckx
T6GavoACgYEAobapDv5p2foH+cG5K07sIFD9r0RD7uKJnlqjYAXzFc8U76wXKgu6
WXup2ac0Co+RnZp7Hsa9G+E+iJ6poI9pOR08XTdPly4yDULNST4PwlfrbSFT9FVh
zkWfpOvAUc8fkQAhZqv/PE6VhFQ8w03Z8GpqXx7b3NvBR+EfIx368KoCFEyfl0vH
Ta7g6mGwIMXrdTQQ8fZs
-----END DSA PRIVATE KEY-----
EOF
      chmod 700 $DEE_DIR/.ascp
    fi
    prefetch -a "/usr/local/bin/ascp|/dee2/.ascp/aspera-license" $SRR \
    || ( echo $SRR failed download with prefetch | tee -a $SRR.log ; sleep 5 ; exit1 ; return 1 )
    mv ~/ncbi/public/sra/${SRR}.sra .
  fi

##########################################################################
  echo $SRR Validate the SRA file
##########################################################################
  echo $SRR SRAfilesize $SRASIZE | tee -a $SRR.log
  md5sum $SRR.sra | tee -a $SRR.log
  VALIDATE_SRA=$(vdb-validate $SRR.sra &> /dev/stdout  | head -4 | awk '{print $NF}' | grep -c ok)
  if [ $VALIDATE_SRA -eq 4 ] ; then
    echo $SRR.sra file validated | tee -a $SRR.log
  else
    echo $SRR.sra md5sums do not match. Deleting and exiting | tee -a $SRR.log
    exit1 ; return 1
  fi

fi

}
export -f main

cd $MYPWD

MEM=$(echo $(free | awk '$1 ~ /Mem:/  {print $2-$3}') \
  $(free | awk '$1 ~ /Swap:/  {print $2-$3}') \
  | awk '{print $1+$2}' )
NUM_CPUS=$(grep -c ^processor /proc/cpuinfo)
CPU_SPEED=$(lscpu | grep MHz | awk '{print $NF}' | sort -k2gr)

ACC_URL="http://dee2.io/acc.html"
ACC_REQUEST="http://dee2.io/cgi-bin/acc.sh"
TMPHTML=/tmp/tmp.$RANDOM.html
wget --no-check-certificate -r -O $TMPHTML "dee2.io/ip"
SFTP_URL=$(cat $TMPHTML )

rm $TMPHTML

if [ ! -z $MY_ORG ] ; then
  ORG_CHECK=$(echo 'athaliana celegans dmelanogaster drerio ecoli hsapiens mmusculus rnorvegicus scerevisiae' \
  | tr ' ' '\n' | grep -wc "$MY_ORG")
  if [ $ORG_CHECK -ne 1 ] ; then
    echo Organism not specified correctly. Check options and try again.
    exit 1
  fi

  MEM_REQD=$(echo 'athaliana        2853904
celegans        2652204
dmelanogaster   3403644
drerio  14616592
ecoli   1576132
hsapiens        28968508
mmusculus       26069664
rnorvegicus     26913880
scerevisiae     1644684' | grep -w $MY_ORG | awk -v f=$MEM_FACTOR '{print $2*f}')

fi

if [ -z $MY_ORG ] ; then
  echo Error, a species name is required
  exit 1
fi

#echo $MY_ORG

myfunc(){
MY_ORG=$1
ACC_REQUEST=$2
TMPHTML=/tmp/tmp.$RANDOM.html
wget --no-check-certificate -r -O $TMPHTML "${ACC_REQUEST}?ORG=${MY_ORG}&Submit"
wget --no-check-certificate -O $TMPHTML $(grep 'frame src=' $TMPHTML | cut -d '"' -f2)
ACCESSION=$(grep 'ACCESSION=' $TMPHTML | cut -d '=' -f2)
echo $ACCESSION
rm $TMPHTML
}
export -f myfunc

key_setup(){
TMPHTML=/tmp/tmp.$RANDOM.html
wget --no-check-certificate -r -O $TMPHTML "dee2.io/ip"
SFTP_URL=$(cat $TMPHTML )
mkdir -p $MYPWD/dee2/.ssh
rm $TMPHTML
touch $MYPWD/dee2/.ssh/known_hosts
chmod 777 $MYPWD/dee2/.ssh/known_hosts

cat << EOF > $MYPWD/dee2/.ssh/guestuser
-----BEGIN RSA PRIVATE KEY-----
MIIEpAIBAAKCAQEAyLJ5TihXnb2ATexgYMIkzpHgopCbctWKH8rrPZNN6PALRYjg
1ozfeMFylSvQilFw+6bCe7HlqUQ3e6pS/jHJukEyzbJOEVR4AwuZxxctI4QH00AL
2eDvWvlEOChlxPg8Er5SjPziUXw8Ov3bNLvFHSQ7qlNb/gbKhKvzl6Lk0n6Yzl9C
/eiwzTKjfEKfXAZ51fjyD2fmSFaVleq+t3zviZaGftFtOLKtDA9wXXiosYrBufEf
zixujQF04Hzv+Eg814bjzgkSpZiDyS735NUzu0PCbnXNjZA6QiymisOkhx0J7w3r
vn/gmlYMmeBa5GZZsnnfRBvj0grQIefkLS30RwIDAQABAoIBAHVdUWzwUJRxPjfT
dGUBA689RaUrdYxI7hY7fyeqHdSLk7vdGMa+6OxgDBbJ4ZERoUW4tmDJnqlGuD98
Uj5OdU6TVBdQHzEpOWlmfk4b8oyjaEQUXxnR3YdQ36ELlsAB/ndjjzjdpafLRBmn
XGpRKCsrhizLxK8f34yIVdImMzQYQ7Enki003AgmEWZ/hZmOJtbXWHq/MIGk67Gq
rD3UJL+w0OVgQMYdD57CNBlpQIVDu4Z2C7NPLW/n2DiatzZ+7wOSWfc3I2Gu1E5o
/YV84Pa0dzwpCnBSNuWtrieSHgF96R2rBk/2slN/q1MV0XAFxFqnup1A/YpmdCI2
04+Q5vkCgYEA/IzR+nGquL/bJevGhtanMMxGMVSZJuCYGQU33R/WrWbDz1AJqPtd
/rQlFWcfkK4hpcZdGNIoVkH3aa/mfhMfGx1DEScxzoFaPEj2vKBDldPYyMW7owVH
ByPP0EiWwmERi7Ds6o/F325b2w0c1+waOTbA1eD9/dUZzExgVeBYes0CgYEAy3BS
TJ+5/wu0XkQi2qqUck4hdB6VrmLujTGcT6MmOyncGL0Y6SR9cvf51UpbPsCQpCEm
bOvRia/3wq9ovKIP3Zx22+SFSeu0bGeYo+i2ofzxl4XzZo3JIMpJtRXahn4BAH5E
PzXd/Hs4AkCgnQB3HXDyp3FSDFxC0V7/jvO+U2MCgYEAioAb47IUg09MOuarsGTl
ucA9Om5/sy92mjofYdhFHkF+XyIwughoivd2Yt90Ex87+rLneWY/ktaIfeBmknug
EnmgvzZ0fSC5QNhu4BEwH2nXuHugJI4PXt4H6Nz2ONGNEsPLmfOQ+7CFFYOCbvPf
icL6TBEgmeUVSdIU/uOTAn0CgYAN6OsnpBAymRlHDL+ZVep6ek8dQm4Xk1oeO1Ml
utEFYJJU+rD2V/Ff6AakB8Z/XulE36Nh9SnJkUeOfzHZG/ebvnP+Cvz2FfCrLNYp
9uJt5v6ZzqXa0Dz9SfeKMylS4tCsuPVvoP5BoictOEADHCII2E0vF7d1cuV6rVUp
8A6GYwKBgQCc8T4sr/wF/BKkk+kCBUpsYibqIxnTw7Rjl/+gUJL5IR3ynmWuJkUt
Qzab+/WnlQMuslmCLxXXOijq5lEDJLJ0m9hZ0sdC+j13jsTCEOnyj/XJ3VgLKifP
8itVEOnDffxs+RKeaXWhPiSll/wp6SlSuIdI2VpYMd15LtmkSkZSYg==
-----END RSA PRIVATE KEY-----
EOF

cat << EOF > $MYPWD/dee2/.ssh/guestuser.pub
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDIsnlOKFedvYBN7GBgwiTOkeCikJty1Yofyus9k03o8AtFiODWjN94wXKVK9CKUXD7psJ7seWpRDd7qlL+Mcm6QTLNsk4RVHgDC5nHFy0jhAfTQAvZ4O9a+UQ4KGXE+DwSvlKM/OJRfDw6/ds0u8UdJDuqU1v+BsqEq/OXouTSfpjOX0L96LDNMqN8Qp9cBnnV+PIPZ+ZIVpWV6r63fO+JloZ+0W04sq0MD3BdeKixisG58R/OLG6NAXTgfO/4SDzXhuPOCRKlmIPJLvfk1TO7Q8Judc2NkDpCLKaKw6SHHQnvDeu+f+CaVgyZ4FrkZlmyed9EG+PSCtAh5+QtLfRH mdz@opti
EOF

chmod -R 700 $MYPWD/dee2/.ssh
}
export -f key_setup

##################################################
# Testing the pipeline with ecoli sample
##################################################
TESTFILE=$MYPWD/dee2/test_pass
if [ ! -r $TESTFILE ] ; then
  echo Initial pipeline test with E. coli dataset
  if [ -d $MYPWD/dee2/data/ecoli/SRR057750 ] ; then
    rm -rf $MYPWD/dee2/data/ecoli/SRR057750
  fi

  #test ssh key setup
  key_setup $SFTP_URL
  mkdir -p $MYPWD/dee2/data/ecoli
  cd $MYPWD/dee2/data/ecoli
  date > date.txt
  sftp -v -i $MYPWD/dee2/.ssh/guestuser -o StrictHostKeyChecking=no guestuser@$SFTP_URL << EOF && KEYTEST="OK"
put date.txt
EOF

  if [ $KEYTEST == "OK" ] ; then
    echo "SSH keys successfully set up"
  else
    echo "SSH keys not set up properly. Quitting now."
    exit 1
  fi

  cd $MYPWD/dee2
  date +"%s" > $TESTFILE

else

##################################################
# Testing whether the user has provided SRR accessions
##################################################
  if [[ $NUMVARS -eq "2" && $OWN_DATA -eq "0" ]] ; then
    TESTACCESSIONS=$(echo $2 | tr ',' '\n' | cut -c2-3 | grep -vc RR)
    if [ $TESTACCESSIONS -eq 0 ] ; then
      for USER_ACCESSION in $(echo $2 | tr ',' ' ') ; do
        DIR=$(pwd)
        echo Starting pipeline with species $1 and accession $USER_ACCESSION
        main $1 $USER_ACCESSION VERBOSE=$VERBOSE
        cd $MYPWD/dee2/data/$MY_ORG
      done
      exit
    else
      echo Looks like accession numbers aren\'t in the correct format ie. SRR123456,ERR789456,DRR321654
      echo Check input parameters and try again.
      exit
    fi
  fi

##################################################
# If no accessions are provided
#################################################

  while true ; do
    count=$(find $MYPWD/data/ecoli/ | grep -c .sra$)

    if [ $count -gt 5 ] ; then
      sleep 10
    else
      echo Run "$count" of 1000
      ACCESSION=$(myfunc $MY_ORG $ACC_REQUEST)
      echo Starting pipeline with species $MY_ORG and accession $ACCESSION
      main "$MY_ORG" "$ACCESSION" VERBOSE=$VERBOSE && COMPLETE=1 || COMPLETE=0
      cd $MYPWD
    fi
  done
fi
