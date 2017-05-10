#!/bin/bash
set -x
#sra2mx
#Copyright Mark Ziemann 2015 to 2017 mark.ziemann@gmail.com

#JOB
##SRR to process
#URL=$1
#SRR=`basename $URL | sed 's/.sra$//'`
SRR_FILE=$1
SRR=$(basename $SRR_FILE .sra)
echo $SRR
##ORGANISM
ORG=$2

#ENVIRONMENT VARS
PIPELINE=$0
CODE_DIR=$(dirname $PIPELINE)
DEE_DIR=$(dirname $CODE_DIR)
SW_DIR=$DEE_DIR/sw
PATH=$PATH:$SW_DIR
DATA_DIR=$DEE_DIR/data/$ORG
REF_DIR=$DEE_DIR/ref
QC_DIR=$DEE_DIR/qc

#echo running $PIPELINE $SRR_FILE $ORG
#exit

#SOFTWARE VARS
VDBVAL=$SW_DIR/vdb-validate
FQDUMP=$SW_DIR/fastq-dump
ABIDUMP=$SW_DIR/abi-dump
SRA_STAT=$SW_DIR/sra-stat
PFQDUMP=$SW_DIR/parallel-fastq-dump
FASTQC=$SW_DIR/FastQC/fastqc
SKEWER=$SW_DIR/skewer
MINION=$SW_DIR/kraken/./minion
BOWTIE2=$SW_DIR/bowtie2-2.3.1/bowtie2
SOLIDTRIMMER=$SW_DIR/solid-trimmer.py
STAR=$SW_DIR/STAR
SUBJUNC=/data/app/bin/subjunc
FEATURECOUNTS=$SW_DIR/featureCounts
KALLISTO=$SW_DIR/kallisto

#REFERENCE SEQ AND ANNOTATIONS
ENS_REFG=$REF_DIR/$ORG/ensembl/star
ENS_GTF=$(readlink -f $(find $REF_DIR/$ORG/ensembl/ -maxdepth 1 | grep .gtf$) )
ENS_REFT=$(readlink -f $(find $REF_DIR/$ORG/ensembl/kallisto/ -maxdepth 1 | grep fa.idx$) )
ENS_REFT_BT2=$(readlink -f $(find $REF_DIR/$ORG/ensembl/bowtie2/ -maxdepth 1 | grep .fa$) )
NCBI_REFG=$REF_DIR/$ORG/ncbi/star
# The ~/bfx/dee2/ref/dmelanogaster/ncbi/fix_gff.sh script should be used to convert gff into something
# that STAR can understand but it must be checked that nonunique genenames dont exist and the
# newly made gtf works properly
NCBI_GTF=$(readlink -f $(find $REF_DIR/$ORG/ncbi/ -maxdepth 1 | grep .gtf$) )
NCBI_REFT=$(readlink -f $(find $REF_DIR/$ORG/ncbi/kallisto/ -maxdepth 1 | grep fna.idx$) )


#LIMITS
DISKLIM=32000000
DLLIM=1
ALNLIM=2
MEMALNLIM=4
THREADS=$(nproc)
DISK=$(df . | awk 'END{print$4}')
MEM=$(free | awk '$1 ~ /Mem:/  {print $2-$3}')

##########################################################################
# Lets test all the input variables
##########################################################################
#check if these are directories
if [ ! -d "$SW_DIR" ] ; then
  echo software directory does not exist! Quitting.
  exit 1
fi

if [ ! -d "$REF_DIR" ] ; then
  echo reference sequence directory does not exist! Quitting.
  exit 1
fi

#check if these are directories
if [ ! -d "$DATA_DIR"  ] ; then
  echo data directory does not exist! Quitting.
  exit 1
elif [ ! -w "$DATA_DIR"  ] ; then
  echo data directory not writable! Quitting.
  exit 1
fi

if [ ! -d "$QC_DIR" ] ; then
  echo QC report directory does not exist! Quitting.
  exit 1
fi

#check that all the software is present and working
FQDUMP_STATUS=$($FQDUMP -h | grep -xc './fastq-dump : 2.8.2')
ABIDUMP_STATUS=$($ABIDUMP -h | grep -xc './abi-dump : 2.8.2')
SKEWER_STATUS=$SW_DIR/skewer
SOLIDTRIMMER_STATUS=$SW_DIR/solid-trimmer.py
STAR_STATUS=$SW_DIR/STAR
SUBJUNC_STATUS=/data/app/bin/subjunc
FEATURECOUNTS_STATUS=$SW_DIR/featureCounts
KALLISTO_STATUS=$SW_DIR/kallisto

#check all the reference sequences exist
#if [[ ! -r $ENS_REFG || ! -r $ENS_GTF || ! -r $ENS_REFT || ! -r $NCBI_REFG || ! -r $NCBI_GTF || ! -r $NCBI_REFT ]] ; then
#  echo One or more of the reference fasta or annotation files is missing or not readable. Quitting | tee -a $SRR.log
#  exit
#fi

##########################################################################
# Lets get started
##########################################################################
cd $DATA_DIR
mkdir $SRR ; cp $0 $SRR ; cd $SRR
if [ -f ../$SRR.sra ] ; then mv ../$SRR.sra . ; fi
ATTEMPTS=$SRR.attempts.txt

echo "Starting $PIPELINE $CFG $URL
  current disk space = $DISK
  free memory = $MEM " | tee $SRR.log

##########################################################################
# Check number of attempts
##########################################################################
if [ -r $SRR.attempts.txt ] ; then
  NUM_ATTEMPTS=$(wc -l < $ATTEMPTS)
  if [ $NUM_ATTEMPTS -gt "2" ] ; then
    echo $SRR has already been tried 3 times, skipping ; exit
  fi
fi

##########################################################################
#Initial disk space check
##########################################################################
DISK=$(df . | awk 'END{print$4}')
if [ $DISK -lt $DISKLIM ] ; then
  echo Error low disk space $DISK available $DISKLIM limit ; exit
fi

##########################################################################
echo $SRR check if SRA file exists and download if neccessary
#might let R do this to maintain constant transfers
##########################################################################
if [ ! -f $SRR.sra ] ; then
  #build URL
  BASEURL=ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra
  PFX1=$(echo $SRR | cut -c-3)
  PFX2=$(echo $SRR | cut -c-6)
  URL=anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/${PFX1}/${PFX2}/${SRR}/${SRR}.sra
  ID=~/.ascp/aspera-license
  ascp -l 500m -O 33001 -T -i $ID $URL . \
  || ( echo $SRR failed ascp download | tee -a $SRR.log ; sleep 5 ; exit)
  SRASIZE=$(du ${SRR}.sra)
fi

##########################################################################
echo $SRR Validate the SRA file
##########################################################################
echo $SRR SRAfilesize $SRASIZE | tee -a $SRR.log
md5sum $SRR.sra | tee -a $SRR.log
VALIDATE_SRA=$($VDBVAL $SRR.sra &> /dev/stdout  | head -4 | awk '{print $NF}' | grep -c ok)
if [ $VALIDATE_SRA -eq 4 ] ; then
  echo $SRR.sra file validated | tee -a $SRR.log
else
  echo $SRR.sra md5sums do not match. Deleting and exiting | tee -a $SRR.log
  rm $SRR.sra ; exit
fi

##########################################################################
echo $SRR diagnose basespace/colorspace, single/paired-end and read length
##########################################################################
$FQDUMP -X 4000 --split-files $SRR.sra
NUM_FQ=$(ls | grep $SRR | grep -v trimmed.fastq | grep -c fastq$)
if [ $NUM_FQ -eq "1" ] ; then
  RDS=SE
  echo $SRR is single end | tee -a $SRR.log
elif [ $NUM_FQ -eq "2" ] ; then
  RDS=PE
  echo $SRR is paired end | tee -a $SRR.log
else
  echo Unable to determine if paired or single end. Quitting. | tee -a $SRR.log
  exit
fi

FQ1=$(ls  | grep $SRR | grep -m1 fastq$)
$FASTQC $FQ1
FQ1BASE=$(basename $FQ1 .fastq)

#diagnose colorspace or conventional
BASECALL_ENCODING=$(unzip -p ${FQ1BASE}_fastqc ${FQ1BASE}_fastqc/fastqc_data.txt \
| grep 'File type' | cut -f2 | awk '{print $1}')

if [ $BASECALL_ENCODING == "Colorspace" ] ; then
  CSPACE=TRUE
  echo $SRR is colorspace | tee -a $SRR.log
elif [ $BASECALL_ENCODING == "Conventional" ] ; then
  CSPACE=FALSE
  echo $SRR is conventional basespace | tee -a $SRR.log
else
  echo Unable to determine if colorspace or basespace. Quitting. | tee -a $SRR.log
fi

#diagnose read length then
#save entire fastq data to log and delete fastqc zip file and html report
FQ1_LEN=$(unzip -p ${FQ1BASE}_fastqc.zip ${FQ1BASE}_fastqc/fastqc_data.txt \
| grep 'Sequence length' | cut -f2)
echo $SRR read1 length is $FQ1_LEN nt | tee -a $SRR.log
unzip -p ${FQ1BASE}_fastqc.zip ${FQ1BASE}_fastqc/fastqc_data.txt | tee -a $SRR.log
rm ${FQ1BASE}_fastqc.zip ${FQ1BASE}_fastqc.html

if [ $RDS == "PE" ] ; then
  FQ2=$(ls  | grep $SRR | grep fastq$ | sed -n 2p)
  $FASTQC $FQ2
  FQ2BASE=$(basename $FQ2 .fastq)
  FQ2_LEN=$(unzip -p ${FQ2BASE}_fastqc.zip ${FQ2BASE}_fastqc/fastqc_data.txt \
  | grep 'Sequence length' | cut -f2)
  echo $SRR read2 length is $FQ2_LEN nt | tee -a $SRR.log
  unzip -p ${FQ2BASE}_fastqc.zip ${FQ2BASE}_fastqc/fastqc_data.txt | tee -a $SRR.log
  rm ${FQ2BASE}_fastqc.zip ${FQ2BASE}_fastqc.html
fi

##########################################################################
echo $SRR if colorspace, then quit
##########################################################################
if [ $CSPACE == "TRUE" ] ; then
  Colorspace data is excluded from analysis for now
  exit
fi

##########################################################################
echo $SRR Dump the fastq file
##########################################################################
rm ${SRR}*fastq
if [ $CSPACE == "FALSE" ] ; then
  #$FQDUMP --split-files --defline-qual '+' ${SRR}.sra
  ##try parallelising with fastq-dump
  $PFQDUMP --threads 8 --outdir . --split-files --defline-qual + -s ${SRR}.sra
fi

NUMRDS=$(sed -n '2~4p' $FQ1 | wc -l)
echo $SRR number of reads is $NUMRDS | tee -a $SRR.log
rm ${SRR}.sra

if [ "$NUMRDS" -eq 0 ] ; then
  echo $SRR has no reads. Aborting | tee -a $ATTEMPTS ; rm $FQ ; exit
fi

echo $SRR completed basic pipeline successfully | tee -a $SRR.log

##########################################################################
echo $SRR Quality trimming
##########################################################################
if [ $RDS == "SE" ] ; then
  skewer -l 18 -q 10 -t $THREADS -o $SRR $FQ1
  rm $FQ1
  FQ1=${SRR}-trimmed.fastq

elif [ $RDS == "PE" ] ; then
  skewer -l 18 -q 10 -t $THREADS -o $SRR $FQ1 $FQ2
  rm $FQ1 $FQ2
  FQ1=${SRR}-trimmed-pair1.fastq
  FQ2=${SRR}-trimmed-pair2.fastq
fi

#append the skewer log and exit if there are no reads passing QC
READCNT=$(wc -l < $FQ1 )
cat ${SRR}-trimmed.log >> $SRR.log && rm ${SRR}-trimmed.log
if [ $READCNT -eq "0" ] ; then
  echo No reads passed QC. Quitting | tee -a $SRR.log
  rm $SRR.sra
  exit
else
  echo $READCNT reads passed initial QC | tee -a $SRR.log
fi


##########################################################################
echo $SRR adapter diagnosis
##########################################################################
# I use minion to diagnose adapters
# http://wwwdev.ebi.ac.uk/enright-dev/kraken/reaper/src/reaper-latest/doc/minion.html
# If adapter is present higher than the threshold then it will be clipped
# using skewer. One problem is that we don't want the arbitrary threshold to kick
# in for some samples and not others. My idea is to clip any reads over this threshold.
# This will mean that the fastq file needs to be broken up into a part that
# won't be clipped and the one which will be clipped

ADAPTER_THRESHOLD=2
if [ $RDS == "SE" ] ; then
  MINION_LOG=$FQ1.minion.log
  $MINION search-adapter -i $FQ1 > $MINION_LOG
  ADAPTER=$(head $MINION_LOG | grep -m1 sequence= | cut -d '=' -f2)
  DENSITY=$(head $MINION_LOG | grep -m1 'sequence-density=' | cut -d '=' -f2 | numround -c)
  cat $MINION_LOG | tee -a $SRR.log && rm $MINION_LOG

  if [[ ! -z $DENSITY ]] ; then
    if [ $DENSITY -gt $ADAPTER_THRESHOLD ] ; then
      echo Potential 3prime adapter identified. Now checking if in reference sequence | tee -a $SRR.log
      # Query to see if adapter sequence present in reference
      ADAPTER_REF_CHECK=$($BOWTIE2 -f -x $ENS_REFT_BT2 -S /dev/stdout <(echo $ADAPTER | sed 's/^/>ADAPTER\n/') 2>>$SRR.log | awk '$1!~/^@/ && $2!=4' | wc -l )

      if [ $ADAPTER_REF_CHECK -eq "0" ] ; then
        echo Adapter seq not found in reference. Now shuffling file before clipping | tee -a $SRR.log
        #FQ1=${SRR}-trimmed.fastq this is the current name

        #shuffling the fq file to remove any tile effects
        paste - - - - < $FQ1 | shuf | tr '\t' '\n' > ${SRR}.fastq
        CLIP_LINE_NUM=$(echo $DENSITY $ADAPTER_THRESHOLD $READCNT | awk '{printf "%.0f\n", ($1-$2)/$1*$3*4}')
        head -$CLIP_LINE_NUM ${SRR}.fastq | skewer -l 18 -t $THREADS -x $ADAPTER -o $SRR -
        cat ${SRR}-trimmed.log >> $SRR.log && rm ${SRR}-trimmed.log
        CLIP_LINE_NUM1=$((CLIP_LINE_NUM+1))
        tail -n+$CLIP_LINE_NUM1 ${SRR}.fastq >> ${SRR}-trimmed.fastq && rm ${SRR}.fastq
        $MINION search-adapter -i ${SRR}-trimmed.fastq | tee -a $SRR.log
      else
        echo Potential adapter found in reference sequence. Continuing without clipping. | tee -a $SRR.log
      fi
    fi
  fi

elif [ $RDS == "PE" ] ; then
  MINION_LOG=$FQ1.minion.log
  $MINION search-adapter -i $FQ1 > $MINION_LOG
  ADAPTER1=$(head $MINION_LOG | grep -m1 sequence= | cut -d '=' -f2)
  DENSITY1=$(head $MINION_LOG | grep -m1 'sequence-density=' | cut -d '=' -f2 | numround -c)
  cat $MINION_LOG | tee -a $SRR.log && rm $MINION_LOG

  MINION_LOG=$FQ2.minion.log
  $MINION search-adapter -i $FQ2 > $MINION_LOG
  ADAPTER1=$(head $MINION_LOG | grep -m1 sequence= | cut -d '=' -f2)
  DENSITY2=$(head $MINION_LOG | grep -m1 'sequence-density=' | cut -d '=' -f2 | numround -c)
  cat $MINION_LOG | tee -a $SRR.log && rm $MINION_LOG

  DENSITY=$(echo $DENSITY1 $DENSITY2 | awk '{print ($1+$2)/2}' | numround)

  if [[ ! -z $DENSITY ]] ; then
    if [ $DENSITY -gt $ADAPTER_THRESHOLD ] ; then
      echo Potential 3prime adapter identified. Now checking if in reference sequence | tee -a $SRR.log
      # Query to see if adapter sequence present in reference
      ADAPTER1_REF_CHECK=$($BOWTIE2 -f -x $ENS_REFT_BT2 -S /dev/stdout <(echo $ADAPTER1 | sed 's/^/>ADAPTER\n/') 2>>$SRR.log | awk '$1!~/^@/ && $2!=4' | wc -l )
      ADAPTER2_REF_CHECK=$($BOWTIE2 -f -x $ENS_REFT_BT2 -S /dev/stdout <(echo $ADAPTER2 | sed 's/^/>ADAPTER\n/') 2>>$SRR.log | awk '$1!~/^@/ && $2!=4' | wc -l )

      if [ $ADAPTER1_REF_CHECK -eq "0" -a $ADAPTER2_REF_CHECK -eq "0" ] ; then
        echo Adapter seq not found in reference. Now shuffling file before clipping | tee -a $SRR.log
        paste $FQ1 $FQ2 | paste - - - - | shuf | awk -F'\t' '{OFS="\n"; print $1,$3,$5,$7 > "R1.fastq"; print $2,$4,$6,$8 > "R2.fastq"}'
        mv R1.fastq ${SRR}_1.tmp.fastq ; mv R2.fastq ${SRR}_2.tmp.fastq
        CLIP_LINE_NUM=$(echo $DENSITY $ADAPTER_THRESHOLD $READCNT | awk '{printf "%.0f\n", ($1-$2)/$1*$3*4}')
        head -$CLIP_LINE_NUM ${SRR}_1.tmp.fastq > ${SRR}_1.fastq &
        head -$CLIP_LINE_NUM ${SRR}_2.tmp.fastq > ${SRR}_2.fastq ; wait
        skewer -l 18 -t $THREADS -x $ADAPTER1 -y $ADAPTER2 -o $SRR ${SRR}_1.fastq ${SRR}_2.fastq
        cat ${SRR}-trimmed.log >> $SRR.log && rm ${SRR}-trimmed.log
        CLIP_LINE_NUM1=$((CLIP_LINE_NUM+1))
        tail -n+$CLIP_LINE_NUM1 ${SRR}_1.tmp.fastq >> ${SRR}-trimmed-pair1.fastq && rm ${SRR}_1.tmp.fastq
        tail -n+$CLIP_LINE_NUM1 ${SRR}_2.tmp.fastq >> ${SRR}-trimmed-pair2.fastq && rm ${SRR}_2.tmp.fastq
        $MINION search-adapter -i ${SRR}-trimmed-pair1.fastq | tee -a $SRR.log
        $MINION search-adapter -i ${SRR}-trimmed-pair2.fastq | tee -a $SRR.log
      else
        echo Potential adapter found in reference sequence. Continuing without clipping. | tee -a $SRR.log
      fi
    fi
  fi
#  skewer -l 18 -q 10 -t $THREADS -o $SRR $FQ1 $FQ2
#  rm $FQ1 $FQ2
#  FQ1=${SRR}-trimmed-pair1.fastq
#  FQ2=${SRR}-trimmed-pair2.fastq
fi

#append the skewer log and exit if there are no reads passing QC
FQSIZE=$(du -s $FQ1 | cut -f1)
cat ${SRR}-trimmed.log >> $SRR.log && rm ${SRR}-trimmed.log
if [ $FQSIZE -eq "0" ] ; then
  echo No reads passed QC. Quitting | tee -a $SRR.log
  rm  $SRR.sra
  exit
fi

##########################################################################
echo $SRR Starting mapping phase
##########################################################################
# need to setup 4 alignments
# -kallisto ensembl
# -kallisto ncbi
# -star/featurecounts ensembl
# -star/featurecounts ncbi

#STAR Ensembl
# runs star in gene-wise quant mode.
# We can then skip samtools and htseq/featurecounts, but
# unfortunately quant mode cannot currently be run with shared memory.
# Also cannot avoid on the fly GTF inclusion as its a requirement of quantMode option
# According to Google groups, quantMode counts only unique mappers but need to confirm this
# quantMode also looks like a nice way to diagnose the mapping strand
# Strand information could then be used for kallisto options
# The ReadsPerGene.out.tab file contains count information
echo $SRR starting STAR mapping to Ensembl genome | tee -a $SRR.log
if [ $RDS == "SE" ] ; then
  $STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad NoSharedMemory \
  --outSAMtype None --genomeDir $ENS_REFG --sjdbGTFfile $ENS_GTF --readFilesIn=$FQ1
elif [ $RDS == "PE" ] ; then
  $STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad NoSharedMemory \
  --outSAMtype None --genomeDir $ENS_REFG --sjdbGTFfile $ENS_GTF --readFilesIn=$FQ1 $FQ2
fi
cat Log.final.out | tee -a $SRR.log && rm Log.final.out Log.out Log.progress.out SJ.out.tab
head -4 ReadsPerGene.out.tab | tee -a $SRR.log
mv ReadsPerGene.out.tab $SRR.se.tsv

# Will need to investigate the compatibility of NCBI gff with STAR quantMode
# the fix_gff.sh script will fix the gff into something STAR can use
echo $SRR starting STAR mapping to NCBI genome | tee -a $SRR.log
if [ $RDS == "SE" ] ; then
  $STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad NoSharedMemory \
  --outSAMtype None --genomeDir $NCBI_REFG --sjdbGTFfile $NCBI_GTF \
  --sjdbGTFtagExonParentTranscript gene_name --readFilesIn=$FQ1
elif [ $RDS == "PE" ] ; then
  $STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad NoSharedMemory \
  --outSAMtype None --genomeDir $NCBI_REFG --sjdbGTFfile $NCBI_GTF \
  --sjdbGTFtagExonParentTranscript gene_name --readFilesIn=$FQ1 $FQ2
fi
cat Log.final.out | tee -a $SRR.log && rm Log.final.out Log.out Log.progress.out SJ.out.tab
head -4 ReadsPerGene.out.tab | tee -a $SRR.log
mv ReadsPerGene.out.tab $SRR.sn.tsv

##########################################################################
echo $SRR diagnose strandedness now
##########################################################################
# 0=untranded 1=posstrant 2=negstrand
UNSTRANDED_CNT=$(cut -f2 $SRR.se.tsv | tail -n +5 | numsum)
POS_STRAND_CNT=$(cut -f3 $SRR.se.tsv | tail -n +5 | numsum)
NEG_STRAND_CNT=$(cut -f4 $SRR.se.tsv | tail -n +5 | numsum)

echo "UnstrandedReadsAssigned:$UNSTRANDED_CNT \
PositiveStrandReadsAssigned:$POS_STRAND_CNT \
NegativeStrandReadsAssigned:$NEG_STRAND_CNT" | tee -a $SRR.log

if [ $POS_STRAND_CNT -ge "$((NEG_STRAND_CNT*5))" ] ; then
  STRAND=1
  KALLISTO_STRAND_PARAMETER='--fr-stranded'
  echo "Dataset is classified positive stranded" | tee -a $SRR.log
elif [ $NEG_STRAND_CNT -ge "$((POS_STRAND_CNT*5))" ] ; then
  STRAND=2
  KALLISTO_STRAND_PARAMETER='--rf-stranded'
  echo "Dataset is classified negative stranded" | tee -a $SRR.log
else
  STRAND=0
  KALLISTO_STRAND_PARAMETER=''
  echo "Dataset is classified unstranded" | tee -a $SRR.log
fi
echo KALLISTO_STRAND_PARAMETER=$KALLISTO_STRAND_PARAMETER
#Now cut out columns to leave us with only the desired strand info
CUTCOL=$((STRAND+2))
cut -f1,$CUTCOL $SRR.sn.tsv | tail -n +5 > $SRR.sn.tsv.tmp && mv $SRR.sn.tsv.tmp $SRR.sn.tsv
cut -f1,$CUTCOL $SRR.se.tsv | tail -n +5 > $SRR.se.tsv.tmp && mv $SRR.se.tsv.tmp $SRR.se.tsv

##########################################################################
echo $SRR checking readlengths now for kmer selection
##########################################################################
## Setting the kallisto kmer correctly is important to getting best accuracy
## Here I measure the median length as well as the 20th percentile
## KMER is set to length at 20th percentile minus 4nt.
MEDIAN_LENGTH=$(sed -n '2~4p' $FQ1 | head -1000000 | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.50 - 0.5)]}')
D20=$(sed -n '2~4p' $FQ1 | head -1000000 | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.20 - 0.5)]}')
KMER=$((D20-4))
ADJUST=$(echo $KMER | awk '{print ($1+1)%2}')
KMER=$((KMER-ADJUST))

echo MeadianReadLen=$MEDIAN_LENGTH 20thPercentileLength=$D20 echo kmer=$KMER | tee -a $SRR.log

if [ $KMER -lt "31" ] ; then
  ENS_REFT=$(echo $ENS_REFT | sed "s#fa.idx#fa.k${KMER}.idx#")
  NCBI_REFT=$(echo $NCBI_REFT | sed "s#fna.idx#fna.k${KMER}.idx#")
fi

##########################################################################
echo $SRR running kallisto now
##########################################################################
#Kallisto Ensembl
if [ $RDS == "SE" ] ; then
  echo $SRR Starting Kallisto single end mapping to ensembl reference transcriptome. kmer=$KMER
############################################
# TODO need intelligent frag size specification
###########################################
  $KALLISTO quant $KALLISTO_STRAND_PARAMETER --single -l 100 -s 20 -t $THREADS -o . -i $ENS_REFT $FQ1 2>&1 \
  | tee -a $SRR.log && mv abundance.tsv $SRR.ke.tsv
  rm abundance.h5
elif [ $RDS == "PE" ] ; then
  echo $SRR Starting Kallisto paired end mapping to ensembl reference transcriptome
  $KALLISTO quant $KALLISTO_STRAND_PARAMETER -t $THREADS -o . -i $ENS_REFT $FQ1 $FQ2 2>&1 \
  | tee -a $SRR.log && mv abundance.tsv $SRR.ke.tsv
  rm abundance.h5
fi

#Kallisto NCBI
if [ $RDS == "SE" ] ; then
  echo $SRR Starting single end mapping to NCBI reference transcriptome
  $KALLISTO quant $KALLISTO_STRAND_PARAMETER --single -l 25 -s 8 -t $THREADS -o . -i $NCBI_REFT $FQ1 2>&1 \
  | tee -a $SRR.log && mv abundance.tsv $SRR.kn.tsv
  rm abundance.h5
elif [ $RDS == "PE" ] ; then
  echo $SRR Starting Kallisto paired end mapping to NCBI reference transcriptome
  $KALLISTO quant $KALLISTO_STRAND_PARAMETER -t $THREADS -o . -i $NCBI_REFT $FQ1 $FQ2 2>&1 \
  | tee -a $SRR.log && mv abundance.tsv $SRR.kn.tsv
  rm abundance.h5
fi

# Tidy up files
rm -rf run_info.json ${SRR}-trimmed.fastq _STARgenome

# Check tsv files
wc -l *tsv | tee -a $SRR.log
# Check that tsv files have the right number of entries

SN_NR=$(wc -l < $SRR.sn.tsv)
SE_NR=$(wc -l < $SRR.se.tsv)
KN_NR=$(wc -l < $SRR.kn.tsv)
KE_NR=$(wc -l < $SRR.ke.tsv)
SN_CNT=$(cat $NCBI_GTF.cnt)
SE_CNT=$(cat $ENS_GTF.cnt)
KN_CNT=$(cat $NCBI_REFT.cnt)
KE_CNT=$(cat $ENS_REFT.cnt)

if [ $SN_NR -eq $SN_CNT -a $SE_NR -eq $SE_CNT -a $KN_NR -eq $((KN_CNT+1)) -a $KE_NR -eq $((KE_CNT+1)) ] ; then
  echo $SRR completed mapping pipeline successfully | tee -a $SRR.log
  touch $SRR.finished
else
  echo "$SRR An error occurred. Count file line numbers don't match the reference." | tee -a $SRR.log
fi
exit
