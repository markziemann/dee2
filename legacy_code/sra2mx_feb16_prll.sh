#!/bin/bash
#sra2mx by mark ziemann mark.ziemann@gmail.com
#2015-02-16

#PARAMETERS
PIPELINE=$0
CFG=$1
ORG=`basename $CFG .cfg`
DATADIR=/data/projects/mziemann/geo2mx_project/v1/data/$ORG
cd $DATADIR
URL=$2
SRR=`basename $URL | sed 's/.sra$//'`
mkdir $SRR ; cp $0 $CFG $SRR ; cd $SRR
FQDUMP=`grep ^FQDUMP= $CFG | cut -d '=' -f2`
ABIDUMP=`grep ^ABIDUMP= $CFG | cut -d '=' -f2`
GENOME=`grep ^GENOME= $CFG | cut -d '=' -f2`
CSGENOME=`grep ^CSGENOME= $CFG | cut -d '=' -f2`
GTF=`grep ^GTF= $CFG | cut -d '=' -f2`
DISKLIM=`grep ^DISKLIM= $CFG | cut -d '=' -f2`
DLLIM=`grep ^DLLIM= $CFG | cut -d '=' -f2`
ALNLIM=`grep ^ALNLIM= $CFG | cut -d '=' -f2`
MEMALNLIM=`grep ^MEMALNLIM= $CFG | cut -d '=' -f2`
ALNLIM=`echo $ALNLIM $MEMALNLIM | tr ' ' '\n' | sort -g | head -1`
THREADS=`grep ^THREADS= $CFG | cut -d '=' -f2`
DISK=`df . | awk 'END{print$4}'`
MEM=`free | awk '/buffers\/cache/{print $4}'`
ATTEMPTS=$SRR.attempts.txt

FQTRIMMER=/data/app/bin/fastq_quality_trimmer
FQ2FA=/data/app/bin/fastq_to_fasta
STAR=/data/app/bin/STAR
SAMTOOLS=/data/app/bin/samtools
SOLIDTRIMMER=/usr/local/bin/solid-trimmer.py
SUBJUNC=/data/app/bin/subjunc
FEATURECOUNTS=/data/app/bin/featureCounts

echo "Starting $PIPELINE $CFG $URL
  current disk space = $DISK
  free memory = $MEM " | tee $SRR.log

##########################################################################
# Check number of attempts
##########################################################################
if [ -r $SRR.attempts.txt ] ; then
	NUM_ATTEMPTS=`wc -l < $ATTEMPTS`
	if [ $NUM_ATTEMPTS -gt "2" ] ; then
	echo $SRR has already been tried 3 times, skipping ; exit
	fi
fi

##########################################################################
#Initial disk space check
##########################################################################
DISK=`df . | awk 'END{print$4}'`
if [ $DISK -lt $DISKLIM ] ; then
echo Error low disk space $DISK available $DISKLIM limit ; exit
fi

##########################################################################
echo $SRR Download the SRA file
##########################################################################
SLEEP=`echo 1 | awk -v R=$RANDOM '{print R/5000}'`
sleep $SLEEP
ASCP=`ps -ax | grep -wv grep | grep -wc ascp `
while test "$ASCP" -ge "$DLLIM" ; do
sleep $SLEEP ; ASCP=`ps -ax | grep -wv grep | grep -c ascp`
done
TIME=`date +"%R" | cut -d ':' -f1`
DAY=`date +"%u"`
if [ $TIME -ge 9 -a $TIME -le 17 -a $DAY -le 5 ] ; then SPEED=100m ; else SPEED=500m ; fi
echo speed $SPEED
#axel -n 8 -q $URL
ID=~/.aspera/connect/etc/asperaweb_id_dsa.openssh
ASCPURL=`echo $URL | sed 's#ftp://ftp-trace.ncbi.nlm.nih.gov/#anonftp@ftp.ncbi.nlm.nih.gov:/#'`
ascp -l $SPEED -O 33001 -T -i $ID $ASCPURL . \
|| ( echo $SRR failed ascp download | tee -a $SRR.log ; sleep 5 ; exit)
SRASIZE=`du ${SRR}.sra`
echo $SRR SRAfilesize $SRASIZE | tee -a $SRR.log
md5sum $SRR.sra | tee -a $SRR.log

##########################################################################
#Test whether $SRR.sra is basespace or colorspace
##########################################################################
SEQRDS=`$FQDUMP -A ${SRR}.sra --split-files -Z -X 1000 \
| $FQTRIMMER -t 1 -l 1 -Q33 | head -4 | wc -l 2>>$SRR.log`
if [ "$SEQRDS" -eq 4 ]
then
CSPACE=FALSE
echo $SRR is basespace | tee -a $SRR.log

##########################################################################
#Check whether $SRR.sra is single or paired end
##########################################################################
$FQDUMP -X 10 --split-files -A ${SRR}.sra 2>>$SRR.log
RDS=`ls ${SRR}*fastq | wc -l`
if [ "$RDS" -eq 1 ] ; then RDS=SE ; else RDS=PE ; fi
echo $SRR is $RDS | tee -a $SRR.log
head -40 *fastq >> $SRR.log
rm *fastq

##########################################################################
#Check read length of ${SRR}.sra and do fastqc
##########################################################################
FQ=${SRR}.fastq
RDLEN=`$FQDUMP --split-files -A ${SRR}*.sra -Z -X 40000 | tee $FQ \
| sed -n '0~4p' | awk '{print length}' | sort -g | tail -1 2>>$SRR.log`
echo $SRR readlength is $RDLEN | tee -a $SRR.log
fastqc --extract -t 10 ${FQ} 2> /dev/null
cat ${SRR}*_fastqc/*txt >> $SRR.log
rm -rf ${SRR}*_fastqc* $FQ

##########################################################################
echo Dump the fastq file from ${SRR}.sra and quality trim
##########################################################################
PE=`$FQDUMP --split-files --defline-qual '+' -A ${SRR}.sra -Z | head -5 \
| sed -n '1~4p' | cut -d ' ' -f1 | paste - - | awk '$1==$2' | wc -l`

if [ $PE -eq 1 ] ; then
echo Start ${SRR}.sra PE
NUMRDS=`$FQDUMP --split-files --defline-qual '+' -A ${SRR}.sra -Z \
| cut -d ' ' -f1 | paste - - - - | sed -n '1~2p'\
| tr '\t' '\n' | $FQTRIMMER -t 20 -l 19 -Q33 \
| cut -d ' ' -f1 | $FQ2FA -Q33 \
| tee $FQ | wc -l 2>>$SRR.log`
else
echo SE
NUMRDS=`$FQDUMP --split-files --defline-qual '+' -A ${SRR}.sra -Z \
| cut -d ' ' -f1 | $FQTRIMMER -t 20 -l 19 -Q33 \
| cut -d ' ' -f1 | $FQ2FA -Q33 \
| tee $FQ | wc -l 2>>$SRR.log`
fi

echo $SRR number of reads is $NUMRDS | tee -a $SRR.log
rm ${SRR}.sra
if [ "$NUMRDS" -eq 0 ] ; then
echo $SRR has no reads - aborting | tee -a $ATTEMPTS ; rm $FQ ; exit
fi

##########################################################################
echo Starting alignment for $FQ
##########################################################################
#STR=`ps -ax | grep -w STAR | grep -c genomeLoad`
#SLEEP=`echo 1 | awk -v R=$RANDOM '{print R/5000}'`
#while test "$STR" -ge "$ALNLIM" ; do
#sleep $SLEEP ; STR=`ps -ax | grep -c STAR | grep -c genomeLoad`
#done
mv $FQ ${SRR}.fasta
FQ=$SRR.fasta
$STAR --genomeLoad LoadAndKeep --genomeDir $GENOME/ --readFilesIn $FQ \
--runThreadN $THREADS --sjdbGTFfile $GTF 2>>$SRR.log
SAM=${FQ}.STAR.sam
rm $FQ Log.progress.out
mv Aligned.out.sam $SAM
cat Log.final.out Log.out >> ${SRR}.log ; rm Log.final.out Log.out
mv SJ.out.tab ${SRR}.STAR.sj ; pbzip2 -qf ${SRR}.STAR.sj

##########################################################################
echo Sorting, index bam file and generate stats file for $SRR
##########################################################################
$SAMTOOLS view -uSh -@$THREADS $SAM \
| $SAMTOOLS sort -@$THREADS - ${SAM}.sort 2>>$SRR.log && rm $SAM
mv ${SAM}.sort.bam ${SRR}.bam ; BAM=${SRR}.bam

else
CSPACE=TRUE
##########################################################################
echo ${SRR} is colorspace | tee -a $SRR.log
##########################################################################
$ABIDUMP -A ${SRR}.sra 2>>$SRR.log
rm ${SRR}.sra
CSFQ=`ls ${SRR}.sra*.csfasta | head -1`
CSQUAL=`ls ${SRR}.sra*qual | head -1`

##########################################################################
echo Quality trimming colorspace data $CSFQ
##########################################################################
FQ=${SRR}.fastq
NUMRDS=`$SOLIDTRIMMER -c $CSFQ -q $CSQUAL --moving-average 7:12 \
 --min-qual 20 --min-read-length 19 | tee $FQ | wc -l`
rm $CSFQ $CSQUAL
RDLEN=`sed -n '0~4p' $FQ | head -100 | awk '{print length}' | sort -g | tail -1`
echo $SRR length is $RDLEN | tee -a $SRR.log
echo $SRR number of reads is $NUMRDS | tee -a $SRR.log
if [ "$NUMRDS" -eq 0 ] ; then
echo $SRR has no reads - aborting | tee -a $ATTEMPTS ; rm $FQ ; exit
fi

##########################################################################
echo FastQC $FQ
##########################################################################
TMPFQ=`echo $FQ | sed 's/.fastq/.fq/'`
head -120000 $FQ > $TMPFQ
fastqc --extract -t 10 $TMPFQ 2> /dev/null
cat ${SRR}*_fastqc/*txt >> $SRR.log
rm -rf ${SRR}*_fastqc* $TMPFQ

##########################################################################
echo Starting colorspace alignment for $FQ when resources become available
##########################################################################
$SUBJUNC -T $THREADS --BAMoutput -i $CSGENOME -r $FQ \
-o ${FQ}.subjunc.unsorted.bam 2>>${SRR}.log
rm $FQ ; rm ${SRR}*fasta ${SRR}*qual

##########################################################################
echo Sorting, index bam file and generate stats file for $SRR
##########################################################################
BAM=${FQ}.subjunc.unsorted.bam ; PFX=`echo $BAM | sed 's/.unsorted.bam//'`
$SAMTOOLS sort -@$THREADS $BAM $PFX 2>>$SRR.log ; rm $BAM ; BAM=${PFX}.bam

##########################################################################
#End the platform specific pipeline part
##########################################################################
fi
$SAMTOOLS index $BAM
NUMRDS=`$SAMTOOLS flagstat $BAM | tee -a ${SRR}.log | awk 'NR==3{print $1}'`
echo $SRR mappedreads $NUMRDS | tee -a $SRR.log
if [ "$NUMRDS" -eq 0 ] ; then
echo $SRR has no mapped reads - aborting | tee -a $ATTEMPTS ; exit
fi

##########################################################################
echo Read counting $SRR using featurecounts
##########################################################################
$FEATURECOUNTS -Q 10 -s 0 -T 8 -a $GTF -o ${SRR}_gene.cnt ${BAM} 2>>${SRR}.log
$FEATURECOUNTS -f -Q 10 -s 0 -T 8 -a $GTF -o ${SRR}_exon.cnt ${BAM} 2>>${SRR}.log
pbzip2 -f ${SRR}_gene.cnt ${SRR}_exon.cnt ; rm $BAM ${BAM}.bai

##########################################################################
#Test whether the pipeline worked
##########################################################################
BAMFAIL=`grep -cx Status ${SRR}_exon.cnt.summary`
if [ "$BAMFAIL" -eq "1" ]
then echo $SRR read counting failed - aborting | tee -a $SRR.log
fi
CNTFAIL=`grep -cxP 'Assigned\t0' ${SRR}_exon.cnt.summary`
if [ "$CNTFAIL" -eq "1" ]
then echo $SRR no assigned reads found - aborting | tee -a $ATTEMPTS | tee -a $SRR.log ; exit
fi
if [ "$BAMFAIL" -eq "0" -a "$CNTFAIL" -eq "0" ]
then echo $SRR pipeline completed successfully | tee -a $ATTEMPTS | tee -a $SRR.log ; exit
fi

##############################################################################
#
#  Run QC analysis of each log file to generate a .qc file
#
##############################################################################

LOG=$SRR.log
CNT_SUMMARY=${1}_gene.cnt.summary

if [ -r $CNT_SUMMARY ]
  then COMPLETE=TRUE
else
  echo $SRR not complete
  echo "SRA_filesize:NA
  SequenceEncoding:NA
  SequenceFormat:NA
  ReadLength:NA
  QualityEncoding:NA
  MeanBaseQual:NA
  NumReads:NA
  UniqMappedReads:NA
  AssignedReads:NA
  UnassignedAmbiguous:NA
  UnassignedNoFeatures:NA
  UniqMappingRate:NA
  AssignmentRate:NA
  QC_SUMMARY:FAIL:incomplete" > $SRR.qc
  exit
fi

echo $LOG start
cat $LOG $CNT_SUMMARY > $LOG.tmp
LOG=$LOG.tmp

FATAL=`grep -c 'FATAL ERROR' $LOG`
if [ "$FATAL" -gt "0" ] ; then
echo $SRR fatal error occurred ; exit
fi

SRA_FILESIZE=`grep SRAfilesize $LOG | cut -d ' ' -f3`
ENCODING=`grep -A2 SRAfilesize $LOG | tail -1 | cut -d ' ' -f3`
FORMAT=`grep -A3 SRAfilesize $LOG | tail -1 | cut -d ' ' -f3`

READLEN=`grep -A100 '>>Sequence Length Distribution' $LOG \
| grep -B100 -m1 '>>END_MODULE' \
| grep ^[0-9] | sort -k2gr | head -1 | cut -f1 | cut -d '-' -f1`

QUALITY_ENCODING=`grep -wm1 ^Encoding $LOG | cut -f2 | tr -d ' '`
MEAN_BASE_Q=`grep -A100 '>>Per base sequence quality' $LOG \
| grep -m1 -B100 '>>END_MODULE' | grep ^[0-9] | cut -f3 | sort -n \
| awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`

if [ $ENCODING = "colorspace" ] ; then
  #COLORSPACE QC
  NUM_READS=`grep Processed $LOG | awk '{print $4}'`
  UNIQ_MAPPED_READS=`grep 'Mapped :' $LOG | awk '{print $4}'`

  else
  #BASESPACE QC
  NUM_READS=`grep -i 'Number of input reads' $LOG | cut -f2`
  UNIQ_MAPPED_READS=`grep 'Uniquely mapped reads number' $LOG | cut -f2`
fi

ASSIGNED_READS=`grep -w ^Assigned $LOG | cut -f2 `
if [ -z "$ASSIGNED_READS" ] ; then ASSIGNED_READS=0 ; fi
UNASSIGNED_AMBIGUITY=`grep -w Unassigned_Ambiguity $LOG | cut -f2`
UNASSIGNED_NOFEATURES=`grep -w Unassigned_NoFeatures $LOG | cut -f2`
UNIQ_MAPPING_RATE=`echo $UNIQ_MAPPED_READS $NUM_READS | awk '{print $1/$2*100}'`
ASSIGNMENT_RATE=`echo $ASSIGNED_READS $NUM_READS | awk '{print $1/$2*100}'`


# Classify QC failures
if [ "$READLEN" -lt 18 ] ; then
  REASON=':Read_length_too_short'
  QC_SUMMARY=FAIL
else
  BASEQ=`echo $MEAN_BASE_Q. | cut -d '.' -f1`
  if [ "$BASEQ" -lt 10 ] ; then
    REASON=':BaseQ_too_poor'
    QC_SUMMARY=FAIL
  else
    MAP_RATE=`echo $UNIQ_MAPPING_RATE. | cut -d '.' -f1`
    if [ "$MAP_RATE" -lt 60 ] ; then
      REASON=':Mapping_rate_too_low'
      QC_SUMMARY=FAIL
    else
      AS_RATE=`echo $ASSIGNMENT_RATE. | cut -d '.' -f1`
      if [ "$AS_RATE" -lt 50 ] ; then
        REASON=':Assignment_rate_too_low'
        QC_SUMMARY=FAIL
      else QC_SUMMARY=PASS
      fi
    fi
  fi
fi

## Output .qc file
echo "SRA_filesize:$SRA_FILESIZE
SequenceEncoding:$ENCODING
SequenceFormat:$FORMAT
ReadLength:$READLEN
QualityEncoding:$QUALITY_ENCODING
MeanBaseQual:$MEAN_BASE_Q
NumReads:$NUM_READS
UniqMappedReads:$UNIQ_MAPPED_READS
AssignedReads:$ASSIGNED_READS
UnassignedAmbiguous:$UNASSIGNED_AMBIGUITY
UnassignedNoFeatures:$UNASSIGNED_NOFEATURES
UniqMappingRate:$UNIQ_MAPPING_RATE%
AssignmentRate:$ASSIGNMENT_RATE%
QC_SUMMARY:${QC_SUMMARY}${REASON}" > $SRR.qc

rm $LOG
echo $LOG successful

