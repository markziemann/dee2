#!/bin/bash
#sra2mx for docker image
#Copyright Mark Ziemann 2015 to 2017 mark.ziemann@gmail.com

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
DEE_DIR=/dee2
cd $DEE_DIR
CODE_DIR=$DEE_DIR/code
PIPELINE=$0
PIPELINE_MD5=$(md5sum $PIPELINE | cut -d ' ' -f1)
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
MEM=$(free | awk '$1 ~ /Mem:/  {print $2-$3}')

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

#check all the reference sequences exist and create if necessary
#REFERENCE SEQ AND ANNOTATIONS
MYREF_DIR=$REF_DIR/$ORG/ensembl/
if [ ! -d $MYREF_DIR ] ; then
  mkdir -p $MYREF_DIR
fi

if [ $ORG == "athaliana" ] ; then
  GTFURL="ftp://ftp.ensemblgenomes.org/pub/release-36/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.36.gtf.gz"
  GDNAURL="ftp://ftp.ensemblgenomes.org/pub/release-36/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz"
  CDNAURL="ftp://ftp.ensemblgenomes.org/pub/release-36/plants/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
  BT2_MD5="b19dc2c188246f6971d1b2386a87299e"
  KAL_MD5="bd508db5be410d08ae9c13f5b4c05353"
  STAR_MD5="1cd7aca6533ceed936d99def34e037e6"
elif [ $ORG == "celegans" ] ; then
  GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.90.gtf.gz"
  GDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz"
  CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz"
  BT2_MD5="53a84e66d63066d0d1aaa98f567f9cb0"
  KAL_MD5="f33aa6faaf2c51ceec84a4b7e66b9388"
  STAR_MD5="23091a11423e7c4b0f234e69934a8281"
elif [ $ORG == "dmelanogaster" ] ; then
  GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.90.gtf.gz"
  GDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz"
  CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz"
  BT2_MD5="6cc6863a80199f7a76bc950a63cca8a4"
  KAL_MD5="086b82dfb67d083d321d8604af4dd5a0"
  STAR_MD5="10808e48408630b2f607a738e37fb2a0"
elif [ $ORG == "drerio" ] ; then
  GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/danio_rerio/Danio_rerio.GRCz10.90.gtf.gz"
  GDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna_sm.toplevel.fa.gz"
  CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/danio_rerio/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz"
  BT2_MD5="8fce1a8287b9f1055825d2557e66f844"
  KAL_MD5="035e2d47bbbf87c41859a200ee185eb5"
  STAR_MD5="2763c8b64543e2ed2381469a5f60b6f8"
elif [ $ORG == "ecoli" ] ; then
  GTFURL="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.36.gtf.gz"
  GDNAURL="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_sm.chromosome.Chromosome.fa.gz"
  CDNAURL="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cdna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa.gz"
  BT2_MD5="9fd53f70df3ba54b851713a514ef3412"
  KAL_MD5="dbc74ab4fa8d55d5f3e88476dc5cc32e"
  STAR_MD5="49dfb0bef4e1c0e34503dc995b9456e5"
elif [ $ORG == "hsapiens" ] ; then
  GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz"
  GDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
  CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
  BT2_MD5="e89b4fc019e93b62c2fca8ec59ed56e2"
  KAL_MD5="3f10ef8e78f2cee6df35d2f679ba1c53"
  STAR_MD5="e362a230513c0807451be52fca93cb8f"
elif [ $ORG == "mmusculus" ] ; then
  GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz"
  GDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
  CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
  BT2_MD5="998fb81262af40202fd0fc053e31f5e2"
  KAL_MD5="6b76f3f4fd724764644128b012808982"
  STAR_MD5="9bd0c660a2d876750733bb6a02f9b7df"
elif [ $ORG == "rnorvegicus" ] ; then
  GTFURL="ftp://ftp.ensembl.org/pub/release-90/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf.gz"
  GDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna_sm.toplevel.fa.gz"
  CDNAURL="ftp://ftp.ensembl.org/pub/release-90/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz"
  BT2_MD5="f79c738fccec2a3febb13b2a5f15dc1b"
  KAL_MD5="6a672fbd88df5c18c5d0622676e8b74c"
  STAR_MD5="6b34404a61b5de697a51fe5f534370f7"
elif [ $ORG == "scerevisiae" ] ; then
  GTFURL="ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.36.gtf.gz"
  GDNAURL="ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
  CDNAURL="ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
  BT2_MD5="d484b1bf9e98c0e93c1f7ec37b5d449e"
  KAL_MD5="eb4dd17423dc9644dbfb6daefb9130d0"
  STAR_MD5="052b2523d3e0912bb19de6739bd1d6ed"
fi

# download the necessary reference files
GTF=$MYREF_DIR/$(basename $GTFURL .gz)
if [ -z $GTF ] || [ ! -r $GTF  ] ; then
  cd $MYREF_DIR
  if [ -r $(basename $GTFURL) ] ; then rm $(basename $GTFURL) ; fi
  wget $GTFURL
  gunzip -f $(basename $GTFURL)
  GTF=$MYREF_DIR/$(basename $GTFURL .gz)
  grep -cw gene $GTF > $GTF.cnt
  cd -
fi

GDNA=$MYREF_DIR/$(basename $GDNAURL .gz)
if [ -z $GDNA ] || [ ! -r $GDNA  ] ; then
  cd $MYREF_DIR
  if [ -r $(basename $GDNAURL) ] ; then rm $(basename $GDNAURL) ; fi
  wget $GDNAURL
  gunzip -f $(basename $GDNAURL)
  GDNA=$MYREF_DIR/$(basename $GDNAURL .gz)
  cd -
fi

CDNA=$MYREF_DIR/$(basename $CDNAURL .gz)
if [ -z $CDNA ] || [ ! -r $CDNA  ] ; then
  cd $MYREF_DIR
  if [ -r $(basename $CDNAURL) ] ; then rm $(basename $CDNAURL) ; fi
  wget $CDNAURL
  gunzip -f $(basename $CDNAURL)
  CDNA=$MYREF_DIR/$(basename $CDNAURL .gz)
  grep -c '>' $CDNA > $CDNA.cnt
  cd -
fi

# setup the necessary genome transcriptome indexes
BT2_DIR=$MYREF_DIR/bowtie2
if [ ! -d $BT2_DIR ] ; then
  mkdir -p $BT2_DIR
fi

BT2_REF=$BT2_DIR/$(basename $CDNA)
if [ -z $BT2_REF ] || [ ! -r $BT2_REF  ] ; then
  cd $BT2_DIR ; ln $CDNA .
  #creating bowtie2 index
  bowtie2-build --quiet --threads $THREADS -f $(basename $CDNA) $(basename $CDNA)
  ENS_REFT_BT2=$BT2_DIR/$(basename $CDNA)
  MY_BT2_MD5=$(md5sum $(ls *bt2 | head -1) | awk '{print $1}')
  if [ $MY_BT2_MD5 != $BT2_MD5 ] ; then
    echo "Error in bowtie2 index found. quitting."
    echo "Solution: Try deleting and reindexing the ref transcriptome."
    exit1 ; return 1
  fi
  cd -
fi

KAL_DIR=$MYREF_DIR/kallisto
if [ ! -d $KAL_DIR ] ; then
  mkdir -p $KAL_DIR
fi

KAL_REF=$KAL_DIR/$(basename $CDNA).idx
if [ -z $KAL_REF ] || [ ! -r $KAL_REF  ] ; then
  cd $KAL_DIR
  #kallisto index here
  ln $CDNA .
#  for KMER in `seq 11 2 29` ; do
#    kallisto index -i $(basename $CDNA).k$KMER.idx -k $KMER $(basename $CDNA)
#  done
  kallisto index -i $(basename $CDNA).idx $(basename $CDNA)
  for IDX in *idx ; do grep -c '>' $(basename $CDNA) > $IDX.cnt ; done
  KAL_REF=$KAL_DIR/$(basename $CDNA).idx
  MY_KAL_MD5=$(md5sum $(ls *idx | head -1) | awk '{print $1}')
  if [ $MY_KAL_MD5 != $KAL_MD5 ] ; then
    echo "Error in kallisto index found. quitting."
    echo "Solution: Try deleting and reindexing the ref transcriptome."
    exit1 ; return 1
  fi
  cd -
fi

STAR_DIR=$MYREF_DIR/star
if [ ! -d $STAR_DIR ] ; then
  mkdir -p $STAR_DIR
fi

if [ ! -r $STAR_DIR/SA ] || [ ! -r $STAR_DIR/SAindex ] ; then
  echo Creating star index
  cd $STAR_DIR
  ln $GDNA $GTF .
  CWD=`pwd`
  STAR --runMode genomeGenerate \
  --sjdbGTFfile $CWD/$(basename $GTF) \
  --genomeDir $CWD  \
  --genomeFastaFiles $CWD/$(basename $GDNA) \
  --runThreadN $THREADS
  MY_STAR_MD5=$(md5sum SAindex | awk '{print $1}')
  if [ $MY_STAR_MD5 != $STAR_MD5 ] ; then
    echo "Error in STAR index found. quitting."
    echo "Solution: Try deleting and reindexing the ref genome."
    exit1 ; return 1
  fi
  cd -
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
    ascp -l 500m -O 33001 -T -i $ID $URL . \
    || ( echo $SRR failed ascp download | tee -a $SRR.log ; sleep 5 ; exit1 ; return 1 )
    SRASIZE=$(du ${SRR}.sra)
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

##########################################################################
  echo $SRR diagnose basespace colorspace, single/paired-end and read length
##########################################################################
  fastq-dump -X 4000 --split-files $SRR.sra
  NUM_FQ=$(ls | grep $SRR | grep -v trimmed.fastq | grep -c fastq$)
  if [ $NUM_FQ -eq "1" ] ; then
    ORIG_RDS=SE
    RDS=SE
    echo $SRR is single end | tee -a $SRR.log
  elif [ $NUM_FQ -eq "2" ] ; then
    ORIG_RDS=PE
    RDS=PE
    echo $SRR is paired end | tee -a $SRR.log
  else
    echo Unable to determine if paired or single end. Quitting. | tee -a $SRR.log
    exit1 ; return 1
  fi

  FQ1=$(ls  | grep $SRR | grep -m1 fastq$)
  echo ; echo Starting FastQC analysis of $FQ1
  fastqc $FQ1
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
    exit1 ; return 1
  fi

  #quality encoding ie Illumina1.9
  QUALITY_ENCODING=$(unzip -p ${FQ1BASE}_fastqc ${FQ1BASE}_fastqc/fastqc_data.txt \
  | grep -wm1 ^Encoding | cut -f2 | tr -d ' ')

  #diagnose read length then
  #save entire fastq data to log and delete fastqc zip file and html report
  FQ1_LEN=$(unzip -p ${FQ1BASE}_fastqc.zip ${FQ1BASE}_fastqc/fastqc_data.txt \
  | grep 'Sequence length' | cut -f2)
  echo $SRR read1 length is $FQ1_LEN nt | tee -a $SRR.log
  unzip -p ${FQ1BASE}_fastqc.zip ${FQ1BASE}_fastqc/fastqc_data.txt >> $SRR.log
  rm ${FQ1BASE}_fastqc.zip ${FQ1BASE}_fastqc.html

  FQ1_MIN_LEN=$(sed -n '2~4p' $FQ1 | awk '{print length($1)}' | sort -g | head -1)
  FQ1_MEDIAN_LEN=$(sed -n '2~4p' $FQ1 | awk '{print length($1)}' | numaverage -M)
  FQ1_MAX_LEN=$(sed -n '2~4p' $FQ1 | awk '{print length($1)}' | sort -gr | head -1)

  FQ2_MIN_LEN=NULL
  FQ2_MEDIAN_LEN=NULL
  FQ2_MAX_LEN=NULL

  if [ $RDS == "PE" ] ; then
    FQ2=$(ls  | grep $SRR | grep fastq$ | sed -n 2p)
    echo ; echo Starting FastQC analysis of $FQ2
    fastqc $FQ2
    FQ2BASE=$(basename $FQ2 .fastq)
    FQ2_LEN=$(unzip -p ${FQ2BASE}_fastqc.zip ${FQ2BASE}_fastqc/fastqc_data.txt \
    | grep 'Sequence length' | cut -f2)
    echo $SRR read2 length is $FQ2_LEN nt | tee -a $SRR.log
    unzip -p ${FQ2BASE}_fastqc.zip ${FQ2BASE}_fastqc/fastqc_data.txt >> $SRR.log
    rm ${FQ2BASE}_fastqc.zip ${FQ2BASE}_fastqc.html

    FQ2_MIN_LEN=$(sed -n '2~4p' $FQ2 | awk '{print length($1)}' | sort -g | head -1)
    FQ2_MEDIAN_LEN=$(sed -n '2~4p' $FQ2 | awk '{print length($1)}' | numaverage -M)
    FQ2_MAX_LEN=$(sed -n '2~4p' $FQ2 | awk '{print length($1)}' | sort -gr | head -1)

    #now checking read lengths and dropping ones too short
    if [[ $FQ1_MAX_LEN -lt 20 && $FQ2_MAX_LEN -lt 20 ]] ; then
      echo Read lengths are too short. Quitting. | tee -a $SRR.log
      exit1 ; return 1
    fi
  fi

##########################################################################
#  If colorspace then quit
##########################################################################
  if [ $CSPACE == "TRUE" ] ; then
    echo Colorspace data is excluded from analysis for now
    exit1 ; return 1
  fi

##########################################################################
  echo $SRR Dump the fastq file
##########################################################################
  rm ${SRR}*fastq
  if [ $CSPACE == "FALSE" ] ; then
    parallel-fastq-dump --threads $THREADS --outdir . --split-files --defline-qual + -s ${SRR}.sra >> $SRR.log 2>&1
  fi

  if [ $RDS == "PE" ] ; then
    if [[ $FQ1_MAX_LEN -ge 20 && $FQ2_MAX_LEN -lt 20 ]] ; then
      echo Read 2 is too short. Omitting from downstream analysis. | tee -a $SRR.log
      rm $FQ2
      FQ1NEW=$(echo $FQ1 | sed 's/_1//')
      mv $FQ1 $FQ1NEW
      FQ1=$FQ1NEW
      RDS=SE
    fi

    if [[ $FQ1_MAX_LEN -lt 20 && $FQ2_MAX_LEN -ge 20 ]] ; then
      echo Read 1 is too short. Omitting from downstream analysis. | tee -a $SRR.log
      rm $FQ1
      FQ1=$(echo $FQ2 | sed 's/_2//')
      mv $FQ2 $FQ1
      RDS=SE
    fi
  fi

  FILESIZE=$(du -s $FQ1 | cut -f1)
  FILESIZE="${FILESIZE:-0}"
  echo $SRR file size $FILESIZE | tee -a $SRR.log
  rm ${SRR}.sra

  if [ "$FILESIZE" -eq 0 ] ; then
    echo $SRR has no reads. Aborting | tee -a $ATTEMPTS ; rm $FQ1 ; exit1 ; return 1
  fi

  echo $SRR completed basic pipeline successfully | tee -a $SRR.log

##########################################################################
# If user has own data
##########################################################################
else
  echo using user data
# then check whether zipped (gzip bzip2) then if unzip if neccessary
# then run fastqc and save it to the log
# if R1 and R2 have a different number of reads, then throw error
# all the above stuff needs to be enclosed in an if statement
# the '-f' switch to specify own data

# The following syntax to support multiple paired end fastq files
# pipeline.sh hsapiens -f sample1_R1.fq.gz,sample2_R1.fq sample1_R2.fq.gz,sample2_R2.fq
# the function is run like this for SE and PE
# main -f $FQ_R1
# main -f $FQ_R1 $FQ_R2

##########################################################################
# OWN data SE
##########################################################################
  if [ $NUMVARS -lt "4" ] ; then
    echo Error, not enough arguments specified. Quitting
    exit1; return 1
  elif [ $NUMVARS -gt "5" ] ; then
    echo Too many arguments. Check syntax and try again. Ensure filenames do not contain spaces. Quitting
    exit1 ; return 1
  elif [ $NUMVARS -eq "4" ] ; then
    echo Processing fastq file $3 in single end mode
    FQ1=$3
    SRR=$(basename $FQ1 | cut -d '.' -f1)
    RDS=SE

    if [ ! -r $FQ1 ] ; then
      Input file $FQ1 does not exist or is not readable. Quitting.
      exit1 ; return 1
    fi

    mkdir $SRR ; cp $PIPELINE $SRR ; cd $SRR
    cp $FQ1 . ; FQ1=$(basename $FQ1)

    echo "Starting $PIPELINE $SRR
      current disk space = $DISK
      free memory = $MEM " | tee -a $SRR.log

    ISGZ=$(echo $FQ1 | grep -c .gz$)
    if [ $ISGZ -eq "1" ] ; then
      pigz -t $FQ1 && GZTEST=OK
      if [ $GZTEST == OK ] ; then
        pigz -d $FQ1
        FQ1=$(basename $FQ1 .gz)
      else
        echo Gzip file is corrupted. Quitting
        exit1 ; return 1
      fi
    fi

    ISBZ=$(echo $FQ1 | grep -c .bz2$)
    if [ $ISBZ -eq "1" ] ; then
      pbzip2 -t $FQ1 && BZTEST=OK
      if [ $BZTEST == OK ] ; then
        pbzip2 -d $FQ1
        FQ1=$(basename $FQ1 .bz2)
      else
        echo Bzip2 file is corrupted. Quitting
        exit1 ; return 1
      fi
    fi

    ISFQ=$(echo $FQ1 | egrep -c '(.fq$|.fastq$)' )
    if [ $ISFQ -ne "1" ] ; then
      echo Error. Unknown input file format. Input file extension should match ".fastq" or ".fq". Quitting
      exit1 ; return 1
    fi

##########################################################################
# OWN data PE
##########################################################################
  elif [ $NUMVARS -eq "5" ] ; then
    echo Processing fastq files $3 and $4 in paired end mode
    FQ1=$3
    FQ2=$4
    SRR=$(basename $FQ1 | cut -d '.' -f1)
    RDS=PE

    if [ ! -r $FQ1 ] ; then
      Input file $FQ1 does not exist or is not readable. Quitting.
      exit1 ; return 1
    fi

    if [ ! -r $FQ2 ] ; then
      Input file $FQ2 does not exist or is not readable. Quitting.
      exit1 ; return 1
    fi

    mkdir $SRR ; cp $PIPELINE $SRR ; cd $SRR
    cp $FQ1 $FQ2 .
    FQ1=$(basename $FQ1)
    FQ2=$(basename $FQ2)

    echo "Starting $PIPELINE $SRR
      current disk space = $DISK
      free memory = $MEM " | tee -a $SRR.log

    ISGZ=$(echo $FQ1 $FQ2 | tr ' ' '\n' | grep -c .gz$)
    if [ $ISGZ -eq "2" ] ; then
      pigz -t $FQ1 && pigz -t $FQ2 && GZTEST=OK
      if [ $GZTEST == OK ] ; then
        pigz -d $FQ1 $FQ2
        FQ1=$(basename $FQ1 .gz)
        FQ2=$(basename $FQ2 .gz)
      else
        echo Gzip file is corrupted. Quitting
        exit1 ; return 1
      fi
    fi

    ISBZ=$(echo $FQ1 $FQ2 | tr ' ' '\n' | grep -c .bz2$)
    if [ $ISBZ -eq "2" ] ; then
      pbzip2 -t $FQ1 && pbzip2 -t $FQ2 && BZTEST=OK
      if [ $BZTEST == OK ] ; then
        pbzip2 -d $FQ1 $FQ2
        FQ1=$(basename $FQ1 .gz)
        FQ2=$(basename $FQ2 .gz)
      else
        echo Bzip2 file is corrupted. Quitting
        exit1 ; return 1
      fi
    fi

    ISFQ=$(echo $FQ1 $FQ2 | tr ' ' '\n' | egrep -c '(.fq$|.fastq$)' )
    if [ $ISFQ -ne "2" ] ; then
      echo Error. Unknown input file format. Input file extension should match ".fastq" or ".fq"$
      exit1 ; return 1
    fi

  fi

  fastqc $FQ1
  FQ1BASE=$(echo $FQ1 | rev | cut -d '.' -f2 | rev)
  SRR=$FQ1BASE

  #diagnose colorspace or conventional
  BASECALL_ENCODING=$(unzip -p ${FQ1BASE}_fastqc ${FQ1BASE}_fastqc/fastqc_data.txt \
  | grep 'File type' | cut -f2 | awk '{print $1}')

  if [ $BASECALL_ENCODING == "Colorspace" ] ; then
    CSPACE=TRUE
    echo $FQ1 is colorspace | tee -a $SRR.log
  elif [ $BASECALL_ENCODING == "Conventional" ] ; then
    CSPACE=FALSE
    echo $FQ1 is conventional basespace | tee -a $SRR.log
  else
    echo Unable to determine if colorspace or basespace. Quitting. | tee -a $SRR.log
    exit1 ; return 1
  fi

  #quality encoding ie Illumina1.9
  QUALITY_ENCODING=$(unzip -p ${FQ1BASE}_fastqc ${FQ1BASE}_fastqc/fastqc_data.txt \
  | grep -wm1 ^Encoding | cut -f2 | tr -d ' ')

  #diagnose read length then
  #save entire fastq data to log and delete fastqc zip file and html report
  FQ1_LEN=$(unzip -p ${FQ1BASE}_fastqc.zip ${FQ1BASE}_fastqc/fastqc_data.txt \
  | grep 'Sequence length' | cut -f2)
  echo $FQ1 read1 length is $FQ1_LEN nt | tee -a $SRR.log
  unzip -p ${FQ1BASE}_fastqc.zip ${FQ1BASE}_fastqc/fastqc_data.txt | tee -a $SRR.log
  rm ${FQ1BASE}_fastqc.zip ${FQ1BASE}_fastqc.html

  FQ1_MIN_LEN=$(sed -n '2~4p' $FQ1 | awk '{print length($1)}' | sort -g | head -1)
  FQ1_MEDIAN_LEN=$(sed -n '2~4p' $FQ1 | awk '{print length($1)}' | numaverage -M)
  FQ1_MAX_LEN=$(sed -n '2~4p' $FQ1 | awk '{print length($1)}' | sort -gr | head -1)

  FQ2_MIN_LEN=NULL
  FQ2_MEDIAN_LEN=NULL
  FQ2_MAX_LEN=NULL


  if [ $RDS == "PE" ] ; then
    fastqc $FQ2
    FQ2BASE=$(echo $FQ2 | rev | cut -d '.' -f2 | rev)
    FQ2_LEN=$(unzip -p ${FQ2BASE}_fastqc.zip ${FQ2BASE}_fastqc/fastqc_data.txt \
    | grep 'Sequence length' | cut -f2)
    echo $SRR read2 length is $FQ2_LEN nt | tee -a $SRR.log
    unzip -p ${FQ2BASE}_fastqc.zip ${FQ2BASE}_fastqc/fastqc_data.txt | tee -a $SRR.log
    rm ${FQ2BASE}_fastqc.zip ${FQ2BASE}_fastqc.html

    FQ2_MIN_LEN=$(sed -n '2~4p' $FQ2 | awk '{print length($1)}' | sort -g | head -1)
    FQ2_MEDIAN_LEN=$(sed -n '2~4p' $FQ2 | awk '{print length($1)}' | numaverage -M)
    FQ2_MAX_LEN=$(sed -n '2~4p' $FQ2 | awk '{print length($1)}' | sort -gr | head -1)

    #now checking read lengths and dropping ones too short
    if [[ $FQ1_MAX_LEN -lt 20 && $FQ2_MAX_LEN -lt 20 ]] ; then
      echo Read lengths are too short. Quitting. | tee -a $SRR.log
      exit1 ; return 1
    fi

##########################################################################
    echo If read 1 and 2 have different number of tags then exit
##########################################################################
    FQ1_NUMRDS=$(sed -n '2~4p' $FQ1 | wc -l)
    FQ2_NUMRDS=$(sed -n '2~4p' $FQ2 | wc -l)

    if [ $FQ1_NUMRDS -ne $FQ2_NUMRDS ] ; then
      echo Number of sequence tags in read 1 and read 2 differ. Quitting.
      exit1 ; return 1
    fi

  fi

##########################################################################
  echo $FQ1 if colorspace, then quit
##########################################################################
  if [ $CSPACE == "TRUE" ] ; then
    echo Colorspace data is excluded from analysis for now
    exit1 ; return 1
  fi

fi

##########################################################################
echo $SRR Quality trimming
##########################################################################
if [ $RDS == "SE" ] ; then
  skewer -l 18 -q 10 -k inf -t $THREADS -o $SRR $FQ1
  rm $FQ1
  FQ1=${SRR}-trimmed.fastq

elif [ $RDS == "PE" ] ; then
  skewer -l 18 -q 10 -k inf -t $THREADS -o $SRR $FQ1 $FQ2
  rm $FQ1 $FQ2
  FQ1=${SRR}-trimmed-pair1.fastq
  FQ2=${SRR}-trimmed-pair2.fastq
fi

# check to see that the skewer log was created - if not then there is a problem
if [ ! -f ${SRR}-trimmed.log ] ; then
  echo Skewer failed. Quitting | tee -a $SRR.log
  exit1 ; return 1
fi

# get read counts and append skewer log and exit if there are no reads passing QC
READ_CNT_TOTAL=$(grep 'processed; of these:' ${SRR}-trimmed.log | awk '{print $1}')
READ_CNT_AVAIL=$(grep 'available; of these:' ${SRR}-trimmed.log | awk '{print $1}')

if [ -z "$READ_CNT_AVAIL" ] ; then READ_CNT_AVAIL=0 ; fi

cat ${SRR}-trimmed.log >> $SRR.log && rm ${SRR}-trimmed.log
if [ $READ_CNT_AVAIL -eq "0" ] ; then
  echo No reads passed QC. Quitting | tee -a $SRR.log
  exit1 ; return 1
else
  echo $READ_CNT_AVAIL reads passed initial QC | tee -a $SRR.log
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
  minion search-adapter -i $FQ1 > $MINION_LOG
  ADAPTER=$(head $MINION_LOG | grep -m1 sequence= | cut -d '=' -f2)
  DENSITY=$(head $MINION_LOG | grep -m1 'sequence-density=' | cut -d '=' -f2 | numround -c)
  cat $MINION_LOG | tee -a $SRR.log && rm $MINION_LOG

  if [[ ! -z $DENSITY ]] ; then
    if [ $DENSITY -gt $ADAPTER_THRESHOLD ] ; then
      echo Potential 3prime adapter identified. Now checking if in reference sequence | tee -a $SRR.log
      # Query to see if adapter sequence present in reference
      ADAPTER_REF_CHECK=$(bowtie2 -f -x $BT2_REF -S /dev/stdout <(echo $ADAPTER | sed 's/^/>ADAPTER\n/') 2>>$SRR.log | awk '$1!~/^@/ && $2!=4' | wc -l )

      if [ $ADAPTER_REF_CHECK -eq "0" ] ; then
        echo Adapter seq not found in reference. Now shuffling file before clipping | tee -a $SRR.log
        paste - - - - < $FQ1 \
        | unsort --seed 42 \
        | tr '\t' '\n' > ${SRR}.fastq
        CLIP_LINE_NUM=$(echo $DENSITY $ADAPTER_THRESHOLD $READ_CNT_AVAIL | awk '{printf "%.0f\n", ($1-$2)/$1*$3*4}' | numround -n 4)
        head -$CLIP_LINE_NUM ${SRR}.fastq | skewer -l 18 -t $THREADS -x $ADAPTER -o $SRR -
        cat ${SRR}-trimmed.log >> $SRR.log && rm ${SRR}-trimmed.log
        CLIP_LINE_NUM1=$((CLIP_LINE_NUM+1))
        tail -n+$CLIP_LINE_NUM1 ${SRR}.fastq >> ${SRR}-trimmed.fastq && rm ${SRR}.fastq
        READ_CNT_AVAIL=$(sed -n '2~4p' ${SRR}-trimmed.fastq | wc -l)
        if [ -z "$READ_CNT_AVAIL" ] ; then READ_CNT_AVAIL=0 ; fi
        minion search-adapter -i ${SRR}-trimmed.fastq | tee -a $SRR.log
      else
        echo Potential adapter found in reference sequence. Continuing without clipping. | tee -a $SRR.log
      fi
    fi
  fi

elif [ $RDS == "PE" ] ; then
  MINION_LOG=$FQ1.minion.log
  minion search-adapter -i $FQ1 > $MINION_LOG
  ADAPTER1=$(head $MINION_LOG | grep -m1 sequence= | cut -d '=' -f2)
  DENSITY1=$(head $MINION_LOG | grep -m1 'sequence-density=' | cut -d '=' -f2 | numround -c)
  cat $MINION_LOG | tee -a $SRR.log && rm $MINION_LOG

  MINION_LOG=$FQ2.minion.log
  minion search-adapter -i $FQ2 > $MINION_LOG
  ADAPTER2=$(head $MINION_LOG | grep -m1 sequence= | cut -d '=' -f2)
  DENSITY2=$(head $MINION_LOG | grep -m1 'sequence-density=' | cut -d '=' -f2 | numround -c)
  cat $MINION_LOG | tee -a $SRR.log && rm $MINION_LOG

  DENSITY=$(echo $DENSITY1 $DENSITY2 | awk '{print ($1+$2)/2}' | numround)

  if [[ ! -z $DENSITY ]] ; then
    if [ $DENSITY -gt $ADAPTER_THRESHOLD ] ; then
      echo Potential 3prime adapter identified. Now checking if in reference sequence | tee -a $SRR.log
      # Query to see if adapter sequence present in reference
      ADAPTER1_REF_CHECK=$($BOWTIE2 -f -x $BT2_REF -S /dev/stdout <(echo $ADAPTER1 | sed 's/^/>ADAPTER\n/') 2>>$SRR.log | awk '$1!~/^@/ && $2!=4' | wc -l )
      ADAPTER2_REF_CHECK=$($BOWTIE2 -f -x $BT2_REF -S /dev/stdout <(echo $ADAPTER2 | sed 's/^/>ADAPTER\n/') 2>>$SRR.log | awk '$1!~/^@/ && $2!=4' | wc -l )

      if [ $ADAPTER1_REF_CHECK -eq "0" -a $ADAPTER2_REF_CHECK -eq "0" ] ; then
        echo Adapter seq not found in reference. Now shuffling file before clipping | tee -a $SRR.log
        paste <(cut -d ' ' -f1 $FQ1) <(cut -d ' ' -f1 $FQ2) | paste - - - - \
        | unsort --seed 42 \
        | awk -F'\t' '{OFS="\n"; print $1,$3,$5,$7 > "R1.fastq"; print $2,$4,$6,$8 > "R2.fastq"}'
        mv R1.fastq ${SRR}_1.tmp.fastq ; mv R2.fastq ${SRR}_2.tmp.fastq
        CLIP_LINE_NUM=$(echo $DENSITY $ADAPTER_THRESHOLD $READ_CNT_AVAIL | awk '{printf "%.0f\n", ($1-$2)/$1*$3*4}' | numround -n 4)
        head -$CLIP_LINE_NUM ${SRR}_1.tmp.fastq > ${SRR}_1.fastq &
        head -$CLIP_LINE_NUM ${SRR}_2.tmp.fastq > ${SRR}_2.fastq ; wait
        skewer -l 18 -t $THREADS -x $ADAPTER1 -y $ADAPTER2 -o $SRR ${SRR}_1.fastq ${SRR}_2.fastq
        READ_CNT_AVAIL=$(grep 'available; of these:' ${SRR}-trimmed.log | cut -d ' ' -f1)
        if [ -z "$READ_CNT_AVAIL" ] ; then READ_CNT_AVAIL=0 ; fi
        cat ${SRR}-trimmed.log >> $SRR.log && rm ${SRR}-trimmed.log
        CLIP_LINE_NUM1=$((CLIP_LINE_NUM+1))
        tail -n+$CLIP_LINE_NUM1 ${SRR}_1.tmp.fastq >> ${SRR}-trimmed-pair1.fastq && rm ${SRR}_1.tmp.fastq
        tail -n+$CLIP_LINE_NUM1 ${SRR}_2.tmp.fastq >> ${SRR}-trimmed-pair2.fastq && rm ${SRR}_2.tmp.fastq
        minion search-adapter -i ${SRR}-trimmed-pair1.fastq | tee -a $SRR.log
        minion search-adapter -i ${SRR}-trimmed-pair2.fastq | tee -a $SRR.log
      else
        echo Potential adapter found in reference sequence. Continuing without clipping. | tee -a $SRR.log
      fi
    fi
  fi
fi

#QcPassRate
QC_PASS_RATE=$(echo $READ_CNT_AVAIL $READ_CNT_TOTAL | awk '{print $1/$2*100"%"}')

FQSIZE=$(du -s $FQ1 | cut -f1)
#cat ${SRR}-trimmed.log >> $SRR.log && rm ${SRR}-trimmed.log
if [ $FQSIZE -eq "0" ] ; then
  echo No reads passed QC. Quitting | tee -a $SRR.log
  exit1 ; return 1
fi

##########################################################################
echo $SRR Starting mapping phase
##########################################################################
# need to setup 2 alignments
# -kallisto ensembl
# -star/featurecounts ensembl

#STAR Ensembl
# runs star in gene-wise quant mode and avoid samtools and htseq/featurecounts.
# quant mode can be used with sheared memory if indexed with GTF file
# quantMode also looks like a nice way to diagnose the mapping strand
# Strand information could then be used for kallisto options
# The ReadsPerGene.out.tab file contains count information

# Before running the full mapping procedure for paired end reads, a sample of
# 10000 forward and reverse reads first.

if [ $RDS == "PE" ] ; then
  echo $SRR testing PE reads STAR mapping to Ensembl genome | tee -a $SRR.log
  head $FQ1 $FQ2
  tail $FQ1 $FQ2

  #test 100k FQ1 and FQ2 reads clipped on the 5 prime end to exclude UMIs and barcodes
  head -10000 $FQ1 > test_R1.fq ; head -1000000 $FQ1 | tail -90000 >> test_R1.fq
  head -10000 $FQ2 > test_R2.fq ; head -1000000 $FQ2 | tail -90000 >> test_R2.fq

#  skewer -m ap --cut 4,4 -l 18 -k inf -t $THREADS test_R1.fq test_R2.fq && mv test_R1-trimmed-pair1.fastq test_R1_clip4.fq && mv test_R1-trimmed-pair2.fastq test_R2_clip4.fq
#  skewer -m ap --cut 8,8 -l 18 -k inf -t $THREADS test_R1.fq test_R2.fq && mv test_R1-trimmed-pair1.fastq test_R1_clip8.fq && mv test_R1-trimmed-pair2.fastq test_R2_clip8.fq
#  skewer -m ap --cut 12,12 -l 18 -k inf -t $THREADS test_R1.fq test_R2.fq && mv test_R1-trimmed-pair1.fastq test_R1_clip12.fq && mv test_R1-trimmed-pair2.fastq test_R2_clip12.fq
#  skewer -m ap --cut 20,20 -l 18 -k inf -t $THREADS test_R1.fq test_R2.fq && mv test_R1-trimmed-pair1.fastq test_R1_clip12.fq && mv test_R1-trimmed-pair2.fastq test_R2_clip12.fq

  fastx_trimmer -f 5 -m 18 -Q 33 -i test_R1.fq > test_R1_clip4.fq &
  fastx_trimmer -f 5 -m 18 -Q 33 -i test_R2.fq > test_R2_clip4.fq &
  fastx_trimmer -f 9 -m 18 -Q 33 -i test_R1.fq > test_R1_clip8.fq &
  fastx_trimmer -f 9 -m 18 -Q 33 -i test_R2.fq > test_R2_clip8.fq &
  fastx_trimmer -f 13 -m 18 -Q 33 -i test_R1.fq > test_R1_clip12.fq &
  fastx_trimmer -f 13 -m 18 -Q 33 -i test_R2.fq > test_R2_clip12.fq &
  fastx_trimmer -f 21 -m 18 -Q 33 -i test_R1.fq > test_R1_clip20.fq &
  fastx_trimmer -f 21 -m 18 -Q 33 -i test_R2.fq > test_R2_clip20.fq &
  wait

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_R1.fq >/dev/null 2>&1

  R1_RD_CNT=$(sed -n '2~4p' < test_R1.fq | wc -l)
  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  UNMAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | head -1)
  R1_MAP_RATE=$(echo $MAPPED_CNT $R1_RD_CNT | awk '{print $1/$2*100}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_R2.fq  >/dev/null 2>&1

  R2_RD_CNT=$(sed -n '2~4p' < test_R2.fq | wc -l)
  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  UNMAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | head -1)
  R2_MAP_RATE=$(echo $MAPPED_CNT $R2_RD_CNT | awk '{print $1/$2*100}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_R1_clip4.fq  >/dev/null 2>&1

  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R1_MAP_RATE_CLIP4=$(echo $MAPPED_CNT $R1_RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_R2_clip4.fq  >/dev/null 2>&1

  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R2_MAP_RATE_CLIP4=$(echo $MAPPED_CNT $R2_RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_R1_clip8.fq  >/dev/null 2>&1

  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R1_MAP_RATE_CLIP8=$(echo $MAPPED_CNT $R1_RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_R2_clip8.fq  >/dev/null 2>&1

  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R2_MAP_RATE_CLIP8=$(echo $MAPPED_CNT $R2_RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_R1_clip12.fq  >/dev/null 2>&1

  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R1_MAP_RATE_CLIP12=$(echo $MAPPED_CNT $R1_RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_R2_clip12.fq  >/dev/null 2>&1

  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R2_MAP_RATE_CLIP12=$(echo $MAPPED_CNT $R2_RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_R1_clip20.fq  >/dev/null 2>&1

  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R1_MAP_RATE_CLIP20=$(echo $MAPPED_CNT $R1_RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_R2_clip20.fq  >/dev/null 2>&1

  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R2_MAP_RATE_CLIP20=$(echo $MAPPED_CNT $R2_RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  rm test_R1.fq test_R2.fq test_R1_clip4.fq test_R2_clip4.fq test_R1_clip8.fq test_R2_clip8.fq test_R1_clip12.fq test_R2_clip12.fq test_R1_clip20.fq test_R2_clip20.fq ReadsPerGene.out.tab

  #NOW logic for clipping entire R1 and R2 dataset
  R1_CLIP_NUM=$(echo $R1_MAP_RATE:0 $R1_MAP_RATE_CLIP4:4 $R1_MAP_RATE_CLIP8:8 $R1_MAP_RATE_CLIP12:12 $R1_MAP_RATE_CLIP20:20\
  | tr ' ' '\n' | sort -gr | head -1 | cut -d ':' -f2)

  R1_MAP_RATE=$(echo $R1_MAP_RATE:0 $R1_MAP_RATE_CLIP4:4 $R1_MAP_RATE_CLIP8:8 $R1_MAP_RATE_CLIP12:12 $R1_MAP_RATE_CLIP20:20 \
  | tr ' ' '\n' | sort -gr | head -1 | cut -d ':' -f1)

  R2_CLIP_NUM=$(echo $R2_MAP_RATE:0 $R2_MAP_RATE_CLIP4:4 $R2_MAP_RATE_CLIP8:8 $R2_MAP_RATE_CLIP12:12 $R2_MAP_RATE_CLIP20:20 \
  | tr ' ' '\n' | sort -gr | head -1 | cut -d ':' -f2)

  R2_MAP_RATE=$(echo $R2_MAP_RATE:0 $R2_MAP_RATE_CLIP4:4 $R2_MAP_RATE_CLIP8:8 $R2_MAP_RATE_CLIP12:12 $R2_MAP_RATE_CLIP20:20 \
  | tr ' ' '\n' | sort -gr | head -1 | cut -d ':' -f1)

  if [[ ( $R1_CLIP_NUM -gt 0 ) || ( $R2_CLIP_NUM -gt 0 ) ]] ; then

    if [[ ( $R1_CLIP_NUM -gt 15 ) || ( $R2_CLIP_NUM -gt 15 ) ]] ; then
      ( fastx_trimmer -f $((R1_CLIP_NUM+1)) -m 18 -Q 33 -i $FQ1 > $FQ1.tmp.fq && mv $FQ1.tmp.fq $FQ1 ) &
      fastx_trimmer -f $((R2_CLIP_NUM+1)) -m 18 -Q 33 -i $FQ2 > $FQ2.tmp.fq && mv $FQ2.tmp.fq $FQ2
      wait
    else
      skewer -m ap --cut $R1_CLIP_NUM,$R2_CLIP_NUM -l 18 -k inf -t $THREADS $FQ1 $FQ2 && \
      mv $(basename $FQ1 .fastq)-trimmed-pair1.fastq $FQ1 && \
      mv $(basename $FQ1 .fastq)-trimmed-pair2.fastq $FQ2 && \
      rm *untrimmed*fastq *trimmed.log
    fi
  fi

  R1R2_DIFF=$((R1_MAP_RATE-R2_MAP_RATE))

  if [ $R2_MAP_RATE -lt "40" -a $R1R2_DIFF -ge "20" ] ; then
    echo Read2 map rate below 40%, dropping it and using Read1 only | tee -a $SRR.log
    DROP_R2=TRUE
    RDS="SE"
    rm $FQ2
  fi

  R2R1_DIFF=$((R2_MAP_RATE-R1_MAP_RATE))

  if [ $R1_MAP_RATE -lt "40" -a $R2R1_DIFF -ge "20" ] ; then
    echo Read1 map rate below 40%, dropping it and using Read2 only | tee -a $SRR.log
    DROP_R1=TRUE
    RDS="SE"
    mv $FQ2 $FQ1
  fi
fi


## Now performing full alignment
if [ $RDS == "SE" ] ; then
  head $FQ1
  tail $FQ1

  #test 100k reads FQ1 and reads clipped on the 5 prime end to exclude UMIs and barcodes
  head -10000 $FQ1 > test.fq
  head -1000000 $FQ1 | tail -90000 >> test.fq
  cp test.fq test2.fq

#  skewer -m ap --cut 4,4 -l 18 -k inf -t $THREADS test.fq test2.fq && mv test-trimmed-pair1.fastq test_clip4.fq
#  skewer -m ap --cut 8,8 -l 18 -k inf -t $THREADS test.fq test2.fq && mv test-trimmed-pair1.fastq test_clip8.fq
#  skewer -m ap --cut 12,12 -l 18 -k inf -t $THREADS test.fq test2.fq && mv test-trimmed-pair1.fastq test_clip12.fq
#  skewer -m ap --cut 20,20 -l 18 -k inf -t $THREADS test.fq test2.fq && mv test-trimmed-pair1.fastq test_clip20.fq
  fastx_trimmer -f 5 -m 18 -Q 33 -i test.fq > test_clip4.fq &
  fastx_trimmer -f 9 -m 18 -Q 33 -i test.fq > test_clip8.fq &
  fastx_trimmer -f 13 -m 18 -Q 33 -i test.fq > test_clip12.fq &
  fastx_trimmer -f 21 -m 18 -Q 33 -i test.fq > test_clip20.fq &
  wait

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test.fq  >/dev/null 2>&1

  RD_CNT=$(sed -n '2~4p' < test.fq | wc -l)
  MAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  UNMAPPED_CNT=$(cut -f2 ReadsPerGene.out.tab | head -1)
  R1_MAP_RATE=$(echo $MAPPED_CNT $RD_CNT | awk '{print $1/$2*100}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_clip4.fq  >/dev/null 2>&1

  MAPPED_CNT_CLIP4=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R1_MAP_RATE_CLIP4=$(echo $MAPPED_CNT_CLIP4 $RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_clip8.fq  >/dev/null 2>&1

  MAPPED_CNT_CLIP8=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R1_MAP_RATE_CLIP8=$(echo $MAPPED_CNT_CLIP8 $RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_clip12.fq  >/dev/null 2>&1

  MAPPED_CNT_CLIP12=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R1_MAP_RATE_CLIP12=$(echo $MAPPED_CNT_CLIP12 $RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=test_clip20.fq  >/dev/null 2>&1

  MAPPED_CNT_CLIP20=$(cut -f2 ReadsPerGene.out.tab | tail -n +3 | numsum)
  R1_MAP_RATE_CLIP20=$(echo $MAPPED_CNT_CLIP12 $RD_CNT | awk '{print ($1/$2*100)-1}' | numround)

  rm test.fq test2.fq test_clip4.fq test_clip8.fq test_clip12.fq test_clip20.fq ReadsPerGene.out.tab

  #NOW logic for clipping entire R1 dataset
  CLIP_NUM=$(echo $R1_MAP_RATE:0 $R1_MAP_RATE_CLIP4:4 $R1_MAP_RATE_CLIP8:8 $R1_MAP_RATE_CLIP12:12 $R1_MAP_RATE_CLIP20:20 \
  | tr ' ' '\n' | sort -gr | head -1 | cut -d ':' -f2)

  R1_MAP_RATE=$(echo $R1_MAP_RATE:0 $R1_MAP_RATE_CLIP4:4 $R1_MAP_RATE_CLIP8:8 $R1_MAP_RATE_CLIP12:12 $R1_MAP_RATE_CLIP20:20 \
  | tr ' ' '\n' | sort -gr | head -1 | cut -d ':' -f1)

  if [[ ( $CLIP_NUM -gt 0 ) && ( $CLIP_NUM -lt 20 ) ]] ; then
    cp $FQ1 $FQ1.tmp.fq
    skewer -m ap --cut $CLIP_NUM,$CLIP_NUM -l 18 -k inf -t $THREADS $FQ1 $FQ1.tmp.fq \
    && mv $(basename $FQ1 .fastq)-trimmed-pair1.fastq $FQ1 \
    && rm $(basename $FQ1 .fastq)-trimmed-pair2.fastq \
    && rm $FQ1.tmp.fq *untrimmed*fastq *trimmed.log
  fi

  if [ $CLIP_NUM -gt 15 ] ; then
    fastx_trimmer -f 21 -m 18 -Q 33 -i $FQ1 > $FQ1.tmp.fq && mv $FQ1.tmp.fq $FQ1
  fi

  # Full SE alignment
  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=$FQ1

elif [ $RDS == "PE" ] ; then
  head $FQ1 $FQ2 ; tail $FQ1 $FQ2
  #proper PE mapping
  STAR --runThreadN $THREADS --quantMode GeneCounts --genomeLoad LoadAndKeep \
  --outSAMtype None --genomeDir $STAR_DIR --readFilesIn=$FQ1 $FQ2
fi

#now grab some qc info from the star alignment for later
UNIQ_MAPPED_READS=$(grep 'Uniquely mapped reads number' Log.final.out | awk '{print $NF}')

cat Log.final.out | tee -a $SRR.log && rm Log.final.out Log.out Log.progress.out SJ.out.tab
head -4 ReadsPerGene.out.tab | tee -a $SRR.log
mv ReadsPerGene.out.tab $SRR.se.tsv

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
  STRANDED=PositiveStrand
  KALLISTO_STRAND_PARAMETER='--fr-stranded'
  echo "Dataset is classified positive stranded" | tee -a $SRR.log
elif [ $NEG_STRAND_CNT -ge "$((POS_STRAND_CNT*5))" ] ; then
  STRAND=2
  STRANDED=NegativeStrand
  KALLISTO_STRAND_PARAMETER='--rf-stranded'
  echo "Dataset is classified negative stranded" | tee -a $SRR.log
else
  STRAND=0
  STRANDED=Unstranded
  KALLISTO_STRAND_PARAMETER=''
  echo "Dataset is classified unstranded" | tee -a $SRR.log
fi
echo KALLISTO_STRAND_PARAMETER=$KALLISTO_STRAND_PARAMETER

#now grab some qc info from the star alignment for later
CUTCOL=$((STRAND+2))
UNMAPPED_CNT=$(cut -f$CUTCOL $SRR.se.tsv | head -1)
MULTIMAPPED_CNT=$(cut -f$CUTCOL $SRR.se.tsv | head -2 | tail -1)
NOFEATURE_CNT=$(cut -f$CUTCOL $SRR.se.tsv | head -3 | tail -1)
AMBIGUOUS_CNT=$(cut -f$CUTCOL $SRR.se.tsv | head -4 | tail -1)
ASSIGNED_CNT=$(cut -f$CUTCOL $SRR.se.tsv | tail -n +5 | numsum)
UNIQ_MAP_RATE=$(echo $UNIQ_MAPPED_READS $READ_CNT_AVAIL | awk '{print $1/$2*100"%"}')
ASSIGNED_RATE=$(echo $ASSIGNED_CNT $READ_CNT_AVAIL | awk '{print $1/$2*100"%"}')

#Now cut out columns to leave us with only the desired strand info
CUTCOL=$((STRAND+2))
cut -f1,$CUTCOL $SRR.se.tsv | tail -n +5 > $SRR.se.tsv.tmp && mv $SRR.se.tsv.tmp $SRR.se.tsv

##########################################################################
echo $SRR checking readlengths now for kmer selection
##########################################################################
## Setting the kallisto kmer correctly is important to getting best accuracy
## Here I measure the median length as well as the 20th percentile
## KMER is set to length at 20th percentile minus 4nt with a lower limit of 19
MEDIAN_LENGTH=$(sed -n '2~4p' $FQ1 | head -1000000 | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.50 - 0.5)]}')
D20=$(sed -n '2~4p' $FQ1 | head -1000000 | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.20 - 0.5)]}')
KMER=$((D20-4))
ADJUST=$(echo $KMER | awk '{print ($1+1)%2}')
KMER=$((KMER-ADJUST))
if [ $KMER -lt 19 ] ; then KMER=19 ; fi
echo MeadianReadLen=$MEDIAN_LENGTH 20thPercentileLength=$D20 echo kmer=$KMER | tee -a $SRR.log

if [ $KMER -lt "31" ] ; then
  KAL_REF=$(echo $KAL_REF | sed "s#fa.idx#fa.k${KMER}.idx#")
  if [ ! -r $KAL_REF ] ; then
    cd $KAL_DIR
    kallisto index -i $(basename $CDNA).k$KMER.idx -k $KMER $(basename $CDNA)
    for IDX in *idx ; do grep -c '>' $(basename $CDNA) > $IDX.cnt ; done
    cd -
  fi
else
  KMER=31
fi

##########################################################################
echo $SRR running kallisto now
##########################################################################
#Kallisto Ensembl
if [ $RDS == "SE" ] ; then
  echo $SRR Starting Kallisto single end mapping to ensembl reference transcriptome. kmer=$KMER | tee -a $SRR.log
############################################
# TODO need intelligent frag size specification
###########################################
  kallisto quant $KALLISTO_STRAND_PARAMETER --single -l 100 -s 20 -t $THREADS -o . -i $KAL_REF $FQ1 2>&1 \
  | tee -a $SRR.log && mv abundance.tsv $SRR.ke.tsv
  rm abundance.h5
elif [ $RDS == "PE" ] ; then
  echo $SRR Starting Kallisto paired end mapping to ensembl reference transcriptome | tee -a $SRR.log
  kallisto quant $KALLISTO_STRAND_PARAMETER -t $THREADS -o . -i $KAL_REF $FQ1 $FQ2 2>&1 \
  | tee -a $SRR.log && mv abundance.tsv $SRR.ke.tsv
  rm abundance.h5
fi

# collect qc data
PSEUDOMAPPED_CNT=$(grep 'reads pseudoaligned' $SRR.log | awk '{print $(NF-2)}' | tr -d ',')
PSEUDOMAP_RATE=$(echo $PSEUDOMAPPED_CNT $READ_CNT_AVAIL | awk '{print $1/$2*100"%"}')

# Tidy up files
rm -rf run_info.json ${SRR}-trimmed*.fastq _STARgenome

# Check tsv files
wc -l *tsv | tee -a $SRR.log
head *tsv | tee -a $SRR.log
# Check that tsv files have the right number of entries

SE_NR=$(wc -l < $SRR.se.tsv)
KE_NR=$(wc -l < $SRR.ke.tsv)
SE_CNT=$(cat $GTF.cnt)
KE_CNT=$(cat $CDNA.cnt)

if [ $SE_NR -eq $SE_CNT -a $KE_NR -eq $((KE_CNT+1)) ] ; then

  #now place header on the file for later
  echo $SRR completed mapping pipeline successfully | tee -a $SRR.log
  sed -e "1i${SRR}" $SRR.se.tsv > tmp ; mv tmp  $SRR.se.tsv
  sed -e "1i${SRR}_target_id\t${SRR}_length\t${SRR}_eff_length\t${SRR}_est_counts\t${SRR}_tpm" $SRR.ke.tsv | sed 2d > tmp ; mv tmp $SRR.ke.tsv

  # gzip *tsv
  touch $SRR.finished
else
  echo "$SRR An error occurred. Count file line numbers don't match the reference." | tee -a $SRR.log
  exit1 ; return 1
fi

## Collect QC information
echo "SequenceFormat:$ORIG_RDS
QualityEncoding:$QUALITY_ENCODING
Read1MinimumLength:$FQ1_MIN_LEN
Read1MedianLength:$FQ1_MEDIAN_LEN
Read1MaxLength:$FQ1_MAX_LEN
Read2MinimumLength:$FQ2_MIN_LEN
Read2MedianLength:$FQ2_MEDIAN_LEN
Read2MaxLength:$FQ2_MAX_LEN
NumReadsTotal:$READ_CNT_TOTAL
NumReadsQcPass:$READ_CNT_AVAIL
QcPassRate:$QC_PASS_RATE
PE_Read1_StarMapRateTest:${R1_MAP_RATE:-NA}
PE_Read2_StarMapRateTest:${R2_MAP_RATE:-NA}
PE_Read1_Excluded:"${DROP_R1:-FALSE}"
PE_Read2_Excluded:"${DROP_R2:-FALSE}"
MappingFormat:$RDS
STAR_UniqMappedReads:$UNIQ_MAPPED_READS
STAR_Strandedness:$STRANDED
STAR_UnmappedReads:$UNMAPPED_CNT
STAR_MultiMappedReads:$MULTIMAPPED_CNT
STAR_NoFeatureReads:$NOFEATURE_CNT
STAR_AmbiguousReads:$AMBIGUOUS_CNT
STAR_AssignedReads:$ASSIGNED_CNT
STAR_UniqMapRate:$UNIQ_MAP_RATE
STAR_AssignRate:$ASSIGNED_RATE
Kallisto_Kmer:$KMER
Kallisto_MappedReads:$PSEUDOMAPPED_CNT
Kallisto_MapRate:$PSEUDOMAP_RATE
QC_SUMMARY:${QC_SUMMARY}${REASON}" > $SRR.qc

rm -rf *fastq
cd ..
#zip -r $SRR.$ORG.zip $SRR
}
export -f main

#TODO
#-allow specific accessions

cd /dee2

#echo Dumping star genomes from memory
for DIR in $(find /dee2/ref/ | grep /ensembl/star$ | sed 's#\/code\/\.\.##' ) ; do
  STAR --genomeLoad Remove --genomeDir $DIR >/dev/null 2>&1
done

MEM=$(free | awk '$1 ~ /Mem:/  {print $2-$3}')
#MEM=$(free | awk 'NR==2{print $4}')
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

  if [ $MEM_REQD -gt $MEM ] ; then
    echo Error, analysis of $ORG data requires at least $(echo $MEM_REQD $MEM_FACTOR | awk '{print $1*$2}') kB in RAM, but there is only $MEM available.
    exit 1
  fi
fi

if [ -z $MY_ORG ] ; then
  ORGS=$(echo 'athaliana	2853904
celegans	2652204
dmelanogaster	3403644
drerio	14616592
ecoli	1576132
hsapiens	28968508
mmusculus	26069664
rnorvegicus	26913880
scerevisiae	1644684' | awk -v M=$MEM -v F=$MEM_FACTOR 'M>($2*F)' | sort -k2gr | awk '{print $1}')

  TMPHTML=/tmp/tmp.$RANDOM.html
  wget --no-check-certificate -O $TMPHTML "$ACC_URL"
  wget --no-check-certificate -O $TMPHTML $(grep 'frame src=' $TMPHTML | cut -d '"' -f2)

  #specify organism if it has not already been specified by user
  MY_ORG=$(join -1 1 -2 1 \
  <(grep ORG $TMPHTML | cut -d '>' -f2 | tr -d ' .' | tr 'A-Z' 'a-z' | tr '()' ' ' | sort -k 1b,1) \
  <(echo $ORGS | tr ' ' '\n' | sort -k 1b,1) | sort -k2gr | awk 'NR==1{print $1}' )
  rm $TMPHTML
fi

#echo $MY_ORG

myfunc(){
MY_ORG=$1
ACC_REQUEST=$2
TMPHTML=/tmp/tmp.$RANDOM.html
wget --no-check-certificate -r -O $TMPHTML "${ACC_REQUEST}?ORG=${MY_ORG}&Submit"
wget --no-check-certificate -O $TMPHTML $(grep 'frame src=' $TMPHTML | cut -d '"' -f2)
ACCESSION=$(grep 'ACCESSION=' $TMPHTML | cut -d '=' -f2)
STAR --genomeLoad LoadAndExit --genomeDir ../ref/$MY_ORG/ensembl/star >/dev/null  2>&1
echo $ACCESSION
rm $TMPHTML
}
export -f myfunc

key_setup(){
TMPHTML=/tmp/tmp.$RANDOM.html
wget --no-check-certificate -r -O $TMPHTML "dee2.io/ip"
SFTP_URL=$(cat $TMPHTML )
mkdir -p /dee2/.ssh
rm $TMPHTML
touch /dee2/.ssh/known_hosts
chmod 777 /dee2/.ssh/known_hosts

cat << EOF > /dee2/.ssh/guestuser
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

cat << EOF > /dee2/.ssh/guestuser.pub
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDIsnlOKFedvYBN7GBgwiTOkeCikJty1Yofyus9k03o8AtFiODWjN94wXKVK9CKUXD7psJ7seWpRDd7qlL+Mcm6QTLNsk4RVHgDC5nHFy0jhAfTQAvZ4O9a+UQ4KGXE+DwSvlKM/OJRfDw6/ds0u8UdJDuqU1v+BsqEq/OXouTSfpjOX0L96LDNMqN8Qp9cBnnV+PIPZ+ZIVpWV6r63fO+JloZ+0W04sq0MD3BdeKixisG58R/OLG6NAXTgfO/4SDzXhuPOCRKlmIPJLvfk1TO7Q8Judc2NkDpCLKaKw6SHHQnvDeu+f+CaVgyZ4FrkZlmyed9EG+PSCtAh5+QtLfRH mdz@opti
EOF

chmod -R 700 /dee2/.ssh
}
export -f key_setup

##################################################
# Testing the pipeline with ecoli sample
##################################################
TESTFILE=test_pass
if [ ! -r $TESTFILE ] ; then
  echo Initial pipeline test with E. coli dataset
  if [ -d /dee2/data/ecoli/SRR057750 ] ; then
    rm -rf /dee2/data/ecoli/SRR057750
  fi

  #test ssh key setup
  key_setup $SFTP_URL
  cd /dee2/data/ecoli
  date > date.txt
  sftp -v -i /dee2/.ssh/guestuser -o StrictHostKeyChecking=no guestuser@$SFTP_URL << EOF && KEYTEST="OK"
put date.txt
EOF

  if [ $KEYTEST == "OK" ] ; then
    echo "SSH keys successfully set up"
  else
    echo "SSH keys not set up properly. Quitting now."
    exit 1
  fi

  #TEST SRA DATASET
  main ecoli SRR057750 VERBOSE=$VERBOSE
  TEST_CHECKSUM=a739998e33947c0a60edbde92e8f0218
  cd /dee2/data/ecoli/
  TEST_DATASET_USER_CHECKSUM=$(cat SRR057750/SRR057750*tsv | md5sum | awk '{print $1}')
  if [ "$TEST_DATASET_USER_CHECKSUM" != "$TEST_CHECKSUM" ] ; then
    echo "Test dataset did not complete properly. Md5sums do not match those provided!"
    echo "Contact the author for help or flag this issue on the GitHub repo"
    exit 1
  fi
  echo "Test SRA dataset completed successfully"

  #TEST OWN PE FQ DATASET
  wget -O /dee2/mnt/SRR5985593_1.fastq.gz "https://github.com/markziemann/dee2/blob/master/misc/example_data/SRR5985593_1.fastq.gz?raw=true"
  wget -O /dee2/mnt/SRR5985593_2.fastq.gz "https://github.com/markziemann/dee2/blob/master/misc/example_data/SRR5985593_2.fastq.gz?raw=true"
  main ecoli -f /dee2/mnt/SRR5985593_1.fastq.gz /dee2/mnt/SRR5985593_2.fastq.gz VERBOSE=$VERBOSE

  TEST_CHECKSUM=68db313456ae8065ff8d0553bd95325f
  TEST_DATASET_USER_CHECKSUM=$(cat SRR5985593_1/SRR5985593_1*tsv | md5sum | awk '{print $1}')
  if [ "$TEST_DATASET_USER_CHECKSUM" != "$TEST_CHECKSUM" ] ; then
    echo "Test dataset did not complete properly. Md5sums do not match those provided!"
    echo "Contact the author for help or flag this issue on the GitHub repo"
    exit 1
  fi
  echo "Test own dataset completed successfully"

  cd /dee2
  date +"%s" > $TESTFILE

else
  echo
##################################################
# Testing whether the user has provided own data
##################################################
  #Putting this bit into a dummy function for now
  OWN_DATA=$(echo $@ | grep -wc '\-f')
  #echo own data? $OWN_DATA
  if [ $OWN_DATA -eq "1" ] ; then
    echo Starting pipeline with own data specified
    NUM_RDS=$(echo $@ | sed 's/\-f/@/' | cut -d '@' -f2 | sed 's/\-V/@/' | cut -d '@' -f1 | wc -w)

    if [ $NUM_RDS -eq 2 ] ; then
      RDS="PE"
      R1_LIST=$(echo $@ | sed 's/\-f/@/' | cut -d '@' -f2 | awk '{print $1}')
      R2_LIST=$(echo $@ | sed 's/\-f/@/' | cut -d '@' -f2 | awk '{print $2}')
      R1_LIST_LEN=$(echo $R1_LIST | sed 's/,/ /g' | wc -w)
      R2_LIST_LEN=$(echo $R2_LIST | sed 's/,/ /g' | wc -w)

      if [ $R1_LIST_LEN -ne $R2_LIST_LEN ] ; then
        echo Number of foward and reverse readsets does not match. Quitting.
      else
        for DATASET_NUM in $(seq $R1_LIST_LEN) ; do
          FQ_R1=/dee2/mnt/$(echo $R1_LIST | cut -d ',' -f$DATASET_NUM)
          FQ_R2=/dee2/mnt/$(echo $R2_LIST | cut -d ',' -f$DATASET_NUM)

          if [ -r $FQ_R1 -a -r $FQ_R2 ] ; then
            echo "running pipeline.sh $MY_ORG -f $FQ_R1 $FQ_R2"
            main $MY_ORG -f $FQ_R1 $FQ_R2 VERBOSE=$VERBOSE
          else
            echo Specified fastq file $FQ_R1 or $FQ_R2 do not exist or not readable. Quitting
          fi
        done
      fi
    fi

    if [ $NUM_RDS -eq 1 ] ; then

      RDS="SE"
      FQ_LIST=$(echo $@ | sed 's/\-f/@/' | cut -d '@' -f2 | awk '{print $1}')
      FQ_LIST_LEN=$(echo $FQ_LIST | sed 's/,/ /' | wc -w)

      for DATASET_NUM in $(seq $FQ_LIST_LEN) ; do

        FQ=/dee2/mnt/$(echo $FQ_LIST | cut -d ',' -f$DATASET_NUM)

        if [ -r $FQ ] ; then
          echo "running pipeline.sh $MY_ORG -f $FQ"
          main $MY_ORG -f $FQ VERBOSE=$VERBOSE
        else
          echo Specified fastq file $FQ does not exist or not readable. Quitting
        fi
      done

    fi

    if [ -z $RDS ] ; then
      echo "Syntax Error. Unknown option to '-f'"
    fi

    exit
  fi


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
        #key_setup
        cd /dee2/data/$MY_ORG
        zip -r $USER_ACCESSION.$MY_ORG.zip $USER_ACCESSION
        sftp -i /dee2/.ssh/guestuser guestuser@$SFTP_URL << EOF
put $USER_ACCESSION.$MY_ORG.zip
EOF
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
  count=0
  while [ $count -lt 1000 ] ; do
  #while true ; do
    (( count++ ))
    cd /dee2
    echo Run "$count" of 1000
    ACCESSION=$(myfunc $MY_ORG $ACC_REQUEST)
    echo Starting pipeline with species $MY_ORG and accession $ACCESSION
    main "$MY_ORG" "$ACCESSION" VERBOSE=$VERBOSE && COMPLETE=1 || COMPLETE=0
    if [ "$COMPLETE" -eq "1" ] ; then
      #key_setup
      cd /dee2/data/$MY_ORG
      zip -r $ACCESSION.$MY_ORG.zip $ACCESSION
      if [ $(du -s $ACCESSION.$MY_ORG.zip | awk '{print $1}') -lt "20000" ] ; then
        sftp -i /dee2/.ssh/guestuser guestuser@$SFTP_URL << EOF
put $ACCESSION.$MY_ORG.zip
EOF
      else
        rm $ACCESSION.$MY_ORG.zip
      fi
    fi
  done
fi
