# Base image
FROM ubuntu:16.04

# Metadata
LABEl base.image="ubuntu:16.04"
LABEL version="1"
LABEL software="Image for DEE2"
LABEL software.version="08252016"
LABEL description="Image for DEE2"
LABEL website=""
LABEL documentation=""
LABEL license=""
LABEL tags="Genomics"

# Maintainer
MAINTAINER Mark Ziemann <mark.ziemann@gmail.com>

ENV DIRPATH /home/data/
WORKDIR $DIRPATH

#numaverage numround numsum
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    num-utils \
    wget \
    git \
    perl \
    unzip \
    bowtie2 \
    python3 \
    python3-pip

########################################
# SRA TOOLKIT WORKING
########################################
ENV VERSION 2.8.2-1
ADD http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-ubuntu64.tar.gz /tmp/
RUN tar zxfv /tmp/sratoolkit.${VERSION}-ubuntu64.tar.gz && \
    cp -r sratoolkit.${VERSION}-ubuntu64/bin/* /usr/local/bin

########################################
# Install parallel-fastq-dump
########################################
#COPY get-pip.py .
RUN pip3 install --upgrade pip
RUN pip3 install parallel-fastq-dump

########################################
# SKEWER WORKING
########################################
RUN \
  wget -O /tmp/skewer-0.2.2-linux-x86_64 "https://downloads.sourceforge.net/project/skewer/Binaries/skewer-0.2.2-linux-x86_64?r=&ts=1504573715&use_mirror=nchc" && \
  mv /tmp/skewer-0.2.2-linux-x86_64 /tmp/skewer && \
  chmod +x /tmp/skewer && \
  cp /tmp/skewer /usr/local/bin/
ENTRYPOINT ["skewer"]

########################################
# MINION from kraken toolkit (ebi)
########################################
RUN \
  wget -c http://wwwdev.ebi.ac.uk/enright-dev/kraken/reaper/binaries/reaper-13-100/linux/minion && \
  chmod +x  minion && \
  cp minion /usr/local/bin/minion
ENTRYPOINT ["minion"]

########################################
# STAR
########################################
RUN \
  wget -c https://github.com/alexdobin/STAR/blob/master/bin/Linux_x86_64/STAR  && \
  chmod +x STAR && \
  cp STAR /usr/local/bin/STAR
ENTRYPOINT ["STAR"]

########################################
# Fastqc
########################################
RUN \
  wget -O fastqc_v0.11.5.zip "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip" && \
  unzip fastqc_v0.11.5.zip && \
  cd FastQC && \
  mv * /usr/local/bin/
ENTRYPOINT ["fastqc"]

########################################
# KALLISTO
########################################
RUN \
  wget -c https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz && \
  tar xf kallisto_linux-v0.43.1.tar.gz && \
  cd kallisto_linux-v0.43.1 && \
  chmod +x kallisto && \
  cp kallisto /usr/local/bin/kallisto  
ENTRYPOINT ["kallisto"]


########################################
# Get the dee2 repo
########################################
RUN git clone https://github.com/markziemann/dee2.git && \
  cd dee2 && \
  mkdir code && \
  cp volunteer_pipeline.sh code && \
  cd code

########################################
# ASCP and the NCBI license WORKING
########################################
ADD http://download.asperasoft.com/download/sw/ascp-client/3.5.4/ascp-install-3.5.4.102989-linux-64.sh /tmp/
## No https, so verify sha1
RUN test $(sha1sum /tmp/ascp-install-3.5.4.102989-linux-64.sh |cut -f1 -d\ ) = a99a63a85fee418d16000a1a51cc70b489755957 && \
   sh /tmp/ascp-install-3.5.4.102989-linux-64.sh
#RUN useradd data
#USER data
ENTRYPOINT ["/usr/local/bin/ascp"]

########################################
# run dee2
########################################
RUN bash volunteer_pipeline.sh ecoli
