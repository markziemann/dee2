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

#numaverage numround numsum
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    num-utils \
    wget \
    fastqc \
    perl \
    unzip \
    bowtie2

########################################
# SRA TOOLKIT WORKING
########################################
ENV VERSION 2.8.2-1
ADD http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-ubuntu64.tar.gz /tmp/
RUN tar zxfv /tmp/sratoolkit.${VERSION}-ubuntu64.tar.gz && \
    cp -r sratoolkit.${VERSION}-ubuntu64/bin/* /usr/

########################################
# ASCP and the NCBI license WORKING
########################################
ADD http://download.asperasoft.com/download/sw/ascp-client/3.5.4/ascp-install-3.5.4.102989-linux-64.sh /tmp/
 No https, so verify sha1
RUN test $(sha1sum /tmp/ascp-install-3.5.4.102989-linux-64.sh |cut -f1 -d\ ) = a99a63a85fee418d16000a1a51cc70b489755957 && \
    sh /tmp/ascp-install-3.5.4.102989-linux-64.sh
RUN useradd data
USER data
ENTRYPOINT ["/usr/local/bin/ascp"]

########################################
# SKEWER WORKING
########################################
RUN \
  wget -c https://downloads.sourceforge.net/project/skewer/Binaries/skewer-0.2.2-linux-x86_64 && \
  chmod +x skewer-0.2.2-linux-x86_64 && \
  cp skewer-0.2.2-linux-x86_64 /usr/local/bin/skewer
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
# KALLISTO
########################################
RUN \
  wget -c https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz && \
  tar xf kallisto_linux-v0.43.1.tar.gz && \
  cd kallisto_linux-v0.43.1 \
  chmod +x kallisto && \
  cp kallisto /usr/local/bin/kallisto  
ENTRYPOINT ["kallisto"]

