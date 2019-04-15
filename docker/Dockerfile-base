FROM debian:stable

LABEL maintainer "Emilio Palumbo <emilio.palumbo@crg.eu>" \
      version "1.0" \
      description "Varian Calling Analysis with RNA-seq data"

# install needed tools
RUN apt-get update --fix-missing -qq && apt-get install -y -q \
    curl \
    locales \
    libncurses5-dev  \
    libncursesw5-dev \
    build-essential \
    pkg-config \
    zlib1g-dev \
    bzip2 \
    r-base \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


RUN curl -fjksSL --header "Cookie: oraclelicense=accept-securebackup-cookie" -L http://download.oracle.com/otn-pub/java/jdk/8u121-b13/e9e7ea248e2c4826b92b3f075a80e441/jre-8u121-linux-x64.tar.gz | tar xz && \
    update-alternatives --install /usr/bin/java java /jre1.8.0_121/bin/java 100

# install Picard Tools
RUN curl -fksSL https://github.com/broadinstitute/picard/releases/download/2.9.0/picard.jar > /usr/local/bin/picard.jar && \
    chmod +x /usr/local/bin/picard.jar

# install SAMtools
RUN curl -fksSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar xj && \
    cd samtools-1.3.1 && \
    make all all-htslib && make install install-htslib

# install VCFtools
RUN curl -fksSL https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz | tar xz && \
    cd vcftools-0.1.14 && \
    ./configure && make && make install

# install STAR
RUN curl -fksSL https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz | tar xz && \
    cp STAR-2.5.2b/bin/Linux_x86_64/* /usr/local/bin

## Install R packages for ggplot2
RUN R -e 'install.packages( c("reshape2","optparse"), repos="http://cloud.r-project.org/");' && \
    apt-get update && apt-get install r-cran-ggplot2 -y -q
