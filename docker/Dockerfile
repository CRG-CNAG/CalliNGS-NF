FROM debian:stable

LABEL maintainer "Emilio Palumbo <emilio.palumbo@crg.eu>" \
      version "1.0" \
      description "Varian Calling Analysis with RNA-seq data"

# install needed tools
RUN apt-get update -qq && apt-get install -y -q --no-install-recommends \
      curl \
      libncurses5-dev  \
      libncursesw5-dev \
      build-essential \
      zlib1g-dev \
      bzip2 && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install Java 8
RUN curl -fjksSL --header "Cookie: oraclelicense=accept-securebackup-cookie" -L http://download.oracle.com/otn-pub/java/jdk/8u121-b13/e9e7ea248e2c4826b92b3f075a80e441/jre-8u121-linux-x64.tar.gz | tar xz && \
    update-alternatives --install /usr/bin/java java /jre1.8.0_121/bin/java 100

# install Picard Tools
RUN curl -fksSL https://github.com/broadinstitute/picard/releases/download/2.9.0/picard.jar > /usr/local/bin/picard.jar && \
    chmod +x /usr/local/bin/picard.jar

# install Samtools
RUN curl -fksSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar xj && \
    cd samtools-1.3.1 && \
    make && make install

# install STAR
RUN curl -fksSL https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz | tar xz && \
    cp STAR-2.5.2b/bin/Linux_x86_64/* /usr/local/bin
