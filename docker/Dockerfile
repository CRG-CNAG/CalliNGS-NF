FROM debian:testing

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
    openjdk-8-jre-headless && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


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
    
    
# install R 

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## Use Debian unstable via pinning -- new style via APT::Default-Release
RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
	&& echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default

ENV R_BASE_VERSION 3.3.3

## Now install R and littler, and create a link for littler in /usr/local/bin
## Also set a default CRAN repo, and make sure littler knows about it too
RUN apt-get update \
	&& apt-get install -t unstable -y --no-install-recommends \
		littler \
                r-cran-littler \
		r-base=${R_BASE_VERSION}* \
		r-base-dev=${R_BASE_VERSION}* \
		r-recommended=${R_BASE_VERSION}* \
        && echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
        && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
	&& ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*    


