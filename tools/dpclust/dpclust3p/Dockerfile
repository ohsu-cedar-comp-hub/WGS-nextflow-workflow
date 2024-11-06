# Use the rocker image with R and Ubuntu as base
FROM rocker/r-ver:4.2.2

USER root
# Install system dependencies for building R packages
RUN apt-get update && apt-get install -y \
    r-base \
    libxml2 \
    libxml2-dev \
    libcurl4-gnutls-dev \
    libssl-dev curl \
    && rm -rf /var/lib/apt/lists/*

# Set CRAN mirror as cloud and use BiocManager to install packages
RUN R -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("BiocManager", lib = "/usr/local/lib/R/site-library")'
RUN R -e 'BiocManager::install(c("VariantAnnotation", "GenomicRanges", "Rsamtools", "IRanges", "S4Vectors"))'
RUN R -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("optparse"); install.packages("ggplot2"); install.packages("reshape2")'

RUN mkdir -p /opt/dpclust3p
COPY . /opt/dpclust3p/
RUN R -q -e 'install.packages("/opt/dpclust3p", repos=NULL, type="source")'

RUN mkdir /tmp/downloads

# Download and configure samtools
RUN curl -sSL -o tmp.tar.gz --retry 10 https://github.com/samtools/htslib/archive/1.7.tar.gz && \
    mkdir /tmp/downloads/htslib && \
    tar -C /tmp/downloads/htslib --strip-components 1 -zxf tmp.tar.gz && \
    make -C /tmp/downloads/htslib && \
    rm -f /tmp/downloads/tmp.tar.gz

ENV HTSLIB=/tmp/downloads/htslib

# Download alleleCount
RUN curl -sSL -o tmp.tar.gz --retry 10 https://github.com/cancerit/alleleCount/archive/v4.0.0.tar.gz && \
    mkdir /tmp/downloads/alleleCount && \
    tar -C /tmp/downloads/alleleCount --strip-components 1 -zxf tmp.tar.gz && \
    cd /tmp/downloads/alleleCount/c && \
    mkdir bin && \
    make && \
    cp /tmp/downloads/alleleCount/c/bin/alleleCounter /usr/local/bin/. && \
    cd /tmp/downloads && \
    rm -rf /tmp/downloads/alleleCount /tmp/downloads/tmp.tar.gz

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
