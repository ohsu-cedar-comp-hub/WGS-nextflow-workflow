FROM ubuntu:20.04 AS builder
RUN apt-get update \
    && apt-get install -y \
       bzip2 \
       gcc \
       g++ \
       libbz2-dev \
       liblzma-dev \
       libncurses5-dev \
       make \
       wget \
       git \
       zlib1g-dev \
    && apt-get clean \
    && apt-get update -y && apt-get upgrade -y \
    && git clone  https://github.com/bwa-mem2/bwa-mem2 \ 
    && cd bwa-mem2 \
    && git submodule init \
    && git submodule update \
    && make \
    && mv bwa-mem2* /usr/local/bin/ \
    && cd ../ \
    && rm -rf bwa-mem2* \
    && wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 \
    && tar xf samtools-1.19.tar.bz2 \
    && cd samtools-1.19 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd ../ \
    && rm -rf samtools-1.19 \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
