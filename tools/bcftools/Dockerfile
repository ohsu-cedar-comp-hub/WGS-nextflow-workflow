FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
  g++ \
  libbz2-dev \
  libcurl4-openssl-dev \
  liblzma-dev \
  make \
  ncurses-dev \
  wget \
  zlib1g-dev

# BCFTOOLS

ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
ENV BCFTOOLS_VERSION=1.12

WORKDIR /tmp

RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
  tar -jxf bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
  cd bcftools-$BCFTOOLS_VERSION && \
  ./configure --prefix=$BCFTOOLS_INSTALL_DIR && \
  make && \
  make install

RUN rm -rf /tmp/bcftools-$BCFTOOLS_VERSION*

# HTSLIB (BGZIP)

ENV HTSLIB_INSTALL_DIR=/opt/htslib
ENV HTSLIB_VERSION=1.20

WORKDIR /tmp

RUN wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 && \
    tar -jxf htslib-$HTSLIB_VERSION.tar.bz2 && \
    cd htslib-$HTSLIB_VERSION && \
    ./configure --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install

RUN rm -rf /tmp/htslib-$HTSLIB_VERSION*

ENV PATH="/opt/bcftools/bin/:${PATH}"
ENV PATH="/opt/htslib/bin/:${PATH}"