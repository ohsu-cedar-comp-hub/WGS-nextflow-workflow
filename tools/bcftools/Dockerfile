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

ENV PATH="$BCFTOOLS_INSTALL_DIR/bin:$PATH"
