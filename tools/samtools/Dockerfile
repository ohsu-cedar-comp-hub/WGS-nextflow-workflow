FROM ubuntu:20.04

RUN apt-get update -qq \
  && apt-get install -qq bzip2 gcc g++ make zlib1g-dev wget libncurses5-dev liblzma-dev libbz2-dev

ENV SAMTOOLS_VERSION=1.19

LABEL \
  version="${SAMTOOLS_VERSION}"
  
RUN cd /opt/ \
	&& wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
	&& cd samtools-${SAMTOOLS_VERSION}/ \
	&& make && make install

ENV PATH="/opt/samtools-${SAMTOOLS_VERSION}/:${PATH}"
