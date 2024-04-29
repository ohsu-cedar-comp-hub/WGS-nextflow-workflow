FROM ubuntu:20.04 AS builder

ARG VERSION=0.38

ENV URL=http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${VERSION}.zip

RUN mkdir -p /usr/share/man/man1 \
    && apt-get update \
	  && apt-get install -y \
	      wget \
	      openjdk-11-jre-headless \
	      unzip \
	  && apt-get clean \
	  && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN wget $URL \
    && unzip Trimmomatic-${VERSION}.zip \
    && ln -s /opt/Trimmomatic-${VERSION}/trimmomatic-${VERSION}.jar /bin/trimmomatic.jar \
	&& chmod 777 /bin/trimmomatic.jar \
    && rm Trimmomatic-${VERSION}.zip

ENV PATH="/usr/bin/:$PATH"
