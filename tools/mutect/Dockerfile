FROM ubuntu:20.04

ENV VERSION=4.2.4.1
ENV URL=https://github.com/broadinstitute/gatk/releases/download/${VERSION}/gatk-${VERSION}.zip
ENV DEBIAN_FRONTEND=noninteractive
ENV GATK_LOCAL_JAR=/usr/local/bin/gatk.jar

RUN apt-get update && apt-get install -y \
        openjdk-8-jre-headless \
        r-base \
        unzip \
        wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && wget $URL \
    && unzip gatk-${VERSION}.zip \
    && mv gatk-${VERSION}/* /usr/local/bin/ \
    && rm -rf gatk-* \
    && ln -s /usr/local/bin/gatk-package-${VERSION}-local.jar $GATK_LOCAL_JAR

VOLUME /data
WORKDIR /data
