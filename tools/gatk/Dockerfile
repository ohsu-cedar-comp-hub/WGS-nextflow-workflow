FROM ubuntu:20.04 AS builder

ENV VERSION=4.4.0.0

ENV URL=https://github.com/broadinstitute/gatk/releases/download/${VERSION}/gatk-${VERSION}.zip
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
        && apt-get install -y \
                build-essential \
                openjdk-17-jre-headless \
                r-base \
                unzip \
                wget \
                python \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

RUN wget $URL \
        && unzip gatk-${VERSION}.zip \
        && mv gatk-${VERSION}/* /usr/local/bin/ \
        && rm -rf gatk-* \
        && ln -s /usr/local/bin/gatk-package-${VERSION}-local.jar /usr/local/bin/gatk.jar \
        && apt-get purge -y \
                build-essential \
                unzip \
                wget \
        && ls -l /usr/local/bin
