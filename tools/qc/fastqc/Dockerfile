FROM ubuntu:bionic-20180426

RUN apt-get update \
    && apt-get install -y \
       libfontconfig1 \
       openjdk-11-jre-headless \
       perl-modules \
       unzip \
       wget \
    && apt-get clean \
    && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip \
    && unzip fastqc_v0.11.7.zip \
    && rm *.zip \
    && mv FastQC /usr/local/ \
    && chmod 755 /usr/local/FastQC/fastqc \
    && sed -i 's/kmer[[:space:]]\+ignore[[:space:]]\+1/kmer ignore 0/' /usr/local/FastQC/Configuration/limits.txt \
    && sed -i 's/assistive_technologies/#assistive_technologies/' /etc/java-11-openjdk/accessibility.properties \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV PATH /usr/local/FastQC/:$PATH