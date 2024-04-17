FROM ubuntu:20.04

RUN apt-get update && apt-get install -y \
    bash \
    coreutils

WORKDIR /data

COPY split_bed_by_chromosome.sh /usr/local/bin

RUN chmod +x /usr/local/bin/split_bed_by_chromosome.sh

CMD ["bash"]

