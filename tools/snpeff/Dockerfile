FROM openjdk:17-jdk-slim

WORKDIR /usr/src/app

RUN apt-get update && \
    apt-get install -y wget unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && \
    rm snpEff_latest_core.zip

# assuming snpeff looks for data in /usr/src/app/snpEff/data (instead maybe we download at runtime and mount as a volume?)
RUN java -jar snpEff/snpEff.jar download -v GRCh38.86

# ENV _JAVA_OPTIONS=-Xmx8g

# ENTRYPOINT ["java", "-jar", "snpEff/snpEff.jar"]
