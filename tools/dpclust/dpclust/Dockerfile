# Use the rocker image with R and Ubuntu as base
FROM rocker/r-ver:4.2.2
USER root
# Install system dependencies for building R packages (like curl, ssl, and libxml)
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Set CRAN mirror for package installation
RUN R -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("BiocManager", lib = "/usr/local/lib/R/site-library")'
RUN R -e 'BiocManager::install(c("KernSmooth", "lattice"))'
RUN R -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("optparse"); install.packages("ks"); install.packages("ggplot2"); install.packages("gridExtra")'

RUN mkdir -p /opt/dpclust
COPY . /opt/dpclust/
RUN R -q -e 'install.packages("/opt/dpclust", repos=NULL, type="source")'

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
