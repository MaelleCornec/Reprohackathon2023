# Spécification de l'image de base
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive

# Installation de TrimGalore, Cutadapt (version 1.11) et FastQC
RUN apt-get update \
    && apt-get install -y wget curl unzip bzip2 build-essential python3.8 python3-dev default-jre \
    && cd /usr/local \
    && curl -fsSL https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip -o fastqc_v0.12.1.zip \
    && unzip fastqc_v0.12.1.zip \
    && rm fastqc_v0.12.1.zip \
    && wget -O cutadapt.tar.gz https://github.com/marcelm/cutadapt/archive/refs/tags/v1.11.tar.gz \
    && tar -xvf cutadapt.tar.gz \
    && rm cutadapt.tar.gz \
    && cd cutadapt-1.11 \
    && python3 setup.py install \
    && cd /usr/local \
    && curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz \
    && tar -xvf trim_galore.tar.gz \
    && rm trim_galore.tar.gz \
    && apt-get clean \
    && ln -s /usr/local/FastQC /usr/local/FastQC_v0.12.1 \
    && ln -s /usr/local/TrimGalore-0.6.10/trim_galore /usr/local/bin/trim_galore

ENV PATH=$PATH:/usr/local/FastQC_v0.12.1:/usr/local/TrimGalore-0.6.10
