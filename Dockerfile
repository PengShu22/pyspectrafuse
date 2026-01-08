################## BASE IMAGE ######################
FROM python:3.10-slim

################## METADATA ######################
LABEL base_image="python:3.10-slim"
LABEL version="1"
LABEL software="pyspectrafuse"
LABEL software.version="0.0.2"
LABEL about.summary="pyspectrafuse - Command-line utilities for spectrum clustering and conversion"
LABEL about.home="https://github.com/bigbio/pyspectrafuse"
LABEL about.documentation="https://github.com/bigbio/pyspectrafuse"
LABEL about.license_file="https://github.com/bigbio/pyspectrafuse/blob/master/LICENSE"
LABEL about.license="SPDX:Apache-2.0"
LABEL about.tags="Proteomics,Multiomics,QuantMS"

################## MAINTAINER ######################
MAINTAINER Yasset Perez-Riverol <ypriverol@gmail.com>

################## INSTALLATION ######################

ENV DEBIAN_FRONTEND=noninteractive

# Disable numba caching to avoid issues in containerized environments
# Numba caching fails in containers because it can't locate source files in site-packages
ENV NUMBA_DISABLE_CACHING=1

## Update and install packages
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    build-essential \
    procps \
    && rm -rf /var/lib/apt/lists/*

## Set working directory
WORKDIR /data/

## Copy requirements first for better caching
COPY requirements.txt .

## Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

## Copy project files
COPY pyproject.toml .
COPY pyspectrafuse/ ./pyspectrafuse/

## Install the package
RUN pip install --no-cache-dir -e .

