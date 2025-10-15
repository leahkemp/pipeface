## build env ##
FROM ubuntu:24.04 AS build

# get necessary libs
RUN apt-get update && \
    apt install -y make wget gcc bzip2 libz-dev libbz2-dev liblzma-dev libcurl4-openssl-dev build-essential python3-dev pip libtool cmake git zlib1g-dev openjdk-11-jdk bash

# samtools
RUN wget -O- "https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2" | tar -xj && \
    cd samtools-* && \
    ./configure --without-curses && \
    make && \
    make install

# htslib
RUN wget -O- "https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2" | tar -xj && \
    cd htslib-* && \
    ./configure && \
    make && \
    make install

# minimap2
RUN wget -O- "https://github.com/lh3/minimap2/archive/refs/tags/v2.28.tar.gz" | tar -zx && \
    cd minimap2-* && \
    make && \
    chmod a+x minimap2 && \
    mv minimap2 /usr/local/bin/

# mosdepth
RUN wget "https://github.com/brentp/mosdepth/releases/download/v0.3.9/mosdepth" && \
    chmod a+x mosdepth && \
    mv mosdepth /usr/local/bin/

# bcftools
RUN wget -O- "https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2" | tar -xj && \
    cd bcftools-* && \
    make && \
    make install && \
    mkdir -p /usr/local/libexec/bcftools/ && \
    mv /bcftools-1.21/plugins/* /usr/local/libexec/bcftools/

# whatshap
RUN pip install whatshap==2.3 --break-system-packages

# minimod
RUN wget -O- "https://github.com/warp9seq/minimod/releases/download/v0.3.0/minimod-v0.3.0-release.tar.gz" | tar -xz && \
    cd minimod-* && \
    ./scripts/install-hts.sh && \
    make && \
    mv minimod /usr/local/bin/

# kentutils (bedGraphToBigWig only)
RUN wget "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig" && \
    chmod a+x bedGraphToBigWig && \
    mv bedGraphToBigWig /usr/local/bin/

# cutesv
RUN pip install cuteSV==2.1.1 --break-system-packages

# sniffles
RUN pip install sniffles==2.6.0 --break-system-packages

# somalier
RUN wget "https://github.com/brentp/somalier/releases/download/v0.2.19/somalier" && \
    chmod a+x somalier && \
    mv somalier /usr/local/bin/

# spectre
RUN pip install spectre-cnv==0.2.1 --break-system-packages

# snf2json
RUN pip install snf2json==0.1.0 --break-system-packages

# bedtools
RUN wget "https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static" && \
    chmod a+x bedtools.static && \
    mv bedtools.static /usr/local/bin/bedtools

# racon
RUN git clone --recursive https://github.com/lbcb-sci/racon.git racon && \
    cd racon && \
    git checkout tags/1.4.3 && \
    git submodule update --init --recursive && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make

# jasmine
RUN git clone --recurse-submodules https://github.com/bioinfomethods/Jasmine.git && \
    cd Jasmine && \
    git checkout tags/1.1.5-r1 && \
    ./build_jar.sh && \
    sed -i '1s|^.*$|#!/bin/bash|' jasmine

## deploy env ##
FROM ubuntu:24.04 AS deploy
LABEL name="pipeface"
LABEL description="docker image containing most software required for pipeface"
LABEL version="0.0.1"
LABEL maintainer.name="Leah Kemp"
LABEL maintainer.email="leahmhkemp@gmail.com"

# get necessary libs
RUN apt-get update && \
    apt install -y python3-dev libcurl4-openssl-dev openjdk-11-jdk

# copy required compiled tools
COPY --from=build \
    /usr/local/bin/samtools \
    /usr/local/bin/bgzip \
    /usr/local/bin/tabix \
    /usr/local/bin/minimap2 \
    /usr/local/bin/mosdepth \
    /usr/local/bin/bcftools \
    /usr/local/bin/whatshap \
    /usr/local/bin/minimod \
    /usr/local/bin/bedGraphToBigWig \
    /usr/local/bin/cuteSV \
    /usr/local/bin/sniffles \
    /usr/local/bin/somalier \
    /usr/local/bin/spectre \
    /usr/local/bin/snf2json \
    /usr/local/bin/bedtools \
    /racon/build/bin/racon \
    /Jasmine/jasmine \
    /Jasmine/jasmine.jar \
    /Jasmine/jasmine_iris.jar \
    /usr/local/bin/

COPY --from=build \
    /usr/local/libexec/bcftools/* \
    /usr/local/libexec/bcftools/

COPY --from=build \
    /usr/local/lib/ \
    /usr/local/lib/

