FROM continuumio/miniconda3
MAINTAINER Anastasia Illarionova <avalillarionova@gmail.com>
	

# Install PROCPS
RUN apt-get update && apt-get install -y procps
RUN yes | apt-get install git make gcc g++ libtool zlib1g-dev libgsl-dev perl

## Install Nanopolish
RUN git clone --recursive https://github.com/jts/nanopolish.git
RUN cd nanopolish; make
ENV PATH $PATH:/nanopolish

## Create conda env for SV detection
COPY SV_env.yml /
RUN conda env create -f /SV_env.yml && conda clean -a


## Install Winnowmap
RUN wget https://github.com/marbl/Winnowmap/archive/refs/tags/v2.0.tar.gz && tar -xf v2.0.tar.gz
RUN cd Winnowmap-2.0; make -j8
ENV PATH $PATH:/Winnowmap-2.0/bin

