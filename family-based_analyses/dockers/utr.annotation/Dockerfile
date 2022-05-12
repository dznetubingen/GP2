FROM r-base:4.1.3

RUN apt-get update && apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev python3 python3-pip
RUN python3 -m pip install h5py tensorflow

# copy requirments R script
COPY dependencies.R /

# run requirments R script
RUN Rscript dependencies.R

# Clean up package registry
RUN rm -rf /var/lib/apt/lists/*
