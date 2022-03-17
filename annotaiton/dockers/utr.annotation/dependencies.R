### Install dependencis

cran_pkgs <- c("parallel", "doParallel", "data.table", "readr", "stringr", "vcfR", "dplyr", "tidyr", "keras", "reticulate", "h5py")
bioc_pkgs <- c("biomaRt", "Biostrings", "AnnotationHub", "ensembldb")

install.packages(c('devtools','curl'))

for (pkg in cran_pkgs) {
  if (!(pkg %in% installed.packages())) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
BiocManager::install(version = "3.14")
}

for (pkg in bioc_pkgs) {
  if (!(pkg %in% installed.packages())) {
    BiocManager::install(pkg)
  }
}

# Install keras package for MRL prediction. 
library(keras)
# Install tensorflow backend
reticulate::install_miniconda()
keras::install_keras(version = "2.2.4", tensorflow = "2.2.0", method = "conda")

# install deep learning model data package

library(devtools)
devtools::install_bitbucket("jdlabteam/mrl.dl.model")


####### Install utr.annotation package

install.packages("utr.annotation")

