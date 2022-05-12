#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

input_csv=args[1]
db = args[2]
version = args[3]

library(keras)
library(readr)
library(utr.annotation)
library(MRL.dl.model)


initUTRAnnotation(variantFile = input_csv, species = "human", dataDir = db, ensemblVersion = version)
