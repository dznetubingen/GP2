#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

input_csv=args[1]
output_csv=args[2]
db = args[3]
tr = args[4]
version = args[5]

library(keras)
library(readr)
library(utr.annotation)
library(MRL.dl.model)

runUTRAnnotation(variantFile = input_csv, annotationResult = output_csv,species = "human", dataDir = db, ensemblVersion = version, mrl_prediction = TRUE, cores = tr)
