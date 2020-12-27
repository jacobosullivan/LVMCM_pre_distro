#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print(args)

### This script will first generate Moran's Eigenvector Maps for an input spatial network
### Required modules to run adespatial on Apocrita HPC: 
# module load gdal/2.3.1 proj/5.2.0 gcc/8.2.0 udunits/2.2.26 geos/3.7.1 R/3.6.1

### The analysis is taken from this vignette:
### https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html#multiscale-analysis-with-mem

### Then it will run the variance partitioning analysis taken from Leibold et al. 2020

### args:
  # 1: ncores for parellel for loop
  # 2: select parameter combination
  # 3: output directory for RDS objects
  # 4: LOAD 
    # 0: build model
    # 1: load pre-built model
  # 5: VARPAR - control algorithm:
    # 1: variance partitioning by species only
    # 2: variance partitioning by species and site
    # 3: variance partitioning by site only

if (length(args) != 5) {
  stop("Lenght of command line arguments not equal to five!")
}

# library(ade4)
# library(adespatial)
# library(adegraphics)
# library(spdep)
# library(maptools)
library(data.table)
library(tidyverse)
library(HMSC)
library(doParallel)
library(ggtern)
library(corrplot)

ncores <- as.numeric(args[1])
# ncores <- 1
LOAD <- as.numeric(args[4])
VARPAR <- as.numeric(args[5])

folderpath <- args[3]
# folderpath <- "~/Desktop/HMSCtestdata/"
if (stringr::str_sub(args[3],start=-1) != "\\/") {
  folderpath <- paste0(folderpath, "/")
}

if(dir.exists(folderpath) == FALSE){
  dir.create(folderpath)
}

print("Loading libraries complete")

# In the past, we tried several combinations of these parameters and there didn't seem to be much of a difference between them. According to the convergence plots it seems to happen pretty fast. We settled on the parameters in the next line, but you are free to try with others, example shown in the commented line after setting the current parameters. 
# hmscPars <- list(niter = 50000, nburn = 25000, thin = 10)
hmscPars <- list(niter = 30000, nburn = 10000, thin = 10) # This is ok, but not enough for raftery.diag
# hmscPars <- list(niter = 10000, nburn = 5000, thin = 5)

# This full process involves running the metacommunity simulation, getting in the hmsc format, doing variation partitioning for species and sites.
source("/data/home/btx188/testingHMSC/manuscript_functions/full_process_fx.R")
# source("/home/jack/gitClones/testingHMSC/manuscript_functions/full_process_fx.R")

# Output processing involves changing structure of data from lists to dataframes and using those to get the figures
source("/data/home/btx188/testingHMSC/manuscript_functions/output_processing_fx.R")
# source("/home/jack/gitClones/testingHMSC/manuscript_functions/output_processing_fx.R")

# Functions to check convergence based on previous conversations with Guillaume. Raftery and Gelman plots
source("/data/home/btx188/testingHMSC/manuscript_functions/convergence_fx.R")
# source("/home/jack/gitClones/testingHMSC/manuscript_functions/convergence_fx.R")

source("/data/home/btx188/testingHMSC/LVMCM/metacom_as_HMSCdataLVMCM.R")
# source("/home/jack/gitClones/testingHMSC/LVMCM/metacom_as_HMSCdataLVMCM.R")

source("/data/home/btx188/testingHMSC/LVMCM/get_species_dataLVMCM.R")
# source("/home/jack/gitClones/testingHMSC/LVMCM/get_species_dataLVMCM.R")

source("/data/home/btx188/testingHMSC/LVMCM/get_site_dataLVMCM.R")
# source("/home/jack/gitClones/testingHMSC/LVMCM/get_site_dataLVMCM.R")

source("/data/home/btx188/testingHMSC/LVMCM/get_VPresults_SITELVMCM.R")
# source("/home/jack/gitClones/testingHMSC/LVMCM/get_VPresults_SITELVMCM.R")

source("/data/home/btx188/testingHMSC/LVMCM/get_VPresults_LVMCM.R")
# source("/home/jack/gitClones/testingHMSC/LVMCM/get_VPresults_LVMCM.R")

print("Loading functions complete")

# Get file names for importing LVMCM output
files <- read.table("/data/home/btx188/LVMCM_src/LVMCM/shFiles/removal/bfileRemoval.txt", as.is = 1) # list of file names
exp <- stringr::str_extract(files$V1, "\\)\\d.*bMat")
exp <- gsub("\\)", "", exp)
exp <- as.numeric(gsub("bMat", "", exp))
files <- data.frame(filename=as.character(files$V1), exp=exp)
files$filename <- as.character(files$filename)
colnames(files)[1] <- "filename" 

# Set presence-absence threshold
thresh <- 1e-4 # simple numerical threshold

# Subset dataset
j = as.numeric(args[2]) # experimental iterator passed as command line argument
#fileSub <- subset(files, exp==unique(exp)[j])
fileSub <- files[c(1:10)+(10*(j-1)),]
namesrds <- paste0("LVMCM_removal_", unique(fileSub$exp), "_", j%%10)
# namesrds <- "LVMCM_removal_test"

print("Subsetting simulations complete")

print(fileSub)

# Import data
print("Importing sim...")
sims <- list()
XY <- list()
E <- list()
MEMsel <- list()

# Test case
# B <- t(as.matrix(read.table("~/Desktop/HMSCtestdata/B10.mat")))
# B <- t(as.matrix(read.table("~/Desktop/HMSCtestdata/2020-8-14_removal(10000)1.1bMat0.mat")))
# P <- B
# P[P<thresh] <- 0
# P[P>=thresh] <- 1 # uncomment for thresholded presence absence data
# sims <- list(as.matrix(P))
# E <- list(t(as.matrix(read.table("~/Desktop/HMSCtestdata/2020-8-14_removal(10000)1.1envMat0.mat"))))
# XY <- list(as.matrix(read.table("~/Desktop/HMSCtestdata/2020-8-14_removal(10000)1.1network0.mat")))
# MEMsel <- list(readRDS("~/Desktop/HMSCtestdata/2020-8-14_removal(10000)1.1MEMsel0.RDS")$MEM.select)

for (i in 1:nrow(fileSub)) {
  # Load data
  print(fileSub$filename[i])
  B <- t(as.matrix(fread(fileSub$filename[i])))
  attr(B, "dimnames") <- NULL
  P <- B
  P[P<thresh] <- 0
  P[P>=thresh] <- 1 # uncomment for thresholded presence absence data
  sims[[i]] <- as.matrix(P)
  XY[[i]] <- as.matrix(fread(gsub("bMat", "network", fileSub$filename[i])))
  E[[i]] <- t(as.matrix(fread(gsub("bMat", "envMat", fileSub$filename[i]))))
  MS <- read_rds(gsub("\\.mat", "\\.RDS", gsub("bMat", "MEMsel", fileSub$filename[i])))$MEM.select
  if (is.null(MS)) {
    MEMsel[[i]] <- NA  
  } else {
    MEMsel[[i]] <- MS
  }
}

write_rds(sims, paste0(folderpath,namesrds,"-metacomSim.RDS"))
print("Sim imported")

if (LOAD==0) {
  print("Building model...")
  system.time(
    model <- metacom_as_HMSCdata(sims, numClusters = ncores, E = E, MEMsel = MEMsel,
                                 hmscPars = hmscPars,
                                 makeRDS = TRUE, whereToSave = folderpath, objName = namesrds)    
  )
  print("Model complete...")
} else {
  print("Importing model...")
  model <- readRDS(paste0(folderpath, namesrds, "-model.RDS"))
  print("Model imported.")
}

if (VARPAR %in% c(1,2)) {
  print("Partitioning variance by species...")
  system.time(
    vpSpp <- get_VPresults(model, MEMsel = MEMsel, numClusters = ncores,
                           makeRDS = TRUE, whereToSave = folderpath, objName = namesrds)
  )
  print("Complete.")
}

if (VARPAR %in% c(2,3)) {
  print("Partitioning variance by site...")
  system.time(
    vpSites <- get_VPresults_SITE(model, MEMsel = MEMsel, E = E, numClusters = ncores,
                                  makeRDS = TRUE, whereToSave = folderpath, objName = namesrds)  
  )
  print("Complete.")
}

# hist(unlist(vpSpp))
# hist(unlist(vpSites))
