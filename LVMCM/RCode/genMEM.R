#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args[1])

### This script will generate Moran's Eigenvector Maps for an input spatial network
### Required modules to run adespatial on Apocrita HPC: 
# module load gdal/2.3.1 proj/5.2.0 gcc/8.2.0 udunits/2.2.26 geos/3.7.1 R/3.6.1

### The following analysis is taken from this vignette:
### https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html#multiscale-analysis-with-mem

library(ade4)
library(adespatial)
library(adegraphics)
library(spdep)
# library(maptools)
library(data.table)
library(tidyverse)

SIMULATE=F
PRINT=F
SELECT=T

if (SIMULATE) {
  ## Gen random network
  # N <- 16
  # S <- 1
  # X <- matrix(runif(2*N), ncol=2)
  # B <- matrix(NA, nrow=S, ncol=N)
  # for (x in 1:N) {
  #   B[,x] <- exp(rnorm(S))
  # }
  ## Load LVMCM output
  X <- as.matrix(fread("~/Desktop/testMEM/2020-8-14_removal(10000)1.1network0.mat"))
  B <- as.matrix(fread("~/Desktop/testMEM/2020-8-14_removal(10000)1.1bMat0.mat"))
} else {
  X <- as.matrix(fread(gsub("bMat", "network",args[1])))
  if (SELECT) {
    B <- as.matrix(fread(args[1]))
  }
}

## Gen Gabriel graph
nbgab <- graph2nb(gabrielneigh(X), sym=T)
if (PRINT) {
  g1 <-  s.label(X, nb = nbgab, pnb.edge.col = "red", main = "Gabriel", plot = FALSE)
  g1 # print call
}

## Gen distance weighting
distgab <- nbdists(nbgab, X)
fdist <- lapply(distgab, function(x) 1 - x/max(dist(X)))
listwgab <- nb2listw(nbgab, glist = fdist)
if (PRINT) {
  listwgab # print call
  print(listw2mat(listwgab)[1:10, 1:10], digits = 3) # print call  
}

if(!SELECT) {
  ## Gen Moran's Eigenvector Map
  mem.x <- mem(listwgab)
  if(PRINT) {
    mem.x # print call
    barplot(attr(mem.x, "values"), 
            main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7) # print call  
  }
} else {
  ## Gen Moran's Eigenvector Map with reduced no. of MEMs to avoid over fitting
  ## This involves assessing which of the MEMs explain most of the variance IN THE RESPONSE VARIABLE(S)
  
  # By default, only MEMs associated to positive eigenvalues are considered
  # (argument MEM.autocor = "positive") and a forward selection (based on R2
  # statistic) is performed after a global test (method = "FWD")
  mem.x <- mem.select(t(B), listw = listwgab, MEM.autocor = "positive", method = "FWD",nperm = 999, alpha = 0.05)
  if (PRINT) {
    mem.x$global.test # print call
    mem.x$summary # print call
    barplot(attr(mem.x$MEM.select, "values"), 
            main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7) # print call    
  }
}

write_rds(mem.x, gsub("\\.mat", "\\.RDS", gsub("bMat", "MEMsel", args[1])))