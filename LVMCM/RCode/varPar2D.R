#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# This script will run a 2D (E|S) variance partitioning on the LVMCM removal experiment database
# Note sure how useful this will be, gives a metacommunity scale metric!

# To run in interactive session:
# args <- c("/data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_data.csv",
#           "/data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_2D_vp.csv",
#           203)

# Example command line run:
# Rscript removalExperiment.R /data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_data.csv /data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_2D_vp.csv 1

## 2D variance partitioning
require(vegan)
require(data.table)
require(tidyverse)

data(mite) # site by species table (abundance)
data(mite.env) # site by environment table
data(mite.pcnm) # reduced Principle Coordinates of Neighbour Matrix (PCNM)

# Two explanatory matrices -- Hellinger-transform Y
# Formula shortcut "~ ." means: use all variables in 'data'.
mod <- varpart(mite, ~ ., mite.pcnm, data=mite.env, transfo="hel")
mod


thresh <- 1e-4 # numerical threshold

dat_complete <- c()

for (sim in as.numeric(args[3]):length(fileName)) {
  # for (sim in 1:2) {
  print(paste("Simulation", sim, " filename ", fileName[sim]))
  
  ### load matrices of unperturbed assembly
  B <- as.matrix(fread(fileName[sim]))
  E <- as.matrix(fread(gsub("bMat", "envMat", fileName[sim])))
  MEM <- readRDS(gsub("mat", "RDS", gsub("bMat", "MEMsel", fileName[sim])))
  
    ### generate unperturbed presence-absence and Bray-Curtis matrices
  P0 <- B0
  P0[P0<thresh] <- 0
  P0[P0>0] <- 1
  B0[B0<0] <- 0
  logBth <- B0
  logBth[logBth<thresh] <- NA
  logBth <- log(logBth)
  BC <- as.matrix(vegdist(t(B0)))
  totalOcc0 <- sum(P0) # total occupied niches
  
  dat <- c() # store diversity metrices etc.
  comm <- c() # store alpha diversity matrix
  alpha0 <- c()
  
  dat$phi <- rep(as.numeric(PARS[which(PARS[,1] == "phi"),2]), nrow(dat))
  dat$dispL <- rep(as.numeric(PARS[which(PARS[,1] == "dispL"),2]), nrow(dat))
  dat$c1 <- rep(as.numeric(PARS[which(PARS[,1] == "c2"),2]), nrow(dat)) # miss-named in C++ program
  dat$N <- rep(as.numeric(PARS[which(PARS[,1] == "no_nodes"),2]), nrow(dat))
  dat$bFile <- rep(fileName[sim], nrow(dat))
  
  if (QUAD) {
    dat$sk <- rep(as.numeric(PARS[which(PARS[,1] == "sk"),2]), nrow(dat))
    dat$env <- E[1,]
  }
  
  dat$alpha0 <- alpha0
  B0_s[B0_s < 0] <- 0
  dat$alphaS <- colSums(B0_s)
  dat$B.tot <- colSums(B0)
  
  ### Network centrality measures
  W <- 1*(D>0) * as.matrix(dist(X)) # weighted adjacency matrix
  g <- graph_from_adjacency_matrix(W, weighted=T)
  dat$cent_degree <- centr_degree(g)$res
  dat$cent_close <- closeness.freeman(g)
  # cent_res <- closeness.currentflow(g)
  dat$cent_eig <- centr_eigen(g)$vector
  dat$cent_betw <- betweenness(g)
  
  comm <- as.data.frame(comm)
  colnames(comm) <- as.character(1:ncol(comm))
  
  dat <- cbind(dat, comm)
  rownames(dat) <- NULL
  
  dat_complete <- rbind(dat_complete, dat)
  
  write.csv(dat_complete, args[2], row.names = F)
}

write.csv(dat_complete, args[2], row.names = F)

