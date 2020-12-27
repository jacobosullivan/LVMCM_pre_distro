#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# To run in interactive session:
# args <- c("/data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_data.csv",
#           "/data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_analysis.csv",
#           1)

# Example command line run:
# Rscript removalExperiment_dist_metrics.R /data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_data.csv /data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_analysis_distance_metrics.csv 1


### This script with analyse the outcome of a simple node removal experiment
require(data.table)
require(vegan)

### Functions for generating presence-absence matrices
anThresh <- function(B,pars) {
  # pars: (c1, c2, emRate, dispL)
  
  Ey2 <- 2.10918; Ey <- 1.28526 # Eq. 17.8 of Rossberg (2013)
  c1 <- pars[1]
  c2 <- pars[2]
  emRate <- pars[3]
  dispL <- pars[4]
  
  # Analytic threshold
  varAij <- c1^2*c2*(1-c2);
  EAij <- c1*c2;
  Spre<-(1-EAij)^2/(Ey2*varAij);
  sig0 <- 1-emRate*(1-exp(-1/dispL))
  std_r <- sig0*(1-EAij)/(Spre*EAij*Ey);
  std_Deta <- sqrt(2/Spre)*std_r;
  std_Dbeta <- std_Deta/(1-EAij);
  threshold <- std_Dbeta
  P <- 1*(B>threshold)
  if (min(rowSums(P))==0) {
    P <- P[-which(rowSums(P)==0),]
  }
  return(P)
}

numThresh <- function(B,thresh) {
  P <- 1*(B>thresh)
  if (min(rowSums(P)==0)) {
    P <- P[-which(rowSums(P)==0),]  
  }
  return(P)
}

### path to list of unperturbed model biomass matrix file names
dat <- read.csv(args[1], as.is=1)
fileName <- as.character(dat$bFile)

thresh <- 1e-4 # numerical threshold

dat_complete <- c()

for (sim in as.numeric(args[3]):length(fileName)) {
# for (sim in 1:2) {
  print(paste("Simulation", sim, " filename ", fileName[sim]))

  ### load matrices of unperturbed assembly
  B0 <- as.matrix(fread(fileName[sim])) # pre-disturbance biomass
  PARS <- as.matrix(fread(gsub("bMat", "params", fileName[sim]))) # dispersal

  ## Remove negative biomasses for distance matrices
  B0[B0<0] <- 0
  
  ### generate unperturbed presence-absence matrices
  # pars: (c1, c2, emRate, dispL)
  pars <- c(c1 = as.numeric(as.character(PARS[which(PARS[,1]=="c1"),2])), 
            c2 = as.numeric(as.character(PARS[which(PARS[,1]=="c2"),2])), 
            emRate = as.numeric(as.character(PARS[which(PARS[,1]=="emRate"),2])),
            dispL = as.numeric(as.character(PARS[which(PARS[,1]=="dispL"),2])))
  P0_a <- anThresh(B0, pars)
  P0_n <- numThresh(B0, 1e-4)
  
  # Transpose matrices
  B0 <- t(B0)
  P0_a <- t(P0_a)
  P0_n <- t(P0_n)
  
  ## Distance matrices
  manhattan <- colMeans(as.matrix(vegdist(floor(B0/thresh), "manhattan")))
  euclidean <- colMeans(as.matrix(vegdist(B0, "euclidean")))
  canberra <- colMeans(as.matrix(vegdist(B0, "canberra")))
  bray <- colMeans(as.matrix(vegdist(B0, "bray")))
  kulczynski <- colMeans(as.matrix(vegdist(B0, "kulczynski")))
  jaccard_n <- colMeans(as.matrix(vegdist(P0_n, "jaccard")))
  jaccard_a <- colMeans(as.matrix(vegdist(P0_a, "jaccard")))
  gower <- colMeans(as.matrix(vegdist(B0, "gower")))
  altGower <- colMeans(as.matrix(vegdist(B0, "altGower")))
  morisita <- colMeans(as.matrix(vegdist(floor(B0/thresh), "morisita")))
  horn <- colMeans(as.matrix(vegdist(B0, "horn")))
  mountford <- colMeans(as.matrix(vegdist(floor(B0/thresh), "mountford")))
  raup_n <- colMeans(as.matrix(vegdist(P0_n, "raup")))
  raup_a <- colMeans(as.matrix(vegdist(P0_a, "raup")))
  binomial_n <- colMeans(as.matrix(vegdist(P0_n, "binomial")))
  binomial_a <- colMeans(as.matrix(vegdist(P0_a, "binomial")))
  chao <- colMeans(as.matrix(vegdist(floor(B0/thresh), "chao")))
  cao <- colMeans(as.matrix(vegdist(floor(B0/thresh), "cao")))
  
  dat <- data.frame(bFile=rep(fileName[sim], nrow(B0)),
                    node=1:nrow(B0),
                    manhattan=manhattan,
                    euclidean=euclidean,
                    canberra=canberra,
                    bray=bray,
                    kulczynski=kulczynski,
                    jaccard_n=jaccard_n,
                    jaccard_a=jaccard_a,
                    gower=gower,
                    altGower=altGower,
                    morisita=morisita,
                    horn=horn,
                    mountford=mountford,
                    raup_n=raup_n,
                    raup_a=raup_a,
                    binomial_n=binomial_n,
                    binomial_a=binomial_a,
                    chao=chao,
                    cao=cao)
                        
  dat_complete <- rbind(dat_complete, dat)
  
  write.csv(dat_complete, args[2], row.names = F)
}

write.csv(dat_complete, args[2], row.names = F)

########### Download data from HPC

downloaded=F
if (downloaded) {

}


