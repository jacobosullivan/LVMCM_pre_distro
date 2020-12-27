#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# To run in interactive session:
# args <- c("/data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_data.csv",
#           "/data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_analysis.csv",
#           1)

# Example command line run:
# Rscript removalExperiment.R /data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_data.csv /data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_analysis2.csv 901


### This script with analyse the outcome of a simple node removal experiment
require(data.table)
require(vegan)
require(spatgraphs)
require(centiserve)
require(igraph)

### path to list of unperturbed model biomass matrix file names
# dat <- read.csv("/data/SBCS-RossbergLab/jack/post_doc/SimulationData/N=20/removalQ_experiment/removal_experiment_data.csv")
dat <- read.csv(args[1], as.is=1)
fileName <- as.character(dat$bFile)

thresh <- 1e-4 # numerical threshold

dat_complete <- c()

for (sim in as.numeric(args[3]):length(fileName)) {
# for (sim in 1:2) {
  print(paste("Simulation", sim, " filename ", fileName[sim]))

  ### load matrices of unperturbed assembly
  B0 <- as.matrix(fread(fileName[sim])) # pre-disturbance biomass
  B0_s <- as.matrix(fread(gsub("bMat", "bMat_src", fileName[sim]))) # source sink matrix
  D <- as.matrix(fread(gsub("bMat", "dMat_n", fileName[sim]))) # dispersal
  X <- as.matrix(fread(gsub("bMat", "network", fileName[sim]))) # dispersal
  diag(D) <- 0
  PARS <- as.matrix(fread(gsub("bMat", "params", fileName[sim]))) # dispersal
  QUAD=T
  if (QUAD) {
    E <- as.matrix(fread(gsub("bMat", "envMat", fileName[sim]))) # dispersal
  }
  
  ### generate path to perturbed biomass matrices
  assembly <- stringr::str_extract(fileName[sim], "\\).*bMat\\d+")
  assembly <- gsub("\\)", "", assembly)
  outputDirectory <- paste0(sub("\\/([^\\/]*)$", "\\/", fileName[sim]), assembly, "/nodeRemoval/")
  fileList <- dir(outputDirectory)
  fileListInit <- grep("init", fileList, value=T) # pre-relaxation
  fileList <- fileList[-grep("init", fileList)] # post-relaxation
  xVec <- as.numeric(stringr::str_extract(fileList, "\\d+"))
  fileList <- fileList[order(xVec)]
  fileListInit <- fileListInit[order(xVec)]
  xVec <- xVec[order(xVec)]
  
  ### generate unperturbed presence-absence and Bray-Curtis matrices
  P0 <- B0
  P0[P0<thresh] <- 0
  P0[P0>0] <- 1
  B0[B0<0] <- 0
  logBth <- B0
  logBth[logBth<thresh] <- NA
  logBth <- log(logBth)
  BC <- as.matrix(vegdist(t(B0)))
  ED <- as.matrix(dist(E[1,], method="euclidean"))
  totalOcc0 <- sum(P0) # total occupied niches
  
  dat <- c() # store diversity metrices etc.
  comm <- c() # store alpha diversity matrix
  alpha0 <- c()
  
  for (x in 1:length(fileList)) { # loop through perturbed matrices
    
    ### Load perturbed matrix and add zero column to removed node for comparison
    B1 <- as.matrix(fread(paste0(outputDirectory,fileListInit[x])))
    B2 <- as.matrix(fread(paste0(outputDirectory,fileList[x])))
    if (xVec[x]==0) {
      B2 <- cbind(numeric(nrow(B2)), B2)
      B1 <- cbind(numeric(nrow(B1)), B1)
    } else if (xVec[x]==(ncol(B0)-1)) {
      B2 <- cbind(B2, numeric(nrow(B2)))
      B1 <- cbind(B1, numeric(nrow(B1)))
    } else {
      B2 <- cbind(B2[,1:xVec[x]], numeric(nrow(B2)), B2[,(xVec[x]+1):ncol(B2)])
      B1 <- cbind(B1[,1:xVec[x]], numeric(nrow(B1)), B1[,(xVec[x]+1):ncol(B1)])
    }
    colnames(B2) <- as.character(1:ncol(B2))
    colnames(B1) <- as.character(1:ncol(B1))
    
    ### 'Remove' site from source sink matrix
    PS1 <- B0_s
    PS1[PS1 < 0] <- 0 # source populations only
    PS1[,xVec[x]] <- numeric(nrow(PS1))
    
    ### generate new presence absence matrix
    P1 <- B1
    P1[P1<thresh] <- 0
    P1[P1>0] <- 1
    totalOcc1 <- sum(P1) # total occupied niches
    
    P2 <- B2
    P2[P2<thresh] <- 0
    P2[P2>0] <- 1
    totalOcc2 <- sum(P2) # total occupied niches
    
    alpha <- colSums(P2) # vector of local diversities
    alpha0 <- c(alpha0, sum(P0[,xVec[x]+1]))
    
    dP <- P0 - P2 # compare pres-abs
    pLost <- dP
    pLost[pLost == -1] <- 0 # populations lost
    pGain <- -dP
    pGain[pGain == -1] <- 0 # populations gained
    
    ### count populations lost which previously occupied focal node
    directPopLoss <- sum(pLost[which(P0[,xVec[x]+1]==1),])
    ### count populations lost which previously did not occupy focal node
    indirectPopLoss <- sum(pLost[which(P0[,xVec[x]+1]==0),])
    ### count populations gained which previously occupied focal node
    directPopGain <- sum(pGain[which(P0[,xVec[x]+1]==1),])
    ### count populations gained which previously did not occupy focal node
    indirectPopGain <- sum(pGain[which(P0[,xVec[x]+1]==0),])
    
    extinctS1 <- length(which(rowSums(PS1)==0))
    extinct1 <- length(which(rowSums(P1)==0))
    extinct2 <- length(which(rowSums(P2)==0))
    
    if (extinct2 > 0) {
      ### count species lost which previously occupied focal node
      directSppLoss <- sum(P0[which(rowSums(P2)==0),xVec[x]+1])
      ### count species lost which previously did not occupy focal node
      indirectSppLoss <- extinct2 - directSppLoss
    } else {
      directSppLoss <- indirectSppLoss <- 0
    }
  
    gamma <- nrow(B0) - extinct2 # new gamma diversity
    gamma.init <- nrow(B0)
    
    unique.sum <- sum(BC[xVec[x]+1,]) # total BC dissimilarity of focal node
    unique.mn <- mean(BC[xVec[x]+1,]) # mean BC dissimilarity of focal node
    envdist.sum <- sum(ED[xVec[x]+1,]) # total Environmental distance of focal node
    envdist.mn <- mean(ED[xVec[x]+1,]) # mean Environmental distance of focal node
    
    ### Compute change in Living Planet index (sum of natural log abundance resulting from patch removal)
    # sum_i [log(abundance of i in all but the focal patch) - log(total abundance of i in all patches)]
    deltaLPI <- sum(rowSums(logBth[,-(xVec[x]+1)], na.rm=T) - rowSums(logBth, na.rm=T))
    
    A <- D
    A[A!=0] <- 1
    diag(A) <- 0
    degree <- sum(A[xVec[x]+1,])
    
    DISP=F
    if (DISP) {
      connectivity <- sum(D[,xVec[x]+1]) # total immigration rate  
    } else {
      DD <- as.matrix(dist(X))
      W <- (DD * A)^-1
      W[is.infinite(W)] <- 0
      connectivity <- sum(W[,xVec[x]+1]) # simple distance weighted adjacency
    }
    
    dat_rep <- data.frame(node=xVec[x],
                          alpha.mn = mean(alpha),
                          gamma = gamma,
                          gamma.init = gamma.init,
                          extinctS1 = extinctS1,
                          extinct1 = extinct1,
                          extinct2 = extinct2,
                          totalOcc0 = totalOcc0,
                          totalOcc1 = totalOcc1,
                          totalOcc2 = totalOcc2,
                          unique.sum = unique.sum,
                          unique.mn = unique.mn,
                          envdist.sum = envdist.sum,
                          envdist.mn = envdist.mn,
                          deltaLPI = deltaLPI,
                          connectivity = connectivity,
                          degree = degree,
                          directSppLoss = directSppLoss,
                          indirectSppLoss = indirectSppLoss,
                          directPopLoss = directPopLoss,
                          indirectPopLoss = indirectPopLoss,
                          directPopGain = directPopGain,
                          indirectPopGain = indirectPopGain)
    dat <- rbind(dat, dat_rep)
    comm <- rbind(comm, alpha)
  }
  
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

########### Download data from HPC

downloaded=F
if (downloaded) {
  ## Some important differences between the high dimensional model
  ## and the abiotic niche model... look into these
  ## in particular connectivity is much less well (anti) correlated with extinctions
  ## suggesting connectivity is an important predictor of regional losses due to local perturbation
  ## in nearly neutral systems
  ## In contrast for systems subject to strong environmental gradients, it is less predictive and instead uniqueness explains most of the variance in regional extinctions
  dat <- read.csv("pCloudDrive/PostDoc/RemovalExperiment/N_20_q_removal_w.csv")
  hist(dat$extinct/dat$gamma.init)
  
  require(ggplot2)
  
  ## Does isolation of removed node predict total species loss? Yes
  ## In the case of the quadratic environmental response, yes, for sk low.
  ## In the case of broad niches, spatial processes are more important than environment
  ## and correspondingly ... extinction declines as a function of focal node connnectivity
  ## I don't quite understand this but it could lead to a concrete recommendation: if space emerges dominant in var par,
  ## then assume the most connected nodes are the least costly to remove...
  ggplot(subset(dat), aes(x=connectivity, y=extinct/gamma.init, col=factor(phi), shape=factor(dispL))) +
    geom_point() +
    scale_shape_manual(values=c(1,19))
  
  # q model:
  p1 <- ggplot(subset(dat), aes(x=degree, y=extinct/gamma.init, col=factor(sk), shape=factor(dispL))) +
    geom_point() +
    scale_shape_manual(values=c(1,19)) +
    labs(x="Degree", y="Percent species loss", col=expression(s[k]), shape="\u2113") +
    theme_bw()
  
  p2 <- ggplot(subset(dat), aes(x=connectivity, y=extinct/gamma.init, col=factor(sk), shape=factor(dispL))) +
    geom_point() +
    scale_shape_manual(values=c(1,19)) +
    labs(x="Closeness (inverse weighted degree)", y="Percent species loss", col=expression(s[k]), shape="\u2113") +
    theme_bw()
    
  ## Does uniqueness of removed node predict total species loss? Yes
  p3 <- ggplot(subset(dat), aes(x=uniqueness, y=extinct/gamma.init, col=factor(sk), shape=factor(dispL))) +
    geom_point() +
    scale_shape_manual(values=c(1,19)) +
    labs(x="Uniqueness", y="Percent species loss", col=expression(s[k]), shape="\u2113") +
    theme_bw()

  require(gridExtra)
  grid.arrange(p1,p2,p3,ncol=3)
  
  ## Does isolation of removed node predict uniqueness? Very parameter dependent
  ggplot(dat, aes(x=connectivity, y=uniqueness, col=factor(sk), shape=factor(dispL))) +
    geom_point() +
    scale_shape_manual(values=c(1,19))
  
  ## Compare the number of metapopulations directly affected by removal which lose and gain grain
  ## Correlated
  ggplot(dat, aes(x=directPopLoss, y=directPopGain)) +
    geom_point()
  
  ## Compare the number of metapopulations directly affected by removal to metapopulations which indirectly benefit
  ## Unimodal
  ggplot(dat, aes(x=directPopLoss, y=indirectPopGain)) +
    geom_point()
  
  ## Compare number of extinctions of species directly and indirectly impacted by removal
  ## Trade-off
  ggplot(dat, aes(x=directSppLoss, y=indirectSppLoss)) +
    geom_point() +
    geom_abline(slope=1, intercept=0)
  
  ## Percent indirect specices losses as a function of parameters and uniqueness
  ## weak negative trend as function of sk - the more heteregenous landscapes are less likely to 
  ## experience indirect species losses
  ## weak negative trend as function of dispL - a shorter dispersal length tends to increase the probability of species losses
  ## these two imply opposing trends with respect to mean range size...
  ggplot(dat, aes(x=uniqueness, y=indirectSppLoss/gamma.init, col=factor(sk), shape=factor(dispL), size=factor(dispL))) +
    geom_point() +
    scale_shape_manual(values=c(19,1))
  ggplot(dat, aes(x=factor(dispL), y=indirectSppLoss/gamma.init)) +
    geom_boxplot()
  
  
  ## Plot percent species loss as function of parameters
  ## Not stronly predicted, perhaps higher for low dispersal, perhaps lower for high homogeneity, perhaps intermediate peak...
  dat$percent.loss <- (1-dat$gamma/dat$gamma.init)*100
  ggplot(dat, aes(x=phi, y=percent.loss, shape=factor(dispL), col=factor(dispL))) +
    geom_point() +
    scale_x_continuous(trans="log10") +
    scale_shape_manual(values=c(19,1))
  
  
  ## Definitely a lot more species loss in quadratic model which is interesting
  ## Up to 20% loss of species for 5% loss of sites...
  dat$percent.loss <- (1-dat$gamma/dat$gamma.init)*100
  ggplot(subset(dat, dispL==1), aes(x=sk, y=percent.loss, shape=factor(dispL), col=factor(dispL))) +
    geom_point() +
    scale_x_continuous(trans="log10") +
    scale_shape_manual(values=c(19,1))
  
  ## Does environmental extremity predict uniqueness and extinctions:
  
  ## definitely 'it depends!'
  ## For sk small (broad niches) the uniqueness and correspondingly the number of extinctions
  ## are maximized at the extremes of the enviromental distribution
  ## However when sk is large and species have narrow niches, the spatial packing is tighter and there is little difference in uniqueness at the extremes
  
  ggplot(subset(dat, sk==1 & dispL == 1), aes(x=extremity, y=uniqueness, col=(factor(sk)), shape=factor(dispL))) +
    geom_point() +
    scale_shape_manual(values=c(1,19))
  
  ggplot(subset(dat, sk==1), aes(x=abs(extremity), y=extinct, col=(factor(sk)))) +
    geom_point()
  
  ## It may be worth keeping phi fixed as in the case of the quadratic environmental response
  ## and widening the scan of dispL since the neutral model focusses on this
  
  ## Seems like alpha scales with uniqueness so maybe uniqueness not well defined
  matplot(dat$uniqueness, (dat[,22:(ncol(dat)-1)]), pch=1)
  
  ## Change in alpha diversity after node removal
  comm <- as.matrix(dat[,22:ncol(dat)])
  for (sim in unique(dat$bFile)) {
    alpha0 <- dat$alpha0[which(dat$bFile == sim)]
    comm[which(dat$bFile == sim),] <- comm[which(dat$bFile == sim),] * matrix(rep(alpha0, length(which(dat$bFile == sim))), ncol=ncol(comm), byrow=T)^-1
  }
  image(comm[1001:1020,])
  image(comm[which(dat$bFile == sim),])
  
}


