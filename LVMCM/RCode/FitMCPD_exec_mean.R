#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## args:
#### 1: Path to list of species by site matrices
#### 2: Transpose - if 1, will transpose matrix after loading
#### 3: Make binary - if 1, will threshold the matrix
#### 4: Threshold options - either: NA, numerical threshold, "D" - dynamic threshold, "A" - analytic threshold
#### 5: Fitting options - either: "ml" - maximum likelihood, "ss" - sum of squared differences

## This script will estimate the mixing rate parameter from the mean patch occupancy distribution of a replicated experiment
## I have programmed it so that either LVMCM output or empirical data can be input
## To switch between data types, different thresholding options must be selected to make the table binary
## For LVMCM, user defined, dynamic, or analytic thresholding is available

## Example command line implementation:
## Rscript FitMCPD_exec.R "/data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_data.csv" 0 1 1e-4 "ml"

## Example R session implementation:
args <- c("/data/home/btx188/LVMCM_src/LVMCM/shFiles/removal/bFile11-20.txt",
          0, 1, "A", "ml")

require(data.table)

regularizedIncompleteGamma <- function(a,z){ # the regularized incomplete gamma function Q(a,z)
  pgamma(1,rate = z,shape = a, lower.tail = FALSE)
}

plainNeutralTheory <- FALSE

nicheNumberDist <- function(lambda,mixing=0.0001,nMax=15){
  S <- rep(NA,nMax)
  S[1] <- 1/(lambda+mixing)
  expFac <- exp((lambda-1)*lambda/mixing)
  mixFac <- lambda/mixing +1
  mixProp <- mixing/(lambda+mixing)
  mixX <- (lambda-1)/mixProp
  if(plainNeutralTheory){
    if(lambda!=1)stop("lambda must be 1")
    for(n in 2:nMax){
      S[n] <- (mixing/(1+mixing))^n/mixing/n
    }
  }else{
    for(n in 2:nMax){
      S[n] <- 1/(n*lambda)*
        ( 1-(mixing/(lambda+mixing))^n -
            regularizedIncompleteGamma(n,lambda-1) +
            expFac * mixFac * mixProp^n * regularizedIncompleteGamma(n,mixX)
        )
    }
  }
  return(S)
}

consistentLambda <- function(mixing=0.0001,nMax=15){
  if(plainNeutralTheory){
    return(1)
  }else{
    l <-
      uniroot(function(l)sum(nicheNumberDist(l,mixing,nMax))*(l-1)-1,
              lower=1,upper=2,extendInt="upX")
    return(l$root)
  }
}

fitMixingSS <- function(data){
  ## simple sum of squared differences
  data <- data/sum(data*(1:length(data)))
  goalFunction <- function(m){ 
    l <- consistentLambda(m,length(data))
    nND <- nicheNumberDist(l,m,length(data))
    nND <- nND/sum(nND*(1:length(data)))
    sum(abs(data - nND)^2)
  }
  o <- optim(3,goalFunction,method="Brent",lower=0.005,upper = 100)
  return(o)
}

fitMixingML <- function(data){  
  ## maximum likelihood fit
  if(any(!(as.integer(as.matrix(data))==as.vector(as.matrix(data))))) stop("Input to fitMixingML must be array of integer counts.")
  patchSum <- sum(data*(1:length(data)))
  goalFunction <- function(m){
    l <- consistentLambda(m,length(data))
    nND <- nicheNumberDist(l,m,length(data))
    nND <- nND*(patchSum/sum(nND*(1:length(data))))
    -sum(-nND + data*log(nND)-lgamma(data+1))
  }
  o <- optim(3,goalFunction,method="Brent",lower=0.006,upper = 100)
  l <- consistentLambda(o$par,length(data))
  nND <- nicheNumberDist(l,o$par,length(data)) * patchSum
  # Compute G2 statistic with grouping of counts according to Wood (2002): 
  # https://doi.org/10.1081/STA-120015014
  c <- 1.6008
  group_size <- 0
  group_mu <- 0
  group_n <- 0
  G2Sum <- 0
  degrees_of_freedom <- -2
  for(i in 1:length(data)){
    mu <- nND[i] 
    minimumGroupSize <- (c - mu + sqrt(c^2+mu^2))/(2*mu) 
    group_size <- group_size + 1
    group_mu <- group_mu + mu
    group_n <- group_n + data[i]
    if(group_size > minimumGroupSize){
      # G2 Sum in case the next group does not finalize:
      if(i < length(data)){
        final_group_n <- group_n + sum(data[i:length(data)])
        final_group_mu <- group_mu + sum(nND[i:length(data)])
      }else{
        final_group_n <- group_n 
        final_group_mu <- group_mu
      }
      finalG2Sum <- G2Sum + 
        2 * ( ifelse( final_group_n>0, final_group_n*log(final_group_n/final_group_mu), 0) 
              - final_group_n + final_group_mu )
      # finalise current group
      G2Sum <- G2Sum +
        2 * ( ifelse( group_n>0, group_n*log(group_n/group_mu), 0) - group_n + group_mu )
      degrees_of_freedom <-
        degrees_of_freedom + 1
      group_size <- 0
      group_mu <- 0
      group_n <- 0
    }
  }
  goodness_of_fit_stats <-
    list(G2=finalG2Sum,dof=degrees_of_freedom,
         p_fit = pchisq(finalG2Sum,degrees_of_freedom, lower.tail = FALSE))  
  return(c(o,goodness_of_fit_stats))
}

makeBinary <- function(comm_cont, thresh=1e-4, P=NULL) {
  # expects a continuous species by site matrix and returns presence absence data 
  # based on various thresholding methods
  # thresh: a) user defined numerical threshold; b) "D" - dynamic; c) "A" - analytic
  # P: LVMCM parameter file
  
  suppressWarnings(
    threshNum <- as.numeric(thresh)  
  )
  threshold <- NA # variable used to threshold matrix
  if (is.na(threshNum)) { # categorical argument passed to function
    if (thresh == "D") { # select dynamic thresholding
      print("Dynamic thresholding...")
      thresholds <- 10^seq(-8,0,length.out = 1000)
      threshold <- 
        thresholds[which.max(mapply(function(t){aa <- colSums(comm_cont > t); mean(aa)/var(aa)}, thresholds))]
    } else if (thresh == "A") { # select analytic thresholding (requires dispersal and competition pars, assumed discrete random matrix)
      print("Analytic thresholding...")
      # Extract params
      Ey2 <- 2.10918; Ey <- 1.28526 # Eq. 17.8 of Rossberg (2013)
      c1 <- as.numeric(as.character(P[which(P[,1]=="c1"),2]))
      c2 <- as.numeric(as.character(P[which(P[,1]=="c2"),2]))
      emRate <- as.numeric(as.character(P[which(P[,1]=="emRate"),2]))
      dispL <- as.numeric(as.character(P[which(P[,1]=="dispL"),2]))
      
      # Analytic threshold
      varAij <- c1^2*c2*(1-c2);
      EAij <- c1*c2;
      Spre<-(1-EAij)^2/(Ey2*varAij);
      sig0 <- 1-emRate*(1-exp(-1/dispL))
      std_r <- sig0*(1-EAij)/(Spre*EAij*Ey);
      std_Deta <- sqrt(2/Spre)*std_r;
      std_Dbeta <- std_Deta/(1-EAij);
      threshold <- std_Dbeta
    }
  } else {
    print("User defined thresholding...")
    threshold <- thresh
  }
  comm <- 1*(comm_cont > threshold)
  return(comm)
}

f <- function(comm=NULL, method="ml", ss=NULL) {
  # wrapper for the fitting procedure, expects a binary species by site matrix
  # This function does not include the dynamic thresholding that can be selected for the case of simulated LVMCM data
  # The vector ss can be passed directly to the function so that the mean patch occupancy for multiple replicate assemblies can be studied
  if (method == "ml") {
    print("Maximum likelihood fitting procedure...")
  } else if (method == "ss") {
    print("Sum of squared differences fitting procedure...")
  } else {
    stop("Requested fitting procedure not recognised, select 'ml' or 'ss'")
  }
  
  if (!is.null(comm)) {
    specs <- data.frame(alpha=NA,alpha2=NA,gamma=NA)
    specs["alpha"] <- mean(colSums(comm==1));
    specs["alpha2"] <- mean(colSums(comm==1)^2);
    specs["gamma"] <- nrow(comm);  
    tt <- table(rowSums(comm==1)); ## rowSums gives a measure of species range size
    if ("0" %in% names(tt)) {
      tt <- tt[-1]  
    }
    ss <- rep(0, ncol(comm))
    ss[as.numeric(names(tt))] <- tt ## add zeros to patch occupancy distribution
  }
  
  names(ss) <- 1:length(ss)
  
  NMAX <- length(ss)
  nNiches <- sum((1:length(ss))*ss)
  if (method == "ss") {
    mixFit <- fitMixingSS(ss/nNiches)
    p_fit <- NA
  } else if (method =="ml") {
    mixFit <- fitMixingML(ss)
    p_fit <- mixFit$p_fit
  }
  mix <- mixFit$par
  lam <- consistentLambda(mix,nMax = NMAX)
  pred <- nicheNumberDist(consistentLambda(mix,nMax = NMAX),mix,nMax = NMAX)
  
  ss <- ss/nNiches
  
  dat <- data.frame(patch.occ = 1:(length(pred)),
                    obs = ss,
                    pred = pred,
                    p_fit = rep(p_fit, length(pred)),
                    mix = rep(mix, length(pred)))
  return(dat)
}

## args:
#### 1: Path to species by site matrix
#### 2: Transpose - if 1, will transpose matrix after loading
#### 3: Make binary - if 1, will threshold the matrix
#### 4: Threshold options - either: NA, numerical threshold, "D" - dynamic threshold, "A" - analytic threshold
#### 5: Fitting options - either: "ml" - maximum likelihood, "ss" - sum of squared differences

## Load matrices
dat <- rbind(read.table("/data/home/btx188/LVMCM_src/LVMCM/shFiles/removal/bFile11-20.txt")) 
dat <- data.frame(bFile=dat[,1])
ss <- c()
sk <- c()
dispL <- c()
emRate <- c()
for (i in 1:nrow(dat)) {
  print(i)
  comm <- as.matrix(fread(as.character(dat$bFile[i])))
  params <- read.table(gsub("bMat", "params", as.character(dat$bFile[i])), as.is=1)
  comm <- makeBinary(comm, thresh=args[4], params)
  tt <- table(rowSums(comm==1)); ## rowSums gives a measure of species range size  
  if ("0" %in% names(tt)) {
    tt <- tt[-1]  
  }
  ss_rep <- rep(0, ncol(comm))
  ss_rep[as.numeric(names(tt))] <- tt ## add zeros to patch occupancy distribution
  ss <- rbind(ss, ss_rep)
  sk <- c(sk, as.numeric(as.character(params[which(params[,1]=="sk"),2])))
  dispL <- c(dispL, as.numeric(as.character(params[which(params[,1]=="dispL"),2])))
  emRate <- c(emRate, as.numeric(as.character(params[which(params[,1]=="emRate"),2])))
}

ss_mn <- aggregate(ss, by=list(dispL, sk, emRate), FUN=mean)
ss_sd <- as.matrix(aggregate(ss, by=list(dispL, sk, emRate), FUN=sd)[,-c(1,2,3)])
dispL <- ss_mn[,1]
sk <- ss_mn[,2]
emRate <- ss_mn[,3]
ss_mn <- floor(as.matrix(ss_mn[,-c(1,2,3)]))

dat_summary <- c()
for (i in 1:nrow(ss_mn)) {
  dat_p <- f(comm=NULL, args[5], ss=ss_mn[i,]) # MCPD fitting procedure
  NMAX <- length(ss_mn[i,])
  nNiches <- sum((1:length(ss_mn[i,]))*ss_mn[i,])
  dat_row <- c(dispL[i], sk[i], emRate[i],
               dat_p$obs,
               ss_sd[i,]/nNiches,
               dat_p$pred,
               dat_p$mix[1],
               dat_p$p_fit[1])
  dat_summary <- rbind(dat_summary, dat_row)
}

dat_summary <- as.data.frame(dat_summary)
names(dat_summary) <- c("dispL", 
                        "sk", 
                        "emRate",
                        paste0("obs.mn",1:20),
                        paste0("obs.sd",1:20),
                        paste0("pred",1:20),
                        "mix",
                        "p_fit")

write.csv(dat_summary, "/data/SBCS-RossbergLab/jack/post_doc/removal/removal_experiment_fitMCPD_ml_A_summary_degree_norm.csv", row.names=F) ## not executable

