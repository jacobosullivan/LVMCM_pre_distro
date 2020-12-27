## Fit Metacommunity Patch Dynamics model to LVMCM output.

regularizedIncompleteGamma <- function(a,z){ # the regularized incomplete gamma function Q(a,z)
  pgamma(q = z, shape = a, lower.tail = FALSE)
}

## go to source file directly and then 3 levels up to match convention:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../../..");getwd()

source("LVMCM_src/LVMCM/RCode/read_results.R")

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
    print(mixing)
    l <- 
      uniroot(function(l)sum(nicheNumberDist(l,mixing,nMax))*(l-1)-1,
              lower=1,upper=2,extendInt="upX" )
    return(l$root)
  }
}

fitMixing <- function(data){
  data <- data/sum(data*(1:length(data)))
  goalFunction <- function(m){ ## simple sum of squared differences
    l <- consistentLambda(m,length(data))
    nND <- nicheNumberDist(l,m,length(data))
    nND <- nND/sum(nND*(1:length(data)))
    sum(abs(data - nND)^2)
  }
  o <- optim(3,goalFunction,method="Brent",lower=0.006,upper = 100)
  return(o)
}

fitMixing2 <- function(data){  ## do a maximum likelihood fit
  if(any(!(as.integer(data)==data))) stop("Input to fitMixing2 must be array of integer counts.")
  patchSum <- sum(data*(1:length(data)))
  goalFunction <- function(m){ ## simple sum of squared differences
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

## Use first matching B file and corresponding other files
#BFile <- system(paste("ls","/tmp/SimulationData/N=48/autoFluct_experiment/*/*.01bMat1.mat"),intern = TRUE)[1]
BFile <- system(paste0("ls ","/tmp/SimulationData/N=16/good_experiment/2020-5-16/2020-5-16_good*1bMat1.mat"),intern = TRUE)[1];cat("'",BFile,"'",sep = "")
#BFile <- system(paste0("ls ","/tmp/SimulationData/N\\=9/traj_experiment/2020-5-18/2020-5-18_traj\\(10000\\)1bMat1.mat"),intern = TRUE)[1];cat("'",BFile,"'",sep = "")
#BFile <- system(paste0("ls ","/tmp/SimulationData/N\\=4/traj_experiment/2020-5-18/2020-5-18_traj\\(10000\\)1bMat1.mat"),intern = TRUE)[1];cat("'",BFile,"'",sep = "")

#BFile <- "/tmp/SimulationData/N=16/good_experiment/2020-5-16/2020-5-16_good(10000)1bMat1.mat"

BFile <- "/home/axel/Projects/NERC-LargeScaleChange/code/messed_up_LVMCM_src/LVMCM/tmp/SimulationData/N=16/axel_experiment/2020-4-11/2020-4-11_axel(30)1bMat1.mat"

## Read the simulation output
fileDir <- sub("^(.*)/[^/]*$","\\1",BFile);fileDir
filePrefix <- sub("^(.*)bMat([0-9]*)[.]mat$","\\1",BFile);filePrefix
fileRep <- as.numeric(sub("^(.*)bMat([0-9]*)[.]mat$","\\2",BFile));fileRep
stepMode <- "a" # "a" or "t"
if(stepMode == "a"){
  stem <- paste0(fileDir,"/1bMat1/assembly/")
}else{
  stem <- paste0(fileDir,"/1bMat1/trajectory/")
}

Ey2 <- 2.10918; Ey <- 1.28526 # Eq. 17.8 of Rossberg (2013)

#repList <- seq(10100,10400)
repList <- 0

specs <- data.frame(alpha=NA,alpha2=NA,gamma=NA);
SS <- matrix(0,nrow=length(repList),ncol=100);
dynamicThreshold <- F
userTreshold <- 0#0.01 # set to 0 to use parameter-derived threshold
o <- readResults(prefix = filePrefix,rep = fileRep)
oldw <- getOption("warn")
options(warn = -1);para <- as.list(mapply(as.numeric,o$para[,2]));options(warn = oldw)
if(dynamicThreshold){
  thresholds <- 10^seq(-8,0,length.out = 1000)
  threshold <- 
    thresholds[which.max(mapply(function(t){aa <- colSums(o$B > t);mean(aa)/var(aa)}, thresholds))]
  if(F){
    plot(thresholds,smooth(mapply(function(t){aa <- colSums(o$B > t);mean(aa)/var(aa)}, thresholds)),type = 'l',
         xlab="Threshold",ylab="Regulations",log="x")
  }
}else{
  if(userTreshold){
    threshold <- userTreshold
  }else{
    threshold <- 
      with(para,{
        varAij <- c1^2*c2*(1-c2);
        EAij <- c1*c2;
        Spre<-(1-EAij)^2/(Ey2*varAij);
        sig0 <- 1-emRate*(1-exp(-1/dispL))
        std_r <- sig0*(1-EAij)/(Spre*EAij*Ey);
        std_Deta <- sqrt(2/Spre)*std_r;
        std_Dbeta <- std_Deta/(1-EAij);
        std_Dbeta
      })
  }
}
cat("Threshold: ",threshold)
for(i in 1:length(repList)){
  if(repList[i]>0){
    B <- readMatrix(paste0(stem,"bMat_",repList[i],".mat"))
  }else{
    B <- o$B
  }
  s <- B > threshold 
  s <- s[rowSums(s>0)>0,];
  specs[i,"alpha"] <- mean(colSums(s==1));
  specs[i,"alpha2"] <- mean(colSums(s==1)^2);
  specs[i,"gamma"] <- nrow(s);
  tt <- table(rowSums(s==1));
  SS[i,as.numeric(names(tt))] <- tt
};
colMeans(SS)[1:20]
colMeans(specs)
NMAX <- ncol(SS)
nNiches <- sum((1:ncol(SS))*colMeans(SS))
mixFit <- fitMixing(colMeans(SS)/nNiches);
mixFit <- fitMixing2(colMeans(SS));
mix <- mixFit$par
lam <- consistentLambda(mix,nMax = NMAX)
plot(colMeans(SS)/nNiches,xlim=c(1,20),col="blue",ylab="Proportion of species",xlab="Number of patches occupied by a species",ylim=c(0,max(colMeans(SS)/nNiches))*1.2)
points(nicheNumberDist(consistentLambda(mix,nMax = NMAX),mix,nMax = NMAX),pch=3,col="red")
legend(x = "topright",legend = c("Simulations","Theory"),pch = c(1,3),col=c("blue","red"))

if(F){
  ## test accuracy of p_fit (should be evenly distributed for samples from the model)
  reps <- 1000
  resampled <- matrix(rpois(n = reps*NMAX,lambda = nicheNumberDist(consistentLambda(mix,nMax = NMAX),mix,nMax = NMAX)*nNiches),nrow=NMAX)
  p_values <-
    apply(resampled,MARGIN = 2, function(v)fitMixing2(v)$p_fit)
  hist(p_values,breaks=sqrt(reps),probability = T) ## not totally flat but reasonalbe: p accurate to within a factor 2 or so
  plot(sort(p_values))
}

## Report results
c(c(mixing = mix, lambda = lam),
  colMeans(specs),betaF = as.numeric(colMeans(specs)["gamma"]/nNiches), 
  regulation = with(as.list(colMeans(specs)),alpha/(alpha2-alpha^2)),
  threshold=threshold)

