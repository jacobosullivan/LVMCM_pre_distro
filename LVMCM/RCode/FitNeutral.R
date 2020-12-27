### Fit neutral theory following https://science.sciencemag.org/content/295/5555/666

## go to source file directly and then 3 levels up to match convention:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../../..");getwd()

source("LVMCM_src/LVMCM/RCode/read_results.R")
source("LVMCM_src/LVMCM/RCode/PrestonPsi.R")

## Use first matching B file and corresponding other files
#BFile <- system(paste("ls","/tmp/SimulationData/N=48/autoFluct_experiment/*/*.01bMat1.mat"),intern = TRUE)[1]
BFile <- system(paste0("ls ","/tmp/N=100/consAreaOM_experiment/*/*_consAreaOM*6bMat0.mat"),intern = TRUE)[1]

## Read the simulation output
filePrefix <- sub("^(.*)bMat([0-9]*)[.]mat$","\\1",BFile);filePrefix
fileRep <- as.numeric(sub("^(.*)bMat([0-9]*)[.]mat$","\\2",BFile));fileRep
o <- readResults(prefix = filePrefix, rep=fileRep)

## Distances and similarity for each patch pair
distMat <- as.matrix(dist(o$X))
diag(distMat) <- 1  ## because we want to check for and eliminate off-diagonal zeros
bad <- colSums(distMat<1e-6)>0
B <- o$B[,!bad]
distMat <- distMat[!bad,!bad]
Btot <- colSums(B)
props <- B/rep(Btot,each=nrow(B))
spatialSimpson <- t(props) %*% props  # The spatial Simpson index is actually defined as 1 - this!
lt <- lower.tri(distMat)  ## Matrices are symmetric, so we need only the lower triangle

## Define a simple moving average filter:
mav <- function(v,n=3*ncol(B)){
  filter(v,filter = rep(1/n,n))
} 
## Check for linear relations
plot(mav(sort(log(distMat[lt]))),mav(spatialSimpson[lt][order(distMat[lt])]),xaxt='n',xlab="Distance",ylab="Conspecific probability F(r)",ylim=c(0,3*median(mav(spatialSimpson[lt][order(distMat[lt])]),na.rm = T)),type='l')
axis(side = 1,at = log(c(10^(-1:3),3*10^(-1:3))),labels = c(10^(-1:3),3*10^(-1:3)))
points(log(distMat[lt]),spatialSimpson[lt],pch=".",col="blue")
sel <- log(distMat[lt])>log(2)
mod <- lm(spatialSimpson[lt][sel]~ log(distMat[lt][sel]))
summary(mod)
abline(mod,col="red")

## Euler's Gamma constant
EulerGamma <- -digamma(1)  # = 0.57721566490153286061

## Name the coefficients of the regression simply A and B:
AB <- mod$coefficients; names(AB) <- c("A","B")

# Compute log(nu) for given sigma, AB, and rho
lognu <- function(sigma,AB=parent.frame()$AB,rho=parent.frame()$rho){
  with(as.list(AB),
         lognu <- 2*(pi * rho * sigma^2 + 1/B)
  )
}

# Goal function for finding sigma:
uniGoal <- function(sigma,AB=parent.frame()$AB,rho=parent.frame()$rho){
  if(is.null(AB)) AB=sys.frame()$AB
  if(is.null(rho)) rho=sys.frame()$rho
  goal <-
    with(as.list(AB),
      EulerGamma  - log(2)/2 + lognu(sigma,AB,rho)/2 - log(sigma) - A/B
  )
  return(goal)
}

rho = 1e4  # density of ind. = number of ind. per patch = 1/(species detection threshold)
## check: round(log10(1/(min(apply(B, 1, max)))))

# sigmas <- exp(seq(log(1e-10),log(10),length.out = 1000))
# plot(sigmas,mapply(uniGoal, sigmas),type = 'l',ylim=c(-1,1)*10)
sigmaEst <- uniroot(f = uniGoal, AB=AB, rho=rho,interval = c(1e-10,10),extendInt = "upX")$root
nuEst <- exp(lognu(sigmaEst))
c(sigmaEst,nuEst)
c(sigmaEst^2*rho,nuEst*rho) # rho-(nearly)-independet values (standard units where rho == 1)
## travel distance:
sigmaEst/sqrt(nuEst)

## Compare richness predicted by Neutral theory with observed richness:
S_contig <- Scontig(A = nrow(o$X),sigma = sigmaEst,nu = nuEst,rho = rho)
c(S_contig, nrow(B) )

sigmaEffRatio30 <- 0.3894111

S_reserve <- Sreserve(A = nrow(o$X),sigmaEff = sigmaEffRatio30*sigmaEst,nu = nuEst,
                      percentageReserve = 30,rho = rho)
S_reserve/S_contig
S_reserve/S_contig * nrow(B)
