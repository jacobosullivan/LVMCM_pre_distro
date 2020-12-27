## go to source file directly and then 3 levels up to match convention:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../../..");getwd()

library(digest)
source("LVMCM_src/LVMCM/RCode/read_results.R")

## Generate outputs after subsequent invasions (one by one), and save them somewhere, the following is just and example.

## There is no deep reason for using indices 2, 3 rather than 1, 2. Anyway, this is highly experimental.

hist(log10(B2),breaks = 400,xlim=c(-6,2)) ## look at abundance distribution
# not sure what the best threshold is, regulation strongest at much larger values
# For comparison with FitMCPD, threshold should be chosen identical in both cases.
threshold <- 0.01 

B2 <- readMatrix("/tmp/N=100/consAreaOM_experiment/2020-4-29/1bMat0/assembly/bavMat_8077.mat")
R2 <- readMatrix("/tmp/N=100/consAreaOM_experiment/2020-4-29/1bMat0/assembly/rMat_8077.mat")
B3 <- readMatrix("/tmp/N=100/consAreaOM_experiment/2020-4-29/1bMat0/assembly/bavMat_8078.mat")
R3 <- readMatrix("/tmp/N=100/consAreaOM_experiment/2020-4-29/1bMat0/assembly/rMat_8078.mat")
speciesNames2 <- apply(R2, 1, digest, algo="md5")  ## Generate species IDs based on R fields
rownames(B2) <- rownames(R2) <- speciesNames2 
speciesNames3 <- apply(R3, 1, digest, algo="md5")  ## Generate species IDs based on R fields
rownames(B3) <- rownames(R3) <- speciesNames3 
commonNames <- intersect(speciesNames2,speciesNames3)
cor(as.vector(B2[commonNames,]),as.vector(B3[commonNames,])) ## how much metacommunity chaged between steps
#image(cor(t(R2),t(R3))>0.9,asp=1)
invaders <- speciesNames3[which(!speciesNames3 %in% speciesNames2)];invaders
which(colSums(cor(t(R2),t(R3))>0.9)==0) ## consistency of this and previous line confirms md5-based naming
extinguished <- speciesNames2[which(!speciesNames2 %in% speciesNames3)];extinguished
which(rowSums(cor(t(R2),t(R3))>0.9)==0)
communities2 <- apply(B2>threshold,2,function(v) sort(speciesNames2[v]))
communities3 <- apply(B3>threshold,2,function(v) sort(speciesNames3[v]))
#mapply(identical,communities2,communities3)
#NumberOfChanges <- function(s1,s2) 2*length(intersect(s1,s2))/(length(s1)+length(s2))
NumberOfChanges <- function(s1,s2) (length(s1)+length(s2))/2-length(intersect(s1,s2))
sum(mapply(NumberOfChanges,communities2,communities3)) # total number of changes... = lambda + mixing ?
hist(mapply(NumberOfChanges,communities2,communities3),breaks = 0:100-0.5)      
image(matrix(mapply(NumberOfChanges,communities2,communities3),nrow=10),asp=1) # ..and where they happened   
#image(matrix(mapply(NumberOfChanges,communities2,communities3) %in% c(1:3),nrow=10),asp=1) # small changes   
#image(matrix(mapply(function(l) sum(invaders %in% l),communities3),nrow=10),asp=1) # invasions here   
#image(matrix(mapply(function(l) -sum(extinguished %in% l),communities2),nrow=10),asp=1) # extinctions here

## Here a threshold-free analysis, to see if NumberOfChanges is representative:
comCor <- rep(NA,ncol(B2))
for(p in seq(comCor)){
  comCor[p] <- cor(as.vector(B2[commonNames,p]),as.vector(B3[commonNames,p]))
}
image(-matrix(comCor,nrow=10),asp=1)   
t(matrix(comCor,nrow=10))

## Check if we have a "Clementian" transition somewhere and where it is:
sort(comCor)
sel <- which.min(comCor);sel
plot(as.vector(B2[commonNames,sel]),as.vector(B3[commonNames,sel]),asp=1)
plot(log10(as.vector(B2[commonNames,sel])),log10(as.vector(B3[commonNames,sel])),xlim=c(-4,2),ylim=c(-4,2))

