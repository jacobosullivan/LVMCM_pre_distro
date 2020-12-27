#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## This script will analyse the scaling relationship between C_ij and A_ij

## Idea: C_ij ~ exp(-eta_ij*N)*A_ij where:
## eta_ij = eta for all i=/=j and eta.prime otherwise

# To run in interactive session:
# args <- c("/data/home/btx188/LVMCM_src/LVMCM/shFiles/scaling/bFileScaling.txt",
#           "/data/SBCS-RossbergLab/jack/post_doc/scaling/scaling_quad.csv")

# Example command line run:
# Rscript scaling_C_ij.R "/data/home/btx188/LVMCM_src/LVMCM/shFiles/scaling/bFileScaling.txt" "/data/SBCS-RossbergLab/jack/post_doc/scaling/scaling_quad.csv"

require(data.table)

bFile <- read.table(args[1], as.is=1)
cFile <- gsub("bMat", "cMat_reg", bFile[,1])

N <- stringr::str_extract(cFile, "N\\=\\d+")
N <- as.numeric(gsub("N=", "", N))

cFile <- cFile[order(N)]

dat <- c()

for (i in 1:length(cFile)) {
#for (i in 1:100) {

  print(cFile[i])
  
  C <- as.matrix(fread(cFile[i]))
  A <- as.matrix(fread(gsub("cMat_reg", "cMat", cFile[i])))
  
  N <- stringr::str_extract(cFile[i], "N\\=\\d+")
  N <- as.numeric(gsub("N=", "", N))
  
  P <- as.matrix(fread(gsub("cMat_reg", "params", cFile[i])))
  
  phi <- as.numeric(P[which(P[,1] == "phi"),2])
  dispL <- as.numeric(P[which(P[,1] == "dispL"),2])
  sk <- as.numeric(P[which(P[,1] == "sk"),2])
  var_e <- as.numeric(P[which(P[,1] == "var_e"),2])
  rp <- stringr::str_extract(cFile[i], "\\d+\\.mat")
  rp <- as.numeric(gsub("\\.mat", "", rp))
  
  gamma <- dim(A)[1]
  
  # filter out unexplained highly negative elements in C  
  C.ij <- c(C[upper.tri(C)], C[lower.tri(C)])
  C.ii <- diag(C)
  rm(C)
  
  q <- 0.99

  q1 <- quantile(C.ij, q)
  q0 <- quantile(C.ij, 1-q)
  C.ij <- C.ij[which(C.ij >= q0 & C.ij <= q1)]

  q1 <- quantile(C.ii, q)
  q0 <- quantile(C.ii, 1-q)
  C.ii <- C.ii[which(C.ii >= q0 & C.ii <= q1)]

  C.ij.mn <- mean(C.ij)
  C.ii.mn <- mean(C.ii)
  C.ij.var <- var(C.ij)
  C.ii.var <- var(C.ii)
  
  A.ij.mn <- mean(c(A[upper.tri(A)],
                    A[lower.tri(A)]))
  A.ij <- c(A[upper.tri(A)],A[lower.tri(A)])
  A.ii.mn <- mean(diag(A))
  
  z_ij <- C.ij.mn / A.ij.mn
  z_ii <- C.ii.mn / A.ii.mn
  
  dat_rep <- data.frame(C_ii=C.ii.mn,
                        C_ij=C.ij.mn,
		                    C_ii.var=C.ii.var,
		                    C_ij.var=C.ij.var,
                        A_ii=A.ii.mn,
                        A_ij=A.ij.mn,
                        gamma=gamma,
                        z_ij=z_ij,
                        z_ii=z_ii,
                        N = N,
                        sk=sk,
                        dispL=dispL,
		                    phi=phi,
		                    var_e=var_e,
                        rep=rp)
  dat <- rbind(dat, dat_rep)
  write.csv(dat, args[2], row.names=F)
}
