require(data.table)

readMatrix <- function(name){
  as.matrix(fread(file = name,header = F))
}

readResults <- function(prefix, rep) {
  cat("Reading output generated", 
      as.character(file.info(paste0(prefix, 'bMat', rep, '.mat'))$ctime), "\n")
  l <-
    list(
      B = NULL,
      B_src = NULL,
      B_c = NULL,
      B_c_src = NULL,
      C = NULL,
      A = NULL,
      D = NULL,
      env = NULL,
      f = NULL,
      p_inv = NULL,
      X = NULL,
      para = NULL,
      R = NULL,
      S.R = NULL,
      S = NULL,
      T = NULL
    )

    if(file.exists(paste0(prefix, 'bMat', rep, '.mat'))) {
      l$B = readMatrix(paste0(prefix, 'bMat', rep, '.mat'))
    } else {
      warning("Biomass matrix not found")
    }
    if(file.exists(paste0(prefix, 'bMat_src', rep, '.mat'))) {
      l$B_src = readMatrix(paste0(prefix, 'bMat_src', rep, '.mat'))
    }
    if(file.exists(paste0(prefix, 'bMat_c', rep, '.mat'))) {
      l$B_c = readMatrix(paste0(prefix, 'bMat_c', rep, '.mat'))
    } 
    if(file.exists(paste0(prefix, 'bMat_c_src', rep, '.mat'))) {
      l$B_c_src = readMatrix(paste0(prefix, 'bMat_c_src', rep, '.mat'))
    }
    if(file.exists(paste0(prefix, 'cMat', rep, '.mat'))) {
      l$C = readMatrix(paste0(prefix, 'cMat', rep, '.mat'))
    }
    if(file.exists(paste0(prefix, 'aMat', rep, '.mat'))) {
      l$A = readMatrix(paste0(prefix, 'aMat', rep, '.mat'))
    }
    if(file.exists(paste0(prefix, 'dMat_n', rep, '.mat'))) {
      l$D = readMatrix(paste0(prefix, 'dMat_n', rep, '.mat'))
    } else {
      warning("Dispersal matrix not found")
    }
    if(file.exists(paste0(prefix, 'environ', rep, '.mat'))) {
      l$env = readMatrix(paste0(prefix, 'environ', rep, '.mat'))
    }
    if(file.exists(paste0(prefix, 'fVec', rep, '.mat'))) {
      l$f = readMatrix(paste0(prefix, 'fVec', rep, '.mat'))
    }
    if(file.exists(paste0(prefix, 'invProb', rep, '.mat'))) {
      l$p_inv = readMatrix(paste0(prefix, 'invProb', rep, '.mat'))
    }
    if(file.exists(paste0(prefix, 'network', rep, '.mat'))) {
      l$X = readMatrix(paste0(prefix, 'network', rep, '.mat'))
    } else {
      warning("Spatial network not found")
    }
    if(file.exists(paste0(prefix, 'params', rep, '.mat'))) {
      l$para = readMatrix(paste0(prefix, 'params', rep, '.mat'))
    } else {
      warning("Model parameters not found")
    }
    if(file.exists(paste0(prefix, 'rMat', rep, '.mat'))) {
      l$R = readMatrix(paste0(prefix, 'rMat', rep, '.mat'))
    } else {
      warning("Growth rate matrix not found")
    }
    if(file.exists(paste0(prefix, 'S', rep, '.mat'))) {
      l$S.R = readMatrix(paste0(prefix, 'S', rep, '.mat'))
    }
    if(file.exists(paste0(prefix, 'sMat', rep, '.mat'))) {
      l$S = readMatrix(paste0(prefix, 'sMat', rep, '.mat'))
    }
    if(file.exists(paste0(prefix, 'tMat', rep, '.mat'))) {
      l$T = readMatrix(paste0(prefix, 'tMat', rep, '.mat'))
    }
  return(l)
}
