## This script tests and demonstrates what I called 
## Probabilistic-Stochastic-Deterministic (PDS) dynamics.

## As populations grow and decline during turnover, the way they are modelled changes.  The result is a better description of the timing of invasions. It s here demonstrated that this can be much different from the timing of invasions predicted by a simple ODE model with low invasion rates, suggesting the ODE models are not correct.

library(deSolve)

if(!exists("orig.par")){
  dev.new()
  orig.par <- par(no.readonly = TRUE)
  dev.off()
}

# n_runs <- 25
# modList <- as.list(rep(NA,n_runs))
# trajectoryList <- as.list(rep(NA,n_runs))
# 
# for(run in 15:n_runs){
# #for(run in c(1,  2,  3,  4,  6,  7,  8,  9, 10, 11, 12, 15, 16, 18, 19, 20, 22, 23, 24, 25)){
#   rseed <- 50+run
# set.seed(rseed)
#   
# print(run)

### defined further below:
# bodyMass <- 1e-6 # invasion pressure for all species
# pars  <- list(
#   C = as.matrix(1),
#   r = 1,
#   log_i = log(bodyMass) 
# )
# BRelaxed <- bodyMass

compModelLog <- function(Time, log_B, Pars) {
  with(as.list(Pars), {
    B <- exp(log_B)
    dB <- as.vector((r - C %*% B) )
    return(list(dB))
  })
}

odeLog <- function(y,times,func,parms,logIO=FALSE){
  out <- ode(y = ifelse(logIO,y,log(y)),times = times,func = func,parms = parms )
  if(! logIO) out[,2:ncol(out)] <- exp(out[,2:ncol(out)])
  return(out)
}

richness <- function(){
  length(pars$r)
}

invade <- function(){
  S <- richness()
  i <- S+1
  pars$r[i] <- rnorm(1,1,0.1)
  pars$C <- pars$C[c(1:S,S),,drop=F][,c(1:S,S),drop=F]
  pars$C[i,] <- 0.4*rbinom(n = i, size = 1,prob = 0.4)
  pars$C[,i] <- 0.4*rbinom(n = i, size = 1,prob = 0.4)
  pars$C[i,i] <- 1
  pars$i[i] <- inv
  BRelaxed[i] <<- bodyMass/10
  return(pars)
}

relax <- function(tRelax = 300){
  times <- c(0,tRelax)
  sol <-
    odeLog(y = c(B=BRelaxed),times = times,
           func = compModelLog,parms = pars )
  return(initfunc = sol[2,-1])
}

relaxX <- function(tRelax = 200){
  times <- c(0,tRelax)
  sol <-
    odeLog(y = c(B=BRelaxed),times = times,
           func = compModel,parms = pars,  atol=1e-14, rtol=1e-14, 
           method='ode45', maxsteps = 1e10 )
  return(sol[2,-1])
}

cleanup <- function(){
  extinct <- BRelaxed < 100*bodyMass
  if(any(extinct)){
    pars$r <- pars$r[!extinct]
    pars$i <- pars$i[!extinct]
    pars$C <- pars$C[!extinct,!extinct,drop=FALSE]
  }
  BRelaxed <<- BRelaxed[!extinct]
  return(pars)
}

bodyMass <- 0.01/ 100 # biomass units
inv <- 1e-5/1000 / 100  # biomass / time
mortalityRate <- 0.2 # should be of order r, otherwise arbitrary
pars  <- list(
  C = as.matrix(1),
  r = 1,
  i = inv
)

BRelaxed <- bodyMass/10
logBRelaxed <- log(BRelaxed)

while(richness() < 400){
  pars <- invade()
}

logBRelaxed <- log(BRelaxed)

mycols <- 
  hsv(runif(richness()),runif(richness(),0.9,1),runif(richness(),0.7,1))

## simulate PSD dynamics:
stepsize <- 0.2
recording_stepsize <- 100
tmax <- 200000
logB <- logBRelaxed
waiting <- rep(TRUE,length(logB))
B <- rep(0,length(logB))
nsteps <- floor(tmax/stepsize)
nrecords <- floor(tmax/recording_stepsize)
PSDtrajectory <- array(data = NA,dim = c(nrecords,length(logB)))
PSDwaiting <- array(data = NA,dim = c(nrecords,length(logB)))
s <- 0
record_number <- 0
PoissonClock <- log(runif(length(logB)))
PSDtime <- system.time(
  while(s < nsteps){
    # update step index
    s <- s + 1 
    #
    ### Situation before the step:
    B <- exp(logB) * (!waiting)
    localGrowthRate <- as.vector((pars$r - pars$C %*% B) )
    
    # Re-classify (in the moment before step)
    was_waiting <- waiting 
    waiting <- (localGrowthRate > 0) & (B < bodyMass)
    new_waiting <- which(waiting & !was_waiting)
    # But new "waiting" species might be extant by chance
    for(nw in new_waiting){
      establishmentProb <- 
        1/(1 + mortalityRate / localGrowthRate[nw])
      PoissonClock[nw] <- 
        PoissonClock[nw] + B[nw]/bodyMass * establishmentProb
      if(PoissonClock[nw] > 0){
        n <- 1
        repeat{
          PoissonClock[nw] <- PoissonClock[nw] + log(runif(1))
          if(PoissonClock[nw] <= 0) break;
          n <- n + 1
        }
        waiting[nw] <- FALSE
        logB[nw] <- log(bodyMass * n / establishmentProb)
      }
    }
    #
    ## Species that ended up not waiting any more because localGrowth 
    ## rate flipped sign now have bogus logB.  
    ## Assume it is zero before the step and compute average amount invaded:
    stopped_waiting <- which(!waiting & was_waiting)
    # Assume they stopped waiting somewhere halfway through *last* step:
    logB[stopped_waiting] <- log(pars$i[stopped_waiting] * stepsize/2)
    ##
    ### Now handle all the different cases during the step:
    ## Deterministic, Probabilistic, or Deterministic -> Probabilistic
    dB <- ( localGrowthRate + exp(log(pars$i) - logB) ) 
    # (note: dB will be overwritten for waiting blow)
    #
    ## Special handlign of Stochastic, Stochastic -> Deterministic  
    waitingList <- which(waiting)
    establishmentProb <- 
      1/(1 + mortalityRate / localGrowthRate[waitingList])
    invasionAndEstablishmentProb <-
      pars$i[waitingList]*(stepsize/bodyMass) * # <- invasion prob
      establishmentProb
    PoissonClock[waitingList] <- PoissonClock[waitingList] + invasionAndEstablishmentProb
    invades <- which(PoissonClock[waitingList] > 0)
    for(i in invades){
      ii <- waitingList[i]
      n <- 1
      repeat{
        PoissonClock[ii] <- PoissonClock[ii] + log(runif(1))
        if(PoissonClock[ii] <= 0) break;
        n <- n + 1
      }
      # We invade more than one individual because we have conditioned on establishment:
      # Suppress rare artifacts for cases where localGrowthRate is close to 0:
      logB[ii] <- min(c(0,log(bodyMass * n/establishmentProb[i])))
      dB[ii] <- dB[ii]/2 # establishment half way through step on average
      waiting[ii] <- FALSE
    }
    dB[waiting] <- 0
    #
    ## Increment step
    logB <- logB + dB * stepsize
    #
    ## Record keeping
    if(ceiling(s/nsteps*nrecords) != record_number){
      record_number <- ceiling(s/nsteps*nrecords)
      PSDtrajectory[record_number,] <- logB
      PSDwaiting[record_number,] <- waiting
      if(record_number %% 10 == 0)
      cat("\r",record_number,"                  ")
    }
  }
)
plot(rowSums(PSDwaiting),type = "l",ylim=c(0,200))
lines(rep(0,nrecords),lty=3)
matplot(exp(PSDtrajectory)*(!PSDwaiting),type = 'l',col=mycols, xlab="time ",
        ylab="Biomass")
#image(PSDwaiting)
#image(exp(PSDtrajectory) > bodyMass)
#matplot(PSDtrajectory*(!PSDwaiting)+PSDwaiting*min(PSDtrajectory-1),type = 'l',col=mycols, xlab="time ",ylab="Biomass")

# median(colSums(PSDwaiting))
# medCols <- which(colSums(PSDwaiting) %in% seq(900,1000));medCols
# focal <- 5
# plot(PSDwaiting[,focal],type = 'l')
# lines( ( PSDtrajectory[,focal] < log(bodyMass) )-0.02,type='l',col="green")
# plot(PSDtrajectory[,focal],type = 'p',pch='.',col=gray(PSDwaiting[,focal]))
# lines(rep(log(bodyMass),nrow(PSDtrajectory)),lty=3)

## simulate PSD2 dynamics:
## The method is different from that proposed by
# Purtan, R.R.P., Udrea, A., 2013. A Modified Stochastic Simulation Algorithm for Time-Dependent Intensity Rates, in: 2013 19th International Conference on Control Systems and Computer Science. Presented at the 2013 19th International Conference on Control Systems and Computer Science, pp. 365–369. https://doi.org/10.1109/CSCS.2013.101
## But the method to handle S states is similar to that proposed by
## Vestergaard, C.L., Génois, M., 2015. Temporal Gillespie Algorithm: Fast Simulation of Contagion Processes on Time-Varying Networks. PLOS Computational Biology 11, e1004579. https://doi.org/10.1371/journal.pcbi.1004579
## The main difference is that here we have a Poisson Clock for each possible event, and use an ODE solver to run it. 

logB <- logBRelaxed
S <- length(logB)
nrecords <- floor(tmax/recording_stepsize)
PSD2trajectory <- array(data = NA,dim = c(nrecords,S))
PSD2waiting <- array(data = NA,dim = c(nrecords,S))
waiting <- rep(TRUE,S)
PoissonClock <- log(runif(S))
last_output_t <- t <- 0
verbose <- FALSE
derivatives <- function(t, vars, pars){
  B <- rep(0,S)
  if(nwaiting < S){
    B[!waiting] <- exp(vars[1:(S-nwaiting)])
    localGrowthRate <- as.vector( (pars$r - pars$C %*% B) )
    establishmentProb <- 
      localGrowthRate[waitingList]/
      (localGrowthRate[waitingList] + mortalityRate)
    #if(any(establishmentProb < 0)) {print(establishmentProb);stop("Negative est. prob! ")}
    invasionAndEstablishmentRate <-
      pars$i[waitingList]*establishmentProb/bodyMass
    dlogB <- as.vector(localGrowthRate[!waiting] + exp(log(pars$i[!waiting]) - vars[1:(S-nwaiting)]))
    if(length(c(dlogB,invasionAndEstablishmentRate))!=S) stop("wrong y vector length")
    return(list(c(dlogB,invasionAndEstablishmentRate)))
  }else{# special case where nwaiting == S, all species are waiting
    localGrowthRate <- as.vector( pars$r )
    establishmentProb <- 
      localGrowthRate/
      (localGrowthRate + mortalityRate)
    invasionAndEstablishmentRate <-
      pars$i*establishmentProb/bodyMass
    return(list(invasionAndEstablishmentRate))
  }
}
rootFun <- function(t, vars, pars){
  B <- rep(0,S)
  if(nwaiting < S){
    B[!waiting] <- exp(vars[1:(S-nwaiting)])
    localGrowthRate <- as.vector( (pars$r - pars$C %*% B) )
    localGrowthRate[B > bodyMass] <- -1
    if(nwaiting){
      return(c(c(localGrowthRate,vars[(S+1-nwaiting):S])))
    }else{
      return(c(c(localGrowthRate)))
    }
  }else{# special case where nwaiting == S, all species are waiting
    localGrowthRate <- as.vector( pars$r )
    if(nwaiting){
      return(c(c(localGrowthRate,vars)))
    }else{
      return(c(c(localGrowthRate)))
    }
  }
}
PSD2time <- system.time(
  while(t < tmax){
    waitingList <- which(waiting)
    nwaiting <- length(waitingList)
    times <- unique(
      c(t, (ceiling(t/recording_stepsize):nrecords)*recording_stepsize) )
    out <- 
      ode(func = derivatives,
          y = c(logB[!waiting],PoissonClock[waitingList]),
          times = times,
          parms = pars,
          events = list(root=TRUE),
          root = rootFun,
          #          hmin = stepsize/100,
          #hini = stepsize,
          # rtol = 0,
          verbose = verbose,
          #atol = 0.01,
          #method = "lsode", mf = 10,
          #hmax = 30, 
          maxsteps = 5000000
      )
    # examine and record output
    lastRow <- nrow(out)
    if(is.null(attributes(out)$troot)){
      dropList <- -1
    }else{
      dropList <- unique(-c(1,lastRow))
    }
    recordTimes <- round(out[dropList,1]/recording_stepsize)
    if(nwaiting < S){
      logB[!waiting] <- out[lastRow, 2:(S-nwaiting+1)]
      PSD2trajectory[recordTimes, !waiting] <-
        exp(out[dropList,2:(S-nwaiting+1)])
    }
    PSD2waiting[recordTimes,] <- rep(waiting,each=length(recordTimes))
    if(nwaiting){
      PoissonClock[waitingList] <- out[lastRow, (S-nwaiting+2):(S+1)]
      PSD2trajectory[round(out[dropList,1]/recording_stepsize),
                     waiting] <- 0
    }
    if(is.null(attributes(out)$troot)){
      t <- out[lastRow,1]
      cat("\r",t,"            ")
      break; # at tmax or some error
    }
    t <- attributes(out)$troot
    ### Current situation:
    B <- exp(logB) * (!waiting)
    localGrowthRate <- as.vector((pars$r - pars$C %*% B) )
    # Recompute waiting to avoid garbling count of growth rate sign flipps
    was_waiting <- waiting 
    waiting <- (localGrowthRate > 0) & (B < bodyMass)
    iroots <- which(attributes(out)$iroot==1)
    for(r in iroots){
      if(r <= S){ # a growth rate flipped
        if(was_waiting[r]){
          waiting[r] <- FALSE
          logB[r] <- log(pars$i[r]/10) # should be log(zero)!
        }else{ # r was not waiting
          waiting[r] <- TRUE
          # set up a PoissonClock
          PoissonClock[r] <- log(runif(1))
        }
      }else{ # a Poisson clock ticked
        rr <- r - S # index of Poisson clock
        i <- waitingList[rr] # index of waiting species
        establishmentProb <- 
          1/(1 + mortalityRate / localGrowthRate[i])
        # We invade more than one individual because we have conditioned on establishment:
        # Suppress rare artifacts for cases where localGrowthRate is close to 0:
        logB[i] <- min(c(0,log(bodyMass/establishmentProb)))
        waiting[i] <- FALSE
      }
    } # end of loop over iroots, usually only one is found.
    ## check if waiting is still correct:
    B <- rep(0,S)
    B[!waiting] <- exp(logB[!waiting])
    localGrowthRate <- as.vector((pars$r - pars$C %*% B) )
    #### FIXEM: The following condition needs to be handled, not ignored!:
    if(any(localGrowthRate[waiting] < 0 )) warning("neg waiting growth rate")
    #
    ## Progress meter:
    if( (!verbose) && t > last_output_t+10*recording_stepsize){
      cat("\r",round(t/recording_stepsize)," ",iroots,"            ")
      last_output_t <- t
    }
  }
)
matplot(PSD2trajectory,type = 'l',col=mycols, xlab="time ",
        ylab="Biomass")
PSD2time
image(PSD2waiting)


ODEtrajectory <- array(data = NA,dim = c(nrecords,length(logB)))
s <- 0
logB <- logBRelaxed
ODEtime <- system.time(
  while(s < nsteps){
    s <- s + 1 
    B <- exp(logB)
    dB <- as.vector((pars$r - pars$C %*% B) )
    dB <- dB + exp(log(pars$i) - logB)
    logB <- logB + dB * stepsize
    if(ceiling(s/nsteps*nrecords) != record_number){
      record_number <- ceiling(s/nsteps*nrecords)
      ODEtrajectory[record_number,] <- logB
      if(record_number %% 10 == 0)
        cat("\r",record_number,"                  ")
    }
  }
)
matplot(exp(ODEtrajectory),type = 'l',col=mycols, xlab="time ",
        ylab="Biomass")
#matplot(ODEtrajectory/log(10),type = 'l',col=mycols, xlab="time ", ylab="Biomass")

IBMtrajectory <- array(data = NA,dim = c(nrecords,length(logB)))
s <- 0
S <- richness()
N <- rep(0,S)#rpois(S,exp(logBRelaxed)/bodyMass)
IBMtime <- system.time(
  while(s < nsteps){
    s <- s + 1 
    B <- N*bodyMass
    localGrowthRate <- as.vector((pars$r - pars$C %*% B) )
    fastDying <- (localGrowthRate < -mortalityRate)
    fullMortality <- rep(mortalityRate,S)
    fullMortality[fastDying] <- -localGrowthRate[fastDying]
    localGrowthRate[fastDying] <- -mortalityRate
    N <- rbinom(S,N,exp(-fullMortality*stepsize)) # death
    N <- N + 
      rpois(S,(exp((localGrowthRate+mortalityRate)*stepsize)-1)*N) + # birth
      rpois(S,pars$i*(stepsize/bodyMass)) # invasion
    ## Record keeping
    if(ceiling(s/nsteps*nrecords) != record_number){
      record_number <- ceiling(s/nsteps*nrecords)
      IBMtrajectory[record_number,] <- N*bodyMass
      if(record_number %% 10 == 0)
        cat("\r",record_number,"                  ")
    }
  }
)
matplot(IBMtrajectory,type = 'l',col=mycols, xlab="time ",
        ylab="Biomass")
#matplot(ODEtrajectory/log(10),type = 'l',col=mycols, xlab="time ", ylab="Biomass")

############## turnover rate analysis ####################

PSD1trajectory <- exp(PSDtrajectory)*(!PSDwaiting)

mav <- function(x,n=100){filter(x,rep(1/n,n), sides=2)}#moving average
bav <- function(x,n=300){ # blockwise average
  M <- matrix(x,nrow = n,ncol = floor(length(x)/n))
  colMeans(M)
}

c1 <- rgb(173,216,230,max = 255, alpha = 60, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 60, names = "lt.pink")
c3 <- rgb(192,255,203, max = 255, alpha = 60, names = "lt.green")
hist(log10(IBMtrajectory[IBMtrajectory>0]), freq = F,col=c1)
hist(PSDtrajectory[!PSDwaiting & PSDtrajectory>=log(bodyMass)]/log(10),freq = F,col=c2,add = T)
hist(log10(PSD2trajectory)[!PSD2waiting & PSD2trajectory>=bodyMass],freq = F,col=c3,add = T)

thresh <- 10*bodyMass

ODEchanges <- 
  rowSums(abs(apply(rbind(0,ODEtrajectory)>log(thresh), 2, diff)))/recording_stepsize
PSDchanges <- 
  rowSums(abs(apply(rbind(0,exp(PSDtrajectory)*(!PSDwaiting))>thresh, 2, diff)))/
  recording_stepsize
PSD2changes <- 
  rowSums(abs(apply(rbind(0,PSD2trajectory)>thresh, 2, diff)))/
  recording_stepsize
IBMchanges <- 
  rowSums(abs(apply(rbind(0,IBMtrajectory)>thresh, 2, diff)))/recording_stepsize

plot(mav(ODEchanges))
plot(mav(PSDchanges))
plot(mav(PSD2changes))
plot(mav(IBMchanges))

ODEalpha <- rowSums(ODEtrajectory>log(thresh))
PSDalpha <- rowSums(exp(PSDtrajectory)*(!PSDwaiting)>thresh)
PSD2alpha <- rowSums(PSD2trajectory>thresh)
IBMalpha <- rowSums(IBMtrajectory>thresh)

plot(ODEalpha,type = 'l')
plot(PSDalpha,type = 'l')
plot(PSD2alpha,type = 'l')
plot(IBMalpha,type = 'l')

mean_se <- function(v){c(mean(v),sqrt(var(v)/length(v)))}

c(median(ODEalpha),median(PSDalpha),median(PSD2alpha),median(IBMalpha))
c(mean_se(bav(ODEalpha,500)),mean_se(bav(PSDalpha,500)),mean_se(bav(PSD2alpha,500)),mean_se(bav(IBMalpha,500)))
c(mean(ODEchanges),mean(PSDchanges),mean(PSD2changes),mean(IBMchanges))
c(mean_se(bav(ODEchanges,500)),mean_se(bav(PSDchanges,500)),mean_se(bav(PSD2changes,500)),mean_se(bav(IBMchanges,500)))
#rbind(ODEtime,PSDtime,PSD2time,IBMtime)
rbind(PSDtime,PSD2time,IBMtime)

par(mfrow=c(1,2))
image(cov(t(IBMtrajectory[seq(1,2000,length.out = 1000),])),asp=1)
image(cov(t(PSD1trajectory[seq(1,2000,length.out = 1000),])),asp=1)
par(orig.par)
