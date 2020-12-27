# based on https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.13398
L <- 512*2
N <- L*L
scaleFac <- 20  # leaves sigmaEff/sigmaVirgin invariant, as long as sigmaVirgin << phiBinary
phiBinary <- 1*scaleFac
percentageReserve <- 30
sigmaVirgin <- 2
nuVirgin <- nuEst/sigmaEst^2*sigmaVirgin^2/scaleFac#10*4e-6/scaleFac

circ <- function(n) return((n-1+L) %% L + 1)
circ2 <- function(n) return((n-1+2*L) %% (2*L) + 1)
circ(257)

distFromOrig <- outer(0:(L-1),0:(L-1),function(x,y) sqrt(x^2+y^2))
AutoCov <- exp(-distFromOrig/(2*phiBinary))
zeros <- rep(0,L)
AutoCov4 <- 
  rbind(cbind(AutoCov,zeros,AutoCov[,rev(2:L)]),
        c(zeros,zeros),
        cbind(AutoCov,zeros,AutoCov[,rev(2:L)])[rev(2:L),] )
#image(AutoCov4[circ2(1:(2*L)-L),circ2(1:(2*L)-L)],asp=1)
AutoCov4Hat <- fft(AutoCov4)
Sig4Hat <- sqrt(AutoCov4Hat/(4*N))
sampleField4 <- Re(fft(Sig4Hat * (rnorm(N)+1i*rnorm(N)),inverse = T))
sampleField <- sampleField4[2*(1:L),2*(1:L)]
#image(sampleField,asp=1)
#var(as.vector(sampleField))
cor(as.vector(sampleField),as.vector(sampleField[circ((1:L)+phiBinary),]))*exp(1)
cor(as.vector(sampleField),as.vector(sampleField[,circ((1:L)+phiBinary)]))*exp(1)
reserveThresh <- quantile(sampleField,probs = (percentageReserve/100))
reserve <- sampleField < reserveThresh
image(reserve,asp=1)

nSteps <- round(1/nuVirgin)
nReps <- 100
travelDist2 <- rep(NA,nReps)
isReserve <- function(x) {y <- circ(x);reserve[y[1],y[2]]}
for(rep in 1:nReps){
  cat("Iteration:",rep,"               \r")
  repeat{
    startPoint <- circ(runif(2,0,L))
    if(isReserve(startPoint)) break;
  }
  currentPoint <- startPoint
  for(s in 1:nSteps){
    repeat{
      nextPoint <- currentPoint + rnorm(2,0,sigmaVirgin)
      if(isReserve(nextPoint)) break;
    }
    currentPoint <- nextPoint
  }
  travelDist2[rep] <- sum((startPoint-currentPoint)^2)   
}
mean_se <- function(x) c(mean(x),sd(x)/sqrt(length(x)))
sigmaEff <- mean_se(sqrt(travelDist2))*sqrt(2/nSteps/pi)
print(sigmaEff/sigmaVirgin)
