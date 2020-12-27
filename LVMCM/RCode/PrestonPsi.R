# Preston function according to O’Dwyer & Cornell (2018) as cited by Thompson et al. (2019, https://doi.org/10.1111/ele.13398):

PrestonPsi <- function(hatA,nu){
  nue <- -nu*log(nu)/(1-nu)
  I0sAp <- besselI(sqrt(hatA/pi),0)
  I1sAp <- besselI(sqrt(hatA/pi),1)
  K0snAp <- besselK(sqrt(nue*hatA/pi),0)
  K1snAp <- besselK(sqrt(nue*hatA/pi),1)
  return(nue*hatA + 
           2*sqrt(hatA*pi)*(1-nue)/ 
           (1/sqrt(nue)*K0snAp/K1snAp+I0sAp/I1sAp) )
}

#
Scontig <- function(A,nu,sigma,rho=1){
  rho*sigma^2*PrestonPsi(A/sigma^2,nu)
}
Sreserve <- function(A,nu,sigmaEff,percentageReserve,rho=1){
  h <- percentageReserve/100
  Scontig(A*h,nu,sqrt(h)*sigmaEff,rho)
}

# test:
c(Scontig(A = 500^2,nu = 1e-4,sigma = 16),1654)


