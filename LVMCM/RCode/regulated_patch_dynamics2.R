alpha_div <- 16
N <- 256/4

occupancy_list <- apply(t(replicate(N,as.array(1:alpha_div))),1,as.list)
occupancy_table <- t(replicate(N,as.array(1:alpha_div)))
last_i <- alpha_div

gamma_from_table <- function(o_table){
  length(unique(as.vector(o_table)))
}
gamma_from_table(occupancy_table)

steps <- 5000 * 2 * 2 * 2 * 4
record <- data.frame(step=1:steps,S=NA)
mixing <- 14.78535
for(step in 1:steps){
  S <- gamma_from_table(occupancy_table)
  if(last_i %% 1000 == 0){
    cat("last_i: ", last_i, ", gamma = ", S, ", oldest species: ", last_i-min(occupancy_table), " \r",sep="")
  }
  record$S[step] <- S
  new_i <- last_i+1
  invaded_patch <- sample.int(ncol(occupancy_table),1)
  invade_other_probability <- alpha_div/S
  invade_patches <- as.logical(rbinom(ncol(occupancy_table),1,invade_other_probability))
  invade_patches[invaded_patch] <- TRUE
  for(patch in which(invade_patches)){
    occupancy_table[patch,sample.int(alpha_div,1)] <- new_i 
  }
  mixers <- occupancy_table[rbinom(length(occupancy_table),1,mixing/length(occupancy_table))==1]
  # for(m in mixers){
  #   candidate_patches <- which(!apply(occupancy_table, 1, function(v) return(m %in% v)))
  #   if(length(candidate_patches) == 0) {cat("Species",m,"is everywhere!                      \n"); next;}
  #   patch <- candidate_patches[sample.int(length(candidate_patches),1)]
  #   occupancy_table[patch,sample.int(alpha_div,1)] <- m
  # }
  for(m in mixers){
    target_patch <- sample.int(nrow(occupancy_table),1)
    if(m %in% occupancy_table[target_patch,]) next;
    occupancy_table[target_patch,sample.int(alpha_div,1)] <- m
  }
  last_i <- new_i
}
plot(record$S,type = "l")
c(mean(record$S[-(1:200)]),0.629084*alpha_div*ncol(occupancy_table),mean(record$S[-(1:200)])/(0.629084*alpha_div*N))

occupation_list_simp  <- table(as.vector(occupancy_table))

hi <- hist(occupation_list_simp,breaks=0.5+seq(0,max(occupation_list_simp)),plot = F)
hi$counts/length(occupancy_table)
mixFit <- fitMixing(hi$counts)
mix <- mixFit$par;mix
NMAX <- length(hi$counts)

plot(hi$counts/length(occupancy_table),xlim=c(1,20),col="blue",ylab="Proportion of species",xlab="Number of patches occupied by a species")
points(nicheNumberDist(consistentLambda(mix,nMax = NMAX),mix,nMax = NMAX),pch=3,col="red")
#points(nicheNumberDist(consistentLambda(100,nMax = NMAX),mix,nMax = NMAX),pch=3,col="red")
legend(x = "topright",legend = c("Simulations","Theory"),pch = c(1,3),col=c("blue","red"))

#occupancy_table <- occupancy_table[1:(N/2),]

#N alpha  alpha2 gamma    lambda   mixing   sumSquares
#1 16     15.65  248.775  61.2     5.055082 14.78535 0.0001268299

Opedal2020 <- read.csv("/tmp/Y.csv")
Ope_hi <- hist(colSums(Opedal2020),breaks = (-1:max(colSums(Opedal2020)))+0.5,plot = F)
Ope_hi$counts <- Ope_hi$counts[-1]
sum(seq(Ope_hi$counts)*Ope_hi$counts)/sum(Ope_hi$counts)
Ope_mix <- fitMixing(Ope_hi$counts)$par;Ope_mix
NMAX=150
plot(Ope_hi$counts/sum(seq(Ope_hi$counts)*Ope_hi$counts),xlim=c(1,NMAX),col="blue",ylab="Proportion of species",xlab="Number of patches occupied by a species",ylim=c(0,6e-4))
Ope_mix <- 40
points(nicheNumberDist(consistentLambda(Ope_mix,nMax = NMAX)/2,Ope_mix,nMax = NMAX)/50,pch=3,col="red")

