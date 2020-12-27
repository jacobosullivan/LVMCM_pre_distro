source("LVMCM_src/LVMCM/RCode/read.results.R")
source("LVMCM_src/LVMCM/RCode/plotting_functions")

## Example assembly - Competitive metacommunity
system("LVMCM_src/LVMCM/build/LVMCM -o 1 testComp 1 -N 16 -Phi 10 -invMax 300 -disp 0.2 1.0 -c_ij 0.3 0.3 -var_e 0.1 -tMax 500 -z F")

## Load output into memory (see Biomass matrix file name in C++ output)
out <- read.results("SimulationData/N=16/testComp_experiment/2020-4-14/2020-4-14_testComp(300)1",1)

## Plot regional diversity as a function of invasion event
pS <- plotAssembly(out)

## Image interation matrix
pI <- plotInteractionMat(out)

## Image growth rate matrix, whole metacommunity and single random species
pR <- plotGrowthDist(out) # (run repeatedly to sample different random species)

## Image biomass matrix, whole metacommunity and single random species
## Source and sink populations (those that do not persist when dispersal is switched off) are coloured independently
pB <- plotBiomassDist(out) # (run repeatedly to sample different random species)



## Example assembly - Bipartite metacommunity
-o 1 test 1 -N 16 -Phi 10 -invMax 1000 -disp 0.0001 0.5 -c_ij 0.3 0.3 -C F -pProd 0.5 -rho 0.1 -F 3 -a 1e-7 -var_e 0.01 -tMax 500 -z F -O T

## Example assembly - Competitive metacommunity
system("LVMCM_src/LVMCM/build/LVMCM -o 1 testComp 1 -N 16 -Phi 10 -invMax 300 -disp 0.2 1.0 -C F -pProd 0.5 -rho 0.1 -F 3 -a 1e-7 -var_e 0.1 -tMax 500 -z F")

## Load output into memory (see Biomass matrix file name in C++ output)
out <- read.results("SimulationData/N=16/testComp_experiment/2020-4-14/2020-4-14_testComp(300)1",1)