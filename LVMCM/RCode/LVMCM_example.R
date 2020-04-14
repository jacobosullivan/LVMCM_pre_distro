source("LVMCM_src/LVMCM/RCode/read_results.R")
source("LVMCM_src/LVMCM/RCode/plotting_functions.R")

## Example assembly - Competitive metacommunity
system("LVMCM_src/LVMCM/build/LVMCM -o 1 testComp 1 -N 16 -Phi 10 -invMax 500 -disp 0.2 1.0 -c_ij 0.3 0.3 -var_e 0.1 -tMax 500 -z F")

## Load output into memory (see Biomass matrix file name in C++ output)
out <- readResults("SimulationData/N=16/testComp_experiment/2020-4-14/2020-4-14_testComp(500)1",1)

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
system("LVMCM_src/LVMCM/build/LVMCM -o 1 testBip 1 -N 16 -Phi 10 -invMax 500 -disp 0.2 1.0 -c_ij 0.3 0.3 -pProd 0.5 -rho 0.1 -sigma 6 -a 1e-4 -var_e 0.1 -tMax 500 -z F")

## Load output into memory (see Biomass matrix file name in C++ output)
out <- readResults("SimulationData/N=16/testBip_experiment/2020-4-14/2020-4-14_testBip(500)1",1)

## Plot regional diversity as a function of invasion event
pS <- plotAssembly(out)

## Image interation matrix
pI <- plotInteractionMat(out)

## Image growth rate matrix, whole metacommunity and single random species
pR <- plotGrowthDist(out) # (run repeatedly to sample different random species)

## Image biomass matrix, whole metacommunity and single random species
## Source and sink populations (those that do not persist when dispersal is switched off) are coloured independently
pB <- plotBiomassDist(out) # (run repeatedly to sample different random species)

