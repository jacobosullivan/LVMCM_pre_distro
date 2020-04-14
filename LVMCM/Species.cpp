////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// The parallelizable Lotka-Volterra Metacommunity assembly Model (pLVMCM) ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Jacob Dinner O'Sullivan -- j.l.dinner@qmul.ac.uk | j.osullivan@zoho.com ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// A metacommunity assembly model parallelised via domain decomposition using MPI ////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
    Copyright (C) 2020  Jacob D. O'Sullivan, Axel G. Rossberg

    This file is part of pLVMCM

    pLVMCM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

//////////////////////////////////// PRE-RELEASE VERSION, PLEASE DO NOT DISTRIBUTE /////////////////////////////////////

/*
 * This class contains the members and methods required for generating and storing the biotic component of the model,
 * the species pool. Methods include invasion (sampling) of new species and extinction (removal) of exluded species.
 */

#include <iostream>
#include <armadillo>
#include <boost/random.hpp>
#include <random>

#include "Species.h"
#include "Topography.h"
#include "LVMCM_rng.h"

using namespace std;
using namespace arma;
using namespace boost;

mat Species::genRVec(rowvec zVecExt) {

    // summary:
        // generate invader spatially autocorrelated growth rate vector
        // modelled using a Gaussian random field generated via spectral decomposition (https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution)
        // environment modelled implicitly 'through the eyes' of species
        // species growth rate distributions are independent

    // arguments:
        // zVecExt - vector of random variables passed from outside, required for modelling time dependent growth rates

    // required members:
        // topo.envVar - number of explicitly modelled environmental variables
        // topo.var_e - variance of the environmental distribution(s)
        // topo.phi - spatial autocorrelation length of the environmental distribution(s)
        // topo.distMat - euclidean distance matrix

    // output:
        // r_i - spatially autocorrelated growth rate vector

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Generate and decompose covariance matrix /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double mu; // mean of growth rate distribution
    if (zVecExt.n_cols == 0) {
        mu = 1.0;
    } else {
        mu = 0.0; // mean set to zero for OU process used to model time dependent growth rates
    }

    vec zVec;
    if (zVecExt.n_cols == 0) {
        zVec.randn(topo.no_nodes); // generate random normal variables
    } else {
        zVec = zVecExt.t();
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Sample from Gaussian random field /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vec r_i;
    r_i = mu + (topo.sigEVec * topo.sigEVal) * zVec; // sample from random field

    if (topo.network.n_rows != topo.no_nodes) { // subset by subdomain
        r_i = r_i.rows(topo.subdom_nodes);
    }

    return(r_i.t()); // return row vector of intrisic growth rates
}

void Species::ouProcess() {

    // summary:
        // simulated abiotic turnover by updating state of matrix of environmental fluctuations
        // efMat updated using spatially autocorrelated ouMat (Ohrstein-Uhlenbeck random process)
        // OU process: Z(t) = (Z(t-1) + sigma_t.z) / sqrt(1 + sigma_t^2)
        // Output of OU random walk passed to genRVec
        // rMat(t) = rMat(0) + efMat(t)
        // if mu != 0 (moving average process), each species responds positively or negatively to abiotic turnover with probability 0.5

    // arguments:
        // sigma_t - temporal autocorrelation (standard deviation of white noise)
        // mu - absolute value of mean of random process, set to zero for stationary OU process

    // output:
        // efMat - current state of environmental fluctuation
        // ouMat process current stat of OU process

    // external function calls:
        // genRVec()

    if (ouMat.n_rows == 0) { // initialize matrices
        ouMat.zeros(bMat_p.n_rows, bMat_p.n_cols);
        efMat.zeros(bMat_p.n_rows, bMat_p.n_cols);
    }

    mat Zmat;
    Zmat.randn(efMat.n_rows, efMat.n_cols); // generate matrix of standard normal random variables
    double mu = 0;
    ouMat = (ouMat + mu + sigma_t * Zmat) / sqrt(1 + sigma_t*sigma_t); // update OU process

    for (int i=0; i<rMat.n_rows; i++) { // use output of OU process to seed random field simulator
        efMat.row(i) = genRVec(ouMat.row(i));
    }
}

mat Species::genRVecERF() {

    // summary:
        // each species is assigned a random environmental response vector which is used to generate r_i by multiplication with the environmental variables
        // tolerances are sample from a spherical distribution so that R is normally distributed

    // required members:
        // topo.envVar - number of explicitly modelled environmental variables
        // topo.envMat - matrix of environmental distribution vectors, dimensions (envVar x no_nodes)

    // output:
        // updates to tMat - matrix of species environmental response functions
        // r_i - spatially autocorrelated growth rate vector

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Sample specific environmental tolerance coefficients ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    rowvec tol;
    tol.randn(topo.envVar);
    rowvec u = pow(tol,2);
    tol = tol / sqrt(sum(u)) * sqrt(topo.envVar);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Generate local internal growth rates ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    rowvec r_i(topo.no_nodes);
    r_i = 1 + (tol * topo.envMat) / sqrt(2 * topo.envVar);

    if (tMat.n_rows == 0) {
        tMat.set_size(1,topo.envVar);
        tMat = tol;
    } else {
        tMat = join_vert(tMat, tol);
    }

    return (r_i);
}

void Species::invade(int trophLev, bool invTest) {

    // summary:
    // generate invader by sampling numerical traits and updating model matrices

    // arguments:
    // port - node to which invader is introduced
    // trophLev - select producer (0) or consumer (1) species
    // invTest - if preparing matrices for invader testing, coupling of invader->residents is suppressed

    // required members:
    // c1, c2 - interspecific competition parameters
    // rho - consumer mortality rate
    // sigma - standard deviation log-normal attack rate distribution
    // alpha - base attack rate
    // pProducer - probability of invading a producer species, set to zero for competitive model
    // emRate - emigration rate
    // dispL - dispersal length
    // prodComp - select coupled (T), uncoupled (F) producer dynamics
    // discr_c_ij - select discrete (T), continuous (F - beta distribution) competition coefficients

    // output:
    // updates to matrices bMat_p, rMat, cMat, bMat_c, aMat

    double inv = 1e-6; // invasion biomass

    if (trophLev == 0) { // generate producer

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// Add row to bMat_p ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (bMat_p.n_rows == 0) {
            bMat_p.set_size(1, topo.network.n_rows);
            bMat_p.fill(inv);
        } else {
            bMat_p.resize(bMat_p.n_rows + 1, topo.network.n_rows);
            bMat_p.row(bMat_p.n_rows - 1).fill(inv);
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Add row to rMat ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (topo.envVar == 0) { // implicitly modelled environment
            if (rMat.n_rows == 0) {
                rMat = genRVec();
            } else {
                rMat = join_vert(rMat, genRVec());
            }
        } else { // explicitly modelled environment
            if (rMat.n_rows == 0) {
                rMat = genRVecERF();
            } else {
                rMat = join_vert(rMat, genRVecERF());
            }
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Add row to emMat ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (emRate < 0) { // sample emigration rate for uniform distribution
            vec emRate_i(1);
            emRate_i.randu();
//            emRate_i(0) = 0.1;
            if (emMat_p.n_rows == 0) {
                emMat_p = emRate_i;
            } else {
                emMat_p = join_vert(emMat_p, emRate_i);
            }
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Add row to aMat ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (pProducer < 1.0) { // only construct aMat in bipartite model
            if (aMat.n_rows == 0) { // a_ij = a_0 * exp(sigma * epsilon_ij)
                rowvec a_pc = randn<rowvec>(bMat_c.n_rows);
                a_pc *= sigma;
                a_pc = alpha * exp(a_pc);

                rowvec a_connect_rand = randu<rowvec>(bMat_c.n_rows);
                a_pc.elem(find(a_connect_rand > a_connectance)).zeros();

                aMat = a_pc;
            } else {
                aMat.resize(aMat.n_rows + 1, aMat.n_cols);
                rowvec a_pc = randn<rowvec>(aMat.n_cols);
                a_pc *= sigma;
                a_pc = alpha * exp(a_pc);

                rowvec a_connect_rand = randu<rowvec>(aMat.n_cols);
                a_pc.elem(find(a_connect_rand > a_connectance)).zeros();

                aMat.row(aMat.n_rows - 1) = a_pc;
            }
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Add row to cMat ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (prodComp) { // generate non-zero interspecific competition terms
            if (discr_c_ij) { // sample from discrete distribution
                if (cMat.n_rows == 0) {
                    cMat.set_size(1, 1);
                    cMat(0, 0) = 1.0;
                } else {
                    boost::random::binomial_distribution<int> distribution(1,c2);
                    cMat.resize(cMat.n_rows + 1, cMat.n_cols + 1);
                    for (int i = 0; i < cMat.n_cols; i++) {
                        cMat(cMat.n_rows - 1, i) = c1 * distribution(LVMCM_rng::boost_rng);
                    }
                    if (invTest) { // set col / diag to zero for invader testing
                        for (int i = 0; i < cMat.n_rows; i++) {
                            cMat(i, cMat.n_cols - 1) = 0;
                        }
                        cMat(cMat.n_rows - 1, cMat.n_cols - 1) = 0.0;
                    } else {
                        for (int i = 0; i < cMat.n_rows; i++) {
                            cMat(i, cMat.n_cols - 1) = c1 * distribution(LVMCM_rng::boost_rng);
                        }
                        cMat(cMat.n_rows - 1, cMat.n_cols - 1) = 1.0;
                    }
                }
            } else { // sample from continuous (beta) distribution
                if (cMat.n_rows == 0) {
                    cMat.set_size(1, 1);
                    cMat(0, 0) = 1.0;
                } else {
                    typedef boost::random::mt19937 RandomNumberGenerator;
                    typedef boost::random::beta_distribution<> BetaDistribution;
                    typedef boost::variate_generator<RandomNumberGenerator&, BetaDistribution> Generator;
                    BetaDistribution distribution(c1,c2);
                    Generator getRandomNumber(LVMCM_rng::boost_rng,distribution);
                    cMat.resize(cMat.n_rows + 1, cMat.n_cols + 1);
                    for (int i = 0; i < cMat.n_cols; i++) {
                        cMat(cMat.n_rows - 1, i) = getRandomNumber();
                    }
                    if (invTest) {
                        for (int i = 0; i < cMat.n_rows; i++) {
                            cMat(i, cMat.n_cols - 1) = 0;
                        }
                        cMat(cMat.n_rows - 1, cMat.n_cols - 1) = 0.0;
                    } else {
                        for (int i = 0; i < cMat.n_rows; i++) {
                            cMat(i, cMat.n_cols - 1) = getRandomNumber();
                        }
                        cMat(cMat.n_rows - 1, cMat.n_cols - 1) = 1.0;
                    }
                }
            }
        }

    } else if (trophLev == 1) { // generate consumer

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Add row to bMat_c /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (bMat_c.n_rows == 0) {
            bMat_c.set_size(1, topo.network.n_rows);
            bMat_c.zeros();
            bMat_c.fill(inv);
        } else {
            bMat_c.resize(bMat_c.n_rows + 1, topo.network.n_rows);
            bMat_c.row(bMat_c.n_rows - 1).fill(inv);
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// Add col to aMat //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (aMat.n_cols == 0) {
            vec a_pc = randn<colvec>(bMat_p.n_rows);
            a_pc *= sigma;
            a_pc = alpha * exp(a_pc);

            vec a_connect_rand = randu<colvec>(bMat_p.n_rows);
            a_pc.elem(find(a_connect_rand > a_connectance)).zeros();

            aMat = a_pc; // dimensional model
        } else {
            aMat.resize(aMat.n_rows, aMat.n_cols + 1);
            vec a_pc = randn<colvec>(aMat.n_rows);
            a_pc *= sigma;
            a_pc = alpha * exp(a_pc);

            vec a_connect_rand = randu<colvec>(bMat_p.n_rows);
            a_pc.elem(find(a_connect_rand > a_connectance)).zeros();

            aMat.col(aMat.n_cols - 1) = a_pc;
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Add row to emMat ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (emRate < 0) { // sample emigration rate for uniform distribution
            vec emRate_i(1);
            emRate_i.randu();
            if (emMat_c.n_rows == 0) {
                emMat_c = emRate_i;
            } else {
                emMat_c = join_vert(emMat_c, emRate_i);
            }
        }
    }
}

field<uvec> Species::extinct(int wholeDom, uvec ind_p, uvec ind_c) {

    // summary:
    // scan biomass matrices for regionally extinct species (root) and remove corresponding vectors from model matrices (all processes)

    //arguments:
    // wholeDom - indicates function called from root process and initiates scan of biomass matrices
    // ind_p - vector of indices of extinct producers - passed to non-root process for vector removal
    // ind_p - vector of indices of extinct consumers - passed to non-root process for vector removal

    // required members:
    // thresh - detection/extinction threshold

    // output:
    // (root) 2D field of producer/consumer extiction indices
    // updates to matrices bMat_p, rMat, cMat, bMat_c, aMat

    field<uvec> indReturn(2); // return object - 2D uvec of prod/cons extinction indices
    int countS=0;
    uvec ind_ext; // temporary extinction index
    int S_tot = bMat_p.n_rows - bMat_c.n_rows; // total species richness

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// Check for extinct producers ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (wholeDom) { // check for extinct species at regional scale
        ind_p.set_size(bMat_p.n_rows);
        for (int i=0; i<bMat_p.n_rows; i++) {
            ind_ext = find(bMat_p.row(i) > thresh); // search for present populations
            if (ind_ext.n_rows == 0) {
                ind_p(countS) = i; // if no populations present, add index to ind_p
                countS++;
            }
        }
        ind_p.resize(countS);
    }
    indReturn(0) = ind_p; // store extinct produce index

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Shed rows/cols model objects ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int i = ind_p.n_rows - 1; i >= 0; i--) {
        bMat_p.shed_row(ind_p(i));
        if (uMat_p.n_rows > 0) {
            uMat_p.shed_row(ind_p(i));
        }
        rMat.shed_row(ind_p(i));
        if (sMat.n_rows > 0) {
            sMat.shed_row(ind_p(i));
        }
        if (pProducer < 1){
            aMat.shed_row(ind_p(i));
        }
        if (prodComp) {
            cMat.shed_row(ind_p(i));
            cMat.shed_col(ind_p(i));
        }
        if (tMat.n_rows != 0) {
            tMat.shed_row(ind_p(i));
        }
        if (emMat_p.n_rows != 0) {
            emMat_p.shed_row(ind_p(i));
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Check for extinct consumers ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (wholeDom) { // as above for consumer species
        ind_c.set_size(bMat_c.n_rows);
        countS=0;
        for (int i=0; i<bMat_c.n_rows; i++) {
            ind_ext = find(bMat_c.row(i) > thresh);
            if (ind_ext.n_rows == 0) {
                ind_c(countS) = i;
                countS++;
            }
        }
        ind_c.resize(countS);
    }
    indReturn(1) = ind_c;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Shed rows/cols model objects ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int i = ind_c.n_rows - 1; i >= 0; i--) {
        bMat_c.shed_row(ind_c(i));
        if (uMat_c.n_rows > 0) {
            uMat_c.shed_row(ind_c(i));
        }
        aMat.shed_col(ind_c(i));
        if (emMat_p.n_rows != 0) {
            emMat_p.shed_row(ind_p(i));
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Update temporal species richness vector and S_p/c ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (sppRichness.n_rows == 0) {
        sppRichness.set_size(1,2);
        sppRichness(sppRichness.n_rows-1,0) = bMat_p.n_rows;
        sppRichness(sppRichness.n_rows-1,1) = bMat_c.n_rows;
    } else { // resize sppRichness vector(s) after each invasion
        sppRichness.resize(sppRichness.n_rows+1,2);
        sppRichness(sppRichness.n_rows-1,0) = bMat_p.n_rows;
        sppRichness(sppRichness.n_rows-1,1) = bMat_c.n_rows;
    }

    return(indReturn);
}

void Species::genDispMat() {

    // summary:
    // generate the regional dispersal operator
    // dMat(x,y) = (e / sum_y(exp(-d_xy/l)) * exp(-d_xy/l)
    // dMat(x,x) = -e
    // if dMat already exists, regenerate for long distance dispersal perturbation

    // required members:
    // dispL - dispersal length
    // emRate - emigration rate
    // topo.distMat
    // topo.adjMat

    // external function calls:
    // genAdjMat()

    // output:
    // dMat_n - regional dispersal operator

    if (topo.no_nodes == 1) {
        dMat_n.set_size(1,1);
        dMat_n.zeros();
    } else {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Generate exponential dispersal kernal /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (topo.distMat.n_rows == 0) {
            topo.genDistMat();
        }
        dMat_n = exp(-topo.distMat / dispL);
        if (topo.adjMat.n_rows == 0) {
            topo.genAdjMat();
        }
        dMat_n %= topo.adjMat;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Normalise immigration rate ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    mat kMat(topo.no_nodes, topo.no_nodes); // matrix encoding degree of each patch
    kMat.zeros();
    for (int i=0; i<topo.no_nodes; i++) { // normalise by weight
        // simple degree normalization
        uvec kVec = find(abs(topo.adjMat.col(i)) == 1); // absolute means that new edges (-1) also included
        for (int j=0; j<topo.no_nodes; j++) {
            if (abs(topo.adjMat(i,j)) == 1) {
                // simple degree normalization
                kMat(i,j) = emRate / kVec.n_rows;
            }
        }
    }

    if (emRate < 0) { // for non-uniform emRate, set to -1.0, emRate included in dynamics
        kMat = abs(kMat);
    }
    dMat_n(find(topo.adjMat == -1)).fill(1); // required for element-wise multiplication
    dMat_n = kMat % dMat_n; // new edges, if any, assigned e/k (distance independent)
    dMat_n.diag().fill(-1*abs(emRate)); // the dispersal operator includes the (negative) emigration terms on the diagonal
}

void Species::subSet(int domain) {

    // summary:
        // subset model spatially resolved model object by domain decomposition using restriction operators and column indices
        // generate the matrix of fixed unknowns representingan approximation of the subdomain interface

    // arguments:
        // domain - index of subdomain, i.e. rank of parallel process, from which function called

    // required members:
        // topo.fVec - indicator vector, elements of which denote subdomain allocation of nodes

    // output:
        // updates to matrices bMat_p, rMat, bMat_c, uMat_p, uMat_c, dMat_n, dMat_m
        // define topo.intIF and topo.adjIF which record the interface nodes

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Infer and store intra- and inter- subdomain edges /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uvec index; // storage indices of specified subdomain
    index = find(topo.fVec == domain);

    std::vector<long long unsigned int> internal, adjacent; // storage objects for internal and adjacent edges - push_back functionality
    for (auto it=index.begin(); it!=index.end(); it++) {
        for (int j=0; j<topo.no_nodes; j++) {
            if (topo.adjMat(*it, j) == 1) {
                uvec in = find(index == j);
                if (in.is_empty()) {
                    internal.push_back(*it);
                    adjacent.push_back(j);
                }
            }
        }
    }

    uvec iIF(&internal[0], internal.size()); // convert to armadillo objects
    uvec aIF(&adjacent[0], internal.size());
    topo.intIF = unique(iIF); // dereplicate
    topo.adjIF = unique(aIF);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Extract row/cols corresponding to specified subdomain /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    topo.network = topo.network.rows(index);
    if (topo.envMat.n_rows != 0) {
        topo.envMat = topo.envMat.cols(index);
    }
    if (bMat_p.n_rows != 0) { // imported objects
        uMat_p = bMat_p.cols(topo.adjIF);
        bMat_p = bMat_p.cols(index);
        rMat = rMat.cols(index);
    }
    if (bMat_c.n_rows != 0) { // imported objects
        uMat_c = bMat_c.cols(topo.adjIF);
        bMat_c = bMat_c.cols(index);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Generate restriction operators for decomposing dispersal operator /////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    umat P_n(index.n_rows, topo.no_nodes), P_m(topo.adjIF.n_rows, topo.no_nodes); // restriction operators
    P_n.zeros();
    P_m.zeros();

    for (int i=0; i<index.n_rows; i++) {
        P_n(i,index(i)) = 1;
    }

    for (int i=0; i<topo.adjIF.n_rows; i++) {
        P_m(i,topo.adjIF(i)) = 1;
    }

    dMat_m = P_m * dMat_n * P_n.t(); // inter-subdomain dispersal operator
    dMat_n = P_n * dMat_n * P_n.t(); // intra-subdomain dispersal operator

    if (uMat_p.n_rows > 0) {
        uMat_p = uMat_p * dMat_m;
    }
    if (uMat_c.n_rows > 0) {
        uMat_c = uMat_c * dMat_m;
    }

    if (topo.distMat.n_rows != 0) {
        topo.distMat = topo.distMat.submat(index,index);
    }

    topo.subdomain = domain;
    topo.subdom_nodes = index;
}