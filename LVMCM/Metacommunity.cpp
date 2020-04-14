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
 * This class contains the members and methods required for simulating metacommunity dynamics,
 * importing and outputing data
 */

#include <iostream>
#include <fstream>
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random.hpp>
#include <vector>

#include "Metacommunity.h"
#include "Topography.h"
#include "Species.h"
#include "ODE.h"
#include "CommunityDynamics.h"
#include "LVMCM_rng.h"

using namespace std;
using namespace arma;
using namespace boost::numeric::ublas;
using namespace boost::filesystem;

void Metacommunity::metaCDynamics(int relaxT, bool dispersal) {

    // summary:
    // numerically solve metacommunity dynamics for trajectory simulation (serial or parallel) or invader testing
    // trajectory approximation uses Sundials CVODE solver (see communityDynamics.h/.cpp)

    // arguments:
    // relaxT - relaxtion time (if != tMax)
    // dispersal - select with/without dispersal term (for invader testing)

    // required memebers:
    // sppPool - object of class Species, where model state is stored
    // tMax - relaxation time
    // storeTraj - selects whether trajectory object, dimensions (N*S)x(tMax), should be stored (large memory cost)
    // if parallel approximation used:
    // iterS - number Schwarz iterations
    // deltaT - size of Schwarz timewindow

    // external function calls:
    // ODE.h::integrate_until()

    // output:
    // updates to matrices bMat_p, and bMat_c

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Initialize community dynamics machinery ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CommunityDynamics dynamics;  // ODE dynamical object

    mat B_p, B_c, R, D; // storage for matrix subsets in case of parallel invader testing only

    if (!dispersal) {
        // invader testing - dispersal switched off for efficiency
        dynamics.bMat_p = &sppPool.bMat_p;
        dynamics.bMat_c = &sppPool.bMat_c;
        dynamics.rMat = &sppPool.rMat;
        if (sppPool.cMat.n_rows != 0) {
            dynamics.cMat = &sppPool.cMat;
        }
        if (sppPool.aMat.size() != 0) {
            dynamics.aMat = &sppPool.aMat;
        }
        dynamics.rho = &sppPool.rho;

    } else {
        // full dynamics - dispersal switched on
        dynamics.bMat_p = &sppPool.bMat_p;
        dynamics.bMat_c = &sppPool.bMat_c;
        dynamics.rMat = &sppPool.rMat;
        if (sppPool.cMat.n_rows != 0) {
            dynamics.cMat = &sppPool.cMat;
        }
        if (sppPool.emMat_p.n_rows != 0) {
            dynamics.emMat_p = &sppPool.emMat_p;
        }
        if (sppPool.emMat_c.n_rows != 0) {
            dynamics.emMat_c = &sppPool.emMat_c;
        }
        if (sppPool.aMat.size() != 0) {
            dynamics.aMat = &sppPool.aMat;
        }
        dynamics.dMat = &sppPool.dMat_n; // if dynamics.dMat is initialized, single domain relaxation selected
        dynamics.rho = &sppPool.rho;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Metacommunity relaxation step /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int res = 10; // set time step resolution for storing trajectories; relaxT must be multiple of res
    if (storeTraj == 0) { // standard relaxation
        {
            ODE_state state(&dynamics);
            state.integrate_until(relaxT);
        }

    } else if (storeTraj == 1) { // relax, storing trajectory object
        {
            // write time point to file in SxN biomass matrix
            // file directory reflects the parOut and rep of the current assembly

            std::size_t pos1 = bMatFileName.find_last_of("/");
            string bMatDir = bMatFileName.substr(0,pos1);
            pos1 = bMatFileName.find_last_of(")");
            std::size_t pos2 = bMatFileName.find_last_of(".");
            string newFolder = bMatFileName.substr(pos1+1, pos2-pos1-1);
            string Bpath = bMatDir + "/" + newFolder + "/trajectory";
            cout << "\nSaving matrices to " << Bpath << endl;

            if (!exists(Bpath)) { // make directory if doesn't currently exist
                create_directories(Bpath);
            }

            ODE_state state(&dynamics);
            for (int t = 0; t < relaxT/res; t++) {
                cout << t << endl;
                state.integrate_until(t*res);
                mat Btofile(&state[0], 1, dynamics.number_of_variables());
                Btofile.reshape(sppPool.S_p, sppPool.topo.no_nodes);
                string Bfile = Bpath + "/bMat" + to_string(t) + ".mat";
                Btofile.save(Bfile, raw_ascii);
            }

            // store parameter file in same directory for convenience
            string filenameP = Bpath + "/pars.mat";
            boost::filesystem::ofstream params;
            params.open(filenameP);
            params << "date " << date << endl;
            params << "jobID " << jobID << endl;
            params << "experiment " << experiment << endl;
            params <<  "invMax " << invMax << endl;
            params << "parOut " << parOut << endl;
            params << "rep " << rep << endl;
            params << "simTime " << simTime << endl;
            params << "no_nodes " << sppPool.topo.no_nodes << endl;
            params << "bisec " << sppPool.topo.bisec << endl;
            params << "phi " << sppPool.topo.phi << endl;
            params << "envVar " << sppPool.topo.envVar << endl;
            params << "c1 " << sppPool.c1 << endl;
            params << "c2 " << sppPool.c2 << endl;
            params << "emRate " << sppPool.emRate << endl;
            params << "dispL " << sppPool.dispL << endl;
            params << "invasion " << sppPool.invasion << endl;
            params << "iterS " << iterS << endl;
            params << "deltaT " << deltaT << endl;
            params << "tMax " << tMax << endl;
            params << "pProducer " << sppPool.pProducer << endl;
            params << "prodComp " << sppPool.prodComp << endl;
            params << "var_e " << sppPool.topo.var_e << endl;
            params << "alpha " << sppPool.alpha << endl;
            params << "sigma " << sppPool.sigma << endl;
            params << "rho " << sppPool.rho << endl;
            params << "discr_c_ij " << sppPool.discr_c_ij << endl;
            params << "T_int " << sppPool.topo.T_int << endl;
            params.close();
        }
    }
}

uvec Metacommunity::invaderSample(int trophLev, int no_invaders) {
    // summary: ...

    // testing parameters
    int no_time_steps = 1; // testing period (run twice) -- seems very short
    double min_b = 1e-7; // minimum biomass for inclusion in model -- set arbitrarily
    double inv = 1e-6; // invasion biomass
    uvec invaderIndex, posGrowth;

    if (trophLev == 0) { // sample producer species
        // size of testing pool
        int invExcess_p = 2; // invade excess species to account for difficulty in finding successful invader -- set arbitrarily
        if (invasionProb_p.n_rows > 0) { // set the invasion excess as equal to local mean invasion probability^-1
            invExcess_p = ceil(1 / invasionProb_p(invasionProb_p.n_rows - 1, 1));
        }

        for (int i = 0; i < invExcess_p * no_invaders; i++) {
            sppPool.invade(0); // invade a producer
        }

        invaderIndex = linspace<uvec>(sppPool.S_p, sppPool.bMat_p.n_rows - 1,
                                      sppPool.bMat_p.n_rows - sppPool.S_p);
        metaCDynamics(no_time_steps, 0); // simulate for no_time_steps without bounday approximation
        mat bInv = -1 * sppPool.bMat_p.rows(invaderIndex);
        metaCDynamics(no_time_steps, 0);
        bInv += sppPool.bMat_p.rows(invaderIndex);
        vec bInv_max(bInv.n_rows);
        for (int i = 0; i < bInv.n_rows; i++) { // store maximum local biomass of each invader
            bInv_max(i) = bInv.row(i).max();
        }

        posGrowth = find(bInv_max >= min_b);

    } else if (trophLev == 1) { // sample consumer species
        // size of testing pool
        int invExcess_c = 1; // invade excess species to account for difficulty in finding successful invader
        if (invasionProb_c.n_rows > 0) {
            invExcess_c = ceil(1 / invasionProb_c(invasionProb_c.n_rows - 1, 1));
        }

        for (int i = 0; i < invExcess_c * no_invaders; i++) {
            sppPool.invade(1); // invade a consumer
        }

        invaderIndex = linspace<uvec>(sppPool.S_c, sppPool.bMat_c.n_rows - 1,
                                      sppPool.bMat_c.n_rows - sppPool.S_c);
        metaCDynamics(no_time_steps+5, 0); // simulate for no_time_steps without bounday approximation
        mat bInv = -1 * sppPool.bMat_c.rows(invaderIndex);
        metaCDynamics(no_time_steps+5, 0);
        bInv += sppPool.bMat_c.rows(invaderIndex);
        vec bInv_max(bInv.n_rows);
        for (int i = 0; i < bInv.n_rows; i++) { // store maximum local biomass of each invader
            bInv_max(i) = bInv.row(i).max();
        }
        posGrowth = find(bInv_max >= min_b);
    }

    return(posGrowth); // return index of successful invaders
}

mat Metacommunity::invaderCleanup(int trophLev, uvec posGrowth) {

    mat bInv_max; // record max b_ix for subdomain selection

    if (trophLev == 0) {// clean up producers

        umat negGrowth(sppPool.bMat_p.n_rows - sppPool.S_p,1);
        negGrowth.col(0) = linspace<uvec>(sppPool.S_p, sppPool.bMat_p.n_rows - 1,
                                          sppPool.bMat_p.n_rows - sppPool.S_p);

        for (int i = posGrowth.n_rows - 1; i>=0; i--) {
            negGrowth.shed_row(posGrowth(i));
        }

        for (int i = negGrowth.n_rows - 1; i >= 0; i--) {
            sppPool.bMat_p.shed_row(negGrowth(i)); // remove prod biomass vec
            if (sppPool.pProducer < 1.0) {
                sppPool.aMat.shed_row(negGrowth(i)); // remove prod trophic link vec
            }
            if (sppPool.rMat.n_rows > 0) {
                sppPool.rMat.shed_row(negGrowth(i)); // remove prod growth vec
            }
            if (sppPool.sMat.n_rows > 0) {
                sppPool.sMat.shed_row(negGrowth(i)); // remove prod growth vec
            }
            if (sppPool.prodComp) {
                sppPool.cMat.shed_row(negGrowth(i)); // remove prod comp term
                sppPool.cMat.shed_col(negGrowth(i));
            }
            if (sppPool.topo.envVar != 0) {
                sppPool.tMat.shed_row(negGrowth(i));  // remove prod env tol vec
            }
        }

        if (posGrowth.n_rows > 0) {
            bInv_max.set_size(posGrowth.n_rows,1);

            for (int i = 0; i < posGrowth.n_rows; i++) { // store maximum local biomass of each invader
                bInv_max.row(i) = sppPool.bMat_p.row(sppPool.S_p + i).max();
            }
        }
    } else if (trophLev == 1) { // clean up consumers

        umat negGrowth(sppPool.bMat_c.n_rows - sppPool.S_c,1);
        negGrowth.col(0) = linspace<uvec>(sppPool.S_c, sppPool.bMat_c.n_rows - 1,
                                          sppPool.bMat_c.n_rows - sppPool.S_c);

        for (int i = posGrowth.n_rows - 1; i>=0; i--) {
            negGrowth.shed_row(posGrowth(i));
        }

        for (int i = negGrowth.n_rows - 1; i >= 0; i--) {
            sppPool.bMat_c.shed_row(negGrowth(i)); // remove cons biomass vec
            sppPool.aMat.shed_col(negGrowth(i)); // remove cons troph link vec
        }

        if (posGrowth.n_rows > 0) {
            bInv_max.set_size(posGrowth.n_rows,1);

            for (int i = 0; i < posGrowth.n_rows; i++) { // store maximum local biomass of each invader
                bInv_max.row(i) = sppPool.bMat_c.row(sppPool.S_c + i).max();
            }
        }
    }
    return(bInv_max);
}

void Metacommunity::invaderPopulate(int trophLev, mat bInv_max) {

    double inv = 1e-6;

    // check which subdomain species are growing fastest
    ucolvec subDom_max = index_max(bInv_max, 1);
    ucolvec bMax_index(1);

    if (trophLev == 0) { // populate producer species
        // reset b_i and invade into favoured node
        for (int i = 0; i < subDom_max.n_rows; i++) {
            if (subDom_max(i) == sppPool.topo.subdomain) {
                bMax_index = index_max(sppPool.bMat_p.row(sppPool.S_p + i));
                sppPool.bMat_p.row(sppPool.S_p + i).zeros();
                sppPool.bMat_p(sppPool.S_p + i, bMax_index(0)) = inv;
            } else {
                sppPool.bMat_p.row(sppPool.S_p + i).zeros();
            }
        }

        if (sppPool.prodComp) { // if producer dynamics coupled, fill cols after checking growth
            if (sppPool.discr_c_ij) { // sample from discrete distribution

                boost::random::binomial_distribution<int> distribution(1, sppPool.c2);
                for (int j = 1; j <= sppPool.bMat_p.n_rows - sppPool.S_p; j++) {
                    for (int i = 0; i < sppPool.cMat.n_cols - j; i++) {
                        sppPool.cMat(i, sppPool.cMat.n_cols - j) = sppPool.c1 * distribution(LVMCM_rng::boost_rng);
                    }
                }
                sppPool.cMat.diag().ones();

            } else { // sample from continuous (beta) distribution

                typedef boost::random::mt19937 RandomNumberGenerator;
                typedef boost::random::beta_distribution<> BetaDistribution;
                typedef boost::variate_generator<RandomNumberGenerator &, BetaDistribution> Generator;
                BetaDistribution distribution(sppPool.c1, sppPool.c2);
                Generator getRandomNumber(LVMCM_rng::boost_rng, distribution);
                for (int j = 1; j <= sppPool.bMat_p.n_rows - sppPool.S_p; j++) {
                    for (int i = 0; i < sppPool.cMat.n_cols - j; i++) {
                        sppPool.cMat(i, sppPool.cMat.n_cols - j) = getRandomNumber();
                    }
                }
                sppPool.cMat.diag().ones();
            }
        }

        if (sppPool.topo.bisec > 0) {
            if (sppPool.uMat_p.n_rows == 0) {
                sppPool.uMat_p.set_size(sppPool.bMat_p.n_rows, sppPool.bMat_p.n_cols);
                sppPool.uMat_p.zeros();
            } else {
                sppPool.uMat_p.resize(sppPool.bMat_p.n_rows, sppPool.bMat_p.n_cols);
                sppPool.uMat_p.rows(sppPool.S_p, sppPool.bMat_p.n_rows - 1).zeros();
            }
        }

    } else if (trophLev == 1) { // populate consumer species
        // reset b_i and invade into favoured node
        for (int i = 0; i < subDom_max.n_rows; i++) {
            if (subDom_max(i) == sppPool.topo.subdomain) {
                bMax_index = index_max(sppPool.bMat_c.row(sppPool.S_c + i));
                sppPool.bMat_c.row(sppPool.S_c + i).zeros();
                sppPool.bMat_c(sppPool.S_c + i, bMax_index(0)) = inv;
            } else {
                sppPool.bMat_c.row(sppPool.S_c + i).zeros();
            }
        }

        // How to test the consumer species invasion fitness at low abundance?

        if (sppPool.topo.bisec > 0) {
            if (sppPool.uMat_c.n_rows == 0) {
                sppPool.uMat_c.set_size(sppPool.bMat_c.n_rows, sppPool.bMat_c.n_cols);
                sppPool.uMat_c.zeros();
            } else {
                sppPool.uMat_c.resize(sppPool.bMat_c.n_rows, sppPool.bMat_c.n_cols);
                sppPool.uMat_c.rows(sppPool.S_c, sppPool.bMat_c.n_rows - 1).zeros();
            }
        }
    }
}

void Metacommunity::genJacobian() {

    // summary:
        // generate the numerical approximation of the Jacobian matrix for computing regional competitive overlap matrix

    // required members:
        // sppPool - object of class Species, where model state is stored

    // external function calls:
        // ODE.h::numerical_jacobian()

    // output:
        // jacobian - object for storing output of numerical_jacobian()

    CommunityDynamics dynamics;

    dynamics.bMat_p = &sppPool.bMat_p;
    dynamics.bMat_c = &sppPool.bMat_c;
    dynamics.rMat = &sppPool.rMat;
    if (sppPool.cMat.n_rows != 0) {
        dynamics.cMat = &sppPool.cMat;
    }
    if (sppPool.aMat.size() != 0) {
        dynamics.aMat = &sppPool.aMat;
    }
    dynamics.dMat = &sppPool.dMat_n; // if dynamics.dMat is initialized, single domain relaxation selected
    dynamics.rho = &sppPool.rho;

    {
        ODE_state state(&dynamics);
        for (double t = 0; t < 0; t++) { // accesses back-end CVode initialization machinery
            state.integrate_until(t);
        }
        jacobian.set_size(sppPool.topo.no_nodes * (sppPool.bMat_p.n_rows + sppPool.bMat_c.n_rows),
                          sppPool.topo.no_nodes * (sppPool.bMat_p.n_rows + sppPool.bMat_c.n_rows));
        dynamics.numerical_Jacobian(jacobian);
    }
}

void Metacommunity::genCMatReg(double h) {

    // summary:
        // generate a numerically approximated regional interaction matrix via a computation harvesting experiment
        // for detail of algorithm see O'Sullivan et al. (2019)

    // arguments:
        // h - harvesting rate

    // required members:
        // jacobian
        // sppPool - object of class Species, where model state is stored

    // external function calls:
        // genJacobian()

    // output:
        // cMat_reg - spatially unresolved approximation of the metacommunity interaction matrix

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Compute numerical Jacobian and invert /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    genJacobian();
    mat J_inv = inv(jacobian);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Initialize population indices and matrix objects /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int S_tot = sppPool.bMat_p.n_rows + sppPool.bMat_c.n_rows;
    uvec index = linspace<uvec>(0, sppPool.topo.no_nodes-1, sppPool.topo.no_nodes);
    index *= S_tot; // indexes all populations of focal species i=1; indices for species j generated by element-wise addition
    cMat_reg.set_size(S_tot,S_tot); //
    vec H_i; // storage for harvesting vector
    vec dB_jx; // storage for perturbed, vectorized state matrix

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Iteratively harvest each species i and compute dB_jx ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int i=0; i<S_tot; i++) {
        H_i.zeros(S_tot*sppPool.topo.no_nodes);
        if (i < sppPool.bMat_p.n_rows) { // harvest producer species
            H_i.elem(index + i) = h * sppPool.bMat_p.row(i);
        } else { // harvest consumer species
            H_i.elem(index + i) = h * sppPool.bMat_c.row(i - sppPool.bMat_p.n_rows);
        }

        // "local shift in biomasses of species j due to harvesting of species i per unit h"
        dB_jx = -1 * (J_inv * H_i) / h; // -1 due to definition of competitive interaction coefficent as (+)c_ij
        mat dB_jx_Mat(&dB_jx(0), S_tot, sppPool.topo.no_nodes); // state vector in matrix form (for row sum operation)
        cMat_reg.col(i) = sum(dB_jx_Mat, 1); // store dB_j = sum_x(dB_jx)
    }


    jacobian.reset();
    cMat_reg = inv(cMat_reg); // invert
    mat norm, sgn;
    norm.zeros(cMat_reg.n_rows, cMat_reg.n_cols);
    sgn.zeros(cMat_reg.n_rows, cMat_reg.n_cols);
    norm.diag() = 1/sqrt(abs(cMat_reg.diag())); // normalization
    sgn.diag() = cMat_reg.diag() / abs(cMat_reg.diag());
    cMat_reg = sgn * norm * cMat_reg * norm; // normalize regional competitive overlap matrix
}

void Metacommunity::genSourceSink(int tFullRelax) {

    // summary:
        // infers which populations are dependent upon immigration for local detectability by switching off dispersal and relaxing to equilibrium

    // arguments:
        // tFullRelax - (long) relaxation time

    // required members:
        // sppPool - object of class Species, where model state is stored

    // external function calls:
        // metaCDynamics()

    // output:
        // matrices of dimensions SxN with 1 indicating source, -1 sink populations

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Set off-diagonal elements of dMat_n to zero and relax ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    mat DStore = sppPool.dMat_n;
    mat BStore_p = sppPool.bMat_p;
    mat BStore_c = sppPool.bMat_c;
    mat Src_p, Src_c, Snk_p, Snk_c;
    sppPool.dMat_n.zeros();
    sppPool.dMat_n.diag().fill(-1*sppPool.emRate);
    metaCDynamics(tFullRelax,0); // relax

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Generate discrete source (1) sink (-1) matrices //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Src_p.zeros(sppPool.bMat_p.n_rows, sppPool.bMat_p.n_cols);
    Src_c.zeros(sppPool.bMat_c.n_rows, sppPool.bMat_c.n_cols);
    Src_p(find(sppPool.bMat_p > sppPool.thresh)).ones(); // 1 allocated to detected populations after dispersal switched off
    Src_c(find(sppPool.bMat_c > sppPool.thresh)).ones();

    sppPool.bMat_p = BStore_p; // reset metacommunity objects
    sppPool.bMat_c = BStore_c;
    sppPool.dMat_n = DStore;

    Snk_p.zeros(sppPool.bMat_p.n_rows,sppPool.bMat_p.n_cols);
    Snk_c.zeros(sppPool.bMat_c.n_rows,sppPool.bMat_c.n_cols);

    Snk_p(find(sppPool.bMat_p > sppPool.thresh)).ones(); // 1 allocated to detected populations before dispersal switched off
    Snk_c(find(sppPool.bMat_c > sppPool.thresh)).ones();
    Snk_p = Src_p - Snk_p; // -1 allocated to sink populations, 1 allocated to pops excluded by immigration
    Snk_c = Src_c - Snk_c;
    Snk_p(find(Snk_p == 1)).zeros(); // weak source populations removed
    Snk_c(find(Snk_c == 1)).zeros();

    sppPool.bMat_p_src = Src_p + Snk_p; // final matrix source (1) sink (-1)
    sppPool.bMat_c_src = Src_c + Snk_c;
}

void Metacommunity::writePars(int repPer) {
    string filenameBPer, filenameBCPer, filenamePPer, filenameDPer, filenameNPer, filenameNOPer;
    ostringstream nameBPer, nameBCPer, namePPer, nameDPer, nameNPer, nameNOPer;
    string directory;

    ostringstream p;
    p << "SimulationData/Assembly_experiment/N=" << sppPool.topo.no_nodes << "/" << date << "/";
    outputDirectory += p.str();

    if (!exists(outputDirectory)) { // make directory if doesn't currently exist
        create_directories(outputDirectory);
    }

    namePPer << outputDirectory << sppPool.topo.bisec << "_paramsP" << repPer << "_" << sppPool.bMat_p.n_rows
             << ".mat";
    filenamePPer = namePPer.str();

    cout << "Saving to " << filenamePPer << endl;

    // generate parameter file
    boost::filesystem::ofstream params;
    params.open(filenamePPer);
    params << "date " << date << endl;
    params << "experiment " << experiment << endl;
    params << "simTime " << simTime << endl;
    params << "no_nodes " << sppPool.topo.no_nodes << endl;
    params << "bisec " << sppPool.topo.bisec << endl;
    params << "phi " << sppPool.topo.phi << endl;
    params << "envVar " << sppPool.topo.envVar << endl;
    params << "c1 " << sppPool.c1 << endl;
    params << "c2 " << sppPool.c2 << endl;
    params << "emRate " << sppPool.emRate << endl;
    params << "dispL " << sppPool.dispL << endl;
    params << "invasion " << sppPool.invasion << endl;
    params << "iterS " << iterS << endl;
    params << "deltaT " << deltaT << endl;
    params << "tMax " << tMax << endl;
    params << "pProducer " << sppPool.pProducer << endl;
    params << "prodComp " << sppPool.prodComp << endl;
    params << "var_e " << sppPool.topo.var_e << endl;
    params << "alpha " << sppPool.alpha << endl;
    params << "sigma " << sppPool.sigma << endl;
    params << "rho " << sppPool.rho << endl;
    params << "discr_c_ij " << sppPool.discr_c_ij << endl;
    params.close();
}

void Metacommunity::printParams() {

    // summary:
        // print model parameterization to console

    printf("\nModel parameters:\n");
    printf("\niterS %d", iterS);
    printf("\ndeltaT %d", deltaT);
    printf("\ntMax %d", tMax);
    printf("\nparOut %f", parOut);
    cout << "\nexperiment " << experiment;
    printf("\nrep %d", rep);
    printf("\nsimTime %f", simTime);
    cout << "\njobID " << jobID;
    cout << "\ndate " << date;
    printf("\n\nc1 %f", sppPool.c1);
    printf("\nc2 %f", sppPool.c2 );
    printf("\nemRate %f", sppPool.emRate);
    printf("\ndispL %f", sppPool.dispL);
    printf("\npProducer %f", sppPool.pProducer) ;
    printf("\nalpha %e", sppPool.alpha);
    printf("\nsigma %f", sppPool.sigma);
    printf("\nrho %f", sppPool.rho);
    printf("\ndiscr_c_ij %d", sppPool.discr_c_ij);
    printf("\n\nno_nodes %d", sppPool.topo.no_nodes);
    printf("\nphi %f", sppPool.topo.phi);
    printf("\nenvVar %d", sppPool.topo.envVar);
    printf("\nT_int %f", sppPool.topo.T_int);
    printf("\nvar_e %f", sppPool.topo.var_e);
    printf("\nrandGraph %d", sppPool.topo.randGraph);
    printf("\ngabriel %d", sppPool.topo.gabriel);
    printf("\nbisec %d", sppPool.topo.bisec);
    printf("\n\ninvMax %d", invMax);
}

void Metacommunity::outputData() {

    // summary:
        // generate file names and save model matrices to file

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Generate strings for file names //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    cout << "\nOutputting data... ";

    string filenameA, filenameR, filenameSR,filenameB, filenameBC, filenameC, filenameF, filenameP, filenameN,
            filenameT, filenameTr, filenameEf, filenameS, filenameD,
            filenameI, filenameE, filenameIP, filenameIPp, filenameIPc, filenameBs,filenameBCs;

    if (bMatFileName.length() == 0) { // generate path from scratch
        ostringstream nameA, nameR, nameSR, nameB, nameBC, nameC, nameF,  nameP, nameN, nameT, nameEf,
                nameTr, nameS, nameD, nameI, nameE, nameIP, nameIPp, nameIPc, nameBs, nameBCs;

        ostringstream p;
        p << "SimulationData/N=" << sppPool.topo.no_nodes << "/" << experiment
          << "_experiment/" << date << "/";
        outputDirectory += p.str();

        if (!exists(outputDirectory)) { // make directory if doesn't currently exist
            create_directories(outputDirectory);
        }

        nameR << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "rMat" << rep << ".mat";
        filenameR = nameR.str();
        nameSR << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "sMat" << rep << ".mat";
        filenameSR = nameSR.str();
        nameA << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "aMat" << rep << ".mat";
        filenameA = nameA.str();
        nameB << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "bMat" << rep << ".mat";
        filenameB = nameB.str();
        nameBC << outputDirectory << date << "_" << experiment << "("
               << invMax << ")" << parOut << "bMat_c" << rep << ".mat";
        filenameBC = nameBC.str();
        nameP << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "params" << rep << ".mat";
        filenameP = nameP.str();
        nameN << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "network" << rep << ".mat";
        filenameN = nameN.str();
        nameTr << outputDirectory << date << "_" << experiment << "("
               << invMax << ")" << parOut << "trajec" << rep << ".mat";
        filenameTr = nameTr.str();
        nameEf << outputDirectory << date << "_" << experiment << "("
               << invMax << ")" << parOut << "envFluct" << rep << ".mat";
        filenameEf = nameEf.str();
        nameS << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "S" << rep << ".mat";
        filenameS = nameS.str();
        nameD << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "dMat_n" << rep << ".mat";
        filenameD = nameD.str();
        nameI << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "cMat" << rep << ".mat";
        filenameI = nameI.str();
        nameC << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "cMat_reg" << rep << ".mat";
        filenameC = nameC.str();
        nameE << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "environ" << rep << ".mat";
        filenameE = nameE.str();
        nameT << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "tMat" << rep << ".mat";
        filenameT = nameT.str();
        nameF << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "fVec" << rep << ".mat";
        filenameF = nameF.str();
        nameIP << outputDirectory << date << "_" << experiment << "("
                << invMax << ")" << parOut << "invProb" << rep << ".mat";
        filenameIP = nameIP.str();
        nameIPp << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "invProb_p" << rep << ".mat";
        filenameIPp = nameIPp.str();
        nameIPc << outputDirectory << date << "_" << experiment << "("
                << invMax << ")" << parOut << "invProb_c" << rep << ".mat";
        filenameIPc = nameIPc.str();
        nameBs << outputDirectory << date << "_" << experiment << "("
                << invMax << ")" << parOut << "bMat_src" << rep << ".mat";
        filenameBs = nameBs.str();
        nameBCs << outputDirectory << date << "_" << experiment << "("
                << invMax << ")" << parOut << "bMat_c_src" << rep << ".mat";
        filenameBCs = nameBCs.str();
        bMatFileName = filenameB;
        cout << "Biomass matrix file name: ";
        cout << bMatFileName << endl;
    } else { // use find and replace to generate paths from bMatFileName

        size_t b = bMatFileName.find("bMat");
        filenameBC = bMatFileName;
        filenameBC.replace(b, string("bMat").length(), "bMat_c");
        filenameD = bMatFileName;
        filenameD.replace(b, string("bMat").length(), "dMat_n");
        filenameI = bMatFileName;
        filenameI.replace(b, string("bMat").length(), "cMat");
        filenameC = bMatFileName;
        filenameC.replace(b, string("bMat").length(), "cMat_reg");
        filenameA = bMatFileName;
        filenameA.replace(b, string("bMat").length(), "aMat");
        filenameN = bMatFileName;
        filenameN.replace(b, string("bMat").length(), "network");
        filenameP = bMatFileName;
        filenameP.replace(b, string("bMat").length(), "params");
        filenameR = bMatFileName;
        filenameR.replace(b, string("bMat").length(), "rMat");
        filenameSR = bMatFileName;
        filenameSR.replace(b, string("bMat").length(), "sMat");
        filenameS = bMatFileName;
        filenameS.replace(b, string("bMat").length(), "S");
        filenameT = bMatFileName;
        filenameT.replace(b, string("bMat").length(), "tMat");
        filenameE = bMatFileName;
        filenameE.replace(b, string("bMat").length(), "environ");
        filenameTr = bMatFileName;
        filenameTr.replace(b, string("bMat").length(), "trajec");
        filenameEf = bMatFileName;
        filenameEf.replace(b, string("bMat").length(), "envFluct");
        filenameF = bMatFileName;
        filenameF.replace(b, string("bMat").length(), "fVec");
        filenameIP = bMatFileName;
        filenameIP.replace(b, string("bMat").length(), "invProb");
        filenameIPp = bMatFileName;
        filenameIPp.replace(b, string("bMat").length(), "invProb_p");
        filenameIPc = bMatFileName;
        filenameIPc.replace(b, string("bMat").length(), "invProb_c");
        filenameBs = bMatFileName;
        filenameBs.replace(b, string("bMat").length(), "bMat_src");
        filenameBCs = bMatFileName;
        filenameBCs.replace(b, string("bMat").length(), "bMat_c_src");
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Create/edit file recording current state ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (jobID != "NA") {
        cout << "\nOutputting current state... ";
        ostringstream st;
        st << jobID << ".txt"; // name of file in which state recorded
        string state = st.str();
        std::ifstream stateFile;
        stateFile.open(state);

        if (!stateFile.is_open()) {
            cout << "generating state file " << state << endl;
            std::ofstream stateFileNew;
            stateFileNew.open(state);
            stateFileNew << sppPool.invasion << endl;
            stateFileNew << invMax << endl;
            stateFileNew << bMatFileName << endl;
            stateFileNew << "1" << endl;
            stateFileNew.close();
        } else {
            cout << "updating state file " << state << endl;
            std::vector<string> state_vec;
            string state_line;
            while ( stateFile.good() ) {
                getline (stateFile,state_line);
                state_vec.push_back(state_line);
            }
            stateFile.close();
            state_vec[0] = to_string(sppPool.invasion);
            std::ofstream stateFileNew;
            stateFileNew.open(state);
            if (stateFileNew.is_open()) {
                for (int i =0; i<state_vec.size()-1; i++) {
                    stateFileNew << state_vec[i] << endl;
                }
                stateFileNew.close();
            }
            stateFileNew.close();
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Store matrix objects ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Note, if specific experiment run, only relevant objects are output
    if (sppPool.trajectories.n_rows != 0) { // output trajectories ONLY
        sppPool.trajectories.save(filenameTr, raw_ascii);
        if (sppPool.efMat.n_rows != 0) {
            mat R = sppPool.rMat;
            R.reshape(1, R.size());
            rowvec N(R.n_cols);
            N.fill(datum::nan);
            R = join_vert(R, N);
            sppPool.fluctuations = join_vert(R, sppPool.fluctuations);
            sppPool.fluctuations.save(filenameEf, raw_ascii);
        }
    } else if (cMat_reg.n_rows != 0) { // output regional competitive overlap matrix ONLY
        cMat_reg.save(filenameC, raw_ascii);
    } else if (sppPool.bMat_p_src.n_rows != 0) { // output source-sink matrices ONLY
        sppPool.bMat_p_src.save(filenameBs, raw_ascii);
        if (sppPool.bMat_c_src.n_rows != 0) {
            sppPool.bMat_c_src.save(filenameBCs, raw_ascii);
        }
    } else { // output standard model objects
        sppPool.topo.network.save(filenameN, raw_ascii);
        sppPool.topo.fVec.save(filenameF, raw_ascii);
        sppPool.rMat.save(filenameR, raw_ascii);
        if (sppPool.sMat.n_rows != 0) {
            sppPool.sMat.save(filenameSR, raw_ascii);
        }
        if (bMatFileName.length() == 0) {
            sppPool.bMat_p.save(filenameB, raw_ascii);
        } else {
            sppPool.bMat_p.save(bMatFileName, raw_ascii);
        }
        if (sppPool.bMat_c.n_rows != 0) {
            sppPool.bMat_c.save(filenameBC, raw_ascii);
        }
        if (sppPool.topo.envMat.n_rows != 0) {
            sppPool.topo.envMat.save(filenameE, raw_ascii);
            sppPool.tMat.save(filenameT, raw_ascii);
        }
        sppPool.sppRichness.save(filenameS, raw_ascii);
        sppPool.dMat_n.save(filenameD, raw_ascii);
        if (sppPool.cMat.n_rows != 0) {
            sppPool.cMat.save(filenameI, raw_ascii);
        }
        if (sppPool.aMat.n_rows != 0) {
            sppPool.aMat.save(filenameA, raw_ascii);
        }
        if (invasionProb.n_rows != 0) {
            invasionProb.save(filenameIP, raw_ascii);
        }
        if (invasionProb_p.n_rows != 0) {
            invasionProb_p.save(filenameIPp, raw_ascii);
        }
        if (invasionProb_c.n_rows != 0) {
            invasionProb_c.save(filenameIPc, raw_ascii);
        }

        boost::filesystem::ofstream params;
        params.open(filenameP);
        params << "date " << date << endl;
        params << "jobID " << jobID << endl;
        params << "experiment " << experiment << endl;
        params <<  "invMax " << invMax << endl;
        params << "parOut " << parOut << endl;
        params << "rep " << rep << endl;
        params << "simTime " << simTime << endl;
        params << "no_nodes " << sppPool.topo.no_nodes << endl;
        params << "bisec " << sppPool.topo.bisec << endl;
        params << "phi " << sppPool.topo.phi << endl;
        params << "envVar " << sppPool.topo.envVar << endl;
        params << "c1 " << sppPool.c1 << endl;
        params << "c2 " << sppPool.c2 << endl;
        params << "emRate " << sppPool.emRate << endl;
        params << "dispL " << sppPool.dispL << endl;
        params << "invasion " << sppPool.invasion << endl;
        params << "iterS " << iterS << endl;
        params << "deltaT " << deltaT << endl;
        params << "tMax " << tMax << endl;
        params << "pProducer " << sppPool.pProducer << endl;
        params << "prodComp " << sppPool.prodComp << endl;
        params << "var_e " << sppPool.topo.var_e << endl;
        params << "alpha " << sppPool.alpha << endl;
        params << "sigma " << sppPool.sigma << endl;
        params << "rho " << sppPool.rho << endl;
        params << "discr_c_ij " << sppPool.discr_c_ij << endl;
        params << "T_int " << sppPool.topo.T_int << endl;
        params << "t_niche " << sppPool.t_niche << endl;
        params.close();
    }
}

bool fexists(const std::string& filename) { // check if file/directory exists
    std::ifstream ifile(filename.c_str());
    return (bool)ifile;
}

void Metacommunity::importData(string bMatFile) {

    // summary:
        // import metacommunity model

    // arguments:
        // bMatFile - path/file name of produce biomass matrix from which the rest of the file names are generated

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Generate strings for file names //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bMatFileName = bMatFile;

    cout << "Importing data, file name " << bMatFileName << endl;

    string filenameA, filenameR, filenameSR, filenameBC, filenameC, filenameF, filenameP, filenameN,
            filenameT, filenameTr, filenameS, filenameD,
            filenameI, filenameE, filenameIPp, filenameIPc;
    size_t b = bMatFileName.find("bMat");
    filenameBC = bMatFileName;
    filenameBC.replace(b, string("bMat").length(), "bMat_c");
    filenameD = bMatFileName;
    filenameD.replace(b, string("bMat").length(), "dMat_n");
    filenameI = bMatFileName;
    filenameI.replace(b, string("bMat").length(), "cMat");
    filenameA = bMatFileName;
    filenameA.replace(b, string("bMat").length(), "aMat");
    filenameN = bMatFileName;
    filenameN.replace(b, string("bMat").length(), "network");
    filenameP = bMatFileName;
    filenameP.replace(b, string("bMat").length(), "params");
    filenameR = bMatFileName;
    filenameR.replace(b, string("bMat").length(), "rMat");
    filenameSR = bMatFileName;
    filenameSR.replace(b, string("bMat").length(), "sMat");
    filenameS = bMatFileName;
    filenameS.replace(b, string("bMat").length(), "S");
    filenameT = bMatFileName;
    filenameT.replace(b, string("bMat").length(), "tMat");
    filenameE = bMatFileName;
    filenameE.replace(b, string("bMat").length(), "environ");
    filenameTr = bMatFileName;
    filenameTr.replace(b, string("bMat").length(), "trajec");
    filenameIPp = bMatFileName;
    filenameIPp.replace(b, string("bMat").length(), "invProb_p");
    filenameIPc = bMatFileName;
    filenameIPc.replace(b, string("bMat").length(), "invProb_c");
    filenameF = bMatFileName;
    filenameF.replace(b, string("bMat").length(), "fVec");

    cout << "Biomass matrix file name:\n" << bMatFileName << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Extract data stored in param file ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int number_of_params = 0;
    std::string paramLine;
    boost::filesystem::ifstream countParams;
    countParams.open(filenameP);

    while (std::getline(countParams, paramLine))
        ++number_of_params;

    boost::filesystem::ifstream param(filenameP);

    string parName[number_of_params];
    string pars[number_of_params];
    if (param.is_open())
    {
        for (int i=0; i<number_of_params; i++) {
            param >> parName[i] >> pars[i];
        }
    }

    for(int i=0; i<number_of_params; i++) {
        if (parName[i] == "date") {
            date = pars[i];
            continue;
        } else if (parName[i] == "jobID") {
            jobID = pars[i];
            continue;
        } else if (parName[i] == "experiment") {
            experiment = pars[i];
            continue;
        } else if (parName[i] == "invMax") {
            invMax = stoi(pars[i]);
            continue;
        } else if (parName[i] == "parOut") {
            parOut = stof(pars[i]);
            continue;
        } else if (parName[i] == "rep") {
            rep = stoi(pars[i]);
            continue;
        } else if (parName[i] == "simTime") {
            simTime += stof(pars[i]);
            continue;
        } else if (parName[i] == "no_nodes"){
            sppPool.topo.no_nodes = stoi(pars[i]);
            continue;
        } else if (parName[i] == "bisec") {
            sppPool.topo.bisec = stoi(pars[i]);
            continue;
        } else if (parName[i] == "phi") {
            sppPool.topo.phi = stof(pars[i]);
            continue;
        } else if (parName[i] == "envVar") {
            sppPool.topo.envVar = stoi(pars[i]);
            continue;
        } else if (parName[i] == "c1") {
            sppPool.c1 = stof(pars[i]);
            continue;
        } else if (parName[i] == "c2") {
            sppPool.c2 = stof(pars[i]);
            continue;
        } else if (parName[i] == "emRate") {
            sppPool.emRate = stof(pars[i]);
            continue;
        } else if (parName[i] == "dispL") {
            sppPool.dispL = stof(pars[i]);
            continue;
        } else if (parName[i] == "invasion") {
            sppPool.invasion = stoi(pars[i]);
            sppPool.invEvent = stoi(pars[i]);
            sppPool.inv0 = sppPool.invasion;
            continue;
        } else if (parName[i] == "iterS") {
            iterS = stoi(pars[i]);
            continue;
        } else if (parName[i] == "deltaT") {
            deltaT = stoi(pars[i]);
            continue;
        } else if (parName[i] == "tMax") {
            tMax = stoi(pars[i]);
            continue;
        } else if (parName[i] == "pProducer") {
            sppPool.pProducer = stof(pars[i]);
            continue;
        } else if (parName[i] == "var_e") {
            sppPool.topo.var_e = stof(pars[i]);
            continue;
        } else if (parName[i] == "alpha") {
            sppPool.alpha = stof(pars[i]);
            continue;
        } else if (parName[i] == "sigma") {
            sppPool.sigma = stof(pars[i]);
            continue;
        } else if (parName[i] == "rho") {
            sppPool.rho = stof(pars[i]);
            continue;
        } else if (parName[i] == "discr_c_ij") {
            sppPool.discr_c_ij = stoi(pars[i]);
            continue;
        } else if (parName[i] == "T_int") {
            sppPool.topo.T_int = stof(pars[i]);
            continue;
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Load matrix objects ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sppPool.bMat_p.load(bMatFileName);
    sppPool.S_p = sppPool.bMat_p.n_rows;
    if(fexists(filenameBC)) {
        sppPool.bMat_c.load(filenameBC);
        sppPool.S_c = sppPool.bMat_c.n_rows;
    }
    sppPool.dMat_n.load(filenameD);
    if(fexists(filenameI)) {
        sppPool.cMat.load(filenameI);
    }
    if(fexists(filenameA)) {
        sppPool.aMat.load(filenameA);
    }
    sppPool.topo.network.load(filenameN);
    if (fexists(filenameT)) {
        sppPool.tMat.load(filenameT);
    }
    if (fexists(filenameE)) {
        sppPool.topo.envMat.load(filenameE);
    }
    if (fexists(filenameF)) {
        sppPool.topo.fVec.load(filenameF);
    }
    sppPool.rMat.load(filenameR);
    if (fexists(filenameSR)) {
        sppPool.sMat.load(filenameSR);
    }
    sppPool.sppRichness.load(filenameS);
    if (fexists(filenameIPp)) {
        invasionProb_p.load(filenameIPp);
    }
    if (fexists(filenameIPc)) {
        invasionProb_c.load(filenameIPc);
    }

    cout << "\nImported network rows 0-4 = " << endl;
    cout << sppPool.topo.network.rows(0,min((int) sppPool.topo.network.n_rows-1, 4));
    sppPool.topo.genDistMat();
    printf("\nSpecies richness, S_p = %d, S_c = %d\n", (int) sppPool.bMat_p.n_rows, (int) sppPool.bMat_c.n_rows);
}
