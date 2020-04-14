////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// The parallelizable Lotka-Volterra Metacommunity assembly Model (pLVMCM) ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Jacob Dinner O'Sullivan -- j.osullivan@qmul.ac.uk | j.osullivan@zoho.com ///////////////////////
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

#include <iostream>
#include <stdio.h>
#include <mpi.h>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <csignal>
#include <csetjmp>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "Metacommunity.h"
#include "LVMCM_rng.h"

using namespace std;
using namespace arma;
using std::vector;

// Assembly flags
bool ASSEMBLE = true; // run assembly algorithm or jump directly to generation of analysis object
bool OUTPUT = true; // select write to file
bool WARMING = false; // select warming experiment
bool LONGDISTDISP = false; // select introduction of long distance dispersal
bool TRAJECTORY = false; // select generate and write trajectory object to file
bool FLUCTUATE = false; // select generate and write to file abiotic fluctuation
bool CMAT_REG = false; // select generate and write to file regional competitive overlap matrix
bool FIX_SEED = false; // set random seeds (all processes) to 1 (true) or generate randomly (false)

// Storage for setjmp/longjmp
jmp_buf jump_buffer;

// Global Metacommunity object and timing variables
Metacommunity meta;
time_t time1, time2;
double mpi_time1, mpi_time2;
int g_seed;

// CVode tolerances
double TolA = 1e-8;
double TolR = 1e-7;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Start of assembly algoritm ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    time_t start;
    time(&start);
    unsigned int sim = 1; // simulation switch

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Store program arguments in variables ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Default Metacommunity parameters
    unsigned int a_init = 1;
    string a_bMat = {};
    unsigned int a_invMax = 0;
    unsigned int a_iterS = 2;
    unsigned int a_deltaT = 100;
    unsigned int a_tMax = 500;
    unsigned int a_perfRep = 0;
    unsigned int a_perfS = 0;
    bool a_autoFluct = false;
    bool a_harvest = false;
    unsigned int a_locRandomize = 0;
    string a_outputDirectory;
    // Default Species parameters
    double a_c1 = 0.5;
    double a_c2 = 0.5;
    double a_emRate = 0.01;
    double a_dispL = 0.1;
    double a_pProducer = 1.0;
    bool a_prodComp = true;
    double a_alpha = 1e-7;
    double a_sigma = 7;
    double a_sigma_t = 0.1;
    double a_rho = 0.1;
    bool a_discr_c_ij = true;
    // Default Topograpy parameters
    unsigned int a_no_nodes = 4;
    double a_phi = 1;
    unsigned int a_envVar = 0;
    double a_var_e = 0.01;
    bool a_randGraph = true;
    bool a_gabriel = true;
    unsigned int a_bisec = 0;
    double a_T_int = -1;
    double a_dTdt = 0.0;
    unsigned int a_edges = 1;
    // Default output variables
    double a_parOut = 0;
    string a_experiment = "DEFAULT";
    unsigned int a_rep = 0;
    string a_jobID = "NA";

    for (int i = 1; i<argc; i++) { // loop through program arguments an allocate to parameters
        char var = argv[i][1];
        i++;
        switch (var) {
            case 'n' : // -new: initialize new model (T/F)
                if (!strcmp(argv[i],"F")) {
                    a_init = 0;
                }
                break;
            case 'b' : // -bFile: path to data for importing initialized model
                a_bMat = argv[i];
                break;
            case 'i' : // -invMax: total number of invasions
                a_invMax = atoi(argv[i]);
                break;
            case 't' : // -tMax: relaxation time
                a_tMax = atoi(argv[i]);
                break;
            case 'c' : // -cij: 2 competition parameters
                a_c1 = atof(argv[i]);
                i++;
                a_c2 = atof(argv[i]);
                break;
            case 'd' : // -disp: emigration rate and dispersal length
                a_emRate = atof(argv[i]);
                i++;
                a_dispL = atof(argv[i]);
                break;
            case 'p' : // -pProd: probabilty of sampling a producer for bipart models
                a_pProducer = atof(argv[i]);
                break;
            case 'C' : // -C: select producer coupling
                if (!strcmp(argv[i],"F")) {
                    a_prodComp = false;
                }
                break;
            case 'a' : // -alpha: base attack rate
                a_alpha = atof(argv[i]);
                break;
            case 's' : // -sigma: trophic link distribution parameter
                a_sigma = atof(argv[i]);
                break;
            case 'F' : // -sigma: trophic link distribution parameter
                a_sigma_t = atof(argv[i]);
                break;
            case 'r' : // -rho: consumer respiration rate
                a_rho = atof(argv[i]);
                break;
            case 'D' : // -Discr: select discrete distribution in competive overlap coefficients
                if (!strcmp(argv[i],"F")) {
                    a_discr_c_ij = false;
                }
                break;
            case 'N' : // -N: number of nodes
                a_no_nodes = atoi(argv[i]);
                break;
            case 'P' : // -P: spatial autocorrelation of the environment
                a_phi = atof(argv[i]);
                break;
            case 'e' : // -envVar: number of explicitly modelled environmental parameters
                a_envVar = atoi(argv[i]);
                break;
            case 'v' : // var_e: variance of environmental/growth rate distribution
                a_var_e = atof(argv[i]);
                break;
            case 'G' : // -Gab: select Gabriel/complete graph
                if (!strcmp(argv[i],"F")) {
                    a_gabriel = false;
                }
                break;
            case 'o' : // -o: output variables - key parameter, experiment name, replicate number
                a_parOut = atof(argv[i]);
                i++;
                a_experiment = argv[i];
                i++;
                a_rep = atoi(argv[i]);
                break;
            case 'O' : // -O: select write to file
                if (!strcmp(argv[i],"F")) {
                    OUTPUT = false;
                }
                break;
            case 'h' : // -harvest: estimate regional scale interaction coefficients using harvesting experiment
                if (!strcmp(argv[i],"F")) {
                    a_harvest = true;
                }
                break;
            case 'f' : // -folder: set output directory
                a_outputDirectory = argv[i];
                break;
            case 'z' : // fix random seed
                if (!strcmp(argv[i],"F")) {
                    FIX_SEED = true;
                }
                break;
            case 'L' : // select lattice
                if (!strcmp(argv[i],"T")) {
                    a_randGraph = false;
                }
                break;
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Initialise MPI communicator ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int rank;
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Set return clause parallel algorithm ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Metacommunity metaSub;

    if (setjmp(jump_buffer)) { // set return call in longjump clause
        for (int s=0; s<size; s++) {
            if (rank == s) {
                if (s == 0) {
                    printf("\n\nMPI:\n");
                    cout << "Process " << rank << " exiting" << endl;
                } else {
                    cout << "Process " << rank << " exiting" << endl;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            time_t finish;
            time(&finish);
            cout << endl << "Simulation stated at: " << ctime(&start);
            cout << "Simulation finished at: " << ctime(&finish);
            printf("\nFinalizing MPI communicator\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return (0);
    }

    if (rank == 0) { // start timer and check for parameter errors
        mpi_time1 = MPI_Wtime();

        if (!a_randGraph) { // check perfect square in case of lattice
            double srN = sqrt(a_no_nodes);
            if (floor(srN) - srN != 0) {
                cout << "\nError: perfect square N expected" << endl;
                longjmp(jump_buffer, 1);
            }
        }
        if (a_pProducer < 1) { // check non-zero trophic parameters
            if ((a_alpha == 0) || (a_sigma == 0) || (a_rho == 0)) {
                cout << "\nError: trophic parameters set to zero!" << endl;
                longjmp(jump_buffer, 1);
            }
        }
        if (a_prodComp == 1) { // check non-zero competitive parameters
            if ((a_c1 == 0) || (a_c2 == 0)) {
                cout << "\nError: competitive parameters set to zero!" << endl;
                longjmp(jump_buffer, 1);
            }
        }
        double no_sub = a_no_nodes / pow(2, a_bisec);
        if (floor(no_sub) - no_sub != 0) { // check bisec/no_nodes correspond
            cout << "\nError: N / 2^bisec non-integer" << endl;
            longjmp(jump_buffer, 1);
        }
        if (a_outputDirectory.length() == 0) { // check output directory set
            cout << a_outputDirectory << endl;
            cout << "\nWarning: no output directory given" << endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Set random seeds  //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (rank == 0) {
        if (FIX_SEED) {
            g_seed = 1;
        } else {
            // generate random seed
            std::random_device rd;
            g_seed = rd();
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (int s = 0; s < size; s++) {
        if (rank == s) {
            cout << "Random seed seen by process " << rank << " " << g_seed << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    LVMCM_rng::boost_rng.seed(g_seed);
    arma_rng::set_seed(g_seed);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Parameterize parallel assembly model ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (rank == 0) { // parameterize root process using program arguments
        printf("\nInitializing LVMCM\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    { // subdomain initialization scope
        // parameterize subdomains
        if (a_init) { // parameterize all processes from program arguments

            metaSub = Metacommunity(
                    // Metacommunity parameters
                    a_init,
                    a_bMat,
                    a_invMax,
                    a_iterS,
                    a_deltaT,
                    a_tMax,
                    a_outputDirectory,
                    // Species parameters
                    a_c1,
                    a_c2,
                    a_emRate,
                    a_dispL,
                    a_pProducer,
                    a_prodComp,
                    a_alpha,
                    a_sigma,
                    a_sigma_t,
                    a_rho,
                    a_discr_c_ij,
                    // Topograpy parameters
                    a_no_nodes,
                    a_phi,
                    a_envVar,
                    a_var_e,
                    a_randGraph,
                    a_gabriel,
                    a_bisec,
                    a_T_int,
                    // output variables
                    a_parOut,
                    a_experiment,
                    a_rep,
                    a_jobID);
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0) {
                cout << "output folder set to " << meta.outputDirectory << endl;
                metaSub.printParams();
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    } // subdomain initialization scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Generate topo - topography/environment (root) ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    { // domain decomposition scope
        if (rank == 0) {

            metaSub.sppPool.topo.genDomainDecomp(); // generate and decompose landscape

            if (meta.sppPool.topo.envVar > 0) {
                meta.sppPool.topo.genEnvironment(); // generate environmental distribution
            }

            metaSub.sppPool.topo.genDistMat();
            metaSub.sppPool.topo.genAdjMat();
            metaSub.sppPool.genDispMat();
            unique(metaSub.sppPool.topo.fVec.t()).print("\nf.unique");
            if (ASSEMBLE) {
                meta = metaSub; // store copy of complete domain at root process for outputting and extinction testing
            }
        }

        meta.rank = rank; // store rank for signal handler
        MPI_Barrier(MPI_COMM_WORLD);

    } // domain decomposition scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Simulation /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("\nStarting metacommunity assembly\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (int s=0; s<size; s++) {
        if (rank == s) {
            if (s == 0) {
              printf("\nMPI:\n");
                cout << "Process " << rank << " active" << endl;
            } else {
                cout << "Process " << rank << " active" << endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (rank == 0) {
        cout << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Initialize MPI buffers //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (ASSEMBLE) {
        { // assembly scope
            // initialise objects
            mat B_p, B_c;
            unsigned int no_invaders_p, no_residents_p, no_extinct_p, no_invaders_c, no_residents_c, no_extinct_c; // counters

            // reset seed for invader sampling (required to ensure all processes synced for invader sampling)
            LVMCM_rng::boost_rng.seed(g_seed);
            arma_rng::set_seed(g_seed);

            do {

                metaSub.sppPool.invEvent++;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Generate invader (all subdomains in parallel) ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                { // invader testing scope

                    // select number of invaders and sample their trophic levels
                    B_p = metaSub.sppPool.bMat_p; // store current state
                    B_c = metaSub.sppPool.bMat_c; // store current state
                    metaSub.sppPool.S_p = metaSub.sppPool.bMat_p.n_rows; // store current diversity
                    metaSub.sppPool.S_c = metaSub.sppPool.bMat_c.n_rows; // store current diversity

                    // select number of invaders and sample their trophic levels
                    int no_trophLev[2] = {0}; // no producers, consumers to invade
                    int no_invaders = 0.05 * (metaSub.sppPool.S_p + metaSub.sppPool.S_c) + 1;

                    if (metaSub.sppPool.pProducer < 1) { // bipartite
                        B_c = metaSub.sppPool.bMat_c; // store current state
                        // random uniform variables generated for allocating trophic level
                        vec trophLev(no_invaders);
                        if (metaSub.sppPool.bMat_p.n_rows > 5) { // seed metacommunity with at least 5 producer
                            trophLev.randu();
                        } else {
                            trophLev.zeros();
                        }
                        uvec invaderIndex = find(trophLev <= metaSub.sppPool.pProducer);
                        no_trophLev[0] = invaderIndex.n_rows; // number of producers to invade
                        no_trophLev[1] = no_invaders - no_trophLev[0]; // number of consumers to invade
                    } else {
                        no_trophLev[0] = no_invaders;
                    }

                    for (int tL = 0; tL < 2; tL++) { // first invade producers, then consumers
                        double spp_tested=0, suc_inv=0;

                        if (no_trophLev[tL] == 0) {
                            continue;
                        } else {

                            do {
                                // sample random invaders and simulate dynamics
                                spp_tested += no_trophLev[tL]*2;
                                uvec posGrowth = metaSub.invaderSample(tL, no_trophLev[tL]);

                                // select desired number of invaders with positive growth rates
                                suc_inv += posGrowth.n_rows;
                                posGrowth.resize(min((int) posGrowth.n_rows, no_trophLev[tL]));

                                // remove unsucessful/excess invaders
                                mat bInv_max;
                                bInv_max = metaSub.invaderCleanup(tL, posGrowth);

                                if (posGrowth.n_rows == 0) {
                                    continue;
                                }

                                // populate interaction coefficients and reset invader biomass
                                metaSub.invaderPopulate(tL, bInv_max);
                                no_trophLev[tL] -= (int) posGrowth.n_rows;

                            } while (no_trophLev[tL] > 0);
                        }

                        // record numerical invasion probability
                        metaSub.invasionProb.resize(metaSub.sppPool.invEvent, 2);
                        if (tL == 0) {
                            metaSub.invasionProb(metaSub.sppPool.invEvent - 1, tL) =
                                    suc_inv / spp_tested;
                        } else if (tL == 1) {
                            metaSub.invasionProb(metaSub.sppPool.invEvent - 1, tL) =
                                    suc_inv / spp_tested;
                        }
                    }

                    // reset resident biomasses and update regional meta object
                    if (metaSub.sppPool.bMat_p.n_rows > B_p.n_rows) {
                        if (B_p.n_rows > 0) {
                            metaSub.sppPool.bMat_p.rows(0, B_p.n_rows - 1) = B_p;
                        }
                    }

                    if (metaSub.sppPool.bMat_c.n_rows > B_c.n_rows) {
                        if (B_c.n_rows > 0) {
                            metaSub.sppPool.bMat_c.rows(0, B_c.n_rows - 1) = B_c;
                        }
                        if (rank == 0) {
                            meta.sppPool.aMat = metaSub.sppPool.aMat;
                        }
                    }

                    metaSub.sppPool.invasion += no_invaders;
                } // invader testing scope


                { // dynamics scope

                    metaSub.metaCDynamics(metaSub.tMax); // simulate metacommunty dynamics

                } // dynamics scope

                { // extinction scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Remove extinct species /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    metaSub.sppPool.extinct(); // remove extinct species

                } // extinction scope

                MPI_Barrier(MPI_COMM_WORLD);

                if (rank == 0) {
                    printf("\rInvasions / S_p / S_c = %d / %d / %d         ",
                           metaSub.sppPool.invasion, (int) metaSub.sppPool.bMat_p.n_rows, (int) metaSub.sppPool.bMat_c.n_rows);
                    fflush(stdout);
                    if (metaSub.sppPool.invasion >= metaSub.invMax) {
                        sim = 0; // switch simulation off
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);

                if (size > 1) {
                    // broadcast switch from root to end simulation if invMax reached
                    MPI_Bcast(&sim, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            } while (sim);
        } // assembly scope

        if (rank == 0) {
            mpi_time2 = MPI_Wtime();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Final book keeping /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (rank == 0) {
        if (OUTPUT) {
            if (ASSEMBLE) {
                if (size > 1) {
                    metaSub = meta; // in serial program meta remains empty during assembly, output from metaSub
                }
                metaSub.metaCDynamics(1000); // final relaxation
                metaSub.sppPool.extinct();
                time(&time2); // stop timing
                metaSub.simTime += time2 - time1; // record assembly time
                metaSub.outputData();

                printf("\n\nDetermining source-sink populations...");
                metaSub.genSourceSink();
                metaSub.outputData();
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    longjmp(jump_buffer, 1); // jump to return
}
