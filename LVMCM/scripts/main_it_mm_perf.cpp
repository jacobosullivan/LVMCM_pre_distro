////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// The parallelizable Lotka-Volterra Metacommunity assembly Model (pLVMCM) ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Jacob Dinner O'Sullivan -- j.l.dinner@qmul.ac.uk | j.osullivan@zoho.com ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// A metacommunity assembly model parallelised via domain decomposition using MPI ////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* This script runs a simple performance test of the invader testing algorithm followed by simulation of metacommunity
 * dynamics and extinction protocols for the serial and parallel numerical approximations of the LVMCM
 * A template spatial network is used that can be spetrally bisected up to 6 times (64 subdomains/CPUs)
 * A template community of 100 species is loaded so that the invasion probability is <1 (as in an empty landscape)
 *
 * The invader testing method is now full parallelized and optimized such that the invasion probabiltity much higher
 * than previously
 *
 * For this single iteration of the assembly algorithm I see a clear improvment in performance on the HPC, particularly
 * when OMP_NUM_THREADS=1. When hyper-threading is permitted in the serial case, I get performance gains, but losses in
 * the MPI case so I need to be careful not to conflate these two effects
*/

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
#include <unistd.h>

#include "Metacommunity.h"
#include "LVMCM_rng.h"

using namespace std;
using namespace arma;
using std::vector;

// Assembly flags
bool PARALLEL = false; // select parallel algorithm
bool OUTPUT = true; // select write to file
bool TRAJECTORY = false; // select generate and write trajectory object to file
bool FLUCTUATE = false; // select generate and write to file abiotic fluctuation
bool CMAT_REG = false; // select generate and write to file regional competitive overlap matrix

// Storage for setjmp/longjmp
jmp_buf jump_buffer;

// Global Metacommunity object and timing variables
Metacommunity meta;
time_t time1, time2;
double mpi_time1, mpi_time2;

// Random seed
int g_seed;

// CVode tolerances
double TolA = 1e-8;
double TolR = 1e-7;

// Signal hander flag
bool sig_recvd = false;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Machinery for MPI message passing of Armadillo objects  ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ArmadilloMPI { // Convert armadillo objects to arrays for MPI send/receive/broadcast
// adapted from https://stackoverflow.com/questions/25663998/pass-an-armadillo-c-matrix-over-mpi
public:
    vector<double> bufferReg, bufferLoc;
    vector<int> bufferLenRcv, bufferLVSnd, bufferVecRcv;
    int nCols;
    int nRows;
    int subCols;
    int totalLen=0;
    vector<int> disp;

    // Constructor
    ArmadilloMPI(int nRows, int nCols=0, int subCols=0, bool vec=false) {
        this->nRows = nRows;
        this->nCols = nCols;
        this->subCols = subCols;
        if (vec) {
            bufferLVSnd.resize(nRows); // 1
            bufferLenRcv.resize(nCols); // size
        } else {
            bufferReg.resize(nRows*nCols);
            bufferLoc.resize(nRows*subCols);
        }
    }

    // Destructor
    ~ArmadilloMPI() {}

    // Store regional scale arma::mat in ArmadilloMPI::bufferReg
    template <class arma_type>
    void armaToArrayReg(arma_type &A) {
        int k=0;
        for(int i = 0; i < A.n_cols; ++i ) {
            for(int j = 0; j < A.n_rows; ++j) {
                bufferReg[k] = A(j,i);
                k++;
            }
        }
    }

    // Store subdomain scale arma::mat in ArmadilloMPI::bufferLoc
    template <class arma_type>
    void armaToArrayLoc(arma_type &A) {
        int k=0;
        for(int i = 0; i < A.n_cols; ++i ) {
            for(int j = 0; j < A.n_rows; ++j) {
                bufferLoc[k] = A(j,i);
                k++;
            }
        }
    }

    // Broadcast data stored in ArmadilloMPI::bufferReg to all processes
    template <class arma_type>
    void armaBcast(arma_type &A, int root=0) {
        MPI_Bcast(&(bufferReg[0]), nRows * nCols, MPI_DOUBLE, root, MPI_COMM_WORLD);
        A.set_size(nRows, nCols);
        int k=0;
        for(int  i = 0; i < nCols; ++i ) {
            for(int j = 0; j < nRows; ++j) {
                A(j,i) = bufferReg[k];
                k++;
            }
        }
    }

    // Scatter data stored in ArmadilloMPI::bufferReg to respective processes
    template <class arma_type>
    void armaScatter(arma_type &A, int root=0) {
        MPI_Scatter(&(bufferReg.at(0)), nRows * subCols, MPI_DOUBLE, &(bufferLoc.at(0)), nRows * subCols, MPI_DOUBLE, root, MPI_COMM_WORLD);
        A.set_size(nRows, subCols);
        int k=0;
        for(int  i = 0; i < subCols; ++i ) {
            for(int j = 0; j < nRows; ++j) {
                A(j,i) = bufferLoc[k];
                k++;
            }
        }
    }

    // Gather data stored in ArmadilloMPI::bufferReg at the root process
    template <class arma_type>
    void armaGather(arma_type &A, int root=0) {
        MPI_Gather(&(bufferLoc.at(0)), nRows * subCols, MPI_DOUBLE, &(bufferReg.at(0)), nRows * subCols, MPI_DOUBLE, root, MPI_COMM_WORLD);
        A.set_size(nRows, nCols);
        int k=0;
        for(int  i = 0; i < nCols; ++i ) {
            for(int j = 0; j < nRows; ++j) {
                A(j,i) = bufferReg[k];
                k++;
            }
        }
    }

    // Allgather data stored in ArmadilloMPI::bufferReg
    template <class arma_type>
    void armaAllgather(arma_type &A) {
        MPI_Allgather(&(bufferLoc.at(0)), nRows * subCols, MPI_DOUBLE, &(bufferReg.at(0)), nRows * subCols, MPI_DOUBLE, MPI_COMM_WORLD);
        A.set_size(nRows, nCols);
        int k=0;
        for(int  i = 0; i < nCols; ++i ) {
            for(int j = 0; j < nRows; ++j) {
                A(j,i) = bufferReg[k];
                k++;
            }
        }
    }

    // Allgatherv data for unequal vector concatenation
    template <class arma_type>
    void armaAllgatherv(arma_type &A) {
        bufferLVSnd[0] = A.n_rows;
        MPI_Allgather(&(bufferLVSnd[0]), 1, MPI_INT, &(bufferLenRcv[0]), 1, MPI_INT, MPI_COMM_WORLD);

        disp.resize(nCols);
        for (int i=0; i<nCols; i++) {
            disp[i] = totalLen;
            totalLen += bufferLenRcv[i];
        }

        bufferLVSnd.resize(A.n_rows);
        bufferVecRcv.resize(totalLen);
        for (int i=0; i<A.n_rows; i++) {
            bufferLVSnd[i] = A(i);
        }

        int* count = &bufferLenRcv[0];
        MPI_Allgatherv(&(bufferLVSnd.at(0)), A.n_rows, MPI_INT, &(bufferVecRcv.at(0)), count, &(disp.at(0)), MPI_INT, MPI_COMM_WORLD);

        A.set_size(totalLen);
        for(int  i = 0; i < totalLen; ++i ) {
            A(i) = bufferVecRcv[i];
        }
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Signal handler for write to file on receipt of interrupt signal  ///////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void sigHandler_term(int signum) {
    if (!sig_recvd) {
        sig_recvd = true; // prevent handling of multiple signals in quick succession
        if (OUTPUT) {
            if (meta.rank == 0) { // root process only outputs data
                printf("\n\nInterrupt signal (%d) received at rank %d. Cleaning up and checkpointing simulation...",
                        signum, meta.rank);
                if (PARALLEL) {
                    mpi_time2 = MPI_Wtime();
                    meta.simTime += mpi_time2 - mpi_time1; // record wall clock
                    mpi_time1 = mpi_time2;
                } else {
                    time(&time2);
                    meta.simTime += time2 - time1; // record wall clock
                    time1 = time2;
                }
                meta.cleanup(); // clean up invader testing algorithm and output
            }
        }
        longjmp(jump_buffer, 1); // jump to return
    } else {
        cout << "Signal already received" << endl;
    }
}

void registerSigHandler(int *sigList) {
    for (int i=0; sigList[i]; i++)
    {
        if (SIG_ERR == signal(sigList[i], sigHandler_term))
        {
            printf("Can not catch signal : %d", sigList[i]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Start of assembly algoritm ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    time_t start;
    time(&start);
    int sim = 1; // simulation switch
    int signalList[] = { //SIGINT, // select interupt signals caugth by handler
            //SIGTERM,
//                          SIGUSR1,
//                          SIGUSR2,
            0};
    registerSigHandler(signalList); // register handlers

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Set random seeds  //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Note to set random seeds uncomment *set_seed(1) AND *boost_rng(1) in the file LVMCM_rng.cpp and recompile
//    arma_rng::set_seed_random();
    arma_rng::set_seed(1);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Store program arguments in variables ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Default Metacommunity parameters
    int a_init = 1;
    string a_bMat = {};
    int a_invMax = 0;
    int a_iterS = 2;
    int a_deltaT = 100;
    int a_tMax = 500;
    int a_perfRep = 0;
    int a_perfS = 0;
    bool a_autoFluct = false;
    bool a_harvest = false;
    int a_locRandomize = 0;
    string a_outputDirectory;
    // Default Species parameters
    double a_c1 = 0.5;
    double a_c2 = 0.5;
    double a_emRate = 0.01;
    double a_dispL = 0.1;
    double a_pProducer = 1.0;
    bool a_prodComp = true;
    double a_alpha = 0;
    double a_sigma = 0;
    double a_sigma_t = 0.1;
    double a_rho = 0;
    bool a_discr_c_ij = true;
    // Default Topograpy parameters
    int a_no_nodes = 1024;
    double a_phi = 1;
    int a_envVar = 0;
    double a_var_e = 1.0;
    bool a_randGraph = true;
    bool a_gabriel = true;
    int a_bisec = 0;
    // Default output variables
    double a_parOut = 0;
    string a_experiment = "Performance";
    int a_rep = 0;
    string a_jobID = "NA";

    for (int i = 1; i < argc; i++) { // loop through program arguments an allocate to parameters
        char var = argv[i][1];
        i++;
        switch (var) {
            case 'n' : // -new: initialize new model (T/F)
                if (!strcmp(argv[i], "F")) {
                    a_init = 0;
                }
                break;
            case 'b' : // -bFile: path to data for importing initialized model
                a_bMat = argv[i];
                break;
            case 'i' : // -invMax: total number of invasions
                a_invMax = atoi(argv[i]);
                break;
            case 'S' : // -Schwartz: Schwartz iteration and time window
                a_iterS = atoi(argv[i]);
                i++;
                a_deltaT = atoi(argv[i]);
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
                if (!strcmp(argv[i], "F")) {
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
                if (!strcmp(argv[i], "F")) {
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
            case 'R' : // -Rand: select random spatial network/lattice
                if (!strcmp(argv[i], "F")) {
                    a_randGraph = false;
                }
                break;
            case 'G' : // -Gab: select Gabriel/complete graph
                if (!strcmp(argv[i], "F")) {
                    a_gabriel = false;
                }
                break;
            case 'B' : // -Bisec: number of recursive spectral bisections
                a_bisec = atoi(argv[i]);
                PARALLEL = true;
                break;
            case 'o' : // -o: output variables - key parameter, experiment name, replicate number
                a_parOut = atof(argv[i]);
                i++;
                a_experiment = argv[i];
                i++;
                a_rep = atoi(argv[i]);
                break;
            case 'O' : // -O: select write to file
                if (!strcmp(argv[i], "F")) {
                    OUTPUT = false;
                }
                break;
            case 'I' : // -ID: job ID used to generate .txt file recording current state of simulation
                a_jobID = argv[i];
                break;
            case 'T' : // -Test: performance testing parameters
                a_perfRep = atoi(argv[i]);
                i++;
                a_perfS = atoi(argv[i]);
                break;
            case 'A' : // -AF: store trajectory
                a_autoFluct = true;
                break;
            case 'h' : // -harvest: estimate regional scale interaction coefficients using harvesting experiment
                if (!strcmp(argv[i], "F")) {
                    a_harvest = true;
                }
                break;
            case 'l' : // -locR: randomize local communities
                a_locRandomize = atoi(argv[i]);
                break;
            case 'M' : // -MPI: select parallel algorithm
                PARALLEL = true;
                break;
            case 'f' : // -folder: set output directory
                a_outputDirectory = argv[i];
                break;
        }
    }

    if (!PARALLEL) { // Run non-parallelized assembly algorithm

        if (setjmp(jump_buffer)) { // set return call in longjump clause
            time_t finish;
            time(&finish);
            cout << "\n\nSimulation stated at: " << ctime(&start);
            cout << "Simulation finished at: " << ctime(&finish);
            return (0);
        }

        printf("\nInitializing LVMCM\n\n");
        a_init = 0;
//        a_bMat = "/home/jack/gitClones/LVMCM_src/LVMCM/SimulationData/N=1024/Performance_experiment/2020-1-23/2020-1-23_Performance(0)0bMat0.mat";
        a_bMat = "/data/home/btx188/SimulationData/N=1024/Performance_experiment/2020-1-23/2020-1-23_Performance(0)0bMat0.mat";


        cout << a_outputDirectory << endl;

        meta = Metacommunity(
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
                // output variables
                a_parOut,
                a_experiment,
                a_rep,
                a_jobID);
        cout << "output folder set to " << meta.outputDirectory << endl;
        meta.printParams();
        cout << endl << endl;
        // generate random seed
        std::random_device rd;
        g_seed = rd();

        LVMCM_rng::boost_rng.seed(g_seed);
        arma_rng::set_seed(g_seed);

        cout << "\nNetwork rows 0-4 = " << endl;
        cout << meta.sppPool.topo.network.rows(0, 4) << endl;
        meta.sppPool.genDispMat(); // generate adjacency/dispersal operator

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Serial metacommunity dynamics  ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        printf("\nStarting non-parallel simulation\n\n");

        time(&time1); // start timing

        meta.sppPool.invEvent++;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Invader testing algorithm  /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        { // invader testing scope

            // select number of invaders and sample their trophic levels
            meta.sppPool.S_p = meta.sppPool.bMat_p.n_rows; // store diversity
            meta.sppPool.S_c = meta.sppPool.bMat_c.n_rows;
            int no_trophLev[2] = {0};
            int no_invaders = 0.05 * (meta.sppPool.S_p + meta.sppPool.S_c) + 1;

            cout << "Generating " << no_invaders << " invaders" << endl;

            mat B_p, B_c;
            B_p = meta.sppPool.bMat_p; // store current state
            B_c = meta.sppPool.bMat_c; // store current state

            if (meta.sppPool.pProducer < 1) {
                B_c = meta.sppPool.bMat_c; // store current state
                // random uniform variables generated for allocating trophic level
                vec trophLev(no_invaders);
                if (meta.sppPool.bMat_p.n_rows > 5) { // seed metacommunity with at least 5 producer
                    trophLev.randu();
                } else {
                    trophLev.zeros();
                }
                uvec invaderIndex = find(trophLev <= meta.sppPool.pProducer);
                no_trophLev[0] = invaderIndex.n_rows; // number of producers to invade
                no_trophLev[1] = no_invaders - no_trophLev[0]; // number of consumers to invade
            } else {
                no_trophLev[0] = no_invaders;
            }

            for (int tL = 0; tL < 2; tL++) { // first invade produces, then consumers
                if (no_trophLev[tL] == 0) {
                    continue;
                } else {

                    do {
                        // sample random invaders and simulate dynamics
                        uvec posGrowth = meta.invaderSample(tL, no_trophLev[tL]);
                        meta.invasionProb.resize(meta.sppPool.invEvent, 2);
                        if (tL == 0) {
                            meta.invasionProb(meta.sppPool.invEvent - 1, tL) =
                                    posGrowth.n_rows / (meta.sppPool.bMat_p.n_rows - meta.sppPool.S_p);
                        } else if (tL == 1) {
                            meta.invasionProb(meta.sppPool.invEvent - 1, tL) =
                                    posGrowth.n_rows / (meta.sppPool.bMat_c.n_rows - meta.sppPool.S_c);
                        }

                        // select desired number of invaders with positive growth rates
                        posGrowth.resize(min((int) posGrowth.n_rows, no_trophLev[tL]));

                        mat bInv_max;

                        // remove unsucessful/excess invaders
                        bInv_max = meta.invaderCleanup(tL, posGrowth);

                        // populate interaction coefficients and reset invader biomass
                        meta.invaderPopulate(tL, bInv_max);

                        no_trophLev[tL] -= (int) posGrowth.n_rows;

                    } while (no_trophLev[tL] > 0); // repeat until desired number of invaders found

                    // reset resident biomasses
                    if (meta.sppPool.bMat_p.n_rows > B_p.n_rows) {
                        meta.sppPool.bMat_p.rows(0, B_p.n_rows - 1) = B_p;
                    }

                    if (meta.sppPool.bMat_c.n_rows > B_c.n_rows) {
                        meta.sppPool.bMat_c.rows(0, B_c.n_rows - 1) = B_c;
                    }

                    meta.sppPool.invasion += no_invaders;
                }
            }
        } // invader testing scope

        { // dynamics scope
            cout << "\nSimulating" << endl;
            meta.storeTraj=1;
            meta.metaCDynamics(meta.tMax, -1); // simulate metacommunty dynamics
        } // dynamics scope

        { // extinction scope
            meta.sppPool.extinct(); // remove extinct species
        } // extinction scope

        time(&time2); // stop timing
        meta.simTime = time2-time1;
        meta.writePars(a_perfRep);
        longjmp(jump_buffer,1); // jump to return

    } else if (PARALLEL) {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Initialise MPI communicator ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        int rank;
        int size;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//        { // machinery for attaching serial debugger to MPI program: https://www.open-mpi.org/faq/?category=debugging
////            if (rank == 0) {
//                volatile int i = 0;
//                char hostname[256];
//                gethostname(hostname, sizeof(hostname));
//                printf("PID %d on %s ready for attach\n", getpid(), hostname);
//                fflush(stdout);
//                while (0 == i) {
//                    sleep(5);
//                    fflush(stdout);
//                }
////            }
////            MPI_Barrier(MPI_COMM_WORLD);
//        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Set return clause parallel algorithm ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        Metacommunity metaSub;

        if (setjmp(jump_buffer)) { // set return call in longjump clause
            for (int s = 0; s < size; s++) {
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
            if (a_bisec > 0) {
                if ((a_iterS == 0) || (a_deltaT == 0)) { // check Schwartz parameters set
                    cout << "\nError: Schwarz iteration parameters not set" << endl;
                    longjmp(jump_buffer, 1);
                }
            }
            if (a_deltaT != 0) {
                double tM_dT = a_tMax / a_deltaT;
                if (floor(tM_dT) - tM_dT != 0) { // check tMax/deltaT correspond
                    cout << "\nError: tMax / deltaT non-integer" << endl;
                    longjmp(jump_buffer, 1);
                }
            }
            if (a_outputDirectory.length() == 0) { // check output directory set
                cout << a_outputDirectory << endl;
                cout << "\nWarning: no output directory given" << endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Parameterize parallel assembly model ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (rank == 0) { // parameterize root process using program arguments

            printf("\nInitializing LVMCM\n\n");
            a_init = 0;
//            a_bMat = "/home/jack/gitClones/LVMCM_src/LVMCM/SimulationData/N=1024/Performance_experiment/2020-1-23/2020-1-23_Performance(0)0bMat0.mat";
            a_bMat = "/data/home/btx188/SimulationData/N=1024/Performance_experiment/2020-1-23/2020-1-23_Performance(0)0bMat0.mat";


            cout << a_outputDirectory << endl;

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
                    // output variables
                    a_parOut,
                    a_experiment,
                    a_rep,
                    a_jobID);
            cout << "output folder set to " << metaSub.outputDirectory << endl;
            metaSub.printParams();
            cout << endl << endl;
            // generate random seed
            std::random_device rd;
            g_seed = rd();
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Broadcast random seed to all processes
        MPI_Bcast(&g_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        for (int s = 0; s < size; s++) {
            if (rank == s) {
                cout << "Random seed seen by process " << rank << " " << g_seed << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Set seeds of boost and armadillo rngs
        LVMCM_rng::boost_rng.seed(g_seed);
        arma_rng::set_seed(g_seed);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Initialize performance test and distribute matrices //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        { // subdomain intialization scope
            int S_p, S_c;
            if (rank == 0) {
                metaSub.sppPool.genDispMat(); // generate adjacency/dispersal operator
                S_p = metaSub.sppPool.bMat_p.n_rows;
                S_c = metaSub.sppPool.bMat_c.n_rows;
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // parameterize subdomains
            // Metacommunity parameters
            MPI_Bcast(&metaSub.invMax, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&metaSub.iterS, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&metaSub.deltaT, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&metaSub.tMax, 1, MPI_INT, 0, MPI_COMM_WORLD);
            // Species parameters
            MPI_Bcast(&metaSub.sppPool.pProducer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&metaSub.sppPool.prodComp, 1, MPI_INT, 0, MPI_COMM_WORLD);
            // Topography parameters
            MPI_Bcast(&metaSub.sppPool.topo.no_nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&metaSub.sppPool.topo.bisec, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&metaSub.sppPool.topo.envVar, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&metaSub.sppPool.topo.var_e, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&metaSub.sppPool.topo.phi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // Species richness
            MPI_Bcast(&S_p, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&S_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
            metaSub.sppPool.S_p = S_p;
            metaSub.sppPool.S_c = S_c;

            MPI_Barrier(MPI_COMM_WORLD);

            // initialize Armadillo buffers
            ArmadilloMPI *bMPI_p, *bMPI_c, *aMPI, *rMPI, *cMPI, *tMPI;

            bMPI_p = new ArmadilloMPI(S_p, metaSub.sppPool.topo.no_nodes);
            if (S_c > 0) {
                bMPI_c = new ArmadilloMPI(S_c, metaSub.sppPool.topo.no_nodes);
                aMPI = new ArmadilloMPI(S_p, S_c);
            }
            rMPI = new ArmadilloMPI(S_p, metaSub.sppPool.topo.no_nodes);
            if (metaSub.sppPool.prodComp) {
                cMPI = new ArmadilloMPI(S_p, S_p);
            }
            if (metaSub.sppPool.topo.envVar > 0) {
                tMPI = new ArmadilloMPI(metaSub.sppPool.topo.envVar, S_p);
            }

            if (rank == 0) { // convert Armadillo matrices to arrays
                bMPI_p->armaToArrayReg(metaSub.sppPool.bMat_p);
                if (S_c > 0) {
                    bMPI_c->armaToArrayReg(metaSub.sppPool.bMat_c);
                    aMPI->armaToArrayReg(metaSub.sppPool.aMat);
                }
                rMPI->armaToArrayReg(metaSub.sppPool.rMat);
                if (metaSub.sppPool.prodComp) {
                    cMPI->armaToArrayReg(metaSub.sppPool.cMat);
                }
                if (metaSub.sppPool.topo.envVar != 0) {
                    tMPI->armaToArrayReg(metaSub.sppPool.tMat);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // broadcast model matrices and store in local Metacommunity objects
            bMPI_p->armaBcast(metaSub.sppPool.bMat_p);
            if (S_c > 0) {
                bMPI_c->armaBcast(metaSub.sppPool.bMat_c);
                aMPI->armaBcast(metaSub.sppPool.aMat);
            }

            rMPI->armaBcast(metaSub.sppPool.rMat);
            if (metaSub.sppPool.prodComp) {
                cMPI->armaBcast(metaSub.sppPool.cMat);
            }
            if (metaSub.sppPool.topo.envVar != 0) {
                tMPI->armaBcast(metaSub.sppPool.tMat);
            }

            MPI_Barrier(MPI_COMM_WORLD);
        } // subdomain intialization scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Generate topo - topography/environment (root) //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        { // domain decomposition scope

            if (rank == 0) {
                metaSub.sppPool.topo.genDomainDecomp(
                        metaSub.sppPool.topo.network); // decompose imported topo
                if (metaSub.sppPool.topo.fVec.n_rows == 0) { // exit if decomposition of imported network fails
                    cout << "\nError: domain decomposition failure" << endl;
                    longjmp(jump_buffer, 1);
                }
                metaSub.sppPool.topo.genDistMat();
                metaSub.sppPool.topo.genAdjMat();
                metaSub.sppPool.genDispMat();
                unique(metaSub.sppPool.topo.fVec.t()).print("\nf.unique");
                meta = metaSub; // store copy of complete domain at root process for extinction testing
            }
            meta.rank = rank; // store rank for signal handler
            MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Broadcast landcape to subdomains //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // initialize Armadillo buffers
            ArmadilloMPI *netMPI, *indMPI, *adjMPI, *dstMPI, *dspMPI, *envMPI;
            netMPI = new ArmadilloMPI(metaSub.sppPool.topo.no_nodes, 2);
            indMPI = new ArmadilloMPI(metaSub.sppPool.topo.no_nodes, 1);
            adjMPI = new ArmadilloMPI(metaSub.sppPool.topo.no_nodes, metaSub.sppPool.topo.no_nodes);
            dstMPI = new ArmadilloMPI(metaSub.sppPool.topo.no_nodes, metaSub.sppPool.topo.no_nodes);
            dspMPI = new ArmadilloMPI(metaSub.sppPool.topo.no_nodes, metaSub.sppPool.topo.no_nodes);
            if (metaSub.sppPool.topo.envVar > 0) {
                envMPI = new ArmadilloMPI(metaSub.sppPool.topo.envVar, metaSub.sppPool.topo.no_nodes);
            }

            if (rank == 0) { // convert Armadillo matrices to arrays
                netMPI->armaToArrayReg(metaSub.sppPool.topo.network);
                indMPI->armaToArrayReg(metaSub.sppPool.topo.fVec);
                adjMPI->armaToArrayReg(metaSub.sppPool.topo.adjMat);
                dstMPI->armaToArrayReg(metaSub.sppPool.topo.distMat);
                dspMPI->armaToArrayReg(metaSub.sppPool.dMat_n);
                if (metaSub.sppPool.topo.envVar != 0) {
                    envMPI->armaToArrayReg(metaSub.sppPool.topo.envMat);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // broadcast model matrices and store in local Metacommunity objects
            netMPI->armaBcast(metaSub.sppPool.topo.network);
            indMPI->armaBcast(metaSub.sppPool.topo.fVec);
            adjMPI->armaBcast(metaSub.sppPool.topo.adjMat);
            dstMPI->armaBcast(metaSub.sppPool.topo.distMat);
            dspMPI->armaBcast(metaSub.sppPool.dMat_n);
            if (metaSub.sppPool.topo.envVar != 0) {
                envMPI->armaBcast(metaSub.sppPool.topo.envMat);
            }
            metaSub.sppPool.subSet(rank); // each process subsets topo objects according to domain decomposition

        } // domain decomposition scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Simulation /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            printf("\nStarting parallel simulation\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);

        for (int s = 0; s < size; s++) {
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

        // initialise objects
        mat B_p, B_c, BStore_p, BStore_c, U_p, U_c, Ub_p, Ub_c, Tr, R, Cr, Cc, Ar, Ac, invIndex_p, invIndex_c; // MPI buffers
        int no_invaders_p, no_residents_p, no_extinct_p, no_invaders_c, no_residents_c, no_extinct_c; // counters
        field<uvec> index_ext(2); // index extinct species

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Generate invader (all subdomains in parallel) ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        metaSub.sppPool.invEvent++;

        { // invader testing scope

            if (rank == 0) {
                mpi_time1 = MPI_Wtime();
            }

            B_p = metaSub.sppPool.bMat_p; // store current state
            B_c = metaSub.sppPool.bMat_c; // store current state

            // select number of invaders and sample their trophic levels
            int no_trophLev[2] = {0};
            int no_invaders = 0.05 * (metaSub.sppPool.S_p + metaSub.sppPool.S_c) + 1;

            cout << "Generating " << no_invaders << " invaders" << endl;

            if (metaSub.sppPool.pProducer < 1) {
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

            ArmadilloMPI *bMPI_p, *rMPI, *crMPI, *ccMPI, *arMPI, *acMPI, *uMPI_p, *uBcastMPI_p, *iMPI_p, *iIMPI_p;
            // MPI buffers consumers
            ArmadilloMPI *bMPI_c, *uMPI_c, *uBcastMPI_c, *iMPI_c, *iIMPI_c;

            for (int tL = 0; tL < 2; tL++) { // first invade produces, then consumers
                if (no_trophLev[tL] == 0) {
                    continue;
                } else {

                    do {

                        // sample random invaders and simulate dynamics
                        uvec posGrowth = metaSub.invaderSample(tL, no_trophLev[tL]);
                        // Based on min and max of r_0 it looks like the seed setting is working but writing to file might be a good idea

                        // gather indices of successful invaders
                        ArmadilloMPI *pgMPI;
                        pgMPI = new ArmadilloMPI(1, size, 0, true);
                        pgMPI->armaAllgatherv(posGrowth);
                        MPI_Barrier(MPI_COMM_WORLD);

                        // select desired number of invaders with positive growth rates
                        posGrowth = unique(posGrowth);

                        // record numerical invasion probability
                        metaSub.invasionProb.resize(metaSub.sppPool.invEvent, 2);
                        if (tL == 0) {
                            metaSub.invasionProb(metaSub.sppPool.invEvent - 1, tL) =
                                    posGrowth.n_rows / (metaSub.sppPool.bMat_p.n_rows - metaSub.sppPool.S_p);
                        } else if (tL == 1) {
                            metaSub.invasionProb(metaSub.sppPool.invEvent - 1, tL) =
                                    posGrowth.n_rows / (metaSub.sppPool.bMat_c.n_rows - metaSub.sppPool.S_c);
                        }
                        posGrowth.resize(min((int) posGrowth.n_rows, no_trophLev[tL]));

                        // remove unsucessful/excess invaders
                        mat bInv_max;
                        bInv_max = metaSub.invaderCleanup(tL, posGrowth);
                        MPI_Barrier(MPI_COMM_WORLD);

                        // gather maximum subdomain biomass for selection of port node
                        ArmadilloMPI *bimMPI;
                        bimMPI = new ArmadilloMPI(no_trophLev[tL], size, 1);
                        bimMPI->armaToArrayLoc(bInv_max);
                        bimMPI->armaAllgather(bInv_max);
                        MPI_Barrier(MPI_COMM_WORLD);

                        // populate interaction coefficients and reset invader biomass
                        metaSub.invaderPopulate(tL, bInv_max);

                        no_trophLev[tL] -= (int) posGrowth.n_rows;

                    } while (no_trophLev[tL] > 0);
                }
            }

            // reset resident biomasses and update meta
            if (metaSub.sppPool.bMat_p.n_rows > B_p.n_rows) {
                if (B_p.n_rows > 0) {
                    metaSub.sppPool.bMat_p.rows(0, B_p.n_rows - 1) = B_p;
                }
                ArmadilloMPI *riMPI;
                riMPI = new ArmadilloMPI(metaSub.sppPool.bMat_p.n_rows - metaSub.sppPool.S_p,
                                         metaSub.sppPool.topo.no_nodes,
                                         metaSub.sppPool.topo.network.n_rows);
                R = metaSub.sppPool.rMat.rows(metaSub.sppPool.S_p, metaSub.sppPool.rMat.n_rows-1);
                riMPI->armaToArrayLoc(R);
                riMPI->armaAllgather(R);
                if (rank == 0) {
                    if (meta.sppPool.rMat.n_rows == 0) {
                        meta.sppPool.rMat = R;
                    } else {
                        meta.sppPool.rMat = join_vert(meta.sppPool.rMat, R);
                    }
                    meta.sppPool.cMat = metaSub.sppPool.cMat;
                    if (meta.sppPool.topo.envVar > 0) {
                        meta.sppPool.tMat = metaSub.sppPool.tMat;
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);
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

            cout << "\nSimulating" << endl;

            ArmadilloMPI *uMPI_p, *uMPI_c;

            int t, it, T = metaSub.tMax / metaSub.deltaT; // number of Schwarz time windows

            cout << "T = " << T << endl;

            uMPI_p = new ArmadilloMPI(metaSub.sppPool.bMat_p.n_rows, metaSub.sppPool.topo.no_nodes,
                                      metaSub.sppPool.topo.network.n_rows);

            if (metaSub.sppPool.pProducer < 1) {
                uMPI_c = new ArmadilloMPI(metaSub.sppPool.bMat_c.n_rows, metaSub.sppPool.topo.no_nodes,
                                          metaSub.sppPool.topo.network.n_rows);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            for (t = 0; t < T; t++) {
                cout << "t " << t << endl;
                mat BStore_p = metaSub.sppPool.bMat_p; // save inital state for Schwartz iteration updates
                mat BStore_c;
                if (metaSub.sppPool.pProducer < 1) {
                    BStore_c = metaSub.sppPool.bMat_c; // save inital state for Schwartz iteration updates
                }
                for (it = 0; it < metaSub.iterS; it++) {
                    cout << "it " << it << endl;
                    metaSub.sppPool.bMat_p = BStore_p; // restore inital state at entry to Schwartz iteration
                    if (metaSub.sppPool.pProducer < 1) {
                        metaSub.sppPool.bMat_c = BStore_c; // restore inital state at entry to Schwartz iteration
                    }
                    metaSub.metaCDynamics(metaSub.deltaT, rank); // simulation subdomain scale dynamics

                    // allgather to update fixed unknowns
                    U_p = (metaSub.sppPool.bMat_p + BStore_p) / 2;
                    if (metaSub.sppPool.pProducer < 1) {
                        U_c = (metaSub.sppPool.bMat_c + BStore_c) / 2;
                    }
                    uMPI_p->armaToArrayLoc(U_p);
                    if (metaSub.sppPool.pProducer < 1) {
                        uMPI_c->armaToArrayLoc(U_c);
                    }
                    uMPI_p->armaAllgather(U_p);
                    if (metaSub.sppPool.pProducer < 1) {
                        uMPI_c->armaAllgather(U_c);
                    }
                    metaSub.sppPool.uMat_p = U_p.cols(metaSub.sppPool.topo.adjIF);

                    // immigration from adjacent subdomains computed once per Schwartz iteration
                    metaSub.sppPool.uMat_p = metaSub.sppPool.uMat_p * metaSub.sppPool.dMat_m;
                    if (metaSub.sppPool.pProducer < 1) {
                        metaSub.sppPool.uMat_c = U_c.cols(metaSub.sppPool.topo.adjIF);
                        metaSub.sppPool.uMat_c = metaSub.sppPool.uMat_c * metaSub.sppPool.dMat_m;
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Gather result at root process //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            ArmadilloMPI *bMPI_p, *bMPI_c;

            bMPI_p = new ArmadilloMPI(metaSub.sppPool.bMat_p.n_rows, metaSub.sppPool.topo.no_nodes,
                                      metaSub.sppPool.topo.network.n_rows);
            bMPI_p->armaToArrayLoc(metaSub.sppPool.bMat_p);
            bMPI_p->armaGather(B_p);

            if (metaSub.sppPool.pProducer < 1) {
                bMPI_c = new ArmadilloMPI(metaSub.sppPool.bMat_c.n_rows, metaSub.sppPool.topo.no_nodes,
                                          metaSub.sppPool.topo.network.n_rows);
                bMPI_c->armaToArrayLoc(metaSub.sppPool.bMat_c);
                bMPI_c->armaGather(B_c);
            }

            // receive current state at root and generate whole domain biomass matrix
            if (rank == 0) {
                meta.sppPool.bMat_p = B_p;
                if (metaSub.sppPool.pProducer < 1) {
                    meta.sppPool.bMat_c = B_c;
                }
            }

        } // dynamics scope

        { // extinction scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Remove extinct species /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            field<uvec> index_ext(2); // index extinct species
            // receive current state at root and generate whole domain biomass matrix
            if (rank == 0) {
                // check for extinctions at the whole domain scale and store indices
                index_ext = meta.sppPool.extinct();
                no_extinct_p = (int) index_ext(0).n_rows;
                no_extinct_c = (int) index_ext(1).n_rows;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&no_extinct_p, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast number of extinctions producers
            MPI_Bcast(&no_extinct_c, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast number of extinctions consumers
            MPI_Barrier(MPI_COMM_WORLD);

            ArmadilloMPI *iMPI_p, *iMPI_c;

            iMPI_p = new ArmadilloMPI(no_extinct_p, 1);
            if (no_extinct_c > 0) {
                iMPI_c = new ArmadilloMPI(no_extinct_c, 1); // set size buffer for extinction index
            }

            if (rank == 0) {
                iMPI_p->armaToArrayReg(index_ext(0));
                if (no_extinct_c > 0) {
                    iMPI_c->armaToArrayReg(index_ext(1)); // set size buffer for extinction index
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            // broadcast index of extinct species
            if (no_extinct_p > 0) {
                iMPI_p->armaBcast(index_ext(0));
            }
            if (no_extinct_c > 0) {
                iMPI_c->armaBcast(index_ext(1));
            }
            if (no_extinct_p + no_extinct_c > 0) {
                metaSub.sppPool.extinct(0, index_ext(0), index_ext(1));
            }

            metaSub.sppPool.S_p = metaSub.sppPool.bMat_p.n_rows;
            metaSub.sppPool.S_c = metaSub.sppPool.bMat_c.n_rows;

        } // extinction scope

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// Final book keeping /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0) {
            mpi_time2 = MPI_Wtime();
            metaSub.simTime = mpi_time2-mpi_time1;
            metaSub.writePars(a_perfRep);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        longjmp(jump_buffer, 1); // all proceses jump to return
    }

//    { // print call for debugging parallel invader test sequence
//        for (int s = 0; s < size; s++) {
//            if (rank == s) {
//                cout << "\nFinal message received at rank " << rank << endl;
//                posGrowth.print("PG total");
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        for (int s = 0; s < size; s++) {
//            if (rank == s) {
//                cout << "\nS_p rank " << rank << " "
//                     << metaSub.sppPool.bMat_p.n_rows << " "
//                     << metaSub.sppPool.rMat.n_rows << " "
//                     << metaSub.sppPool.cMat.n_rows << " "
//                     << metaSub.sppPool.cMat.n_cols << " " << endl;
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        for (int s = 0; s < size; s++) {
//            if (rank == s) {
//                cout << "\n1bInv_max from rank " << rank << endl;
//                bInv_max.print("BIM");
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        for (int s = 0; s < size; s++) {
//            if (rank == s) {
//                cout << "\n2bInv_max from rank " << rank << endl;
//                bInv_max.print("BIM");
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        for (int s = 0; s < size; s++) {
//            if (rank == s) {
//                cout << "\nno_trophLev at rank " << rank << " " << no_trophLev[tL] << endl;
//
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
}
