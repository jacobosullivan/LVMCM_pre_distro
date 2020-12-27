////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// The parallelizable Lotka-Volterra Metacommunity assembly Model (pLVMCM) ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Jacob Dinner O'Sullivan -- j.l.dinner@qmul.ac.uk | j.osullivan@zoho.com ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// A metacommunity assembly model parallelised via domain decomposition using MPI ////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* This script runs a simple performance test of the invader testing algorithm (only) for the serial and
 * parallel numerical approximations of the LVMCM
 * A template spatial network is used that can be spetrally bisected up to 6 times (64 subdomains/CPUs)
 * A template community of 100 species is loaded so that the invasion probability is <1 (as in an empty landscape)
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

// CVode tolerances
double TolA = 1e-8;
double TolR = 1e-7;

// Signal hander flag
bool sig_recvd = false;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Machinery for MPI message passing of Armadillo objects  ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ArmadilloMPI { // Convert armadillo objects to arrays for MPI send/receive/broadcast
public:
    double *bufferReg, *bufferLoc;
    int nCols;
    int nRows;
    int subCols;

    // Constructor
    ArmadilloMPI(int nRows, int nCols, int subCols=0) {
        this->nRows = nRows;
        this->nCols = nCols;
        this->subCols = subCols;
        bufferReg = new double[nRows*nCols];
        bufferLoc = new double[nRows*subCols];
    }

    // Deconstructor
    ~ArmadilloMPI() {
        delete [] bufferReg;
        delete [] bufferLoc;
    }

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
        MPI_Scatter(&(bufferReg[0]), nRows * subCols, MPI_DOUBLE, &(bufferLoc[0]),
                nRows * subCols, MPI_DOUBLE, root, MPI_COMM_WORLD);
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
        MPI_Gather(&(bufferLoc[0]), nRows * subCols, MPI_DOUBLE, &(bufferReg[0]),
                nRows * subCols, MPI_DOUBLE, root, MPI_COMM_WORLD);
        A.set_size(nRows, nCols);
        int k=0;
        for(int  i = 0; i < nCols; ++i ) {
            for(int j = 0; j < nRows; ++j) {
                A(j,i) = bufferReg[k];
                k++;
            }
        }
    }

    // Allgather data stored in ArmadilloMPI::bufferReg at the root process
    template <class arma_type>
    void armaAllgather(arma_type &A) {
        MPI_Allgather(&(bufferLoc[0]), nRows * subCols, MPI_DOUBLE, &(bufferReg[0]),
                nRows * subCols, MPI_DOUBLE, MPI_COMM_WORLD);
        A.set_size(nRows, nCols);
        int k=0;
        for(int  i = 0; i < nCols; ++i ) {
            for(int j = 0; j < nRows; ++j) {
                A(j,i) = bufferReg[k];
                k++;
            }
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
    int signalList[] =  { //SIGINT, // select interupt signals caugth by handler
                          //SIGTERM,
//                          SIGUSR1,
//                          SIGUSR2,
                          0 };
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
    int a_iterS = 0;
    int a_deltaT = 0;
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
            case 'R' : // -Rand: select random spatial network/lattice
                if (!strcmp(argv[i],"F")) {
                    a_randGraph = false;
                }
                break;
            case 'G' : // -Gab: select Gabriel/complete graph
                if (!strcmp(argv[i],"F")) {
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
                if (!strcmp(argv[i],"F")) {
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
                if (!strcmp(argv[i],"F")) {
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Set return clause and check for parameter errors ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (setjmp(jump_buffer)) { // set return call in longjump clause
            time_t finish;
            time(&finish);
            cout << "\n\nSimulation stated at: " << ctime(&start);
            cout << "Simulation finished at: " << ctime(&finish);
            return (0);
        }

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
        if (OUTPUT) {
            if (a_outputDirectory.length() == 0) { // check output directory set
                cout << a_outputDirectory << endl;
                cout << "\nWarning: no output directory given" << endl;
            }
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Initialize serial assembly model ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        a_init = 0;
        a_bMat = "/home/jack/gitClones/LVMCM_src/LVMCM/SimulationData/N=1024/Performance_experiment/2020-1-23/2020-1-23_Performance(0)0bMat0.mat";
        printf("\nInitializing LVMCM\n\n");
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
        meta.printParams();
//        meta.sppPool.topo.network.load("../PerformanceTemplate/X.mat");
//        meta.sppPool.topo.network.load("/data/home/btx188/LVMCM_src/LVMCM/PerformanceTemplate/X.mat");
        cout << "\nNetwork rows 0-4 = " << endl;
        cout << meta.sppPool.topo.network.rows(0,4) << endl;
        meta.sppPool.genDispMat(); // generate adjacency/dispersal operator

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Serial metacommunity dynamics  ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        printf("\nStarting non-parallel simulation\n\n");
        // Create function for assembly performance testing - all it needs to do is to output the parameter file
        // Also add to the parallel branch and try out the MPI profiling - how much time spent on invader testing?!

        time(&time1); // start timing
        meta.invaderTest(-1); // search for successful invader
        cout << "Invader found" << endl;
//        meta.metaCDynamics(meta.tMax, -1); // simulate metacommunty dynamics
//        meta.sppPool.extinct(); // remove extinct species
        time(&time2); // start timing
        meta.simTime = time2-time1;
        meta.writePars(a_perfRep);
        longjmp(jump_buffer,1); // jump to return

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Final book keeping /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    } else if (PARALLEL) { // Run parallelized assembly algorithm

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Initialise MPI communicator ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        int rank;
        int size;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//        { // machinery for attaching serial debugger to MPI program: https://www.open-mpi.org/faq/?category=debugging
//            if (rank == 0) {
//                volatile int i = 0;
//                char hostname[256];
//                gethostname(hostname, sizeof(hostname));
//                printf("PID %d on %s ready for attach\n", getpid(), hostname);
//                fflush(stdout);
//                while (0 == i) {
//                    sleep(5);
//                    fflush(stdout);
//                }
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }

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
            a_bMat = "/home/jack/gitClones/LVMCM_src/LVMCM/SimulationData/N=1024/Performance_experiment/2020-1-23/2020-1-23_Performance(0)0bMat0.mat";

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
        }
        MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Initialize performance test and distribute matrices //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            int S_p, S_c;
            if (rank == 0) {
//                metaSub.sppPool.topo.network.load("/home/jack/gitClones/LVMCM_src/LVMCM/PerformanceTemplate/X.mat");
//                metaSub.sppPool.topo.network.load("/data/home/btx188/LVMCM_src/LVMCM/PerformanceTemplate/X.mat");
                cout << "\nNetwork rows 0-4 = " << endl;
                cout << metaSub.sppPool.topo.network.rows(0,4) << endl;
                metaSub.sppPool.genDispMat(); // generate adjacency/dispersal operator
//                metaSub.performanceTest(true, a_perfRep, a_perfS, floor(a_perfS) / 2);
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
            // Species richness
            MPI_Bcast(&S_p, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&S_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
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

            // delete Armadillo buffers
            delete bMPI_p, bMPI_c, aMPI, rMPI, cMPI, tMPI;
            MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Generate topo - topography/environment (root) //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (rank == 0) {
            metaSub.sppPool.topo.genDomainDecomp(
                    metaSub.sppPool.topo.network); // decompose imported topo
            if (metaSub.sppPool.topo.fVec.n_rows == 0) { // exit if decomposition of imported network fails
                cout << "\nError: domain decomposition failure" << endl;
                longjmp(jump_buffer,1);
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
        ArmadilloMPI* netMPI, *indMPI, *adjMPI, *dspMPI, *envMPI;
        netMPI = new ArmadilloMPI(metaSub.sppPool.topo.no_nodes, 2);
        indMPI = new ArmadilloMPI(metaSub.sppPool.topo.no_nodes, 1);
        adjMPI = new ArmadilloMPI(metaSub.sppPool.topo.no_nodes, metaSub.sppPool.topo.no_nodes);
        dspMPI = new ArmadilloMPI(metaSub.sppPool.topo.no_nodes, metaSub.sppPool.topo.no_nodes);
        if (metaSub.sppPool.topo.envVar > 0) {
            envMPI = new ArmadilloMPI(metaSub.sppPool.topo.envVar, metaSub.sppPool.topo.no_nodes);
        }

        if (rank == 0) { // convert Armadillo matrices to arrays
            netMPI->armaToArrayReg(metaSub.sppPool.topo.network);
            indMPI->armaToArrayReg(metaSub.sppPool.topo.fVec);
            adjMPI->armaToArrayReg(metaSub.sppPool.topo.adjMat);
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
        dspMPI->armaBcast(metaSub.sppPool.dMat_n);
        if (metaSub.sppPool.topo.envVar != 0) {
            envMPI->armaBcast(metaSub.sppPool.topo.envMat);
        }
        metaSub.sppPool.subSet(rank); // each process subsets topo objects according to domain decomposition
        // delete Armadillo buffers
        delete netMPI, envMPI, indMPI, adjMPI, dspMPI;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Simulation /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            printf("\nStarting parallel simulation\n");
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

        // initialise objects
        mat B_p, B_c, BStore_p, BStore_c, U_p, U_c, Ub_p, Ub_c, Tr, R, Cr, Cc, Ar, Ac; // MPI buffers
        int t, it, T = metaSub.tMax/metaSub.deltaT; // number of Schwarz time windows
        int no_invaders_p, no_residents_p, no_extinct_p, no_invaders_c, no_residents_c, no_extinct_c; // counters
        field<uvec> index_ext(2); // index extinct species

        {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Generate invader (random subdomain) ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if (rank == 0) {
                mpi_time1 = MPI_Wtime();
            }

            // MPI buffers producers
            ArmadilloMPI *bMPI_p, *rMPI, *crMPI, *ccMPI, *arMPI, *acMPI, *uMPI_p, *uBcastMPI_p, *iMPI_p;
            // MPI buffers consumers
            ArmadilloMPI *bMPI_c, *uMPI_c, *uBcastMPI_c, *iMPI_c;

            if (rank == 0) { // sample invader and scatter ecology
                no_residents_p = (int) meta.sppPool.bMat_p.n_rows;
                no_residents_c = (int) meta.sppPool.bMat_c.n_rows;
                ivec randSubDom = randi<ivec>(1, distr_param(0, size - 1));
                meta.invaderTest(randSubDom(0)); // search for successful invader
                no_invaders_p = (int) meta.sppPool.bMat_p.n_rows - no_residents_p;
                no_invaders_c = (int) meta.sppPool.bMat_c.n_rows - no_residents_c;
            }
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Bcast(&no_invaders_p, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast number of invaders
            MPI_Bcast(&no_invaders_c, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast number of invaders
            if (no_invaders_p > 0) {
                bMPI_p = new ArmadilloMPI(no_invaders_p, metaSub.sppPool.topo.no_nodes);
                rMPI = new ArmadilloMPI(no_invaders_p, metaSub.sppPool.topo.no_nodes,
                                        metaSub.sppPool.topo.network.n_rows);
                if (metaSub.sppPool.prodComp) {
                    crMPI = new ArmadilloMPI(no_invaders_p, metaSub.sppPool.cMat.n_cols + no_invaders_p);
                    if (metaSub.sppPool.cMat.n_rows != 0) {
                        ccMPI = new ArmadilloMPI(metaSub.sppPool.bMat_p.n_rows, no_invaders_p);
                    }
                }
                if (metaSub.sppPool.pProducer < 1) {
                    arMPI = new ArmadilloMPI(no_invaders_p, metaSub.sppPool.aMat.n_cols + no_invaders_c);
                }
            }

            if (no_invaders_c > 0) {
                bMPI_c = new ArmadilloMPI(no_invaders_c, metaSub.sppPool.topo.no_nodes);
                acMPI = new ArmadilloMPI(metaSub.sppPool.bMat_p.n_rows, no_invaders_c);
            }

            if (rank == 0) {
                if (no_invaders_p > 0) {
                    B_p = meta.sppPool.bMat_p.rows(no_residents_p, meta.sppPool.bMat_p.n_rows - 1);
                    bMPI_p->armaToArrayReg(B_p);
                    if (meta.sppPool.aMat.n_rows > 0) {
                        Ar = meta.sppPool.aMat.rows(no_residents_p, meta.sppPool.aMat.n_rows - 1);
                        arMPI->armaToArrayReg(Ar);
                    }
                    R = meta.sppPool.rMat.rows(no_residents_p, meta.sppPool.rMat.n_rows - 1);
                    rMPI->armaToArrayReg(R);
                    if (metaSub.sppPool.prodComp) {
                        Cr = meta.sppPool.cMat.rows(no_residents_p, meta.sppPool.cMat.n_rows - 1);
                        crMPI->armaToArrayReg(Cr);
                        if (no_residents_p > 0) {
                            Cc = meta.sppPool.cMat.submat(0, no_residents_p, no_residents_p - 1,
                                                          meta.sppPool.cMat.n_cols - 1);
                            ccMPI->armaToArrayReg(Cc);
                        }
                    }
                }

                if (no_invaders_c > 0) {
                    B_c = meta.sppPool.bMat_c.rows(no_residents_c, meta.sppPool.bMat_c.n_rows-1);
                    bMPI_c->armaToArrayReg(B_c);
                    Ac = meta.sppPool.aMat.submat(0,no_residents_c, no_residents_p-1, meta.sppPool.aMat.n_cols-1);
                    acMPI->armaToArrayReg(Ac);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Broadcast invader ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if (no_invaders_c > 0) {
                bMPI_c->armaBcast(B_c); // broadcast B_c
                if (metaSub.sppPool.bMat_c.n_rows == 0) { // add invader rows to biomass matrix consumers
                    metaSub.sppPool.bMat_c = B_c.cols(find(metaSub.sppPool.topo.fVec == rank));
                } else {
                    metaSub.sppPool.bMat_c = join_vert(metaSub.sppPool.bMat_c,
                                                       B_c.cols(find(metaSub.sppPool.topo.fVec == rank)));
                }
                if (metaSub.sppPool.uMat_c.n_rows == 0) { // seed subdomain fixed unknowns
                    metaSub.sppPool.uMat_c = B_c.cols(metaSub.sppPool.topo.adjIF) * metaSub.sppPool.dMat_m;
                } else {
                    metaSub.sppPool.uMat_c = join_vert(metaSub.sppPool.uMat_c,
                                                       B_c.cols(metaSub.sppPool.topo.adjIF) * metaSub.sppPool.dMat_m);
                }
                acMPI->armaBcast(Ac); // col local interaction matrix
                metaSub.sppPool.aMat = join_horiz(metaSub.sppPool.aMat, Ac);
            }

            if (no_invaders_p > 0) {
                bMPI_p->armaBcast(B_p); // broadcast B_p
                if (metaSub.sppPool.bMat_p.n_rows == 0) { // add invader rows to biomass matrix producers
                    metaSub.sppPool.bMat_p = B_p.cols(find(metaSub.sppPool.topo.fVec == rank));
                } else {
                    metaSub.sppPool.bMat_p = join_vert(metaSub.sppPool.bMat_p,
                                                       B_p.cols(find(metaSub.sppPool.topo.fVec == rank)));
                }

                if (metaSub.sppPool.uMat_p.n_rows == 0) { // seed subdomain fixed unknowns
                    metaSub.sppPool.uMat_p = B_p.cols(metaSub.sppPool.topo.adjIF) * metaSub.sppPool.dMat_m;
                } else {
                    metaSub.sppPool.uMat_p = join_vert(metaSub.sppPool.uMat_p,
                                                       B_p.cols(metaSub.sppPool.topo.adjIF) * metaSub.sppPool.dMat_m);
                }

                if (metaSub.sppPool.pProducer < 1) {
                    arMPI->armaBcast(Ar); // broadcast row/
                    if (metaSub.sppPool.aMat.n_rows == 0) { // add invader rows/cols to competitive overlap matrix
                        metaSub.sppPool.aMat = Ar;
                    } else {
                        metaSub.sppPool.aMat = join_vert(metaSub.sppPool.aMat, Ar);
                    }
                }

                rMPI->armaScatter(R); // scatter r_i
                if (metaSub.sppPool.rMat.n_rows == 0) { // add invader rows to growth matrix
                    metaSub.sppPool.rMat = R;
                } else {
                    metaSub.sppPool.rMat = join_vert(metaSub.sppPool.rMat, R);
                }

                if (metaSub.sppPool.prodComp) {
                    crMPI->armaBcast(Cr); // broadcast row/
                    if (metaSub.sppPool.cMat.n_rows != 0) {
                        // segmentation fault occurs here when Release configuration (-O3) selected
                        ccMPI->armaBcast(Cc); // col local interaction matrix
                    }
                    if (metaSub.sppPool.cMat.n_rows == 0) { // add invader rows/cols to competitive overlap matrix
                        metaSub.sppPool.cMat = Cr;
                    } else {
                        metaSub.sppPool.cMat = join_horiz(metaSub.sppPool.cMat, Cc);
                        metaSub.sppPool.cMat = join_vert(metaSub.sppPool.cMat, Cr);
                    }
                }
            }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Domain decomposed numerical solution ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//            if (rank == 0) {
//                mpi_time1 = MPI_Wtime();
//            }
//
//            MPI_Barrier(MPI_COMM_WORLD);
//
//            uMPI_p = new ArmadilloMPI(metaSub.sppPool.bMat_p.n_rows, metaSub.sppPool.topo.no_nodes,
//                                      metaSub.sppPool.topo.network.n_rows);
//            if (metaSub.sppPool.pProducer < 1) {
//                uMPI_c = new ArmadilloMPI(metaSub.sppPool.bMat_c.n_rows, metaSub.sppPool.topo.no_nodes,
//                                          metaSub.sppPool.topo.network.n_rows);
//            }
//
//            MPI_Barrier(MPI_COMM_WORLD);
//            for (t = 0; t < T; t++) {
//                mat BStore_p = metaSub.sppPool.bMat_p; // save inital state for Schwartz iteration updates
//                mat BStore_c;
//                if (metaSub.sppPool.pProducer < 1) {
//                    BStore_c = metaSub.sppPool.bMat_c; // save inital state for Schwartz iteration updates
//                }
//                for (it = 0; it < metaSub.iterS; it++) {
//                    metaSub.sppPool.bMat_p = BStore_p; // restore inital state at entry to Schwartz iteration
//                    if (metaSub.sppPool.pProducer < 1) {
//                        metaSub.sppPool.bMat_c = BStore_c; // restore inital state at entry to Schwartz iteration
//                    }
//                    metaSub.metaCDynamics(metaSub.deltaT); // simulation subdomain scale dynamics
//
//                    // allgather to update fixed unknowns
//                    U_p = (metaSub.sppPool.bMat_p + BStore_p) / 2;
//                    if (metaSub.sppPool.pProducer < 1) {
//                        U_c = (metaSub.sppPool.bMat_c + BStore_c) / 2;
//                    }
//                    uMPI_p->armaToArrayLoc(U_p);
//                    if (metaSub.sppPool.pProducer < 1) {
//                        uMPI_c->armaToArrayLoc(U_c);
//                    }
//                    uMPI_p->armaAllgather(U_p);
//                    if (metaSub.sppPool.pProducer < 1) {
//                        uMPI_c->armaAllgather(U_c);
//                    }
//                    metaSub.sppPool.uMat_p = U_p.cols(metaSub.sppPool.topo.adjIF);
//
//                    // immigration from adjacent subdomains computed once per Schwartz iteration
//                    metaSub.sppPool.uMat_p = metaSub.sppPool.uMat_p * metaSub.sppPool.dMat_m;
//                    if (metaSub.sppPool.pProducer < 1) {
//                        metaSub.sppPool.uMat_c = U_c.cols(metaSub.sppPool.topo.adjIF);
//                        metaSub.sppPool.uMat_c = metaSub.sppPool.uMat_c * metaSub.sppPool.dMat_m;
//                    }
//                    MPI_Barrier(MPI_COMM_WORLD);
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//            }
//            MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Gather result at root process //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//            bMPI_p = new ArmadilloMPI(metaSub.sppPool.bMat_p.n_rows, metaSub.sppPool.topo.no_nodes,
//                                      metaSub.sppPool.topo.network.n_rows);
//            bMPI_p->armaToArrayLoc(metaSub.sppPool.bMat_p);
//            bMPI_p->armaGather(B_p);
//
//            if (metaSub.sppPool.pProducer < 1) {
//                bMPI_c = new ArmadilloMPI(metaSub.sppPool.bMat_c.n_rows, metaSub.sppPool.topo.no_nodes,
//                                          metaSub.sppPool.topo.network.n_rows);
//                bMPI_c->armaToArrayLoc(metaSub.sppPool.bMat_c);
//                bMPI_c->armaGather(B_c);
//            }
//
//            // receive current state at root and generate whole domain biomass matrix
//            if (rank == 0) {
//                meta.sppPool.bMat_p = B_p;
//                if (metaSub.sppPool.pProducer < 1) {
//                    meta.sppPool.bMat_c = B_c;
//                }
//            }

            if (rank == 0) {
                mpi_time2 = MPI_Wtime();
            }

            // deallocate memory assigned to MPI buffers
            delete bMPI_p, rMPI, crMPI, ccMPI, arMPI, acMPI, uMPI_p, uBcastMPI_p, iMPI_p;
            delete bMPI_c, uMPI_c, uBcastMPI_c, iMPI_c;
            MPI_Barrier(MPI_COMM_WORLD);

        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Final book keeping /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (rank == 0) {
            meta.simTime = mpi_time2-mpi_time1;
//             meta.performanceTest(false, a_perfRep);
            meta.writePars(a_perfRep);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        longjmp(jump_buffer,1); // all proceses jump to return
    }
}
