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

/*
 * NEED TO UPDATE THE COMPUTION OF U*D_m IN PARALLEL ALGORITHM TO ACCOMODATE THE DISTRIBUTION IN EMIGRATION RATES
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

#include "Metacommunity.h"
#include "LVMCM_rng.h"

using namespace std;
using namespace arma;
using std::vector;

// Assembly flags
bool ASSEMBLE = true; // run assembly algorithm or jump directly to generation of analysis object
bool SNAPSHOT = false; // store snapshot biomass and growth rate matrices
bool OUTPUT = true; // select write to file
bool WARMING = false; // select warming experiment
bool LONGDISTDISP = false; // select introduction of long distance dispersal
bool CONSAREA = false; // select conservation area experiment
bool TRAJECTORY = false; // select generate and write trajectory object to file
bool FLUCTUATE = false; // select generate and write to file abiotic fluctuation
bool CMAT_REG = false; // select generate and write to file regional competitive overlap matrix
int FIX_SEED = 0; // if 0 , all seeds random, if 1 landscape seed fixed, if 2 all seeds fixed
bool NODE_REMOVAL = false; // select node removal experiment
string DD_SOL_PATH = {}; // path to location for dd solution output

// Storage for setjmp/longjmp
jmp_buf jump_buffer;
jmp_buf jump_buffer_continue_assembly;

// Global Metacommunity object and timing variables
Metacommunity meta;
time_t time1, time2;
double mpi_time1, mpi_time2;
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
    const int nCols;
    const int nRows;
    const int subCols;
    int totalLen=0;
    vector<int> disp;

    // Constructor
    ArmadilloMPI(int nRows, int nCols, int subCols=0, bool vec=false):
        nRows(nRows),
        nCols(nCols),
        subCols(subCols) {
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

        if (totalLen != 0) {
            if (A.n_rows != 0) {
                bufferLVSnd.resize(A.n_rows);
            }
            bufferVecRcv.resize(totalLen);
            for (int i=0; i<A.n_rows; i++) {
                bufferLVSnd[i] = A(i);
            }
            int* count = &bufferLenRcv[0];
            MPI_Allgatherv(&(bufferLVSnd.at(0)), A.n_rows, MPI_INT, &(bufferVecRcv.at(0)), count, &(disp.at(0)), MPI_INT, MPI_COMM_WORLD);
        }

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
                mpi_time2 = MPI_Wtime();
                meta.simTime += mpi_time2 - mpi_time1; // record wall clock
                mpi_time1 = mpi_time2;
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
/////////////////////////////////////// Store program arguments in variables ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Default Metacommunity parameters
    int a_init = 1;
    string a_bMat = {};
    string a_xMat = {};
    string a_scMat = {};
    int a_invMax = 0;
    int a_iterS = 2;
    int a_deltaT = 100;
    int a_tMax = 500;
    int a_perfRep = 0;
    int a_perfS = 0;
    double a_autoFluct = 0; // if nonzero, gives trajectory recording time
    bool a_harvest = false;
    string a_outputDirectory;
    bool a_continue_assembly = false;
    bool a_write_continue_assembly = false;
    double a_invasionSize = 0.05;
    int a_g_max = 0;
    int a_nodeRemoval = 0;

    // Default Species parameters
    double a_c1 = 0.5;
    double a_c2 = 0.5;
    double a_emRate = 0.1;
    double a_dispL = 0.1;
    double a_pProducer = 1.0;
    bool a_prodComp = true;
    bool a_symComp = false;
    double a_alpha = 0;
    double a_sigma = 0;
    double a_sigma_t = 0.1;
    double a_rho = 0;
    bool a_discr_c_ij = true;
    double a_omega = 0.5;
    int a_dispNorm = 0;

    // Default Topograpy parameters
    int a_no_nodes = 4;
    double a_phi = 1;
    int a_envVar = 0;
    vec a_skVec = {0.0};
    double a_var_e = 0.01;
    bool a_randGraph = true;
    bool a_gabriel = true;
    int a_bisec = 0;
    double a_T_int = -1;
    double a_dTdt = 0.0;
    int a_edges = 1;
    int a_cons_percent = 60;
    double a_cons_phi = 1e-4;
    double a_cons_rep = 10;

    // Default output variables
    double a_parOut = 0;
    string a_experiment = "DEFAULT";
    int a_rep = 0;
    string a_jobID = "NA";

    for (int i = 1; i<argc; i++) { // loop through program arguments an allocate to parameters
        char var1 = argv[i][1];
        char var2 = argv[i][2];
        if (!isalpha(var2)) {
            var2 = '0';
        }
        i++;

        switch (var1) {
            // input parameters
            case 'a' : // set alpha - base attack rate (double)
                a_alpha = atof(argv[i]);
                break;

            case 'b' : // bFile: path to data for importing initialized model (string)
                a_bMat = argv[i];
                a_init = 0;
                break;

            case 'c' : // set c1, c2 - competition parameters (2x double)
                a_c1 = atof(argv[i]);
                i++;
                a_c2 = atof(argv[i]);
                break;

            case 'd' :
                switch (var2) {
                    case '0' : // set emRate, dispL - emigration rate and dispersal length (2x double)
                        a_emRate = atof(argv[i]);
                        i++;
                        a_dispL = atof(argv[i]);
                        break;
                    case 'd' : // set jobID - job ID used for automatic checkpointing (string)
                        DD_SOL_PATH = argv[i];
                        break;
                    case 'n' : // set normalization of dispersal model
                        a_dispNorm = atoi(argv[i]);
                        break;
 
                }
                break;

            case 'e' : // set envVar - number of explicitly modelled environmental variables (int)
                a_envVar = atoi(argv[i]);
                break;

            case 'f' : // set outputDirectory - location for write to file (string)
                a_outputDirectory = argv[i];
                break;

            case 'g' : // select maximum gamma diversity
                       // CAUTION: if regional limit is less than requested diversity, model will not converge!
                a_g_max = atoi(argv[i]);
                a_invasionSize = 0; // This will ensure gamma is not exceeded but if gamma set high, will slow assembly
                break;

            case 'i' :
                switch (var2) {
                    case '0' : // set invMax - total number of invasions (int)
                        a_invMax = atoi(argv[i]);
                        break;
                    case 'd' : // set jobID - job ID used for automatic checkpointing (string)
                        a_jobID = argv[i];
                        break;
                    case 's' : // set proportion of extra new invaders in each iteration (double)
                        a_invasionSize = atof(argv[i]);
                        break;
                }
                break;

            case 'n' : // set N - number of nodes (int)
                a_no_nodes = atoi(argv[i]);
                break;

            case 'o' : // set parOut, experiment, rep - key parameter, experiment name, replicate number for output filenames (double, string, int)
                a_parOut = atof(argv[i]);
                i++;
                a_experiment = argv[i];
                i++;
                a_rep = atoi(argv[i]);
                break;

            case 'p' :
                switch (var2) {
                    case '0' : // set phi - spatial autocorrelation length of the environment (double)
                        a_phi = atof(argv[i]);
                        break;
                    case 'p' : // set pProducer - probabilty of sampling a producer, bipartite models (double)
                        a_pProducer = atof(argv[i]);
                        break;
                }
                break;

            case 'r' : // set rho - consumer respiration rate (double)
                a_rho = atof(argv[i]);
                break;

            case 's' :
                switch(var2) {
                    case '0' : // set sigma - trophic link distribution parameter (double)
                        a_sigma = atof(argv[i]);
                        break;
                    case 'c' : // scFile: path to data for scaling of local interaction matrix
                        a_scMat = argv[i];
                        break;
                    case 'i' : // set iterS, deltaT - Schwartz iteration and time window (2x int)
                        a_iterS = atoi(argv[i]);
                        i++;
                        a_deltaT = atoi(argv[i]);
                        break;
                    case 'k' : // sk: environmental sensitivity shape parameter
                        a_envVar = atoi(argv[i]);
                        a_skVec.set_size(a_envVar);
                        i++;
                        for (int s=0; s<a_skVec.n_rows; s++) {
                            a_skVec(s) =  atof(argv[i]);
                            if (s < a_skVec.n_rows-1) {
                                i++;
                            }
                        }
                        break;
                    case 't' : // set sigma_t - standard deviation of random environmental fluctuation (double)
                        a_sigma_t = atof(argv[i]);
                        break;
                }
                break;

            case 't' : // set tMax - relaxation time (int)
                switch(var2) {
                    case '0' :
                        a_tMax = atoi(argv[i]);
                        break;
                    case 'n' : // set omega - temperature niche width in units 1/sqrt(N) (double)
                        a_omega = atof(argv[i]);
                        a_T_int = 1.0;
                        a_envVar = 1;
                        break;
                }
                break;

            case 'v' : // set var_e - variance of base environmental distribution (double)
                a_var_e = atof(argv[i]);
                break;

            case 'x' : // xFile: path to data for importing spatial network
                a_xMat = argv[i];
                break;

            // Switches
            case 'C' : // set prodComp - select producer coupling on/off (bool)
                if (!strcmp(argv[i],"F")) {
                    a_prodComp = false;
                }
                break;

            case 'D' : // set discr_c_ij - select discrete distribution in competive overlap coefficients (bool)
                if (!strcmp(argv[i],"F")) {
                    a_discr_c_ij = false;
                }
                break;

            case 'G' : // set gabriel - select Gabriel/complete graph (bool)
                if (!strcmp(argv[i],"F")) {
                    a_gabriel = false;
                }
                break;

            case 'O' : // set OUTPUT - select write to file (bool)
                if (!strcmp(argv[i],"F")) {
                    OUTPUT = false;
                }
                break;

            case 'R' : // set randGraph - select random spatial network/lattice (bool)
                if (!strcmp(argv[i],"F")) {
                    a_randGraph = false;
                }
                break;

            case 'S' : // set SNAPSHOT - select regular write to file (bool)
                switch(var2) {
                    case '0' :
                        if (!strcmp(argv[i],"T")) {
                            SNAPSHOT = true;
                        }
                        break;
                    case 'C' : // select symmetric competition model
                        if (!strcmp(argv[i],"T")) {
                            a_symComp = true;
                        }
                        break;
                }
                break;


            case 'Z' : // set FIX_SEED: 0 - random seed generated randomly; 1 - network/environment fixed; 2 - all sampling fixed
                FIX_SEED = atoi(argv[i]);
                break;

            // Perturbation experiments/analysis objects

            case 'F' : // set tMax, cons_percent, cons_phi, cons_rep - select conservation area/fragmentation experiment and set
                       // relaxation time, percent landscape conserved, correlation length of binary conservation area, and replicate number
                       // (int, int, double, int)
                a_tMax = atoi(argv[i]);
                i++;
                a_cons_percent = atoi(argv[i]);
                i++;
                a_cons_phi = atof(argv[i]);
                i++;
                if (atoi(argv[i]) > 0) {
                    a_continue_assembly = true;
                    a_invMax = atoi(argv[i]);
                }
                i++;
                a_cons_rep = atoi(argv[i]);
                CONSAREA = true;
                ASSEMBLE = false;
                break;

            case 'H' : // set harvest - estimate regional scale interaction coefficients using harvesting experiment (bool)
                if (!strcmp(argv[i],"T")) {
                    a_harvest = true;
                    ASSEMBLE = false;
                    CMAT_REG = true;
                }
                break;

            case 'K' : // set nodeRemoval
                a_nodeRemoval = atoi(argv[i]);
                ASSEMBLE = false;
                NODE_REMOVAL = true;
                break;

            case 'L' : // set tMax, edges - select long distance dispersal experment and set
                       // relaxation time and number of edges per perturbation event (2x int)
                a_tMax = atoi(argv[i]);
                i++;
                a_edges = atoi(argv[i]);
                LONGDISTDISP = true;
                ASSEMBLE = false;
                break;

            case 'T' : // set autoFluct - store trajectory for studying autonomous fluctuations (bool)
                a_autoFluct = atof(argv[i]);
                break;

            case 'W' : // set dTdt - select temperature warming experment and set rate of warming (double)
                a_dTdt = atof(argv[i]);
                WARMING = true;
                ASSEMBLE = false;
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

        // infer number of bisections from size of MPI communicator
        a_bisec = (int) log2(size);

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
            if (a_emRate < 0) { // check abs(emRate) = 1, required for non-uniform emRate;
                if (a_emRate != -1.0) {
                    cout << "\nError: for non-uniform emigration rate model, emRate = -1.0 required" << endl;
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
///////////////////////////////////////////// Set random seeds 1  //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (rank == 0) {
            if (FIX_SEED > 0) {
                g_seed = 1;
            } else {
                // generate random seed
                std::random_device rd;
                g_seed = rd();
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if (size > 1) {
            MPI_Bcast(&g_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }

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
                        a_xMat,
                        a_scMat,
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
                        a_symComp,
                        a_alpha,
                        a_sigma,
                        a_sigma_t,
                        a_rho,
                        a_discr_c_ij,
                        a_omega,
                        a_dispNorm,
                        // Topograpy parameters
                        a_no_nodes,
                        a_phi,
                        a_envVar,
                        a_skVec,
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Broadcast imported metacommunity objects   //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if (!a_init) { // broadcast imported metacommunity model objects and parameters to non-root processes
                // Is there a risk in loading data simultaneously to multiple processes? If not don't bother with MPI parameterization of subdomains
                int S_p, S_c, prodComp, envVar;
                if (rank == 0) {
                    metaSub = Metacommunity(
                            // Metacommunity parameters
                            a_init,
                            a_bMat,
                            a_xMat,
                            a_scMat,
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
                            a_symComp,
                            a_alpha,
                            a_sigma,
                            a_sigma_t,
                            a_rho,
                            a_discr_c_ij,
                            a_omega,
                            a_dispNorm,
                            // Topograpy parameters
                            a_no_nodes,
                            a_phi,
                            a_envVar,
                            a_skVec,
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
                    S_p = metaSub.sppPool.rMat.n_rows;
                    S_c = metaSub.sppPool.bMat_p.n_rows - S_p;
                    prodComp = metaSub.sppPool.prodComp;
                    envVar = metaSub.sppPool.topo.envVar;
                }
                MPI_Barrier(MPI_COMM_WORLD);

                if (size > 1) { // THIS NEEDS UPDATING
                    // parameterize subdomains
                    // Metacommunity parameters
                    MPI_Bcast(&metaSub.invMax, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.iterS, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.deltaT, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.tMax, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    // Species parameters
                    MPI_Bcast(&metaSub.sppPool.pProducer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.prodComp, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.c1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.c2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.rho, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.sigma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.discr_c_ij, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    // Topography parameters
                    MPI_Bcast(&metaSub.sppPool.topo.no_nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.topo.bisec, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.topo.envVar, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.topo.var_e, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.topo.phi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&metaSub.sppPool.topo.T_int, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    // Species richness
                    MPI_Bcast(&S_p, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&S_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    metaSub.sppPool.S_p = S_p;
                    metaSub.sppPool.S_c = S_c;
                    MPI_Barrier(MPI_COMM_WORLD);

                    // initialize Armadillo buffers
                    ArmadilloMPI bMPI_p(S_p, metaSub.sppPool.topo.no_nodes); // producer biomass MPI
                    ArmadilloMPI rMPI(S_p, metaSub.sppPool.topo.no_nodes); // producer growth rate MPI
                    ArmadilloMPI sMPI(S_p, metaSub.sppPool.topo.no_nodes); // producer growth rate MPI (w/o temp grad)
		    ArmadilloMPI cMPI(S_p, S_p); // producer competition MPI
		    ArmadilloMPI tMPI(metaSub.sppPool.topo.envVar, S_p); // producer environmental tolerance MPI

                    if (rank == 0) { // convert Armadillo matrices to arrays
                        bMPI_p.armaToArrayReg(metaSub.sppPool.bMat_p);
                        rMPI.armaToArrayReg(metaSub.sppPool.rMat);
                        rMPI.armaToArrayReg(metaSub.sppPool.sMat);
                        if (prodComp) {
                            cMPI.armaToArrayReg(metaSub.sppPool.cMat);
                        }
                        if (metaSub.sppPool.topo.envVar != 0) {
                            tMPI.armaToArrayReg(metaSub.sppPool.tMat);
                        }
                    }
                    MPI_Barrier(MPI_COMM_WORLD);

                    // broadcast model matrices and store in local Metacommunity objects
                    bMPI_p.armaBcast(metaSub.sppPool.bMat_p);
                    rMPI.armaBcast(metaSub.sppPool.rMat);
                    sMPI.armaBcast(metaSub.sppPool.sMat);
                    if (metaSub.sppPool.prodComp) {
                        cMPI.armaBcast(metaSub.sppPool.cMat);
                    }
                    if (metaSub.sppPool.topo.envVar != 0) {
                        tMPI.armaBcast(metaSub.sppPool.tMat);
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        } // subdomain initialization scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Generate topo - topography/environment (root) ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        { // domain decomposition scope
            if (rank == 0) {
                if (a_init) {
                    metaSub.sppPool.topo.genDomainDecomp(); // generate and decompose landscape

                    if (meta.sppPool.topo.envVar > 0) {
                        meta.sppPool.topo.genEnvironment(); // generate environmental distribution
                    }

                    if (metaSub.sppPool.topo.T_int != -1.0) {
                        metaSub.sppPool.topo.genTempGrad();
                    }
                } else {
                    metaSub.sppPool.topo.genDomainDecomp(
                        metaSub.sppPool.topo.network); // decompose imported topo
                    if (size > 1) {
                        if (metaSub.sppPool.topo.fVec.n_rows == 0) { // exit if decomposition of imported network fails
                            cout << "\nError: domain decomposition failure" << endl;
                            longjmp(jump_buffer,1);
                        }
                    }
                }

                metaSub.sppPool.genDispMat();

                unique(metaSub.sppPool.topo.fVec.t()).print("\nf.unique");
                if (ASSEMBLE) {
                    meta = metaSub; // store copy of complete domain at root process for outputting and extinction testing
                }
            }

            meta.rank = rank; // store rank for signal handler
            MPI_Barrier(MPI_COMM_WORLD);

            if (size > 1) {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Broadcast landcape to subdomains //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                // initialize Armadillo buffers
                ArmadilloMPI netMPI(metaSub.sppPool.topo.no_nodes, 2); // spatial network MPI
                ArmadilloMPI indMPI(metaSub.sppPool.topo.no_nodes, 1); // indicator vector MPI
                ArmadilloMPI adjMPI(metaSub.sppPool.topo.no_nodes,
				    metaSub.sppPool.topo.no_nodes); // adjacency MPI
                ArmadilloMPI dstMPI(metaSub.sppPool.topo.no_nodes,
				    metaSub.sppPool.topo.no_nodes); // distance mat MPI
                ArmadilloMPI seVeMPI(metaSub.sppPool.topo.no_nodes,
				     metaSub.sppPool.topo.no_nodes); // eigen vec cov mat MPI
                ArmadilloMPI seVaMPI(metaSub.sppPool.topo.no_nodes,
				     metaSub.sppPool.topo.no_nodes); // eigen val cov mat MPI
                ArmadilloMPI dspMPI(metaSub.sppPool.topo.no_nodes,
				    metaSub.sppPool.topo.no_nodes); // dispersal mat MPI
		ArmadilloMPI envMPI(metaSub.sppPool.topo.envVar,
				    metaSub.sppPool.topo.no_nodes); // environment MPI

                if (rank == 0) { // convert Armadillo matrices to arrays
                    netMPI.armaToArrayReg(metaSub.sppPool.topo.network);
                    indMPI.armaToArrayReg(metaSub.sppPool.topo.fVec);
                    adjMPI.armaToArrayReg(metaSub.sppPool.topo.adjMat);
                    dstMPI.armaToArrayReg(metaSub.sppPool.topo.distMat);
                    seVeMPI.armaToArrayReg(metaSub.sppPool.topo.sigEVec);
                    seVaMPI.armaToArrayReg(metaSub.sppPool.topo.sigEVal);
                    dspMPI.armaToArrayReg(metaSub.sppPool.dMat_n);
                    if (metaSub.sppPool.topo.envVar != 0) {
                        envMPI.armaToArrayReg(metaSub.sppPool.topo.envMat);
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);

                // broadcast model matrices and store in local Metacommunity objects
                netMPI.armaBcast(metaSub.sppPool.topo.network);
                indMPI.armaBcast(metaSub.sppPool.topo.fVec);
                adjMPI.armaBcast(metaSub.sppPool.topo.adjMat);
                dstMPI.armaBcast(metaSub.sppPool.topo.distMat);
                seVeMPI.armaBcast(metaSub.sppPool.topo.sigEVec);
                seVaMPI.armaBcast(metaSub.sppPool.topo.sigEVal);
                dspMPI.armaBcast(metaSub.sppPool.dMat_n);
                if (metaSub.sppPool.topo.envVar != 0) {
                    envMPI.armaBcast(metaSub.sppPool.topo.envMat);
                }
                metaSub.sppPool.subSet(rank); // each process subsets topo objects according to domain decomposition

                // gather at root to check objects properly ordered - required for MPI gather calls
                mat N;
		ArmadilloMPI netMPI2(metaSub.sppPool.topo.network.n_rows, 2 * size, 2);
                netMPI2.armaToArrayLoc(metaSub.sppPool.topo.network);
                netMPI2.armaGather(N);

                mat E;
                if (metaSub.sppPool.topo.envVar > 0) {
                    ArmadilloMPI eMPI(metaSub.sppPool.topo.envMat.n_rows, size * metaSub.sppPool.topo.network.n_rows,
                                            metaSub.sppPool.topo.network.n_rows);
                    eMPI.armaToArrayLoc(metaSub.sppPool.topo.envMat);
                    eMPI.armaGather(E);
                    MPI_Barrier(MPI_COMM_WORLD);
                }

                if (rank == 0) { // overwrite meta with node order after domain decomposition algorithm
                    mat netTemp = N.submat(0, 0, N.n_rows - 1, 1);
                    for (int s = 1; s < size; s++) {
                        netTemp = join_vert(netTemp, N.submat(0, 2 * s, N.n_rows - 1, 2 * s + 1));
                    }
                    N = netTemp;
                    meta.sppPool.topo.network = N;
                    meta.sppPool.topo.envMat = E;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }

        } // domain decomposition scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Simulation /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            printf("\nStarting parallel assembly\n");
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
///////////////////////////////////////////// Set random seeds 2  //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (FIX_SEED == 1) {
        if (rank == 0) {
            // generate random seed for invaders
            std::random_device rd;
            g_seed = rd();
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if (size > 1) {
            MPI_Bcast(&g_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        for (int s = 0; s < size; s++) {
            if (rank == s) {
                cout << "Random seed (species only) seen by process " << rank << " " << g_seed << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        LVMCM_rng::boost_rng.seed(g_seed);
        arma_rng::set_seed(g_seed);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Initialize MPI buffers //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        setjmp(jump_buffer_continue_assembly);

        if (ASSEMBLE) {

            { // assembly scope

                // initialise objects
                mat B_p, BStore_p, BStore_c, U_p, Ub_p, Tr, R, S, Cr, Cc; // MPI buffers
                mat B_dd; // for storage of domain decomposed solution if required
                int t, it, T = metaSub.tMax/metaSub.deltaT; // number of Schwarz time windows
                int no_invaders_p, no_residents_p, no_extinct_p, no_invaders_c, no_residents_c, no_extinct_c; // counters

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
                        int S_p0 = metaSub.sppPool.S_p; // record diversity prior to invader testing
                        int S_c0 = metaSub.sppPool.S_c;

                        // select number of invaders and sample their trophic levels
                        int no_trophLev[2] = {0}; // no producers, consumers to invade
                        int no_invaders = a_invasionSize * (metaSub.sppPool.S_p + metaSub.sppPool.S_c) + 1;

                        if (metaSub.sppPool.pProducer < 1) { // bipartite
                            // random uniform variables generated for allocating trophic level
                            vec trophLev(no_invaders);
                            int seedProd = 5;
                            if (metaSub.sppPool.bMat_p.n_rows > seedProd) { // seed metacommunity with at least seedProd producers
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

                                    if (size > 1) {
                                        // gather indices of successful invaders
                                        ArmadilloMPI pgMPI(1, size, 0, true);
                                        pgMPI.armaAllgatherv(posGrowth);
                                        posGrowth = unique(posGrowth);
                                        MPI_Barrier(MPI_COMM_WORLD);
                                    }

                                    // select desired number of invaders with positive growth rates
                                    suc_inv += posGrowth.n_rows;
                                    posGrowth.resize(min((int) posGrowth.n_rows, no_trophLev[tL]));

                                    // remove unsucessful/excess invaders
                                    mat bInv_max;
                                    bInv_max = metaSub.invaderCleanup(tL, posGrowth);

                                    if (posGrowth.n_rows == 0) { // in this case all processes restart invader testing
                                        metaSub.invaderPopulate(tL, bInv_max);
                                        continue;
                                    }

                                    if (size > 1) {
                                        // gather maximum subdomain biomass for selection of port node
                                        ArmadilloMPI bimMPI(bInv_max.n_rows, size, 1);
                                        bimMPI.armaToArrayLoc(bInv_max);
                                        bimMPI.armaAllgather(bInv_max);
                                        MPI_Barrier(MPI_COMM_WORLD);
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

                        if (size > 1) {
                            if (metaSub.sppPool.S_p > S_p0) { // collect R/S/T at root
                                ArmadilloMPI riMPI(metaSub.sppPool.S_p - S_p0,
                                                   metaSub.sppPool.topo.no_nodes,
                                                   metaSub.sppPool.topo.network.n_rows);
                                R = metaSub.sppPool.rMat.rows(S_p0, metaSub.sppPool.rMat.n_rows-1);

                                // replace with gather?
                                riMPI.armaToArrayLoc(R);
                                riMPI.armaAllgather(R);

                                if (metaSub.sppPool.topo.T_int > -1) {
                                    ArmadilloMPI siMPI(metaSub.sppPool.S_p - S_p0,
						                metaSub.sppPool.topo.no_nodes,
						                metaSub.sppPool.topo.network.n_rows);
                                    S = metaSub.sppPool.sMat.rows(S_p0, metaSub.sppPool.sMat.n_rows-1);

                                    siMPI.armaToArrayLoc(S);
                                    siMPI.armaAllgather(S);
                                }
                                MPI_Barrier(MPI_COMM_WORLD);

                                if (rank == 0) {
                                    if (meta.sppPool.rMat.n_rows == 0) {
                                        meta.sppPool.rMat = R;
                                    } else {
                                        meta.sppPool.rMat = join_vert(meta.sppPool.rMat, R);
                                    }
                                    if (S.n_rows > 0) {
                                        if (meta.sppPool.sMat.n_rows == 0) {
                                            meta.sppPool.sMat = S;
                                        } else {
                                            meta.sppPool.sMat = join_vert(meta.sppPool.sMat, S);
                                        }
                                    }
                                    meta.sppPool.cMat = metaSub.sppPool.cMat;
                                    if (meta.sppPool.topo.envVar > 0) {
                                        meta.sppPool.tMat = metaSub.sppPool.tMat;
                                    }
                                }
                                MPI_Barrier(MPI_COMM_WORLD);
                            }
                        }

                        if (rank == 0) { // update cMat for both producer and consumer invasions
                            meta.sppPool.cMat = metaSub.sppPool.cMat; // straight copy since C is not spatially decomposed
                        }
                        MPI_Barrier(MPI_COMM_WORLD);

                        metaSub.sppPool.invasion += no_invaders;
                    } // invader testing scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Domain decomposed numerical solution ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    { // dynamics scope

                        if (size > 1) {
                            // intialize MPI objects for gathering fixed unknowns at boundary
                            ArmadilloMPI uMPI_p(metaSub.sppPool.bMat_p.n_rows, metaSub.sppPool.topo.no_nodes,
                                metaSub.sppPool.topo.network.n_rows);

                            if ((metaSub.sppPool.invasion == metaSub.invMax) && (DD_SOL_PATH.length() > 0)) {
                                if (rank==0) {
                                    B_dd.set_size(0,metaSub.sppPool.bMat_p.n_cols*metaSub.sppPool.bMat_p.n_rows);
                                }
                            }

			                MPI_Barrier(MPI_COMM_WORLD);
                            for (t = 0; t < T; t++) {
                                mat BStore_p = metaSub.sppPool.bMat_p; // save inital state for Schwartz iteration updates
                                MPI_Barrier(MPI_COMM_WORLD);
                                
                                for (it = 0; it < metaSub.iterS; it++) {
                                    metaSub.sppPool.bMat_p = BStore_p; // restore inital state at entry to Schwartz iteration
                                    metaSub.metaCDynamics(metaSub.deltaT); // simulation subdomain scale dynamics
                                    // allgather to update fixed unknowns
#if 1  // set to 1 for old version using mid-points
                                    U_p = (metaSub.sppPool.bMat_p + BStore_p) / 2; // store MID-POINT
#else
                                    U_p = metaSub.sppPool.bavMat_p; // store AVERAGE
#endif
                                    uMPI_p.armaToArrayLoc(U_p);
                                    uMPI_p.armaAllgather(U_p);
                                    MPI_Barrier(MPI_COMM_WORLD);
                                    metaSub.sppPool.uMat_p = U_p.cols(metaSub.sppPool.topo.adjIF);
                                    MPI_Barrier(MPI_COMM_WORLD);

                                    // fixed unknowns multiplied by dispersal matrix
                                    sp_mat dMat_m_sp(metaSub.sppPool.dMat_m); // cast as sparse matrix
                                    metaSub.sppPool.uMat_p = metaSub.sppPool.uMat_p * dMat_m_sp;
                                    MPI_Barrier(MPI_COMM_WORLD);
                                    if ((metaSub.sppPool.invasion == metaSub.invMax) && (DD_SOL_PATH.length() > 0)) {
                                        if (rank == 0) { // record solution at end of Schwartz iteration for testing
                                            rowvec bVec = vectorise(metaSub.sppPool.bMat_p, 1);
                                            B_dd.resize(B_dd.n_rows + 1, B_dd.n_cols);
                                            B_dd.row(B_dd.n_rows - 1) = bVec;;
                                        }
                                    }
                                }
                                MPI_Barrier(MPI_COMM_WORLD);
                            }
                            MPI_Barrier(MPI_COMM_WORLD);

                            if ((metaSub.sppPool.invasion == metaSub.invMax) && (DD_SOL_PATH.length() > 0)) {
                                // save dd solution for testing
                                if (rank == 0) {
                                    B_dd.save(DD_SOL_PATH, raw_ascii);
                                }
                            }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Gather result at root process //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                            ArmadilloMPI bMPI_p(metaSub.sppPool.bMat_p.n_rows, metaSub.sppPool.topo.no_nodes,
						metaSub.sppPool.topo.network.n_rows);
                            bMPI_p.armaToArrayLoc(metaSub.sppPool.bMat_p);
                            bMPI_p.armaGather(B_p);

                            // receive current state at root and generate whole domain biomass matrix (required for extinction)
                            if (rank == 0) {
                                meta.sppPool.bMat_p = B_p;
                                meta.sppPool.S_p = metaSub.sppPool.S_p;
                                meta.sppPool.S_c = metaSub.sppPool.S_c;
                            }
                        } else {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Undecomposed numerical solution /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            metaSub.metaCDynamics(metaSub.tMax); // simulate metacommunty dynamics
                        }
                    } // dynamics scope

                    { // extinction scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Remove extinct species /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        if (size > 1) {
                            field<uvec> index_ext(2); // extinct species index
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

                            ArmadilloMPI iMPI_p(no_extinct_p, 1); // set size buffer for extinction index
                            ArmadilloMPI iMPI_c(no_extinct_c, 1); // set size buffer for extinction index

                            // broadcast indices of extinct species
                            if (rank == 0) {
                                iMPI_p.armaToArrayReg(index_ext(0));
                                if (no_extinct_c > 0) {
                                    iMPI_c.armaToArrayReg(index_ext(1)); // set size buffer for extinction index
                                }
                            }
                            MPI_Barrier(MPI_COMM_WORLD);
                            // broadcast index of extinct species
                            if (no_extinct_p > 0) {
                                iMPI_p.armaBcast(index_ext(0));
                            }
                            if (no_extinct_c > 0) {
                                iMPI_c.armaBcast(index_ext(1));
                            }
                            if (no_extinct_p + no_extinct_c > 0) { // remove regionally extict species at subdomain scale
                                metaSub.sppPool.extinct(0, index_ext(0), index_ext(1));
                            }

                        } else {
                            metaSub.sppPool.extinct(); // remove extinct species

                        }

                        metaSub.sppPool.S_p = metaSub.sppPool.rMat.n_rows;
                        metaSub.sppPool.S_c = metaSub.sppPool.bMat_p.n_rows - metaSub.sppPool.rMat.n_rows;
                    } // extinction scope

                    MPI_Barrier(MPI_COMM_WORLD);

                    if (rank == 0) {
                        mat presAbs;
                        presAbs.zeros(metaSub.sppPool.bMat_p.n_rows, metaSub.sppPool.bMat_p.n_cols);
                        double thresh = 1e-4;
                        presAbs.elem(find(metaSub.sppPool.bMat_p > thresh)).ones();
                        rowvec alpha = sum(presAbs,0);

                        printf("\rInvasions / S_p / S_c = %d / %d / %d         ",
                               metaSub.sppPool.invasion, (int) metaSub.sppPool.rMat.n_rows, (int) metaSub.sppPool.bMat_p.n_rows - (int) metaSub.sppPool.rMat.n_rows);
                        fflush(stdout);

                        if (a_g_max > 0) { // select regional number of species
                            if (metaSub.sppPool.bMat_p.n_rows == a_g_max) {
                                sim = 0; // switch simulation off
                            }
                        } else if (metaSub.sppPool.invasion >= metaSub.invMax) {
                            sim = 0; // switch simulation off
                        }
                    }
                    MPI_Barrier(MPI_COMM_WORLD);

//                    { // machinery for attaching serial debugger to MPI program: https://www.open-mpi.org/faq/?category=debugging
//                        if (rank == 0) {
//                            if (metaSub.sppPool.invasion == 28) {
//                                volatile int i = 0;
//                                char hostname[256];
//                                gethostname(hostname, sizeof(hostname));
//                                printf("\nPID %d on %s ready for attach\n", getpid(), hostname);
//                                fflush(stdout);
//                                while (0 == i) {
//                                    sleep(5);
//                                    fflush(stdout);
//                                }
//                            }
//                        }
//                        MPI_Barrier(MPI_COMM_WORLD);
//                    }

                    if (size > 1) {
                        // broadcast switch from root to end simulation if invMax reached
                        MPI_Bcast(&sim, 1, MPI_INT, 0, MPI_COMM_WORLD);
                        MPI_Barrier(MPI_COMM_WORLD);
                    }

                    if (SNAPSHOT) {
                        if (size > 1) {
//                            meta.snapShot();
                            meta.outputData();
                        } else {
//                            metaSub.snapShot();
                            metaSub.outputData();
                        }
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
                if (metaSub.sppPool.topo.consArea_bin.n_rows > 0) {
                    ASSEMBLE = false; // write to file handled by perturbation clause in case of continued assemblies
                    a_write_continue_assembly = true;
                }
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

                if (CMAT_REG) {
                    printf("\n\nGenerating regional scale interaction matrix...");
                    metaSub.genCMatReg();
                    metaSub.outputData();
                }

                if (a_autoFluct) {
                    printf("\n\nGenerating static environment trajectory...");
                    metaSub.storeTraj = 1; // store trajectories in (NxS)xtRelax object
                    metaSub.metaCDynamics(a_autoFluct);
                }

                if (FLUCTUATE) {
                    printf(" Generating dynamic environment trajectory...");
                    metaSub.storeTraj = 2; // concatenate static/dynamic environment trajectories
                    int tRelax = 200; // relaxtion time for trajectories object
                    for (int t = 0; t < tRelax; t++) { // simulate temporally fluctuating environment
                        metaSub.envFluct();
                    }
                }

                if (WARMING) {
                    // Begin with static environment step to demonstrate degree of autonomous fluctuations
                    printf("\n\nGenerating static environment trajectory...");
                    metaSub.storeTraj = 1; // store trajectories in (NxS)xtRelax object
                    int tRelax = 100; // relaxtion time for trajectories object
//                    metaSub.metaCDynamics(tRelax);

                    // Dial up intercept of temperature gradient, driving temperature optima to the right
                    printf("\n\nSimulating regional warming...");
                    metaSub.storeTraj = 1; // concatenate static/dynamic environment trajectories
                    metaSub.sppPool.topo.T_int = 1.0; // required for importing models prior to adding T_int to output command
                    tRelax = 10000; // relaxtion time for warming experiment object
                    int res = 100;
                    for (int t = 0; t < (tRelax / res); t++) { // simulate temporally fluctuating environment
                        cout << "\nYear " << t << endl;
                        metaSub.warming(a_dTdt, res, t);
                    }
                }

                if (LONGDISTDISP) {
                    printf("Randomly allocating long distance spatial coupling...\n");
                    metaSub.longDistDisp(a_tMax, a_edges);
                }

                if (NODE_REMOVAL) {
//                    metaSub.metaCDynamics(5000); // final relaxation
                    metaSub.outputData();
                    printf("Removing nodes and simulating dynamics...\n");
                    metaSub.nodeRemoval(a_tMax, a_nodeRemoval);
                }
            }
        }
    MPI_Barrier(MPI_COMM_WORLD);
    longjmp(jump_buffer, 1); // jump to return
}

// Local Variables:
// c-file-style: "stroustrup"
// End:
