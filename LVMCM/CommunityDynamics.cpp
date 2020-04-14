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
 * This class contains the members and methods required interfacing with the SUNDIALS numerical ODE solver.
 */

#include "CommunityDynamics.h"

CommunityDynamics::CommunityDynamics(const CommunityDynamics & C) {

    bMat_p = C.bMat_p;
    bMat_c = C.bMat_c;
    emMat_p = C.emMat_p;
    emMat_c = C.emMat_c;
    rMat = C.rMat;
    efMat = C.efMat;
    uMat_p = C.uMat_p;
    uMat_c = C.uMat_c;
    cMat = C.cMat;
    aMat = C.aMat;
    dMat = C.dMat;
    dMat_n = C.dMat_n;
    dMat_m = C.dMat_m;
    rho = C.rho;

}

CommunityDynamics & CommunityDynamics::operator=(const CommunityDynamics &C) {

    if (this == &C) {
        return *this;
    } else {
        delete bMat_p;
        delete bMat_c;
        delete emMat_p;
        delete emMat_c;
        delete rMat;
        delete efMat;
        delete uMat_p;
        delete uMat_c;
        delete cMat;
        delete aMat;
        delete dMat;
        delete dMat_n;
        delete dMat_m;
        delete rho;

        bMat_p = C.bMat_p;
        bMat_c = C.bMat_c;
        emMat_p = C.emMat_p;
        emMat_c = C.emMat_c;
        rMat = C.rMat;
        efMat = C.efMat;
        uMat_p = C.uMat_p;
        uMat_c = C.uMat_c;
        cMat = C.cMat;
        aMat = C.aMat;
        dMat = C.dMat;
        dMat_n = C.dMat_n;
        dMat_m = C.dMat_m;
        rho = C.rho;

        return *this;
    }
}

int CommunityDynamics::number_of_variables() const {
    // define number of state variables
    return (bMat_p->n_rows*bMat_p->n_cols + bMat_c->n_rows*bMat_c->n_cols);
}


void CommunityDynamics::read_state_from(const ODE_vector & state) {
    // converts 1D vector to N dimensional biomassMat
    int k = 0;
    for (int i=0; i<bMat_p->n_cols; i++) {
        for (int j=0; j<bMat_p->n_rows; j++) {
            (*bMat_p)(j,i) = state[k];
            k++;
        }
    }

    for (int i=0; i<bMat_c->n_cols; i++) {
        for (int j=0; j<bMat_c->n_rows; j++) {
            (*bMat_c)(j,i) = state[k];
            k++;
        }
    }
}

void CommunityDynamics::write_state_to(ODE_vector & state) const {
    // converts N dimensional biomassMat to 1D vector
    // note B_p and B_c are concatenated into single ODE_vector state
    int k = 0;
    for (int i=0; i<bMat_p->n_cols; i++) {
        for (int j=0; j<bMat_p->n_rows; j++) {
            state[k] = (*bMat_p)(j,i);
            k++;
        }
    }

    for (int i=0; i<bMat_c->n_cols; i++) {
        for (int j=0; j<bMat_c->n_rows; j++) {
            state[k] = (*bMat_c)(j,i);
            k++;
        }
    }
}

int CommunityDynamics::dynamics(ODE_vector const & state, ODE_vector & time_derivative) {

    // summary:
        // numerically approximate (serial/parallel) the solution to the LVMCM (competitive/bipartite)
        // dynamics equations are constructed in a modular way for flexibility
        // interfaces with SUNDIALs ODE solver

    // required members:
        // members of CommunityDynamics object are pointers to matrix objects stored in a object of class Species

    // arguments:
        // state - ODE_vector representing the current state (biomass) of the metacommunity
        // time_derivative - ODE_vector representing the dynamics of the system constructed in a modular way by successive matrix operations

    // output:
        // updates to sppPool matrices - no return value since memory is shared

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Initialize temporary matrix objects ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    mat B_p(&state[0], bMat_p->n_rows, bMat_p->n_cols);
    mat B_c(&state[bMat_p->n_rows*bMat_p->n_cols], bMat_c->n_rows, bMat_c->n_cols);

    mat dBdt_p(&time_derivative[0], bMat_p->n_rows, bMat_p->n_cols, false, true);
    mat dBdt_c(&time_derivative[bMat_p->n_rows*bMat_p->n_cols], bMat_c->n_rows, bMat_c->n_cols, false, true);

    // growth and dispersal objects used to construct model dynamics in _modular_ way
    mat prodGrowth(B_p.n_rows, B_p.n_cols);
    mat consGrowth(B_c.n_rows, B_c.n_cols);
    mat prodDisp;
    mat consDisp;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Modular construction of model objects ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // producer internal growth rate
    prodGrowth = *rMat; // note R = ones in case var_e = 0 (uniform topo)

    // producer density dependence
    if (cMat) { // intra-/inter-specific
        prodGrowth -= (*cMat) * B_p;
    } else {
        prodGrowth -= B_p; // intra-specific only
    }

    // trophic interactions
    if (aMat) {
        prodGrowth -= (*aMat) * B_c;
        consGrowth = *rho * (aMat->t() * B_p - 1);
    }

    if (dMat) {
        // non-parallel dispersal operator
        sp_mat dMat_sp(*dMat); // cast as sparse matrix
        if (bMat_p) {
            if (!emMat_p) {
                prodDisp.set_size(B_p.n_rows, B_p.n_cols);
//            prodDisp = B_p * (*dMat);
                prodDisp = B_p * dMat_sp;
            } else {
                mat emMat_p_N(B_p.n_rows, B_p.n_cols);
                for (int i=0; i<B_p.n_rows; i++) {
                    emMat_p_N.row(i).fill(emMat_p->at(i));
                }
//            prodDisp = (emMat_p_N % B_p) * (*dMat);
                prodDisp = (emMat_p_N % B_p) * dMat_sp;
            }
        }
        if (bMat_c->n_rows != 0) {
            consDisp.set_size(B_c.n_rows, B_c.n_cols);
            consDisp = B_c * (*dMat);
            consDisp = B_c * dMat_sp;
        }
    }

    dBdt_p = B_p % prodGrowth;
    dBdt_c = B_c % consGrowth;

    if (prodDisp.n_rows > 0) {
        dBdt_p += prodDisp;
    }

    if (consDisp.n_rows > 0) {
        dBdt_c += consDisp;
    }

    // set negative elements to absolute - accounts for small numerical errors but produces large 'oscillations' in pathological case
    uvec negativeB_p = find(B_p < 0);
    dBdt_p(negativeB_p) = - B_p(negativeB_p);
    uvec negativeB_c = find(B_c < 0);
    dBdt_c(negativeB_c) = - B_c(negativeB_c);

    return 0;
}