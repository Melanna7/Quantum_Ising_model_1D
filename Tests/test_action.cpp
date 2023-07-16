/*******************************************************************************
*
* Test program for the action of the hamiltonian
*
*******************************************************************************/

// g++ test_action.cpp -I ../include/ -o action.out

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <chrono>
#include <cmath>

// Import the Class hamiltonian
#include "class_Eigen.h"

using namespace std;

/*******************************************************************************
* PARAMETERS OF THE LATTICE
*
* SIDE = size of the lattice's side.
*
* SPARSE_FLAG = condition on storing or not the Hamiltonian.
*               0 build the Hamiltonian and LAPACK diagonalization.
*               else number of desired eigenvalues with Lancsoz method.
*
* PBC_FLAG = flag for the boundary conditions.
*            0 for open boundary conditions.
*            else periodic boundary conditions.
*
* G_FIELD = amplitude of the spin coupling field
*
* H_FIELD = amplitude of the single spin field
*
*******************************************************************************/

#define DIM_HILBERT 2

#define SIDE 4
#define SPARSE_FLAG 0
#define PBC_FLAG 0

//--- Main Test ----------------------------------------------------------------

int main(){
    /* Test the methods of the hamiltonian Class. */

    HamiltParameters param;
    param.sparse_flag = SPARSE_FLAG;
    param.pbc_flag = PBC_FLAG;
    param.num_sites = SIDE;
    param.num_states = pow(DIM_HILBERT, SIDE);
    param.gz_field = -1.;
    param.gy_field = 0.;
    param.gx_field = 0.;
    param.hz_field = 0.;
    param.hy_field = 0.;
    param.hx_field = 0.;
    // Build the computational basis
    param.comp_basis.resize(param.num_states, vector<int>(param.num_sites, 0.));
    for (int n = 0; n < param.num_states; n++) {
        int index = n;
        for (int i = 0; i < param.num_sites; i++) {
            param.comp_basis[n][param.num_sites - i - 1 ] = index % 2;
            index = index / 2;
        }
    }

    double lambda;
    VectorXcd v(param.num_states), w(param.num_states), z(param.num_states);
    hamiltonian HamOp(param);
    HamOp.diagonalize();

    v[0] = 1.;
    w[param.num_states - 1] = 1.;

    cout << "The average energy of the state below is ";
    cout << HamOp.average_energy(v) << endl;
    cout << v.transpose() << endl << endl;
    cout << "It is normalized to " << (v.transpose() * v).real();
    cout << endl << endl;

    cout << "The average energy of the state below is ";
    cout << HamOp.average_energy(w) << endl;
    cout << w.transpose() << endl << endl;
    cout << "It is normalized to " << (w.transpose() * w).real();
    cout << endl << endl;

    z = (v + w) / sqrt(2);
    cout << "The average energy of the state below is ";
    cout << HamOp.average_energy(z) << endl;
    cout << z.transpose() << endl << endl;
    cout << "It is normalized to " << (z.transpose() * z).real();
    cout << endl << endl;

    z = (v - w) / sqrt(2);
    cout << "The average energy of the state below is ";
    cout << HamOp.average_energy(z) << endl;
    cout << z.transpose() << endl << endl;
    cout << "It is normalized to " << (z.transpose() * z).real();
    cout << endl << endl;

    z = HamOp.get_eigenvector(0);
    cout << "The average energy of the state below is ";
    cout << HamOp.average_energy(z) << endl;
    cout << z.transpose() << endl << endl;
    cout << "It is normalized to " << (z.transpose() * z).real();
    cout << endl << endl;
}
