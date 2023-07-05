/*******************************************************************************
*
* Test program for the Class hamiltonian
*
*******************************************************************************/

// g++ ../class_Eigen.h test_Eigen.cpp -I ../include/ -o e_eigen.out

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <chrono>
#include <cmath>

// Import the Class hamiltonian
#include "../class_Eigen.h"

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

#define SIDE 2
#define SPARSE_FLAG 3
#define PBC_FLAG 1

//--- Main Test ----------------------------------------------------------------

int main(){
    /* Test the methods of the hamiltonian Class. */

    HamiltParameters param;
    param.sparse_flag = SPARSE_FLAG;
    param.pbc_flag = PBC_FLAG;
    param.num_sites = SIDE;
    param.num_states = pow(DIM_HILBERT, SIDE);
    param.gz_field = 1.;
    param.gy_field = 1.;
    param.gx_field = 0.5;
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
    VectorXcd v, w;
    hamiltonian HamOp(param);

    // Test if the basis and the Hamiltonian are correct
    HamOp.show_comput_basis();
    HamOp.show_hamiltonian();
    param.hz_field = 1.;
    HamOp.set_fields(param, 1);
    HamOp.show_hamiltonian();

    // Testing and timing the diagonalization process
    cout << "Start diagonalization..." << endl;
    auto start = chrono::steady_clock::now();
    HamOp.diagonalize();
    auto end = chrono::steady_clock::now();
    // Print elapsed time for diagonalization
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Elapsed time : " << elapsed_seconds.count() << "s" << endl << endl;
    // Show eigen stuff and test show procedures
    cout << "The eigenvalues and eigenstates of H are:" << endl << endl;
    HamOp.show_eigenvalues();
    HamOp.show_eigenvectors();
    HamOp.show_eigen();

    // Test the Hamiltonian action on the ground state vector
    lambda = HamOp.get_eigenvalue(0);
    v = HamOp.get_eigenvector(0);
    w = HamOp.action(v);
    cout << "Consider first eigenvalue, lambda = " << lambda << endl;
    cout << "If v is the corresponding eigenvector" << endl << v.transpose();
    cout << endl << "then lambda * v = " << endl << lambda * v.transpose();
    cout << endl << "... action(v) = " << endl << w.transpose() << endl;
    // More test on diagonalization
    if (param.sparse_flag) {
        cout << "... and H * v = " << endl;
        cout << (HamOp.spars_hamilt * v).transpose() << endl << endl;
    } else {
        cout << "... and H * v = " << endl;
        cout << (HamOp.dense_hamilt * v).transpose() << endl << endl;

        cout << "Finally |V * D * V^(-1)  - H | = ";
        cout << norm((HamOp.eigenvectors
                * HamOp.eigenvalues.asDiagonal()
                * HamOp.eigenvectors.inverse()
                - HamOp.dense_hamilt).sum()) << endl << endl;
    }

    // Test average_energy
    // Let us prepare a superposition of the two lowest
    // energy state and evaluate the average energy
    v = HamOp.get_eigenvector(0);
    w = HamOp.get_eigenvector(2);
    v = (v + w) / sqrt(2);
    cout << "The average energy of the state below is ";
    cout << HamOp.average_energy(v) << endl;
    cout << v.transpose() << endl << endl;
    cout << "It is normalized to " << v.transpose() * v << endl << endl;
}
