/*******************************************************************************
*
* Test program for the Class hamiltonian
*
*******************************************************************************/

// g++ ../class_Eigen.h test_Eigen.cpp -o e_eigen

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

//--- Main Test ----------------------------------------------------------------

int main(){
    /* Test the methods of the hamiltonian Class. */

    HamiltParameters param;
    param.sparse_flag = 3;
    param.pbc_flag = 1;
    param.num_sites = 3;
    param.gz_field = 0.;
    param.gy_field = 0.;
    param.gx_field = 0.;
    param.hz_field = 0.5;
    param.hy_field = 1.;
    param.hx_field = 1.;

    double lambda;
    VectorXcd v, w;
    hamiltonian HamOp(param);

    // Test if the basis and the Hamiltonian are correct
    HamOp.show_comput_basis();
    HamOp.show_hamiltonian();

    // Test set_fields
    param.gz_field = -1.;
    param.gx_field = 0.;
    param.hy_field = 0.;
    param.hz_field = -0.3;
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
