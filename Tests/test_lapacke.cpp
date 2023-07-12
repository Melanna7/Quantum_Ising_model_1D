/*******************************************************************************
*
* Test program for the Class hamilt_L
*
*******************************************************************************/

// g++ test_lapacke.cpp -llapacke -I ../include/ -o lapacke.out

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <chrono>

// Import the Class hamilt_L
#include "class_lapacke.h"

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

    vector<complex<double>> state, ground;
    hamilt_L HamOp(param);

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
    ground = HamOp.get_eigenvector(0);
    state = HamOp.action(ground);
    cout << "The state |Psi> = H|GS> is :" << endl;
    for (auto val: state) {
        cout << fixed << setprecision(2) << val << " ";
    }
    cout << endl << endl;

    // Let us prepare a superposition of the two lowest
    // energy state and evaluate the average energy
    ground = HamOp.get_eigenvector(0);
    state = HamOp.get_eigenvector(2);
    for (int i = 0; i < state.size(); i++){
        state[i] = (ground[i] + state[i]) / sqrt(2);
    }
    cout << "The average energy of the state below is ";
    cout << HamOp.average_energy(state) << endl;
    for (auto val: state) {
        cout << fixed << setprecision(2) << val << " ";
    }
    cout << endl << endl;

}
