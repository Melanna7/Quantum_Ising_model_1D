/*******************************************************************************
*
* Test program for the Ising simulation
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <chrono>
#include <cmath>

// Import the Class lattice
#include "../class_sparse_hamilt.h"

using namespace std;

/*******************************************************************************
* PARAMETERS OF THE LATTICE
*
* SIDE = size of the lattice's side.
*
* PBC_FLAG = flag for the boundary conditions.
*            0 for open boundary conditions.
*            else periodic boundary conditions.
*
* G_FIELD = amplitude of the spin coupling field
*
* H_FIELD = amplitude of the longitudinal field
*
* T_FIELD = amplitude of the transverse field.
*
* NUM_STATES = number of eigenvalue to compute.
*
*******************************************************************************/

#define SIDE 3      // Max 8 for LAPACK complete diagonalization
#define PBC_FLAG 1
#define G_FIELD 0.5
#define H_FIELD 0.1
#define T_FIELD 1.

#define NUM_STATES 3

//--- Main Test ----------------------------------------------------------------

int main(){
    /* Test the methods of the hamiltonian Class. */

    HamiltParameters param;
    param.pbc_flag = PBC_FLAG;
    param.num_sites = SIDE;
    param.g_field = G_FIELD;
    param.h_field = H_FIELD;
    param.t_field = T_FIELD;

    vector<complex<double>> state, ground;
    sparse_hamilt HamOp(param);

    HamOp.show_comput_basis();

    // Testing and timing the diagonalization process
    cout << "Start diagonalization..." << endl;
    auto start = chrono::steady_clock::now();
    HamOp.jacobi_davidson(NUM_STATES);
    auto end = chrono::steady_clock::now();

    // Print elapsed time for diagonalization
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Elapsed time : " << elapsed_seconds.count() << "s" << endl << endl;

    HamOp.show_eigen();

    ground = HamOp.get_eigenstate(0);

    state = HamOp.action(ground);
    cout << "The state |Psi> = H|GS> is :" << endl;
    for (auto val: state) {
        cout << fixed << setprecision(2) << val << " ";
    }
    cout << endl << endl;

    // Let us prepare a superposition of the two lowest
    // energy state and evaluate the average energy
    state = HamOp.get_eigenstate(2);
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
