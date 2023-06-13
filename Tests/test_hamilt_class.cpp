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
#include "../class_hamiltonian.h"

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
* H_FIELD = amplitude of the longitudinal field
*
* T_FIELD = amplitude of the transverse field
*
*******************************************************************************/

#define SIDE 3         // Max 8 for LAPACK complete diagonalization
#define SPARSE_FLAG 0
#define PBC_FLAG 0
#define GZ_FIELD -1.
#define HZ_FIELD 0.
#define HX_FIELD 1.

//--- Main Test ----------------------------------------------------------------

int main(){
    /* Test the methods of the hamiltonian Class. */

    HamiltParameters param;
    param.sparse_flag = SPARSE_FLAG;
    param.pbc_flag = PBC_FLAG;
    param.num_sites = SIDE;
    param.gz_field = GZ_FIELD;
    param.hz_field = HZ_FIELD;
    param.hx_field = HX_FIELD;

    vector<complex<double>> state, ground;
    hamiltonian HamOp(param);

    HamOp.show_comput_basis();

    HamOp.show_hamiltonian();

    // HamOp.set_gz_field(-1.);
    // HamOp.show_hamiltonian();
    //
    // HamOp.set_hz_field(0.);
    // HamOp.show_hamiltonian();
    //
    // HamOp.set_hx_field(-1.);
    // HamOp.show_hamiltonian();

    // Testing and timing the diagonalization process
    cout << "Start diagonalization..." << endl;
    auto start = chrono::steady_clock::now();
    HamOp.diagonalize();
    auto end = chrono::steady_clock::now();

    // Print elapsed time for diagonalization
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Elapsed time : " << elapsed_seconds.count() << "s" << endl << endl;

    // HamOp.show_eigenvalues();
    // HamOp.show_eigenvectors();
    HamOp.show_eigen();

    ground = HamOp.compute_GS();

    state = HamOp.action(ground);
    cout << "The state |Psi> = H|GS> is :" << endl;
    for (auto val: state) {
        cout << fixed << setprecision(2) << val << " ";
    }
    cout << endl << endl;

    // Let us prepare a superposition of the two lowest
    // energy state and evaluate the average energy
    ground = HamOp.get_eigenstate(0);
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
