/*******************************************************************************
*
* Test program for the Ising simulation
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <chrono>
#include <cmath>

#include "../class_hamiltonian.h"
#include "../magnetization.h"

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

#define SIDE 2         // Max 8 for LAPACK complete diagonalization
#define SPARSE_FLAG 0
#define PBC_FLAG 1
#define G_FIELD 0.
#define H_FIELD 0.
#define T_FIELD 1.

//--- Main Test ----------------------------------------------------------------

int main(){
    /* Test the methods of the hamiltonian Class. */

    HamiltParameters param;
    param.sparse_flag = SPARSE_FLAG;
    param.pbc_flag = PBC_FLAG;
    param.num_sites = SIDE;
    param.g_field = G_FIELD;
    param.h_field = H_FIELD;
    param.t_field = T_FIELD;

    vector<vector<double>> magZX;
    vector<complex<double>> state, ground;
    hamiltonian HamOp(param);

    HamOp.show_hamiltonian();

    ground = HamOp.compute_GS();

    magZX = magnetization(ground);

    cout << "Average of sigma_Z " << endl;
    for (int i = 0; i < param.num_sites; i++) {
        cout << magZX[0][i] << " ";
    }
    cout << endl;
    cout << "Average of sigma_X " << endl;
    for (int i = 0; i < param.num_sites; i++) {
        cout << magZX[1][i] << " ";
    }
    cout << endl;

    // Let us prepare a superposition of the two lowest
    // energy state and evaluate the average energy
    ground = HamOp.get_eigenstate(0);
    state = HamOp.get_eigenstate(2);
    for (int i = 0; i < state.size(); i++){
        state[i] = (ground[i] + state[i]) / sqrt(2);
    }

    magZX = magnetization(state);

    cout << "Average of sigma_Z " << endl;
    for (int i = 0; i < param.num_sites; i++) {
        cout << magZX[0][i] << " ";
    }
    cout << endl;
    cout << "Average of sigma_X " << endl;
    for (int i = 0; i < param.num_sites; i++) {
        cout << magZX[1][i] << " ";
    }
    cout << endl;


}
