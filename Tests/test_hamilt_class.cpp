/*******************************************************************************
*
* Test program for the Ising simulation
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <chrono>

// Import the Class lattice
#include "../hamiltonian_class.h"

using namespace std;

/*******************************************************************************
* PARAMETERS OF THE LATTICE
*
* SIDE = size of the lattice's side.
*
* EIG_FLAG = flag for the diagonalization procedure.
*            0 for complete LAPACK diagonalization.
*            else for Lancsoz method.................
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

#define SIDE 2      // Max 8 for LAPACK complete diagonalization 
#define EIG_FLAG 0
#define PBC_FLAG 1
#define G_FIELD 1.
#define H_FIELD 0.1
#define T_FIELD 0.33

//--- Main Test ----------------------------------------------------------------

int main(){
    /* Test the methods of the hamiltonian Class. */

    vector<complex<double>> ground_state;
    hamiltonian HamOp(SIDE, PBC_FLAG, G_FIELD, H_FIELD, T_FIELD);

    HamOp.show_comput_basis();

    HamOp.show_hamiltonian();

    HamOp.show_eigenvalues();

    auto start = chrono::steady_clock::now();
    HamOp.diagonalize();
    auto end = chrono::steady_clock::now();

    // Print elapsed time
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Elapsed time : " << elapsed_seconds.count() << "s" << endl << endl;

    HamOp.show_eigenvalues();

    HamOp.show_eigenvectors();

    ground_state = HamOp.get_GS(); // Assicurati sia quello corretto

}
