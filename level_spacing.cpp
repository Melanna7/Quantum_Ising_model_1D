/*******************************************************************************
*
* Full diagonalization of the Ising Hamiltonian
*
*******************************************************************************/

// g++ class_Eigen.h level_spacing.cpp -I ./include/ -o levels.out 12550.2s = 3.5 h

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <fstream>
#include <chrono>

// Import the Class hamiltonian
#include "class_Eigen.h"

using namespace std;

/*******************************************************************************
* PARAMETERS OF THE LATTICE
*
* SIDE = size of the lattice's side.
*
* G_FIELD = amplitude of the spin coupling field
*
* H_FIELD = amplitude of the single spin field
*
*******************************************************************************/

#define DIM_HILBERT 2
#define SIDE 12

//--- Main Test ----------------------------------------------------------------

int main(){
    /* Test the methods of the hamiltonian Class. */

    ofstream file;
    string directory, file_name;

    HamiltParameters param;
    param.sparse_flag = 0;
    param.pbc_flag = 0;
    param.num_sites = SIDE;
    param.num_states = pow(DIM_HILBERT, SIDE);
    param.gz_field = -1.;
    param.gy_field = 0.;
    param.gx_field = 0.;
    param.hz_field = 0.;
    param.hy_field = 0.;
    param.hx_field = -1.;
    // Build the computational basis
    param.comp_basis.resize(param.num_states, vector<int>(param.num_sites, 0.));
    for (int n = 0; n < param.num_states; n++) {
        int index = n;
        for (int i = 0; i < param.num_sites; i++) {
            param.comp_basis[n][param.num_sites - i - 1 ] = index % 2;
            index = index / 2;
        }
    }

    hamiltonian HamOp(param);

    // Testing and timing the diagonalization process
    cout << "Start diagonalization..." << endl;
    auto start = chrono::steady_clock::now();
    HamOp.diagonalize();
    auto end = chrono::steady_clock::now();
    // Print elapsed time for diagonalization
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Elapsed time : " << elapsed_seconds.count() << "s" << endl << endl;

    // Define path data directory
    directory = "Data_Ising/";
    // Define name file last configuration of the lattice
    file_name = "levels_hz_" + to_string(param.hz_field) + ".dat";
    // Open file and save eigenvalues
    file.open(directory + file_name);
    for(int i = 0; i < param.num_states; i++){
        file << HamOp.get_eigenvalue(i) << endl;
    }
    file.close();
    cout << "The work is done." << endl << endl;

}
