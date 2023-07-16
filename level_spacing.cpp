/*******************************************************************************
*
* Full diagonalization of the Ising Hamiltonian
*
*******************************************************************************/

// g++ level_spacing.cpp -I ./include/ -o levels.out

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

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
#define SIDE_MAX 12
#define SIDE_MIN 12
#define HZ -1.

//--- Subroutine ---------------------------------------------------------------

void run_diagonalization(HamiltParameters param){
    /* Diagonalize the Hamiltonian and store the eigenvalues into a file */

    ofstream file;
    string directory, file_name;
    hamiltonian HamOp(param);

    // Define path data directory
    directory = "Data_Ising/";
    // Define name file last configuration of the lattice
    file_name = "levels_" + to_string(param.num_sites);
    file_name += "_hz_" + to_string(abs(param.hz_field)) + ".dat";

    // Timing the diagonalization process
    cout << "Start diagonalization..." << endl << file_name << endl << endl;
    auto start = chrono::steady_clock::now();
    HamOp.diagonalize();
    auto end = chrono::steady_clock::now();
    // Print elapsed time for diagonalization
    chrono::duration<double> elapsed_seconds = end - start;
    cout << file_name << endl << "Elapsed time : ";
    cout << elapsed_seconds.count() / 60 << " m" << endl <<  endl;

    // Open file and save eigenvalues
    file.open(directory + file_name);
    for(int i = 0; i < param.num_states; i++){
        file << HamOp.get_eigenvalue(i) << endl;
    }
    file.close();
}

//--- Main Test ----------------------------------------------------------------

int main(){
    /* Full diagonalization for level spacing stats */

    vector<thread> threadPool;

    HamiltParameters param;
    param.sparse_flag = 0;
    param.pbc_flag = 0;
    param.gz_field = -1.;
    param.gy_field = 0.;
    param.gx_field = 0.;
    param.hy_field = 0.;
    param.hx_field = -1.;

    for(int side = SIDE_MAX; side >= SIDE_MIN ; side--){

        // Set parameters
        param.num_sites = side;
        param.num_states = pow(DIM_HILBERT, side);
        // Build the computational basis
        param.comp_basis.resize(param.num_states, vector<int>(param.num_sites, 0.));
        for (int n = 0; n < param.num_states; n++) {
            int index = n;
            for (int i = 0; i < param.num_sites; i++) {
                param.comp_basis[n][param.num_sites - i - 1 ] = index % 2;
                index = index / 2;
            }
        }

        param.hz_field = 0.;
        cout << "Run for side = " + to_string(side);
        cout << " , hz = " + to_string(param.hz_field) + " ..." << endl;
        threadPool.emplace_back([param]() {run_diagonalization(param);});

        // param.hz_field = HZ;
        // cout << "Run for side = " + to_string(side);
        // cout << " , hz = " + to_string(param.hz_field) + " ..." << endl;
        // threadPool.emplace_back([param]() {run_diagonalization(param);});
    }

    for (auto& thr : threadPool) thr.join();

    cout << "The work is done." << endl << endl;
    return 0;
}
