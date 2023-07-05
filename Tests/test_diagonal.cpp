/*******************************************************************************
*
* Test program for confronting eigenvalues from different methods
*
*******************************************************************************/

// g++ ../class_Eigen.h ../class_lapacke.h test_diagonal.cpp -llapacke -I ../include/ -o e_differences.out

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <chrono>

// Import the Class hamiltonian
#include "../class_Eigen.h"
// Import the Class hamilt_L
#include "../class_lapacke.h"

using namespace std;

#define DIM_HILBERT 2
#define SIDE 8
#define GET_VECTORS 0

//--- Main Test ----------------------------------------------------------------

int main(){
    HamiltParameters param;
    param.sparse_flag = 0;
    param.pbc_flag = 1;
    param.num_sites = SIDE;
    param.num_states = pow(DIM_HILBERT, SIDE);
    param.gz_field = 1.;
    param.gy_field = 1.;
    param.gx_field = 0.5;
    param.hz_field = 0.6;
    param.hy_field = 0.;
    param.hx_field = 0.1;
    // Build the computational basis
    param.comp_basis.resize(param.num_states, vector<int>(param.num_sites, 0.));
    for (int n = 0; n < param.num_states; n++) {
        int index = n;
        for (int i = 0; i < param.num_sites; i++) {
            param.comp_basis[n][param.num_sites - i - 1 ] = index % 2;
            index = index / 2;
        }
    }

    vector<complex<double>> state;
    VectorXcd v, w;

    hamilt_L Ham_LD(param);       // Class Lapacke Dense
    hamiltonian Ham_ED(param);    // Class Eigen Dense
    param.sparse_flag = param.num_states;
    hamilt_L Ham_LS(param);       // Class Lapacke Sparse
    hamiltonian Ham_ES(param);    // Class Eigen Sparse

    // Testing and timing the diagonalization process
    cout << "Start diagonalization..." << endl;

    auto start = chrono::steady_clock::now();
    Ham_LD.diagonalize();
    auto end = chrono::steady_clock::now();
    // Print elapsed time for diagonalization
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Elapsed time LD: " << elapsed_seconds.count() << "s" << endl;

    start = chrono::steady_clock::now();
    Ham_ED.diagonalize();
    end = chrono::steady_clock::now();
    // Print elapsed time for diagonalization
    elapsed_seconds = end - start;
    cout << "Elapsed time ED: " << elapsed_seconds.count() << "s" << endl;

    start = chrono::steady_clock::now();
    Ham_LS.diagonalize();
    end = chrono::steady_clock::now();
    // Print elapsed time for diagonalization
    elapsed_seconds = end - start;
    cout << "Elapsed time LS: " << elapsed_seconds.count() << "s" << endl;

    start = chrono::steady_clock::now();
    Ham_ES.diagonalize();
    end = chrono::steady_clock::now();
    // Print elapsed time for diagonalization
    elapsed_seconds = end - start;
    cout << "Elapsed time ES: " << elapsed_seconds.count() << "s" << endl;
    cout << "End diagonalization." << endl << endl;

    // Test difference between eigenvalues
    cout << "Difference between lapack, eigen and sparse eigenvalues:" << endl;
    for (int i = 0; i < param.sparse_flag; i++) {
        cout << i << ") " << Ham_LD.get_eigenvalue(i) << " | "
             << Ham_ED.get_eigenvalue(i) << " | "
             << Ham_ES.get_eigenvalue(i) << endl;

        if(GET_VECTORS){
            state = Ham_LD.get_eigenvector(i);
            cout << "lapck = ";
            for (auto val : state) cout << val.real() << " ";
            cout << endl;
            v = Ham_ED.get_eigenvector(i);
            cout << "dense = " << v.transpose().real() << endl;
            w = Ham_ES.get_eigenvector(i);
            cout << "spars = " << w.transpose().real() << endl;
            v = (v - w) / 2 ;
            cout << "diff = " << v.transpose() * v << endl << endl;
        }
    }

}
