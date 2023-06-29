/*******************************************************************************
*
* Test program for confronting eigenvalues from different methods
*
*******************************************************************************/

// g++ ../class_Eigen.h ../class_lapacke.h test_diagonal.cpp -llapacke -o e_diff

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <chrono>

// Import the Class hamiltonian
#include "../class_Eigen.h"
// Import the Class hamilt_L
#include "../class_lapacke.h"

using namespace std;

//--- Main Test ----------------------------------------------------------------

int main(){
    HamiltParameters param;
    param.sparse_flag = 0;
    param.pbc_flag = 1;
    param.num_sites = 3;
    param.gz_field = -3.2;
    param.gy_field = -0.2;
    param.gx_field = -1.;
    param.hz_field = 0.7;
    param.hy_field = 0.;
    param.hx_field = 0.3;

    vector<complex<double>> state;
    VectorXcd v, w;

    hamilt_L Ham_L(param);
    hamiltonian Ham_E(param);
    param.sparse_flag = pow(2, param.num_sites);
    hamiltonian Ham_S(param);

    // Testing and timing the diagonalization process
    cout << "Start diagonalization..." << endl;
    auto start = chrono::steady_clock::now();
    Ham_E.diagonalize();
    Ham_S.diagonalize();
    Ham_L.diagonalize();
    auto end = chrono::steady_clock::now();
    // Print elapsed time for diagonalization
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Elapsed time : " << elapsed_seconds.count() << "s" << endl << endl;

    // Test difference between eigenvalues
    cout << "Difference between lapack, dense and sparse eigenvalues:" << endl;
    for (int i = 0; i < param.sparse_flag; i++) {
        cout << i << ") " << Ham_L.get_eigenvalue(i) << " | "
             << Ham_E.get_eigenvalue(i) << " | "
             << Ham_S.get_eigenvalue(i) << endl;
        state = Ham_L.get_eigenstate(i);
        cout << "lapck = ";
        for (auto val : state) cout << val.real() << " ";
        cout << endl;
        v = Ham_E.get_eigenvector(i);
        cout << "dense = " << v.transpose().real() << endl;
        w = Ham_S.get_eigenvector(i);
        cout << "spars = " << w.transpose().real() << endl;
        v = (v - w) / 2 ;
        cout << "diff = " << v.transpose() * v << endl << endl;
    }

}
