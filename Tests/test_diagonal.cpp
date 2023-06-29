/*******************************************************************************
*
* Test program for confronting eigenvalues from different methods
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <chrono>
#include <cmath>

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

    Parameters_L para_L;
    para_L.sparse_flag = param.sparse_flag;
    para_L.pbc_flag = param.pbc_flag;
    para_L.num_sites = param.num_sites;
    para_L.gz_field = param.gz_field;
    para_L.gy_field = param.gy_field;
    para_L.gx_field = param.gx_field;
    para_L.hz_field = param.hz_field;
    para_L.hy_field = param.hy_field;
    para_L.hx_field = param.hx_field;

    vector<complex<double>> state, ground;
    complex<double> lambda;
    VectorXcd v, w;

    hamilt_L Ham_L(para_L);
    hamiltonian Ham_E(param);
    param.sparse_flag = 8;
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
