/*******************************************************************************
*
* Hamiltonian class definition
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

// include Eigen modules
//#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
// include Lanczos module
#include <lambda_lanczos/lambda_lanczos.hpp>

using namespace std;
using namespace Eigen;
using lambda_lanczos::LambdaLanczos;

typedef Triplet<complex<double>> Tcd;

//--- Contents -----------------------------------------------------------------

#ifndef PARAM_CLASS_H
#define PARAM_CLASS_H

struct HamiltParameters {
    int pbc_flag=0, sparse_flag=0;
    int num_sites=2, num_states=4;
    double gz_field=0., gx_field=0., gy_field=0.;
    double hz_field=0., hx_field=0., hy_field=0.;
    vector<vector<int>> comp_basis;
};

#endif

#ifndef HAMILT_CLASS_H
#define HAMILT_CLASS_H

class hamiltonian {
    /* Hamiltonian operator class */
    /* Class constructor ******************************************

    Inputs:

    LENGHT = total number of sites.

    SPARSE_FLAG = condition on storing or not the Hamiltonian.
                  0 build the complete Hamiltonian.
                  else number of desired eigenvalues.

    PBC_FLAG = flag for the boundary conditions.
               0 open boundary conditions.
               else periodic boundary conditions.

    G_FIELD = amplitude of the spin coupling field

    H_FIELD = amplitude of the single spin field

    *************************************************************/
private:
    int spr_flag_, pbc_flag_;
    int tot_length_, tot_states_;
    double gx_field, gy_field, gz_field;
    double hx_field, hy_field, hz_field;
    MatrixXi basis;

public:
    VectorXd eigenvalues;
    MatrixXcd dense_hamilt, eigenvectors;
    SparseMatrix<complex<double>> spars_hamilt;

    hamiltonian(const HamiltParameters& param):
        tot_states_(param.num_states),
        tot_length_(param.num_sites),
        spr_flag_(param.sparse_flag),
        pbc_flag_(param.pbc_flag),
        gx_field(param.gx_field),
        gy_field(param.gy_field),
        gz_field(param.gz_field),
        hx_field(param.hx_field),
        hy_field(param.hy_field),
        hz_field(param.hz_field)
    {// CONSTRUCTION BEGIN

        // Build the computational basis
        basis.resize(tot_states_, tot_length_);
        for (int n = 0; n < tot_states_; n++) {
            for (int i = 0; i < tot_length_; i++) {
                basis(n, i) = param.comp_basis[n][i];
            }
        }

        if (spr_flag_) {
            // Initialize partial eigen-stuff
            eigenvalues = VectorXd::Zero(spr_flag_);
            eigenvectors= MatrixXcd::Zero(tot_states_, spr_flag_);
            // Build sparse Hamiltonian
            buildSparse();
        } else {
            // Initialize complete eigen-stuff
            eigenvalues = VectorXd::Zero(tot_states_);
            eigenvectors = MatrixXcd::Zero(tot_states_, tot_states_);
            // Build dense Hamiltonian
            buildHamilt();
        }

    }// CONSTRUCTION END

    vector<Tcd> build_coefficients() {
        /* Generate a list of the Hamiltonian coefficients */

        int index = 0;
        complex<double> value;
        vector<Tcd> tripletList;

        // Sigma_z [longitudinal field]
        if (hz_field != 0.) for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == 1)
                    tripletList.push_back(Tcd(n, n, -hz_field));
                if (basis(n, i) == 0)
                    tripletList.push_back(Tcd(n, n, hz_field));
            }
        }
        // Sigma_x [transverse field]
        if (hx_field != 0.) for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == 1)
                    index = n - pow(2, tot_length_ - 1 - i);
                if (basis(n, i) == 0)
                    index = n + pow(2, tot_length_ - 1 - i);
                tripletList.push_back(Tcd(index, n, hx_field));
            }
        }
        // Sigma_y [transverse field]
        if (hy_field != 0.) for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == 1) {
                    index = n - pow(2, tot_length_ - 1 - i);
                    value = - hy_field * complex(0., 1.);
                }
                if (basis(n, i) == 0) {
                    index = n + pow(2, tot_length_ - 1 - i);
                    value = hy_field * complex(0., 1.);
                }
                tripletList.push_back(Tcd(index, n, value));
            }
        }
        // Sigma_z Sigma_z [coupling]
        if (gz_field != 0.) for (int i = 0; i < tot_length_ - 1; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == basis(n, i+1))
                    tripletList.push_back(Tcd(n, n, gz_field));
                if (basis(n, i) != basis(n, i+1))
                    tripletList.push_back(Tcd(n, n, -gz_field));
            }
        }
        // Sigma_z Sigma_z [boundary conditions]
        if (pbc_flag_ & (gz_field != 0.)) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, tot_length_ - 1) == basis(n, 0))
                    tripletList.push_back(Tcd(n, n, gz_field));
                if (basis(n, tot_length_ - 1) != basis(n, 0))
                    tripletList.push_back(Tcd(n, n, -gz_field));
            }
        }
        // Sigma_x Sigma_x [coupling]
        if (gx_field != 0.) for (int i = 0; i < tot_length_ - 1; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == 1)
                    index = n - pow(2, tot_length_ - 1 - i);
                if (basis(n, i) == 0)
                    index = n + pow(2, tot_length_ - 1 - i);
                if (basis(n, i+1) == 1)
                    index -= pow(2, tot_length_ - 2 - i);
                if (basis(n, i+1) == 0)
                    index += pow(2, tot_length_ - 2 - i);
                tripletList.push_back(Tcd(index, n, gx_field));
            }
        }
        // Sigma_x Sigma_x [boundary conditions]
        if (pbc_flag_ & (gx_field != 0.)) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, 0) == 1)
                    index = n - pow(2, tot_length_ - 1);
                if (basis(n, 0) == 0)
                    index = n + pow(2, tot_length_ - 1);
                if (basis(n, tot_length_ - 1) == 1)
                    index -= 1;
                if (basis(n, tot_length_ - 1) == 0)
                    index += 1;
                tripletList.push_back(Tcd(index, n, gx_field));
            }
        }
        // Sigma_y Sigma_y [coupling]
        if (gy_field != 0.) for (int i = 0; i < tot_length_ - 1; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == 1)
                    index = n - pow(2, tot_length_ - 1 - i);
                if (basis(n, i) == 0)
                    index = n + pow(2, tot_length_ - 1 - i);
                if (basis(n, i+1) == 1)
                    index -= pow(2, tot_length_ - 2 - i);
                if (basis(n, i+1) == 0)
                    index += pow(2, tot_length_ - 2 - i);

                if (basis(n, i) == basis(n, i+1))
                    tripletList.push_back(Tcd(index, n, -gy_field));
                if (basis(n, i) != basis(n, i+1))
                    tripletList.push_back(Tcd(index, n, gy_field));
            }
        }
        // Sigma_y Sigma_y [boundary conditions]
        if (pbc_flag_ & (gy_field != 0.)) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, 0) == 1)
                    index = n - pow(2, tot_length_ - 1);
                if (basis(n, 0) == 0)
                    index = n + pow(2, tot_length_ - 1);
                if (basis(n, tot_length_ - 1) == 1)
                    index -= 1;
                if (basis(n, tot_length_ - 1) == 0)
                    index += 1;

                if (basis(n, tot_length_ - 1) == basis(n, 0))
                    tripletList.push_back(Tcd(index, n, -gy_field));
                if (basis(n, tot_length_ - 1) != basis(n, 0))
                    tripletList.push_back(Tcd(index, n, gy_field));
            }
        }

        return tripletList;
    }

    void buildHamilt() {
        /* Construct the dense Hamiltonian */

        // Get the list of coefficients
        vector<Tcd> tripletList = build_coefficients();
        // Initialize the Hamiltonian
        dense_hamilt = MatrixXcd::Zero(tot_states_, tot_states_);
        // Fill the Hamiltonian with the coefficients
        for (const auto& tripl : tripletList)
            dense_hamilt.coeffRef(tripl.row(),tripl.col()) += tripl.value();

        // Check that the Hamiltonian is hermitian
        complex<double> rest;
        rest = (dense_hamilt - dense_hamilt.transpose().adjoint()).sum();
        if (rest != 0.) {
          cout << "Error occurred building the matrix: ";
          cout << "it is not hermitian!" << endl;
          show_hamiltonian();
          exit(1);
        }

        cout << "Building complete!" << endl;
    }

    void buildSparse() {
        /* Construct the Sparse Hamiltonian */

        // Get the list of coefficients
        vector<Tcd> tripletList = build_coefficients();
        // Initialize the sparse Hamiltonian
        spars_hamilt.resize(tot_states_, tot_states_);
        // Fill the Hamiltonian with the coefficients
        for (const auto& tripl : tripletList)
            spars_hamilt.coeffRef(tripl.row(),tripl.col()) += tripl.value();
        spars_hamilt.makeCompressed();

        // Check that the Hamiltonian is hermitian
        complex<double> rest;
        rest = (spars_hamilt - spars_hamilt.transpose().adjoint()).sum();
        if (rest != 0.) {
          cout << "Error occurred building the matrix: ";
          cout << "it is not hermitian!" << endl;
          show_hamiltonian();
          exit(1);
        }

        cout << "Building sparse complete!" << endl;
    }

    VectorXcd action(const VectorXcd& psi_in) {
        /* Action of the Hamiltonian on a given input state */

        int index;
        VectorXcd psi_out = VectorXcd::Zero(tot_states_);

        // Sigma_z [longitudinal field]
        for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == 1) psi_out[n] -= hz_field * psi_in[n];
                if (basis(n, i) == 0) psi_out[n] += hz_field * psi_in[n];
            }
        }
        // Sigma_x [transverse field]
        for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == 1) index = n - pow(2, tot_length_ - 1 - i);
                if (basis(n, i) == 0) index = n + pow(2, tot_length_ - 1 - i);
                psi_out[index] += hx_field * psi_in[n];
            }
        }
        // Sigma_y [transverse field]
        for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == 1) {
                    index = n - pow(2, tot_length_ - 1 - i);
                    psi_out[index] += hy_field * complex(0., -1.) * psi_in[n];
                }
                if (basis(n, i) == 0) {
                    index = n + pow(2, tot_length_ - 1 - i);
                    psi_out[index] += hy_field * complex(0., 1.) * psi_in[n];
                }
            }
        }
        // Sigma_z Sigma_z [coupling]
        for (int i = 0; i < tot_length_ - 1; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == basis(n, i + 1))
                    psi_out[n] += gz_field * psi_in[n];
                if (basis(n, i) != basis(n, i + 1))
                    psi_out[n] -= gz_field * psi_in[n];
            }
        }
        // Sigma_z Sigma_z [boundary conditions]
        if (pbc_flag_) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, tot_length_ - 1) == basis(n, 0))
                    psi_out[n] += gz_field * psi_in[n];
                if (basis(n, tot_length_ - 1) != basis(n, 0))
                    psi_out[n] -= gz_field * psi_in[n];
          }
        }
        // Sigma_x Sigma_x [coupling]
        for (int i = 0; i < tot_length_ - 1; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == 1)
                    index = n - pow(2, tot_length_ - 1 - i);
                if (basis(n, i) == 0)
                    index = n + pow(2, tot_length_ - 1 - i);
                if (basis(n, i+1) == 1)
                    index -= pow(2, tot_length_ - 2 - i);
                if (basis(n, i+1) == 0)
                    index += pow(2, tot_length_ - 2 - i);
                psi_out[index] += gx_field * psi_in[n];
            }
        }
        // Sigma_x Sigma_x [boundary conditions]
        if (pbc_flag_) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, 0) == 1)
                    index = n - pow(2, tot_length_ - 1);
                if (basis(n, 0) == 0)
                    index = n + pow(2, tot_length_ - 1);
                if (basis(n, tot_length_ - 1) == 1)
                    index -= 1;
                if (basis(n, tot_length_ - 1) == 0)
                    index += 1;
                psi_out[index] += gx_field * psi_in[n];
            }
        }
        // Sigma_y Sigma_y [coupling]
        for (int i = 0; i < tot_length_ - 1; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, i) == 1)
                    index = n - pow(2, tot_length_ - 1 - i);
                if (basis(n, i) == 0)
                    index = n + pow(2, tot_length_ - 1 - i);
                if (basis(n, i+1) == 1)
                    index -= pow(2, tot_length_ - 2 - i);
                if (basis(n, i+1) == 0)
                    index += pow(2, tot_length_ - 2 - i);
                if (basis(n, i) == basis(n, i+1))
                    psi_out[index] -= gy_field * psi_in[n];
                if (basis(n, i) != basis(n, i+1))
                    psi_out[index] += gy_field * psi_in[n];
            }
        }
        // Sigma_y Sigma_y [boundary conditions]
        if (pbc_flag_) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis(n, 0) == 1)
                    index = n - pow(2, tot_length_ - 1);
                if (basis(n, 0) == 0)
                    index = n + pow(2, tot_length_ - 1);
                if (basis(n, tot_length_ - 1) == 1)
                    index -= 1;
                if (basis(n, tot_length_ - 1) == 0)
                    index += 1;

                if (basis(n, tot_length_ - 1) == basis(n, 0))
                    psi_out[index] -= gy_field * psi_in[n];
                if (basis(n, tot_length_ - 1) != basis(n, 0))
                    psi_out[index] += gy_field * psi_in[n];
            }
        }

        return psi_out;
    }

    double average_energy(const VectorXcd& psi_in){
        /* Braket of the Hamiltonian on a given input state */

        complex<double> energy = 0.;
        VectorXcd psi_out = action(psi_in);
        energy = psi_in.adjoint() * psi_out;
        return energy.real();
    }

    void diagonalize(){
        /* Diagonalization subroutine with Eigen and Lanczos */

        if (!spr_flag_) {
            // Complete diagonalization with Eigen -----------------------------
            // The matrix is first reduced to Schur form. The Schur decomposition
            // is then used to compute the eigenvalues and the corresponding eigenvectors.
            // https://eigen.tuxfamily.org/dox/classEigen_1_1ComplexEigenSolver.html

            // Construct solver object and compute
            //ComplexEigenSolver<MatrixXcd> solver;
            SelfAdjointEigenSolver<MatrixXcd> solver;
            solver.compute(dense_hamilt);
            // Retrieve results
            if (solver.info() == 0) {
                eigenvalues = solver.eigenvalues();
                eigenvectors = solver.eigenvectors();
            } else {
                cerr << "An error occurred in dense diagonalization!" << endl;
            }

            //------------------------------------------------------------------
        } else {
            // Partial diagonalization with Lanczos ----------------------------
            // Lambda Lanczos calculates the smallest or largest eigenvalue and
            // the corresponding eigenvector of a symmetric (Hermitian) matrix.
            // https://github.com/Dario-Maglio/lambda-lanczos/blob/master/src/samples/sample4_use_Eigen_library.cpp

            // Define matrix-vector multiplication routine
            auto amul = [&](const vector<double>& in, vector<double>& out) {
                auto psi_in = Map<const VectorXd>(&in[0], in.size());
                auto psi_out = Map<VectorXd>(&out[0], out.size());
                VectorXcd psi = psi_in;
                psi_out = action(psi_in).real();
            };
            // Construct solver object and compute
            vector<double> eigvalues; vector<vector<double>> eigvectors;
            LambdaLanczos<double> solver(amul, tot_states_, false, spr_flag_);
            solver.run(eigvalues, eigvectors);
            // Retrieve results
            for (int i = 0; i < eigenvalues.size(); i++) {
                eigenvalues[i] = eigvalues[i];
                for (int n = 0; n < tot_states_; n++)
                    eigenvectors(n, i) = eigvectors[i][n];
            }

            //------------------------------------------------------------------
        }
    }

    //--- Set and Get methods --------------------------------------------------

    void set_fields(const HamiltParameters& param, bool debug=0) {
        /* Rebuild the Hamiltonian changing the field values */

        // Print changes
        if (debug) cout << "New couplings are "<< endl << setprecision(4)
            << "gz: " << gz_field << " --> " << param.gz_field << endl
            << "gx: " << gx_field << " --> " << param.gx_field << endl
            << "gy: " << gy_field << " --> " << param.gy_field << endl
            << "hz: " << hz_field << " --> " << param.hz_field << endl
            << "hx: " << hx_field << " --> " << param.hx_field << endl
            << "hy: " << hy_field << " --> " << param.hy_field << endl << endl;
        // Save the new couplings
        gx_field = param.gx_field;
        gy_field = param.gy_field;
        gz_field = param.gz_field;
        hx_field = param.hx_field;
        hy_field = param.hy_field;
        hz_field = param.hz_field;
        // Re-construct the Hamiltonian
        (spr_flag_)? buildSparse() : buildHamilt();

    }

    double get_eigenvalue(int k){
        /* After Diagonalization returns the k-th eigenvalue */

        int limit = tot_states_;
        if (spr_flag_) limit = spr_flag_;

        if (k >= limit)
            cerr << "Index exceeds the amount of eigenvectors." << endl;

        return eigenvalues[k];
    }

    VectorXcd get_eigenvector(int k){
        /* After Diagonalization returns the k-th eigenvector */

        int limit = tot_states_;
        if (spr_flag_) limit = spr_flag_;

        if (k >= limit)
            cerr << "Index exceeds the amount of eigenvectors." << endl;

        return eigenvectors.col(k);
    }

    //--- Show methods ---------------------------------------------------------

    void show_comput_basis() {
        cout << "Total number of states: " << tot_states_ << endl;
        cout << "Lattice computational basis: " << endl;
        for (int n = 0; n < tot_states_; n++) {
            cout << "|";
            for (int i = 0; i < tot_length_; i++) cout << basis(n, i);
            cout << ">  ";
        }
        cout << endl << endl;
    }

    void show_hamiltonian() {
        /* Print the Hamiltonian matrix */

        if (spr_flag_) {
            cout << "--- Sparse Hamiltonian" << endl << endl;
            cout << setprecision(2) << spars_hamilt << endl << endl;
        } else {
            cout << "--- Dense Hamiltonian" << endl;
            cout << setprecision(2) << dense_hamilt << endl << endl;
        }
    }

    void show_eigenvalues() {
        /* Print stored eigenvalues */

        cout << "Eigenvalues" << endl;
        for (const auto& val : eigenvalues)
            cout << setprecision(4) << val << " ";
        cout << endl << endl;
    }

    void show_eigenvectors() {
        /* Print stored eigenvectors */

        cout << "Eigenvectors" << endl;
        for (int i = 0; i < eigenvalues.size(); i++) {
            cout << i << ": ";
            for (int n = 0; n < tot_states_; n++)
                cout << fixed << setprecision(3) << eigenvectors(n, i) << " ";
            cout << endl;
        }
        cout << endl;
    }

    void show_eigen() {
        /* Print stored eigenvalues and corresponding eigenvectors */

        cout << "List of eigenvalues and corresponding eigenvectors" << endl;
        for (int i = 0; i < eigenvalues.size(); i++) {
            cout << "Eigenvalue " << i << " -> " << get_eigenvalue(i) << endl;
            for (int n = 0; n < tot_states_; n++) {
                cout << fixed << setprecision(3) << eigenvectors(n, i) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

};

#endif
