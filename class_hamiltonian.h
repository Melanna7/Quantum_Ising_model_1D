/*******************************************************************************
*
* Hamiltonian class definition
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#ifndef ISING_HAMILT_CLASS_H
#define ISING_HAMILT_CLASS_H

#include <iostream>
#include <iomanip>

#include <vector>
#include <complex>
#include <cmath>

#include <lapacke.h>

using namespace std;

#define DIM_HILBERT 2

struct HamiltParameters {
    int pbc_flag, num_sites, sparse_flag;
    double gz_field, hz_field, hx_field;
};

//--- Contents -----------------------------------------------------------------

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

    H_FIELD = amplitude of the longitudinal field

    T_FIELD = amplitude of the transverse field

    *************************************************************/
private:
    int spr_flag_, pbc_flag_;
    int tot_length_, tot_states_;
    double gz_field, hz_field, hx_field;

public:
    vector<double> eigval;
    vector<vector<int>> basis;
    vector<vector<complex<double>>> hamilt, eigvec;

    hamiltonian(const HamiltParameters& param):
        tot_length_(param.num_sites),
        spr_flag_(param.sparse_flag),
        pbc_flag_(param.pbc_flag),
        gz_field(param.gz_field),
        hz_field(param.hz_field),
        hx_field(param.hx_field)
    {// CONSTRUCTION BEGIN
        int index = 0;

        // Compute the number of states
        tot_states_ = pow(DIM_HILBERT, tot_length_);

        // Build the computational basis
        basis.resize(tot_states_, vector<int>(tot_length_, 0.));
        for (int n = 0; n < tot_states_; n++) {
            index = n;
            for (int i = 0; i < tot_length_; i++) {
                basis[n][tot_length_ - i - 1 ] = index % 2;
                index = index / 2;
            }
        }

        if (spr_flag_) {
            // Initialize partial eigen-stuff
            eigval.resize(spr_flag_);
            eigvec.resize(spr_flag_, vector<complex<double>>(tot_states_, 0.));
        } else {
            // Initialize complete eigen-stuff
            eigval.resize(tot_states_);
            eigvec.resize(tot_states_, vector<complex<double>>(tot_states_, 0.));
            hamilt.resize(tot_states_, vector<complex<double>>(tot_states_, 0.));
            // Build the Hamiltonian
            buildHamilt();
        }

    }// CONSTRUCTION END

    void buildHamilt() {
        /* Subroutine to build the Hamiltonian */

        int index = 0;
        // Initialize Hamiltonian
        for (int i = 0; i < tot_states_; i++)
            for (int j = 0; j < tot_states_; j++)
                hamilt[i][j] = 0.;
        // Sigma_z Sigma_z [coupling]
        for (int i = 0; i < tot_length_ - 1; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][i] == basis[n][i+1]) hamilt[n][n] += gz_field;
                if (basis[n][i] != basis[n][i+1]) hamilt[n][n] -= gz_field;
            }
        }
        // Sigma_z Sigma_z boundary conditions
        if (pbc_flag_) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][tot_length_ - 1] == basis[n][0])
                    hamilt[n][n] += gz_field;
                if (basis[n][tot_length_ - 1] != basis[n][0])
                    hamilt[n][n] -= gz_field;
            }
        }
        // Sigma_z [longitudinal field]
        for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][i] == 1) hamilt[n][n] -= hz_field;
                if (basis[n][i] == 0) hamilt[n][n] += hz_field;
            }
        }
        // Sigma_x [transverse field]
        for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][i] == 1) index = n - pow(2, tot_length_ - 1 - i);
                if (basis[n][i] == 0) index = n + pow(2, tot_length_ - 1 - i);
                hamilt[index][n] += hx_field;
            }
        }
        // Check if the Hamiltonian is hermitian
        for (int i = 0; i < tot_states_; i++) {
            for (int j = 0; j < tot_states_; j++) {
                if (hamilt[i][j] != conj(hamilt[j][i])) {
                  cout << "Error occurred building the matrix: ";
                  cout << "it is not hermitian!" << endl;
                  show_hamiltonian();
                  exit(1);
                }
            }
        }

        cout << "Building complete!" << endl;
    }

    vector<complex<double>> action(const vector<complex<double>>& psi_in){
        /* Action of the Hamiltonian on a given input state */

        int index;
        vector<complex<double>> psi_out(tot_states_, 0.);

        // Sigma_z Sigma_z [coupling]
        for (int i = 0; i < tot_length_ - 1; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][i] == basis[n][i + 1])
                    psi_out[n] += gz_field * psi_in[n];
                if (basis[n][i] != basis[n][i + 1])
                    psi_out[n] -= gz_field * psi_in[n];
            }
        }
        // Sigma_z Sigma_z boundary conditions
        if (pbc_flag_) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][tot_length_ - 1] == basis[n][0])
                    psi_out[n] += gz_field * psi_in[n];
                if (basis[n][tot_length_ - 1] != basis[n][0])
                    psi_out[n] -= gz_field * psi_in[n];
          }
        }
        // Sigma_z [longitudinal field]
        for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][i] == 1) psi_out[n] -= hz_field * psi_in[n];
                if (basis[n][i] == 0) psi_out[n] += hz_field * psi_in[n];
            }
        }
        // Sigma_x [transverse field]
        for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][i] == 1) index = n - pow(2, tot_length_ - 1 - i);
                if (basis[n][i] == 0) index = n + pow(2, tot_length_ - 1 - i);
                psi_out[index] += hx_field * psi_in[n];
            }
        }

        return psi_out;
    }

    vector<complex<double>> evolution(vector<complex<double>> state_in, double t){
        /* Hamiltonian evolution U(t) of a given input state */

        vector<complex<double>> state_out(tot_states_);
        cout << "Unitary evolution not implemented yet!" << endl << endl;
        return state_out;
    }

    void set_fields(const HamiltParameters& param) {
        /* Rebuild the Hamiltonian changing the value of gz_field */

        gz_field = param.gz_field;
        hz_field = param.hz_field;
        hx_field = param.hx_field;
        if (!spr_flag_) buildHamilt();
        cout << "Coupling setted to: "<< endl;
        cout << "gz -> " << setprecision(4) << param.gz_field << endl;
        cout << "hz -> " << setprecision(4) << param.hz_field << endl;
        cout << "hx -> " << setprecision(4) << param.hx_field << endl << endl;
    }

    vector<complex<double>> get_eigenstate(int k){
        /* After Diagonalization returns the k-th eigenvector */

        if (k >= tot_states_)
            cerr << "Index exceeds the amount of eigenvectors." << endl;

        return eigvec[k];
    }

    double get_eigenvalue(int k){
        /* After Diagonalization returns the k-th eigenvector */

        if (k >= tot_states_)
            cerr << "Index exceeds the amount of eigenvectors." << endl;

        return eigval[k];
    }

    vector<complex<double>> compute_GS(){
        /* Diagonalize and return the ground state eigenvector */

        vector<complex<double>> state(tot_states_);

        diagonalize();

        cout << "The ground state energy is: " << eigval[0] << endl;
        cout << "It is associated to the eigenstate:" << endl;
        for (int n = 0; n < tot_states_; n++) {
            cout << fixed << setprecision(2) << eigvec[0][n] << " ";
            state[n] = eigvec[0][n];
        }
        cout << endl << endl;

        return state;
    }

    double average_energy(const vector<complex<double>>& psi_in){
        /* Braket of the Hamiltonian on a given input state */

        complex<double> energy = 0.;
        vector<complex<double>> psi_out(tot_states_, complex<double>(0.0, 0.0));

        psi_out = action(psi_in);
        for(int n = 0; n < tot_states_; n++){
            energy += conj(psi_in[n]) * psi_out[n];
        }

        return energy.real();
    }

    void diagonalize(){
        /* Diagonalization subroutine with LAPACK */

        if (!spr_flag_) {
            // Complete diagonalization with LAPACKE ---------------------------

            lapack_int info;
            lapack_int n = tot_states_;
            lapack_int lda = n; // leading dimension of the matrix
            double eigenval[n]; // eigenvalues vector for zheev
            double val_1, val_2;
            lapack_complex_double matr[lda*n]; // hamiltonian

            // Cast the Hamiltonian into a LAPACK object
            for (int i = 0; i < n; i++) {
                eigenval[i] = 0.;
                for (int j = 0; j < n; j++) {
                    val_1 = hamilt[i][j].real();
                    val_2 = hamilt[i][j].imag();
                    matr[i*lda + j] = lapack_make_complex_double(val_1, val_2);
                }
            }
            // Solve eigenproblem
            // https://www.netlib.org/lapack/explore-html/df/d9a/
            //group__complex16_h_eeigen_gaf23fb5b3ae38072ef4890ba
            //43d5cfea2.html#gaf23fb5b3ae38072ef4890ba43d5cfea2
            info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', n, matr, lda, eigenval);
            // Check for convergence
            if (info > 0) cerr << "Algorithm zheev failed to compute eigenvalues." << endl;
            // Update eigenvalues and eigenvectors
            for (int i = 0; i < n; i++) {
                eigval[i] = eigenval[i];
                for (int j = 0; j < lda; j++) {
                    eigvec[i][j] = matr[j*lda + i];
                }
            }
            //------------------------------------------------------------------
        } else {
            // Partial diagonalization with ARPACK++

            //------------------------------------------------------------------
        }

    }

    //--- Show methods ---------------------------------------------------------

    void show_comput_basis() {
        cout << "Total number of states: " << tot_states_ << endl;
        cout << "Lattice computational basis: " << endl;
        for (const auto& state : basis) {
            cout << "|";
            for (int i = 0; i < tot_length_; i++) cout << state[i];
            cout << ">  ";
        }
        cout << endl << endl;
    }

    void show_hamiltonian() {
        /* Print the Hamiltonian matrix */

        if (spr_flag_) {
            cout << "Not allowed to print sparse Hamiltonian!" << endl << endl;
        } else {
            cout << "Hamiltonian" << endl;
            for (const auto& row : hamilt) {
                for (const auto& element : row) {
                    cout << fixed << setprecision(2) << element << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }

    void show_eigenvectors() {
        /* Print stored eigenvectors */

        cout << "Eigenvectors" << endl;
        for (const auto& row : eigvec) {
            for (const auto& element : row) {
                cout << fixed << setprecision(2) << element << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    void show_eigenvalues() {
        /* Print stored eigenvalues */

        cout << "Eigenvalues" << endl;
        for (const auto& val : eigval) cout << val << " ";
        cout << endl << endl;
    }

    void show_eigen() {
        /* Print stored eigenvalues and corresponding eigenvectors */

        for (int i = 0; i < eigval.size(); i++) {
            cout << "Eigenvalue " << i << " -> " << get_eigenvalue(i) << endl;
            for (const auto& element : eigvec[i]) {
                cout << fixed << setprecision(2) << element << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

};

#endif
