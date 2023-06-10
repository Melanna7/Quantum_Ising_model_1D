/*******************************************************************************
*
* Hamiltonian class definition
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#ifndef ISING_HAMILT_CLASS_H
#define ISING_HAMILT_CLASS_H

#include <iostream>
#include <iomanip> // for setting output precision

#include <vector>
#include <complex>
#include <cmath>

#include <lapacke.h>

using namespace std;

#define DIM_HILBERT 2

//--- Contents -----------------------------------------------------------------

class hamiltonian {
private:
    int pbc_flag_;
    int tot_length_, tot_states_;

public:
    const double g_coupling, h_field, t_field;
    vector<vector<int>> basis;
    vector<double> eigenvalues;
    vector<vector<complex<double>>> hamilt, eigenvectors;

    hamiltonian(const int &LENGTH, const int &PBC_FLAG, const double &G_FIELD, const double &H_FIELD, const double &T_FIELD):
        /* Class lattice constructor **********************************

        Inputs:

        LENGHT = total number of sites.

        PBC_FLAG = flag for the boundary conditions.
                   0 open boundary conditions.
                   else periodic boundary conditions.

        G_FIELD = amplitude of the spin coupling field

        H_FIELD = amplitude of the longitudinal field

        T_FIELD = amplitude of the transverse field

        *************************************************************/
        tot_length_(LENGTH),
        pbc_flag_(PBC_FLAG),
        g_coupling(G_FIELD),
        h_field(H_FIELD),
        t_field(T_FIELD)
    {// CONSTRUCTION BEGIN
        int index;
        complex<double> zero_compl(0.0, 0.0);

        // Compute the number of states
        tot_states_ = pow(DIM_HILBERT, tot_length_);

        // Initialize eigen-stuff
        eigenvalues.resize(tot_states_);
        eigenvectors.resize(tot_states_, vector<complex<double>>(tot_states_, zero_compl));

        // Build the computational basis
        basis.resize(tot_states_, vector<int>(tot_length_, 0.));
        for (int i = 0; i < tot_states_; i++) {
            index = i;
            for (int j = 0; j < tot_length_; j++) {
                basis[i][tot_length_ - j - 1 ] = index % 2;
                index = index / 2;
            }
        }

        // Build the Hamiltonian
        hamilt.resize(tot_states_, vector<complex<double>>(tot_states_, zero_compl));
        // Sigma_z Sigma_z [coupling]
        for (int i = 0; i < tot_length_ - 1; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][i] == basis[n][i+1]) hamilt[n][n] += g_coupling;
                if (basis[n][i] != basis[n][i+1]) hamilt[n][n] -= g_coupling;
            }
        }
        // Sigma_z Sigma_z boundary conditions
        if (pbc_flag_) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][tot_length_ - 1] == basis[n][0])
                    hamilt[n][n] += g_coupling;
                if (basis[n][tot_length_ - 1] != basis[n][0])
                    hamilt[n][n] -= g_coupling;
            }
        }
        // Sigma_z [longitudinal field]
        for (int i = 0; i < tot_length_; i++) {
            for (int n = 0; n < tot_states_; n++) {
                if (basis[n][i] == 1) hamilt[n][n] += h_field;
                if (basis[n][i] == 0) hamilt[n][n] -= h_field;
            }
        }
        // Sigma_x [transverse field]
        // for (int i = 0; i < tot_length_; i++) {
        //     for (int n = 0; n < tot_states_; n++) {
        //         if (basis[n][i] == 0) index = n + pow(2, i + 1);
        //         if (basis[n][i] == 1) index = n - pow(2, i + 1);
        //         hamilt[index][n] += t_field;
        //     }
        // }

    }// CONSTRUCTION END

    void show_comput_basis() {
        cout << "Total number of states: " << tot_states_ << endl;
        cout << "Lattice computational basis: " << endl;
        for (auto state : basis) {
            cout << "|";
            for (int i = 0; i < tot_length_; i++) cout << state[i];
            cout << ">  ";
        }
        cout << endl << endl;
    }

    void show_hamiltonian() {
        cout << "Hamiltonian" << endl;
        for (const auto& row : hamilt) {
            for (const auto& element : row) {
                cout << fixed << setprecision(2) << element << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    void show_eigenvectors() {
        cout << "Eigenvectors" << endl;
        for (const auto& row : eigenvectors) {
            for (const auto& element : row) {
                cout << fixed << setprecision(2) << element << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    void show_eigenvalues() {
        cout << "Eigenvalues" << endl;
        for (auto val : eigenvalues) cout << val << " ";
        cout << endl << endl;
    }

    void diagonalize(){
        /* Diagonalization subroutine with LAPACK */

        lapack_int n = tot_states_, lda = n, info;
        double val_1, val_2, eigenval[n];
        lapack_complex_double lap_matr[n*lda] = {0};
        lapack_complex_double eigenvec[lda*n] = {0};

        // Cast the Hamiltonian in a LAPACK object
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                val_1 = hamilt[i][j].real();
                val_2 = hamilt[i][j].imag();
                lap_matr[i*lda + j] = lapack_make_complex_double(val_1, val_2);
            }
        }

        // Solve eigenproblem
        info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', n, lap_matr, lda, eigenval);
        //info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'N', 'V', n, lap_matr, lda, eigenval, 0, lda, eigenvec, lda);

        // Check for convergence
        if (info > 0) {
            cout << "ZHEEV algorithm failed to compute eigenvalues." << endl;
            exit(1);
        }

        // Update eigenvalues and eigenvectors
        for (int i = 0; i < n; i++) {
            eigenvalues[i] = eigenval[i];
            for (int j = 0; j < lda; j++) {
                eigenvectors[j][i] = lap_matr[j*lda + i];
            }
        }
    }

    vector<complex<double>> get_GS(){
        /* After Diagonalization returns the ground state eigenvector */

        vector<complex<double>> state(tot_states_);

        cout << "The ground state energy is: " << eigenvalues[0] << endl;
        cout << "It is associated to the eigenstate:" << endl;
        for (int i = 0; i < tot_states_; i++) {
            cout << fixed << setprecision(2) << eigenvectors[i][0] << " ";
            state.push_back(eigenvectors[i][0]);
        }
        cout << endl;

        return state;
    }

    vector<complex<double>> action(vector<complex<double>> state_in){
        /* Action of the Hamiltonian on a given input state */

        vector<complex<double>> state_out(tot_states_);

        return state_out;
    }

    vector<complex<double>> evolution(vector<complex<double>> state_in, double t){
        /* Hamiltonian evolution U(t) of a given input state */

        vector<complex<double>> state_out(tot_states_);

        return state_out;
    }

};

#endif
