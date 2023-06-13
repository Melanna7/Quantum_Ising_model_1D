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
#include <dsaupd.h>

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

    vector<complex<double>> evolution(vector<complex<double>> state_in, double t){
        /* Hamiltonian evolution U(t) of a given input state */

        vector<complex<double>> state_out(tot_states_);
        cout << "Unitary evolution not implemented yet!" << endl << endl;
        return state_out;
    }

    void set_gz_field(double gz) {
        /* Rebuild the Hamiltonian changing the value of gz_field */

        gz_field = gz;
        if (!spr_flag_) buildHamilt();
        cout << "Coupling g setted to " << gz << endl << endl;
    }

    void set_hz_field(double hz) {
        /* Rebuild the Hamiltonian changing the value of hz_field */

        hz_field = hz;
        if (!spr_flag_) buildHamilt();
        cout << "Coupling h setted to " << hz << endl << endl;
    }

    void set_hx_field(double hx) {
        /* Rebuild the Hamiltonian changing the value of gz_field */

        hx_field = hx;
        if (!spr_flag_) buildHamilt();
        cout << "Coupling t setted to " << hx << endl << endl;
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

    void diagonalize(){
        /* Diagonalization subroutine with LAPACK */

        if (!spr_flag_) {
            // Complete diagonalization with LAPACKE ---------------------------

            lapack_int info;
            lapack_int n = tot_states_;
            lapack_int lda = n; // leading dimension of the matrix
            double eigenval[n]; // eigenvalues vector for zheev
            double val_1, val_2;
            lapack_complex_double val_c;
            lapack_complex_double lap_matr[lda*n]; // hamiltonian

            // Cast the Hamiltonian into a LAPACK object
            for (int i = 0; i < n; i++) {
                eigenval[i] = 0.;
                for (int j = 0; j < n; j++) {
                    val_1 = hamilt[i][j].real();
                    val_2 = hamilt[i][j].imag();
                    val_c = lapack_make_complex_double(val_1, val_2);
                    lap_matr[i*lda + j] = val_c;
                }
            }
            // Solve eigenproblem
            // https://www.netlib.org/lapack/explore-html/df/d9a/
            //group__complex16_h_eeigen_gaf23fb5b3ae38072ef4890ba
            //43d5cfea2.html#gaf23fb5b3ae38072ef4890ba43d5cfea2
            info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', n, lap_matr, lda, eigenval);
            // Check for convergence
            if (info > 0) cerr << "Algorithm zheev failed to compute eigenvalues." << endl;
            // Update eigenvalues and eigenvectors
            for (int i = 0; i < n; i++) {
                eigval[i] = eigenval[i];
                for (int j = 0; j < lda; j++) {
                    eigvec[i][j] = lap_matr[j*lda + i];
                }
            }
            //------------------------------------------------------------------
        } else {
            // Partial diagonalization with ARPACK++
            // https://gist.github.com/wpoely86/043d0eb92f7e42563e7f -----------

            int info = 0;
            int n = tot_states_; // dimension of the matrix
            int nev = spr_flag_; // number of eigenvalues to calculate
            int ido = 0; // reverse communication param, zero on first iter
            char bmat = 'I'; // standard eigenvalue problem A*x=lambda*x
            char which[] = {'S','A'}; // smallest algebraic eigenvalue
            double tol = 0; // calculate until machine precision
            // array used for reverse communication
            double *workd = new double[3*n];
            for(int i=0;i<3*n;i++) workd[i] = 0;
            // the residual vector
            double *resid = new double[n];
            // lanczos vectors generated at each iteration
            int ncv = 42;
            if (n < ncv) ncv = n;
            double *v = new double[n*ncv];
            // length of the workl array
            int lworkl = ncv*(ncv+8);
            double *workl = new double[lworkl];
            // parameters for mode of dsaupd
            int *iparam = new int[11];
            /* Sets the mode of dsaupd.   *
            *  1 is exact shifting,       *
            *  2 is user-supplied shifts, *
            *  3 is shift-invert mode,    *
            *  4 is buckling mode,        *
            *  5 is Cayley mode.          */
            iparam[6] = 1;
            iparam[0] = 1;   // specifies the shift strategy (1->exact)
            iparam[2] = 3*n; // maximum number of iterations
            // other parameters
            int *ipntr = new int[11];
            /* Indicates the locations in the work array *
            *  workd where the input and output vectors  *
            *  in thecallback routine are located.       */
            double sigma; // not used if iparam[6] == 1

            // This vector will return the eigenvalues from dseupd
            double *d = new double[nev];
            // rvec == 0 : calculate only eigenvalue
            // rvec > 0 : calculate eigenvalue and eigenvector
            int rvec = 1;
            double *z = 0;
            // This vector will return the eigenvectors from dseupd
            if (rvec) z = new double[n*nev];
            char howmny = 'A'; // 'A' => nev eigenvectors
            int *select; // if howmny == 'A', workspace to reorder eigenvec
            if (howmny == 'A') select = new int[ncv];

            // Solve eigenproblem: first iteration
            dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
            // Solve eigenproblem: loop...
            while (ido != 99) {
                // matrix-vector multiplication
                mvprod(workd+ipntr[0]-1, workd+ipntr[1]-1,0);
                // following iterations
                dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
            }
            // Check for errors in the loop
            if (info < 0)
                cerr << "Error with dsaupd, info = " << info << endl;
            else if (info == 1)
                cerr << "Maximum number of Lanczos iterations reached." << endl;
            else if (info == 3)
                cerr << "No shifts could be applied during implicit Arnoldi update, try increasing NCV." << endl;

            // Solve eigenproblem: last step
            dseupd_(&rvec, &howmny, select, d, z, &n, &sigma, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
            // Check for convergence
            if ( info != 0 ) cerr << "Error with dseupd, info = " << info << endl;
            // Update eigenvalues and eigenvectors
            for (int i = 0; i < nev; i++) {
                eigval[i] = d[i];
                for (int j = 0; j < n; j++) {
                    eigvec[i][j] = z[i*n + j];
                }
            }

            // Remember to free up memory
            delete [] resid;
            delete [] v;
            delete [] iparam;
            delete [] ipntr;
            delete [] workd;
            delete [] workl;
            delete [] d;
            if (rvec) delete [] z;
            if (howmny == 'A') delete [] select;
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
