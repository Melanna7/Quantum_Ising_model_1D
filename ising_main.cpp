/*******************************************************************************
*
* Main program for the Ising simulation
*
*******************************************************************************/

// g++ ising_main.cpp -llapacke -I ./include/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

// #include "class_Eigen.h"
#include "class_lapacke.h"

using namespace std;

/*******************************************************************************
* PARAMETERS OF THE LATTICE
*
* SIDE = size of the lattice's side.
*
* SPARSE_FLAG = condition on storing or not the Hamiltonian.
*               0 build the Hamiltonian and LAPACK diagonalization.
*               else number of desired eigenvalues with Lancsoz method.
*
* PBC_FLAG = flag for the boundary conditions.
*            0 for open boundary conditions.
*            else periodic boundary conditions.
*
* G_FIELD = amplitude of the spin coupling field
*
* H_FIELD = amplitude of the external field
*
*******************************************************************************/

#define DIM_HILBERT 2

// Range sides
#define MIN_SIDE 4
#define MAX_SIDE 14
#define LIMIT_SIDE 8
// General settings
#define SPARSE_FLAG 3
#define PBC_FLAG 1
#define GZ_FIELD -1.
#define HZ_FIELD 0.
// Range of HX
#define MIN_HX_FIELD -2.00
#define MAX_HX_FIELD 0.00
#define SEP_HX_FIELD 0.02
// For susceptibility
#define DELTA_HZ -0.001

//--- Contents -----------------------------------------------------------------

vector<double> magnetization(const vector<complex<double>>& psi, const HamiltParameters& param) {

    int index, tot_states = param.num_states;
    int sites = param.num_sites;
    double mag_z_tilde, mag_x, mag_y, mag_z;
    complex<double> val_x, val_y;
    vector<double> sigmaX(sites, 0.), sigmaY(sites, 0.), sigmaZ(sites, 0.);
    vector<double> output;

    // Compute Symmetry-broken magnetization
    mag_z_tilde = 0.;
    for (int n = 0; n < tot_states; n++) {
        index = 0;
        for (int i = 0; i < sites; i++) {
            if (param.comp_basis[n][i] == 1) index += 1;
            if (param.comp_basis[n][i] == 0) index -= 1;
        }
        mag_z_tilde += abs(1.0 * index) * norm(psi[n]);
    }
    mag_z_tilde = mag_z_tilde / sites;
    output.push_back(mag_z_tilde);
    cout << "Symmetry-broken magnet along Z: " << mag_z_tilde << endl;

    // Compute average of sigmaX, sigmaY and sigmaZ over psi
    for (int i = 0; i < sites; i++) {
        val_x = 0.;
        val_y = 0.;
        for (int n = 0; n < tot_states; n++) {
            // Compute sigmaX_i
            if (param.comp_basis[n][i] == 1) index = n - pow(2, sites - 1 - i);
            if (param.comp_basis[n][i] == 0) index = n + pow(2, sites - 1 - i);
            val_x += conj(psi[n]) * psi[index];
            // Compute sigmaY_i
            if (param.comp_basis[n][i] == 1) {
                index = n - pow(2, sites - 1 - i);
                val_y += conj(psi[n]) * complex(0., -1.) * psi[index];
            }
            if (param.comp_basis[n][i] == 0) {
                index = n + pow(2, sites - 1 - i);
                val_y += conj(psi[n]) * complex(0., 1.) * psi[index];
            }
            // Compute sigmaZ_i
            if (param.comp_basis[n][i] == 1) sigmaZ[i] += norm(psi[n]);
            if (param.comp_basis[n][i] == 0) sigmaZ[i] -= norm(psi[n]);
        }

        if (abs(imag(val_x)) > 1.0e-10) cerr << "Non-real X magnetiz!" << endl;
        if (abs(imag(val_y)) > 1.0e-10) cerr << "Non-real Y magnetiz!" << endl;
        sigmaX[i] = real(val_x);
        sigmaY[i] = real(val_y);
    }

    mag_x = 0.;
    mag_y = 0.;
    mag_z = 0.;
    for (int i = 0; i < sites; i++) {
        mag_x += sigmaX[i];
        mag_y += sigmaY[i];
        mag_z += sigmaZ[i];
    }
    mag_x = mag_x / sites;
    mag_y = mag_y / sites;
    mag_z = mag_z / sites;

    output.push_back(mag_x);
    output.push_back(mag_y);
    output.push_back(mag_z);
    return output;
}

double susceptibility(const vector<complex<double>>& psi, const HamiltParameters& param) {
    int index, tot_states = param.num_states;
    int sites = param.num_sites;
    double mag_z, chi;
    vector<double> sigmaZ(sites, 0.);

    // Compute average of sigmaZ over psi
    for (int i = 0; i < sites; i++) {
        for (int n = 0; n < tot_states; n++) {
            // Compute sigmaZ_i
            if (param.comp_basis[n][i] == 1) sigmaZ[i] += norm(psi[n]);
            if (param.comp_basis[n][i] == 0) sigmaZ[i] -= norm(psi[n]);
        }
    }

    mag_z = 0.;
    for (int i = 0; i < sites; i++) mag_z += sigmaZ[i];
    mag_z = mag_z / sites;
    chi = mag_z / DELTA_HZ;

    return chi;
}

void run_simulation(HamiltParameters param){
    /* Simulation for a given side */

    ofstream file;
    string directory, file_name;
    double ener_gs, ener_gap, chi;
    vector<double> magZXYZ;
    vector<complex<double>> psi;
    hamilt_L HamOp(param);
    // VectorXcd psi_gs;
    // hamiltonian HamOp(param);

    // Define path data directory
    directory = "Data_Ising/";
    // Define name file last configuration of the lattice
    file_name = "side_" + to_string(param.num_sites) + ".dat";

    file.open(directory + file_name);
    // Update transverse field and take measures
    for(double hx = MIN_HX_FIELD; hx < MAX_HX_FIELD; hx += SEP_HX_FIELD){
        // Store hx field value
        file << -hx << " ";
        // Set field and diagonalize
        param.hx_field = hx;
        param.hz_field = HZ_FIELD;
        HamOp.set_fields(param);
        HamOp.diagonalize();
        // Get and store GS energy and gaps
        cout << "Side " << param.num_sites;
        cout << " - Field " << setprecision(4) << -hx << endl;
        ener_gs = HamOp.get_eigenvalue(0);
        file << ener_gs << " ";
        cout << "Ground state energy: " << ener_gs << endl;
        ener_gap = HamOp.get_eigenvalue(1) - ener_gs;
        file << ener_gap << " ";
        cout << "First energy gap: " << ener_gap << endl;
        ener_gap = HamOp.get_eigenvalue(2) - ener_gs;
        file << ener_gap << " ";
        cout << "Second energy gap: " << ener_gap << endl;
        // Get and store GS magnetizations
        psi = HamOp.get_eigenvector(0);
        // psi_gs = HamOp.get_eigenvector(0);
        // for(int i = 0; i < psi_gs.size(); i++) psi[i] = psi_gs[i];
        magZXYZ = magnetization(psi, param);
        file << magZXYZ[0] << " ";
        file << magZXYZ[1] << " ";
        file << magZXYZ[2] << " ";
        file << magZXYZ[3] << " ";
        // Set field and diagonalize to find chi
        param.hz_field = DELTA_HZ;
        HamOp.set_fields(param);
        HamOp.diagonalize();
        psi = HamOp.get_eigenvector(0);
        // psi_gs = HamOp.get_eigenvector(0);
        // for(int i = 0; i < psi_gs.size(); i++) psi[i] = psi_gs[i];
        chi = susceptibility(psi, param);
        file << chi << endl;
        cout << "Finished." << endl << endl;
    }

    file.close();
}

//--- Main sim -----------------------------------------------------------------

int main(){

    HamiltParameters param;
    param.sparse_flag = SPARSE_FLAG;
    param.pbc_flag = PBC_FLAG;
    param.num_sites = MIN_SIDE;
    param.gz_field = GZ_FIELD;
    param.hz_field = HZ_FIELD;
    param.hx_field = MIN_HX_FIELD;

    for (int side = MAX_SIDE; side >= MIN_SIDE; side--){

        param.num_sites = side;
        param.num_states = pow(DIM_HILBERT, side);
        // Check maximum side for complete diagonalization
        if(side == LIMIT_SIDE) param.sparse_flag = 0;

        // Store computational basis
        param.comp_basis.resize(param.num_states, vector<int>(side, 0.));
        for(int n = 0; n < param.num_states; n++) {
            int index = n;
            for (int i = 0; i < side; i++) {
                param.comp_basis[n][side - i - 1 ] = index % 2;
                index = index / 2;
            }
        }
        // Start simulation
        cout << "Start new side" << endl;
        auto start = chrono::steady_clock::now();
        run_simulation(param);
        auto end = chrono::steady_clock::now();
        chrono::duration<double> elapsed_sec = end - start;
        cout << "Elapsed time: " << elapsed_sec.count() << "s" << endl << endl;
    }
}
