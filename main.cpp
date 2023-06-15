/*******************************************************************************
*
* Main program for the Ising simulation
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <cmath>

#include "class_hamiltonian.h"

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
* H_FIELD = amplitude of the longitudinal field
*
* T_FIELD = amplitude of the transverse field
*
*******************************************************************************/

#define MIN_SIDE 4
#define MAX_SIDE 10

#define SPARSE_FLAG 0
#define PBC_FLAG 1
#define GZ_FIELD -1.
#define HZ_FIELD 0.
#define HX_FIELD 0.01
// Range of HX (we use -HX)
#define MIN_HX_FIELD HX_FIELD
#define MAX_HX_FIELD 2.01
#define SEP_HX_FIELD 0.02

#define DELTA_HZ -0.001

//--- Contents -----------------------------------------------------------------

double susceptibility(const vector<complex<double>>& psi) {
    int index, tot_states = psi.size();
    int sites = int(log2(tot_states));
    double mag_z, chi;
    vector<double> sigmaZ(sites, 0.);
    vector<vector<int>> basis;

    // Store computational basis
    basis.resize(tot_states, vector<int>(sites, 0.));
    for (int n = 0; n < tot_states; n++) {
        index = n;
        for (int i = 0; i < sites; i++) {
            basis[n][sites - i - 1] = index % 2;
            index = index / 2;
        }
    }

    // Compute average of sigmaZ over psi
    for (int i = 0; i < sites; i++) {
        for (int n = 0; n < tot_states; n++) {
            // Compute sigmaZ_i
            if (basis[n][i] == 1) sigmaZ[i] += norm(psi[n]);
            if (basis[n][i] == 0) sigmaZ[i] -= norm(psi[n]);
        }
    }

    mag_z = 0.;
    for (int i = 0; i < sites; i++) mag_z += sigmaZ[i];
    mag_z = mag_z / sites;
    chi = mag_z / DELTA_HZ;

    return chi;
}

vector<double> magnetization(const vector<complex<double>>& psi) {

    int index, tot_states = psi.size();
    int sites = int(log2(tot_states));
    double mag_z_tilde, mag_z, mag_x;
    complex<double> val_c;
    vector<double> sigmaX(sites, 0.), sigmaZ(sites, 0.), output;
    vector<vector<int>> basis;

    // Store computational basis
    basis.resize(tot_states, vector<int>(sites, 0.));
    for (int n = 0; n < tot_states; n++) {
        index = n;
        for (int i = 0; i < sites; i++) {
            basis[n][sites - i - 1] = index % 2;
            index = index / 2;
        }
    }

    // Compute Symmetry-broken magnetization
    mag_z_tilde = 0.;
    for (int n = 0; n < tot_states; n++) {
        index = 0;
        for (int i = 0; i < sites; i++) {
            if (basis[n][i] == 1) index += 1;
            if (basis[n][i] == 0) index -= 1;
        }
        mag_z_tilde += abs(1.0 * index) * norm(psi[n]);
    }
    mag_z_tilde = mag_z_tilde / sites;
    output.push_back(mag_z_tilde);
    cout << "Symmetry-broken magnet along Z: " << mag_z_tilde << endl << endl;

    // Compute average of sigmaX and sigmaZ over psi
    for (int i = 0; i < sites; i++) {
        val_c = 0.;
        for (int n = 0; n < tot_states; n++) {
            // Compute sigmaZ_i
            if (basis[n][i] == 1) sigmaZ[i] += norm(psi[n]);
            if (basis[n][i] == 0) sigmaZ[i] -= norm(psi[n]);
            // Compute sigmaX_i
            if (basis[n][i] == 1) index = n - pow(2, sites - 1 - i);
            if (basis[n][i] == 0) index = n + pow(2, sites - 1 - i);
            val_c += conj(psi[n]) * psi[index];
        }

        if (abs(imag(val_c)) > 1.0e-10) cerr << "Non-real X magnetiz!" << endl;
        sigmaX[i] = real(val_c);
    }

    mag_z = 0.;
    mag_x = 0.;
    for (int i = 0; i < sites; i++) {
        mag_z += sigmaZ[i];
        mag_x += sigmaX[i];
    }
    mag_z = mag_z / sites;
    mag_x = mag_x / sites;

    output.push_back(mag_z);
    output.push_back(mag_x);
    return output;
}

void run_simulation(HamiltParameters param){
    /* Simulation for a given side */

    ofstream file;
    string directory, file_name, message;
    double ener_gs, ener_gap, chi;
    vector<double> magZZX;
    vector<complex<double>> psi_gs;
    hamiltonian HamOp(param);

    // Define path data directory
    directory = "Data_simulations/";
    // Define name file last configuration of the lattice
    file_name = "side_" + to_string(param.num_sites) + ".dat";

    file.open(directory + file_name);
    // Update transverse field and take measures
    for(double hx = MIN_HX_FIELD; hx < MAX_HX_FIELD; hx += SEP_HX_FIELD){
        // Store hx field value
        file << hx << " ";
        // Set field and diagonalize
        param.hx_field = -hx;
        param.hz_field = HZ_FIELD;
        HamOp.set_fields(param);
        HamOp.diagonalize();
        // Get and store GS energy and gaps
        cout << "Side " << param.num_sites << " - Field " << hx << endl;
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
        psi_gs = HamOp.get_eigenstate(0);
        magZZX = magnetization(psi_gs);
        file << magZZX[0] << " ";
        file << magZZX[1] << " ";
        file << magZZX[2] << " ";
        // Set field and diagonalize to find chi
        param.hz_field = DELTA_HZ;
        HamOp.set_fields(param);
        HamOp.diagonalize();
        psi_gs = HamOp.get_eigenstate(0);
        chi = susceptibility(psi_gs);
        file << chi << endl;
    }

    file.close();
}

//--- Main sim -----------------------------------------------------------------

int main(){
    /* Test the methods of the hamiltonian Class. */

    HamiltParameters param;
    param.sparse_flag = SPARSE_FLAG;
    param.pbc_flag = PBC_FLAG;
    param.num_sites = MIN_SIDE;
    param.gz_field = GZ_FIELD;
    param.hz_field = HZ_FIELD;
    param.hx_field = HX_FIELD;

    for (int side = MIN_SIDE; side < MAX_SIDE; side++){
        param.num_sites = side;
        if (side == 10) param.sparse_flag = 3;
        run_simulation(param);
    }

}
