/*******************************************************************************
*
* Main program for the Heisenberg simulation
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "class_Eigen.h"

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

// Range sides
#define MIN_SIDE 4
#define MAX_SIDE 12
#define SEP_SIDE 2
// General settings
#define SPARSE_FLAG 3
#define PBC_FLAG 0
#define GX_FIELD 1.
#define GY_FIELD 1.
#define HZ_FIELD 0.0001
// Range of GZ
#define MIN_GZ_FIELD -10.0
#define MAX_GZ_FIELD 10.0
#define SEP_GZ_FIELD 0.05

vector<vector<int>> basis;

//--- Contents -----------------------------------------------------------------

vector<double> magnetization(const VectorXcd& psi) {

    int index, tot_states = psi.size();
    int sites = int(log2(tot_states));
    double mag_x, mag_y, mag_z, val_o;
    complex<double> val_x, val_y;
    vector<double> sigmaX(sites, 0.), sigmaY(sites, 0.), sigmaZ(sites, 0.);
    vector<double> output;

    // Compute average of sigmaX, sigmaY and sigmaZ over psi
    for (int i = 0; i < sites; i++) {
        val_x = 0.;
        val_y = 0.;
        for (int n = 0; n < tot_states; n++) {
            // Compute sigmaX_i
            if (basis[n][i] == 1) index = n - pow(2, sites - 1 - i);
            if (basis[n][i] == 0) index = n + pow(2, sites - 1 - i);
            val_x += conj(psi[n]) * psi[index];
            // Compute sigmaY_i
            if (basis[n][i] == 1) {
                index = n - pow(2, sites - 1 - i);
                val_y += conj(psi[n]) * complex(0., -1.) * psi[index];
            }
            if (basis[n][i] == 0) {
                index = n + pow(2, sites - 1 - i);
                val_y += conj(psi[n]) * complex(0., 1.) * psi[index];
            }
            // Compute sigmaZ_i
            if (basis[n][i] == 1) sigmaZ[i] += norm(psi[n]);
            if (basis[n][i] == 0) sigmaZ[i] -= norm(psi[n]);
        }

        if (abs(imag(val_x)) > 1.0e-10) cerr << "Non-real X magnetiz!" << endl;
        if (abs(imag(val_y)) > 1.0e-10) cerr << "Non-real Y magnetiz!" << endl;
        sigmaX[i] = real(val_x);
        sigmaY[i] = real(val_y);
    }

    val_o = 0.;
    mag_x = 0.;
    mag_y = 0.;
    mag_z = 0.;
    for (int i = 0; i < sites; i++) {
        mag_x += sigmaX[i];
        mag_y += sigmaY[i];
        mag_z += sigmaZ[i];
        val_o += pow(-1, i) * sigmaZ[i];
    }
    mag_x = mag_x / sites;
    mag_y = mag_y / sites;
    mag_z = mag_z / sites;
    val_o = val_o / sites;

    output.push_back(val_o);
    output.push_back(mag_x);
    output.push_back(mag_y);
    output.push_back(mag_z);
    return output;
}

void run_simulation(HamiltParameters param){
    /* Simulation for a given side */

    ofstream file;
    string directory, file_name, message;
    double ener_gs, ener_gap, chi;
    vector<double> magZXYZ;
    VectorXcd psi_gs;
    hamiltonian HamOp(param);

    // Define path data directory
    directory = "Data_Heise/";
    // Define name file last configuration of the lattice
    file_name = "side_" + to_string(param.num_sites) + ".dat";

    file.open(directory + file_name);
    // Update transverse field and take measures
    for(double gz = MIN_GZ_FIELD; gz < MAX_GZ_FIELD; gz += SEP_GZ_FIELD){
        // Store gx field value
        file << gz << " ";
        // Set field and diagonalize
        param.gz_field = gz;
        HamOp.set_fields(param);
        HamOp.diagonalize();
        // Get and store GS energy and gaps
        cout << "Side " << param.num_sites;
        cout << " - Field " << setprecision(4) << gz << endl;
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
        psi_gs = HamOp.get_eigenvector(0);
        magZXYZ = magnetization(psi_gs);
        cout << "Symmetry-broken parameter: " << magZXYZ[0] << endl;
        file << magZXYZ[0] << " ";
        file << magZXYZ[1] << " ";
        file << magZXYZ[2] << " ";
        file << magZXYZ[3] << " ";
        file << 0. << endl; // Eventually for susceptibility
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
    param.gx_field = GX_FIELD;
    param.gy_field = GY_FIELD;
    param.gz_field = MIN_GZ_FIELD;
    param.hz_field = HZ_FIELD;



    for (int side = MAX_SIDE; side >= MIN_SIDE; side -= SEP_SIDE){

        int index, tot_states = pow(2, side);

        // Store computational basis
        basis.resize(tot_states, vector<int>(side, 0.));
        for (int n = 0; n < tot_states; n++) {
            index = n;
            for (int i = 0; i < side; i++) {
                basis[n][side - i - 1] = index % 2;
                index = index / 2;
            }
        }

        // Start simulation
        cout << "Start new side" << endl;
        auto start = chrono::steady_clock::now();

        param.num_sites = side;
        run_simulation(param);

        auto end = chrono::steady_clock::now();
        chrono::duration<double> elapsed_sec = end - start;
        cout << "Elapsed time: " << elapsed_sec.count() << "s" << endl << endl;
    }

}
