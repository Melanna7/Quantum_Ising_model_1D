#include <complex>
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

vector<double> magnetization(const vector<complex<double>>& psi) {

    int index;
    int tot_states = psi.size();
    int sites = int(log2(tot_states));
    double mag_z_tilde, mag_z, mag_x;
    complex<double> val_c;
    vector<double> sigmaX(sites, 0.);
    vector<double> sigmaZ(sites, 0.);
    vector<vector<int>> basis;
    vector<double> output;

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
    cout << "Symmetry-broken magnetization along Z: " << mag_z_tilde << endl;

    // Compute average of sigmaX and sigmaZ over psi
    for (int i = 0; i < sites; i++) {
        for (int n = 0; n < tot_states; n++) {
            val_c = 0.;
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
