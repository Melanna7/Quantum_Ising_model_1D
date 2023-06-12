#include <complex>
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

vector<vector<double>> magnetization(const vector<complex<double>>& psi) {

    int index;
    int tot_states = psi.size();
    int sites = int(log2(tot_states));
    double val = 0.;
    complex<double> val_c(0., 0.);
    vector<double> sigmaX(sites, 0.);
    vector<double> sigmaZ(sites, 0.);
    vector<vector<int>> basis;
    vector<vector<double>> output;

    // Store computational basis
    basis.resize(tot_states, vector<int>(sites, 0.));
    for (int n = 0; n < tot_states; n++) {
        index = n;
        for (int i = 0; i < sites; i++) {
            basis[n][sites - i - 1 ] = index % 2;
            index = index / 2;
        }
    }

    // Compute Symmetry-broken magnetization
    index = 0;
    for (int n = 0; n < tot_states; n++) {
        for (int i = 0; i < sites; i++) {
            if (basis[n][i] == 1) index += 1;
            if (basis[n][i] == 0) index -= 1;
        }
        val += abs(1.0 * index) * norm(psi[n]);
    }
    cout << "Symmetry-broken magnetization along Z: " << val / sites << endl;

    // Compute average of sigmaX and sigmaZ over psi
    for (int i = 0; i < sites; i++) {
        for (int n = 0; n < tot_states; n++) {
            // Compute sigmaZ_i
            if (basis[n][i] == 1) sigmaZ[i] += norm(psi[n]);
            if (basis[n][i] == 0) sigmaZ[i] -= norm(psi[n]);
            // Compute sigmaX_i
            if (basis[n][i] == 1) index = n - pow(2, sites - 1 - i);
            if (basis[n][i] == 0) index = n + pow(2, sites - 1 - i);
            val_c += conj(psi[n]) * psi[index];
        }

        if (abs(imag(val_c)) > 1.0e-10) cout << "Err: non-real mag" << endl;
        sigmaX[i] = real(val_c);
    }

    output.push_back(sigmaZ);
    output.push_back(sigmaX);
    return output;
}
