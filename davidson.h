#include <complex>
#include <vector>
#include <cmath>

using namespace std;

void Davidson(int n, int Kmax, vector<double>& Eigenvalue, vector<complex<double>>& EigenVector) {
    bool Useguess = false;  // Flag for using a guess
    int kmax, jmax, jmin, maxstep, method, m, l, maxnmv, order, testspace, j, lwork, istate, ii;
    double tol, lock, targetEn, Norm, emin, etemp;
    vector<double> alpha, beta;  // Eigenvalue parameters
    vector<complex<double>> tmp, residu;  // Temporary variables
    vector<vector<complex<double>>> eivec, zwork;  // Eigenvalue and workspace variables

    // INIZIALIZATION OF PARAMETERS
    Useguess = false;  // Set Useguess flag to false
    kmax = Kmax;  // Set maximum number of eigenvalues
    targetEn = -5.0 * ell;  // Set target energy value
    tol = 1.0e-9;  // Set tolerance for eigensolutions
    maxnmv = 100;  // Set maximum number of matvecs in cgstab or gmres
    bool wanted = true;  // Flag to compute converged eigenvectors
    order = -1;  // Selection criterion for Ritz values
    if (order == 0)
        testspace = 3;  // Use 3 as the test space value if a reasonable target is known
    else
        testspace = 2;  // Otherwise, use 2 as the test space value

    if (3 * kmax <= 20)
        jmax = 20;  // Set maximum size of the search space
    else
        jmax = 3 * kmax;
    jmin = 2 * kmax;  // Set minimum size of the search space
    maxstep = 1000;  // Set maximum number of Jacobi-Davidson iterations
    lock = 1.0e-12;  // Set tracking parameter
    method = 2;  // Set method for linear equation solver (2: cgstab(l))
    m = 30;  // Set maximum dimension of search space for gmres(m)
    l = 2;  // Set degree of gmres-polynomial in bi-cgstab(l)
    if (method == 1)
        lwork = 4 + m + 5 * jmax + 3 * kmax;  // Calculate size of workspace for gmres(m)
    else if (method == 2)
        lwork = 10 + 6 * l + 5 * jmax + 3 * kmax;  // Calculate size of workspace for cgstab(l)

    // END OF INIZIALIZATION

    alpha.resize(jmax);  // Resize alpha vector
    beta.resize(jmax);  // Resize beta vector
    eivec.resize(n, vector<complex<double>>(Kmax));  // Resize eivec matrix
    zwork.resize(n, vector<complex<double>>(lwork));  // Resize zwork matrix

    for (int i = 0; i < jmax; i++) {
        alpha[i] = 0.0;  // Initialize alpha values to 0
        beta[i] = 0.0;  // Initialize beta values to 0
    }

    tmp.resize(n, 0.0);  // Initialize tmp vector to 0
    residu.resize(n, 0.0);  // Initialize residu vector to 0
    for (int i = 0; i < n; i++) {
        zwork[i].resize(lwork, 0.0);  // Initialize zwork matrix to 0
    }

    JDQZ(alpha, beta, eivec, wanted, n, targetEn, tol, kmax, jmax, jmin,
      method, m, l, maxnmv, maxstep, lock, order, testspace, zwork, lwork, Useguess);

    // Computes the norms of the residuals
    for (int j = 0; j < Kmax; j++) {
       AMUL(n, eivec[0][j], residu);
       ZSCAL(n, beta[j], residu, 1);
       BMUL(n, eivec[0][j], tmp);
       ZAXPY(n, -alpha[j], tmp, 1, residu, 1);
    }

    Eigenvalue.resize(Kmax);  // Resize Eigenvalue vector
    for (int i = 0; i < Kmax; i++) {
        Eigenvalue[i] = alpha[i] / beta[i];  // Calculate Eigenvalue values
    }

    emin = Eigenvalue[0];  // Calculate the smallest eigenvalue (ground state)
    istate = 0;
    for (int ii = 1; ii < Kmax; ii++) {
        if (Eigenvalue[ii] < emin) {
            emin = Eigenvalue[ii];
            istate = ii;
        }
    }

    if (istate != 0) {
        etemp = Eigenvalue[0];
        Eigenvalue[0] = Eigenvalue[istate];  // Swap Eigenvalue values
        Eigenvalue[istate] = etemp;
    }

    EigenVector = eivec[istate];  // Choose the eigenvector corresponding to the selected eigenvalue

    Norm = 0.0;
    for (int i = 0; i < n; i++) {
        Norm += conj(EigenVector[i]) * EigenVector[i];  // Calculate the norm of the eigenvector
    }

    Norm = sqrt(Norm);
    for (int i = 0; i < n; i++) {
        EigenVector[i] /= Norm;  // Normalize the eigenvector
    }
}
