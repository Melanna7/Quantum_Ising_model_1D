/*******************************************************************************
*
* Test program for the lapacke library
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <stdlib.h>
#include <stdio.h>

#include <lapacke.h>

/*******************************************************************************
* PARAMETERS OF THE SIMULATION
*
* SIDE_SEP = separation between the sides of different simulations.
*
* BETA_SEP = separation between the betas of different simulations.
*
*******************************************************************************/

#define N 4
#define LDA N
#define lint lapack_int
#define ldcmplex lapack_complex_double

typedef struct double2 {
  double v[2];
} double2_t;

//Auxiliary routines prototypes
extern void print_matrix( char* desc, lint m, lint n, ldcmplex* a, lint lda);
extern void print_rmatrix( char* desc, lint m, lint n, double* a, lint lda);
extern void set_matrix(lint n, ldcmplex* a, lint lda, double2_t *a2);

//--- Contents -----------------------------------------------------------------

//Auxiliary routine: printing a matrix
void print_matrix( char* desc, lint m, lint n, ldcmplex* a, lint lda ) {
        lint i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", creal(a[i*lda+j]), cimag(a[i*lda+j]) );
                printf( "\n" );
        }
}

//Auxiliary routine: printing a real matrix
void print_rmatrix( char* desc, lint m, lint n, double* a, lint lda ) {
        lint i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}

//Auxiliary routine: set a complex matrix from a double[2] type matrix
void set_matrix(lint n, lapack_complex_double* a, lint lda, double2_t *a2) {
        lapack_int i, j;

        for( i = 0; i < n; i++ ) {
                for( j = 0; j < n; j++ )
                        a[i*lda+j] = lapack_make_complex_double(a2[i*lda+j].v[0], a2[i*lda+j].v[1]);
        }
}

//--- Main ---------------------------------------------------------------------

int main()
  {
    //Locals
    lint n = N, lda = LDA, info;

    //Local arrays
    double wr[N];
    ldcmplex ah[LDA*N] = {0};
    double2_t ah2[LDA*N] = {
           { 9.14,  0.00}, { 0.00,  0.00}, { 0.00,  0.00}, { 0.00,  0.00},
           {-4.37,  9.22}, {-3.35,  0.00}, { 0.00,  0.00}, { 0.00,  0.00},
           {-1.98,  1.72}, { 2.25,  9.51}, {-4.82,  0.00}, { 0.00,  0.00},
           {-8.96,  9.50}, { 2.57, -2.40}, {-3.24, -2.04}, { 8.44,  0.00}
    };

    //Executable statements
    set_matrix(n, ah, lda, ah2);
    printf( "LAPACKE_zheev (row-major, high-level) Example Program Results\n" )    ;

    //Print martix
    print_matrix( "Input Matrix", n, n, ah, lda );

    //Solve eigenproblem
    info = LAPACKE_zheev( LAPACK_ROW_MAJOR, 'V', 'L', n, ah, lda, wr );

    //Check for convergence
    if( info > 0 ) {
            printf( "zheev algorithm failed to compute eigenvalues.\n" );
            exit( 1 );
    }

    //Print eigenvalues
    print_rmatrix( "zheev Eigenvalues", 1, n, wr, 1 );

    //Print eigenvectors
    print_matrix( "Eigenvectors (stored columnwise)", n, n, ah, lda );

    //Local arrays
    ldcmplex wc[N];
    ldcmplex ag[LDA*N] = {0};
    ldcmplex wgr[LDA*N] = {0};
    double2_t ag2[LDA*N] = {
      { 9.14,  0.00}, {-4.37, -9.22}, {-1.98, -1.72}, {-8.96, -9.50},
      {-4.37,  9.22}, {-3.35,  0.00}, { 2.25, -9.51}, { 2.57,  2.40},
      {-1.98,  1.72}, { 2.25,  9.51}, {-4.82,  0.00}, {-3.24,  2.04},
      {-8.96,  9.50}, { 2.57, -2.40}, {-3.24, -2.04}, { 8.44,  0.00},
    };

    printf("\n\n");

    //Executable statements
    set_matrix(n, ag, lda, ag2);
    printf( "LAPACKE_zgeev (row-major, high-level) Example Program Results\n" );

    //Print martix
    print_matrix( "Input Matrix", n, n, ag, lda );

    //Solve eigenproblem
    info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, ag, lda, wc, 0, lda, wgr, lda);

    //Check for convergence
    if( info > 0 ) {
            printf( "zgeev algorithm failed to compute eigenvalues.\n" );
            exit( 1 );
    }

    //Print eigenvalues
    print_matrix( "zgeev Eigenvalues", 1, n, wc, 1);

    //Print eigenvectors
    print_matrix( "Eigenvectors (stored columnwise)", n, n, wgr, lda );
    exit( 0 );
  }
