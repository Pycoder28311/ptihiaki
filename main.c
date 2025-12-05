#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "head.h"
#include "globals.h"

// ---------- Helper Functions ----------

// Function to find and print roots
double find_and_print_roots(double xmin, double xmax, double dx, double tol, double tol_zero, int max_roots, int n_max) {
    double *roots = malloc(max_roots * sizeof(double));
    if (!roots) {
        fprintf(stderr, "Error: Could not allocate memory for roots.\n");
        return 0.0 / 0.0;
    }

    int n = find_roots(xmin, xmax, dx, tol, tol_zero, roots, max_roots, n_max);

    printf("Found %d root(s):\n", n);
    for (int i = 0; i < n; ++i) {
        printf("  root[%d] = %.12f  f(root)=%.12e\n", i, roots[i], f(roots[i], n_max));
    }

    double first_root = (n > 0) ? roots[0] : 0.0 / 0.0; // NaN if no root
    free(roots);
    return first_root;
}

// Function to create matrix, compute determinant, and free memory
void compute_matrix_determinant(int n_max, double k_perp) {
    int N_mtx = 2 * n_max + 1;

    // --- CREATE MATRIX ---
    double complex **M = create_matrix(N_mtx, k_perp);
    if (!M) {
        fprintf(stderr, "Error: Could not allocate matrix.\n");
        return;
    }

    printf("Matrix values (row x column):\n");
    for (int i = 0; i < N_mtx; i++) {
        for (int j = 0; j < N_mtx; j++) {
            printf("(% .6lf + % .6lfi)  ", creal(M[i][j]), cimag(M[i][j]));
        }
        printf("\n");
    }

    // Convert to pointer array for determinant function
    double complex *M_ptrs[N_mtx];
    for (int i = 0; i < N_mtx; i++)
        M_ptrs[i] = M[i];

    // Compute determinant
    double complex det = determinant_complex(M_ptrs, N_mtx);
    printf("Determinant for check = %.10lf + %.10lf i\n", creal(det), cimag(det));

    // Free memory
    for (int i = 0; i < N_mtx; i++)
        free(M[i]);
    free(M);
}


#define MAX_ITER 10000
#define TOL 1e-10

// Multiply matrix by vector: y = A * x
void mat_vec_mul(int N, double **A, double *x, double *y) {
    for (int i = 0; i < N; i++) {
        y[i] = 0.0;
        for (int j = 0; j < N; j++)
            y[i] += A[i][j] * x[j];
    }
}
// Compute vector norm
double norm(int N, double *v) {
    double s = 0.0;
    for (int i = 0; i < N; i++)
        s += v[i] * v[i];
    return sqrt(s);
}
// Power iteration for dominant eigenvalue/vector
void power_iteration(int N, double **A, double *eigvec, double *eigval) {
    double *b = malloc(N * sizeof(double));
    double *b_new = malloc(N * sizeof(double));

    // initial guess
    for (int i = 0; i < N; i++)
        b[i] = 1.0;

    for (int iter = 0; iter < MAX_ITER; iter++) {
        mat_vec_mul(N, A, b, b_new);
        double b_new_norm = norm(N, b_new);

        for (int i = 0; i < N; i++)
            b_new[i] /= b_new_norm;

        double diff = 0.0;
        for (int i = 0; i < N; i++)
            diff += fabs(b_new[i] - b[i]);

        if (diff < TOL) break;

        for (int i = 0; i < N; i++)
            b[i] = b_new[i];
    }

    // Rayleigh quotient for eigenvalue
    double numerator = 0.0, denominator = 0.0;
    double *Ab = malloc(N * sizeof(double));
    mat_vec_mul(N, A, b, Ab);
    for (int i = 0; i < N; i++) {
        numerator += b[i] * Ab[i];
        denominator += b[i] * b[i];
        eigvec[i] = b[i];
    }
    *eigval = numerator / denominator;

    free(b);
    free(b_new);
    free(Ab);
}
// Subtract outer product λ*v*v^T from matrix for deflation
void deflate_matrix(int N, double **A, double *eigvec, double eigval) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] -= eigval * eigvec[i] * eigvec[j];
}
// Find all eigenvalues and eigenvectors
void find_all_eigenvalues(int N, double **A) {
    // Copy of A to modify during deflation
    double **B = malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        B[i] = malloc(N * sizeof(double));
        for (int j = 0; j < N; j++)
            B[i][j] = A[i][j];
    }

    double *eigvec = malloc(N * sizeof(double));
    double eigval;

    printf("Eigenvalues and eigenvectors:\n");
    for (int k = 0; k < N; k++) {
        power_iteration(N, B, eigvec, &eigval);

        printf("Eigenvalue %d: %.6lf\n", k+1, eigval);
        printf("Eigenvector %d: [ ", k+1);
        for (int i = 0; i < N; i++)
            printf("%.6lf ", eigvec[i]);
        printf("]\n\n");

        deflate_matrix(N, B, eigvec, eigval);
    }

    for (int i = 0; i < N; i++)
        free(B[i]);
    free(B);
    free(eigvec);
}



// ---------- Main Function ----------
int main() {
    // ----- Root Finding -----
    int max_roots = 1;
    int n_max = 1;

    double xmin = 0.0, xmax = 5.0; // search interval
    double dx = 0.25;         // scanning step (smaller -> more robust, slower)
    double tol = 1e-8;        // brent tolerance
    double tol_zero = 1e-12;  // consider f(x) ~ 0

    double k_perp = find_and_print_roots(xmin, xmax, dx, tol, tol_zero, max_roots, n_max);
    printf("\n");

    if (isnan(k_perp)) {
        fprintf(stderr, "No roots found! Exiting.\n");
        return 1;
    }

    printf("Using k_perp = %.12f\n", k_perp);
    compute_matrix_determinant(n_max, k_perp);
    printf("\n");
    


    int N = 2;
    double **A = malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++)
        A[i] = malloc(N * sizeof(double));

    A[0][0] = 4; A[0][1] = 1;
    A[1][0] = 2; A[1][1] = 3;

    find_all_eigenvalues(N, A);

    for (int i = 0; i < N; i++)
        free(A[i]);
    free(A);

    return 0;
}
