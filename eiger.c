#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "head.h"

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
void find_all_eigenvalues(int N, double complex **A, FILE *fp) {
    // Copy of A to modify during deflation
    double **B = malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        B[i] = malloc(N * sizeof(double));
        for (int j = 0; j < N; j++)
            B[i][j] = A[i][j];
    }

    double *eigvec = malloc(N * sizeof(double));
    double eigval;

    fprintf(fp, "\nEIGENVALUES AND EIGENVECTORS:\n");
    for (int k = 0; k < N; k++) {
        power_iteration(N, B, eigvec, &eigval);

        fprintf(fp, "Eigenvalue %d: %.6lf\n", k+1, eigval);
        fprintf(fp, "Eigenvector %d: [ ", k+1);
        for (int i = 0; i < N; i++)
            fprintf(fp, "%.6lf ", eigvec[i]);
        fprintf(fp, "]\n\n");

        deflate_matrix(N, B, eigvec, eigval);
    }

    for (int i = 0; i < N; i++)
        free(B[i]);
    free(B);
    free(eigvec);
}
void compute_matrix_determinant(int n_max, double k_perp, FILE *fp) {
    int N_mtx = 2 * n_max + 1;
    int sub_N = N_mtx - 1;

    // --- CREATE MATRIX ---
    double complex **M_full = create_matrix(N_mtx, k_perp);
    if (!M_full) {
        fprintf(stderr, "Error: Could not allocate matrix.\n");
        return;
    }
    double complex *M_ptrs[N_mtx];
    for (int i = 0; i < N_mtx; i++)
        M_ptrs[i] = M_full[i];

    // Compute determinant
    double complex det = determinant_complex(M_ptrs, N_mtx);
    fprintf(fp, "Determinant for check = %.10lf + %.10lf i\n", creal(det), cimag(det));



    // --- CREATE SUBMATRIX (excluding last row and column) ---
    double complex **M = malloc(sub_N * sizeof(double complex *));
    for (int i = 0; i < sub_N; i++) {
        M[i] = malloc(sub_N * sizeof(double complex));
        for (int j = 0; j < sub_N; j++)
            M[i][j] = M_full[i][j];  // copy only top-left submatrix
    }

    fprintf(fp, "Matrix values (row x column):\n");
    for (int i = 0; i < sub_N; i++) {
        for (int j = 0; j < sub_N; j++) {
            fprintf(fp, "(% .6lf + % .6lfi)  ", creal(M[i][j]), cimag(M[i][j]));
        }
        fprintf(fp, "\n");
    }

    find_all_eigenvalues(sub_N, M_ptrs, fp);

    // Free memory
    for (int i = 0; i < sub_N; i++)
        free(M[i]);
    free(M);
}