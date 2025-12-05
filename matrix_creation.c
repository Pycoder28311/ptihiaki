#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "head.h"
#include "globals.h"

// Υπολογισμός ορίζουσας πραγματικού πίνακα
void getCofactor(double **mat, double **temp, int p, int q, int n) {
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                temp[i][j++] = mat[row][col];
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

double determinant(double **mat, int n) {
    if (n == 1) return mat[0][0];

    double det = 0.0;
    double **temp = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        temp[i] = (double *)malloc(n * sizeof(double));

    int sign = 1;

    for (int f = 0; f < n; f++) {
        getCofactor(mat, temp, 0, f, n);
        det += sign * mat[0][f] * determinant(temp, n - 1);
        sign = -sign;
    }

    for (int i = 0; i < n; i++) free(temp[i]);
    free(temp);

    return det;
}

// Υπολογισμός ορίζουσας μιγαδικού πίνακα με Gauss elimination
double complex determinant_complex(double complex **mat, int N_mtx) {
    double complex **A = malloc(N_mtx * sizeof(double complex *));
    for (int i = 0; i < N_mtx; i++) {
        A[i] = malloc(N_mtx * sizeof(double complex));
        for (int j = 0; j < N_mtx; j++)
            A[i][j] = mat[i][j];
    }

    double complex det = 1.0 + 0.0*I;

    for (int i = 0; i < N_mtx; i++) {
        int pivot = i;
        for (int j = i+1; j < N_mtx; j++)
            if (cabs(A[j][i]) > cabs(A[pivot][i]))
                pivot = j;

        if (cabs(A[pivot][i]) < 1e-12) {
            det = 0.0 + 0.0*I;
            break;
        }

        if (pivot != i) {
            double complex *tmp = A[i];
            A[i] = A[pivot];
            A[pivot] = tmp;
            det *= -1.0;
        }

        det *= A[i][i];

        for (int j = i+1; j < N_mtx; j++) {
            double complex factor = A[j][i] / A[i][i];
            for (int k = i; k < N_mtx; k++)
                A[j][k] -= factor * A[i][k];
        }
    }

    for (int i = 0; i < N_mtx; i++)
        free(A[i]);
    free(A);

    return det;
}

double complex **create_matrix(int N_mtx, double k_perp) {
    double complex **M = malloc(N_mtx * sizeof(double complex *));
    if (!M) {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(1);
    }

    int m = 0;

    for (int i = 0; i < N_mtx; i++)
        M[i] = malloc(N_mtx * sizeof(double complex));
    for (int n = 0; n < N_mtx; n++) {

        for (int q = 0; q < N_mtx; q++) {
        
            double complex Znq_val = Z_nq(q, n, m, k_perp);
            double complex Zn_val = Z_diff(q, n, m, k_perp);

            if (n == q)
                M[n][q] = Znq_val - Zn_val;
            else
                M[n][q] = Znq_val;
        }
    }

    return M;
}

