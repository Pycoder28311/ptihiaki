#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "globals.h"
#include "utils.h"

#define NM (2*n_max + 1)

static inline int IDX(int n) {
    return n + n_max;   // Χάρτης από [-n_max, +n_max] σε [0, 2*n_max]
}

void build_Z_matrix(double matrix[NM][NM], double k_perp) {
    for (int i = 0; i < NM; i++) {       
        for (int j = 0; j < NM; j++) {   
            int n = i - n_max;  
            int q = j - n_max;  

            double Z_n  = Zn(phi_c, k_perp, n, a_radius);
            double Z_nq = Znq(n, q, k_perp, a_radius, phi_c, phi_i, lmax);

            matrix[i][j] = Z_nq - delta(n, q) * Z_n;
        }
    }
}

double compute_det_matrix(int N, double mat[N][N]) {
    double A[N][N];
    int i, j, k;
    double det = 1.0;

    // Αντιγραφή πίνακα
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            A[i][j] = mat[i][j];

    // LU αποσύνθεση
    for (k = 0; k < N; k++) {
        if (fabs(A[k][k]) < 1e-12) return 0.0; 
        for (i = k + 1; i < N; i++) {
            double factor = A[i][k] / A[k][k];
            for (j = k; j < N; j++)
                A[i][j] -= factor * A[k][j];
        }
    }

    // Υπολογισμός ορίζουσας
    for (i = 0; i < N; i++)
        det *= A[i][i];
    return det;
}

double compute_det_for_k(double k_perp) {
    double matrix[NM][NM];
    build_Z_matrix(matrix, k_perp);

    return compute_det_matrix(NM, matrix);
}

double** create_mat_reduced(double k_perp) {
    double mat[NM][NM];
    int i, j;

    build_Z_matrix(mat, k_perp);

    /*printf("Initial matrix (full):\n");
    for (int n = -n_max; n <= n_max; n++) {
        for (int q = -n_max; q <= n_max; q++) {
            printf("%lf\t", mat[IDX(n)][IDX(q)]);
        }
        printf("\n");
    }*/

    // Πίνακα b (τελευταία στήλη χωρίς την πρώτη τιμή)
    double *b = malloc((NM-1) * sizeof(double));
    for (int i = 1; i <= 2*n_max; i++) { // παράβλεψη πρώτης σειράς
        b[i-1] = - mat[i][0]; // τελευταία στήλη
    }
    
    printf("\nVector b:\n");
    for (i = 0; i < NM-1; i++) {
        printf("%lf\n", b[i]);
    }
    // Δημιουργία μειωμένου πίνακα (παράβλεψη πρώτης σειράς και τελευταίας στήλης)
    double **reduced = malloc((NM-1) * sizeof(double*));
    for (i = 0; i < NM-1; i++)
        reduced[i] = malloc((NM-1) * sizeof(double));

    // Αντιγραφή δεδομένων στον μειωμένο πίνακα
    for (i = 1; i < NM; i++) {        
        for (j = 1; j < NM; j++) {  
            reduced[i-1][j-1] = mat[i][j];
        }
    }

    // ΕΚτυπωση μειωμένου πίνακα
    printf("\nReduced matrix:\n");
    for (i = 0; i < NM - 1; i++) {
        for (j = 0; j < NM - 1; j++) {
            printf("%lf\t", reduced[i][j]);
        }
        printf("\n");
    }

    // Λύση συστήματος Ax = b με μέθοδο Gaussian elimination
    double *x = malloc((NM-1) * sizeof(double));

    double *b_copy = malloc((NM-1) * sizeof(double));
    for (i = 0; i < NM-1; i++) 
        b_copy[i] = b[i];

    for (i = 0; i < NM-1; i++) {
        int max_row = i;
        for (int k = i+1; k < NM-1; k++)
            if (fabs(reduced[k][i]) > fabs(reduced[max_row][i]))
                max_row = k;

        double *temp_row = reduced[i];
        reduced[i] = reduced[max_row];
        reduced[max_row] = temp_row;

        double temp_b = b_copy[i];
        b_copy[i] = b_copy[max_row];
        b_copy[max_row] = temp_b;

        for (int k = i+1; k < NM-1; k++) {
            double factor = reduced[k][i] / reduced[i][i];
            for (j = i; j < NM-1; j++)
                reduced[k][j] -= factor * reduced[i][j];
            b_copy[k] -= factor * b_copy[i];
        }
    }

    // Back-substitution
    for (i = NM-2; i >= 0; i--) {
        x[i] = b_copy[i];
        for (j = i+1; j < NM-1; j++)
            x[i] -= reduced[i][j] * x[j];
        x[i] /= reduced[i][i];
    }

    // Find the maximum absolute value
    double max_val = fabs(x[0]);
    for (i = 1; i < NM-1; i++) {
        if (fabs(x[i]) > max_val)
            max_val = fabs(x[i]);
    }

    // Normalize so that the largest becomes 1
    for (i = 0; i < NM-1; i++) {
        x[i] /= max_val;
    }

    // Print normalized solution
    printf("\nNormalized solution x:\n");
    printf("%lf\n", 1.0 / max_val);
    for (i = 0; i < NM-1; i++)
        printf("%lf\n", x[i]);
    printf("\n");


    free(b);
    free(b_copy);
    free(x);

    return reduced;
}