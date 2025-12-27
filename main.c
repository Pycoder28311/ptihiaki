#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
//#include "head.h"

#include "final/globals.h"
#include "final/utils.h"
#include "functions/root.h"

#define NM (2*n_max + 1)

// ---------- Main Function ----------
int main() {

    int nr;
    double roots_initial[100];

    // Example: scan for real k from 0.0 to 20.0
    riza(compute_det_for_k, down, up, step, tol_riza, &nr, roots_initial);

    int valid_nr = 0;
    double roots[100];

    for (int i = 0; i < nr; i++) {
        if (fabs(roots_initial[i]) > 1e-15) {  // skip k_perp = 0
            roots[valid_nr++] = roots_initial[i];
        }
    }

    printf("\nFound %d roots for k_perp\n\n", valid_nr);

    for (int i = 0; i < valid_nr; i++) { 
        double k_real = roots[i]; 
        double complex k_perp = k_real + 0.0*I; 
        printf("Root %d: %g", i+1, k_real);

        double matrix[NM][NM];

        build_Z_matrix(matrix, k_perp);

        // Optional: print the matrix
        /*printf("Matrix (Z_nq - delta_nq * Z_n):\n");
        for (int n = 0; n < N_MAX; n++) {
            for (int q = 0; q < N_MAX; q++) {
                printf("%10.6f ", matrix[n][q]);
            }
            printf("\n");
        }*/

        // Compute determinant
        double det = compute_det_matrix(NM, matrix);
        printf(" | Determinant: %g\n", det); 
    }

    if (nr > 0) {
        double k_real = roots[0];
        double complex k_perp = k_real + 0.0*I;
        printf("\nMatrix for first root (k_perp = %g):\n", k_real);

        double matrix[NM][NM];
        build_Z_matrix(matrix, k_perp);

        for (int n = 0; n < NM; n++) {
            for (int q = 0; q < NM; q++) {
                printf("%10.6f ", matrix[n][q]);
            }
            printf("\n");
        }

        create_mat_reduced(k_perp);
    }

    /*
        // ----- Root Finding -----
        int max_roots = 10;
        int n_max = 1;

        double xmin = 0.0, xmax = 10.0; // search interval
        double dx = 1.0;         // scanning step (smaller -> more robust, slower)
        double tol = 1e-10;        // brent tolerance
        double tol_zero = 1e-12;  // consider f(x) ~ 0
        double k_perp[max_roots]; 

        FILE *fp = fopen("output.txt", "w");
        if (!fp) {
            fprintf(stderr, "Error: Could not open file for writing.\n");
            return 1;
        }
        
        double N_roots = find_roots(xmin, xmax, dx, tol, tol_zero, k_perp, max_roots, n_max);

        for (int i = 0; i < N_roots; i++) {
            compute_matrix_determinant(n_max, k_perp[i], fp);
            printf("Found root k_perp = %.12f\n", k_perp[i]);
        }

        if (N_roots == 0) {
            fprintf(stderr, "No roots found! Exiting.\n");
            return 1;
        }
    */

    return 0;
}
