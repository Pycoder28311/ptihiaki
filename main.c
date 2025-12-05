#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "head.h"
#include "globals.h"

// ---------- Main Function ----------
int main() {
    // ----- Root Finding -----
    int max_roots = 50;
    int n_max = 1;

    double xmin = 0.0, xmax = 20.0; // search interval
    double dx = 0.1;         // scanning step (smaller -> more robust, slower)
    double tol = 1e-8;        // brent tolerance
    double tol_zero = 1e-12;  // consider f(x) ~ 0
    double k_perp[max_roots]; 

    FILE *fp = fopen("output.txt", "w");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file for writing.\n");
        return 1;
    }

    double N_roots = find_roots(xmin, xmax, dx, tol, tol_zero, k_perp, max_roots, n_max);

    for (int i = 0; i < N_roots; i++) {
        fprintf(fp, "\nUSING K_PERP = %.12f. ITERATION %d\n\n", k_perp[i], i+1);
        compute_matrix_determinant(n_max, k_perp[i], fp);
    }

    if (N_roots == 0) {
        fprintf(stderr, "No roots found! Exiting.\n");
        return 1;
    }

    return 0;
}
