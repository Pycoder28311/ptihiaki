#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "head.h"
#include "globals.h"

// ---------- Helper Functions ----------

// Function to find and print roots
void find_and_print_roots(double xmin, double xmax, double dx, double tol, double tol_zero, int max_roots, int n_max) {
    double *roots = malloc(max_roots * sizeof(double));
    if (!roots) {
        fprintf(stderr, "Error: Could not allocate memory for roots.\n");
        return;
    }

    int n = find_roots(xmin, xmax, dx, tol, tol_zero, roots, max_roots, n_max);

    printf("Found %d root(s):\n", n);
    for (int i = 0; i < n; ++i) {
        printf("  root[%d] = %.12f  f(root)=%.12e\n", i, roots[i], f(roots[i], n_max));
    }

    free(roots);
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

    // Convert to pointer array for determinant function
    double complex *M_ptrs[N_mtx];
    for (int i = 0; i < N_mtx; i++)
        M_ptrs[i] = M[i];

    // Compute determinant
    double complex det = determinant_complex(M_ptrs, N_mtx);
    printf("Determinant = %.10lf + %.10lf i\n", creal(det), cimag(det));

    // Free memory
    for (int i = 0; i < N_mtx; i++)
        free(M[i]);
    free(M);
}

// ---------- Main Function ----------
int main() {
    // ----- Root Finding -----
    double xmin = -1.0, xmax = 5.0;
    double dx = 0.25;         // scanning step (smaller -> more robust, slower)
    double tol = 1e-8;        // brent tolerance
    double tol_zero = 1e-12;  // consider f(x) ~ 0
    int max_roots = 20;
    int n_max = 1;

    find_and_print_roots(xmin, xmax, dx, tol, tol_zero, max_roots, n_max);

    // ----- Matrix & Determinant -----
    int n_max = 1;
    double k_perp = 2.0; // <-- replace with actual value if needed
    compute_matrix_determinant(n_max, k_perp);

    return 0;
}
