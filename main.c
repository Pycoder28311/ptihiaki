#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "head.h"
#include "globals.h"

int bessik01a(double x,
              double *i0, double *i1,
              double *k0, double *k1,
              double *i0p, double *i1p,
              double *k0p, double *k1p);

int bessik01b(double x,
              double *i0, double *i1,
              double *k0, double *k1,
              double *i0p, double *i1p,
              double *k0p, double *k1p);

int bessikna(int n, double x, int *nm,
             double *in, double *kn,
             double *inp, double *knp);

int bessiknb(int n, double x, int *nm,
             double *in, double *kn,
             double *inp, double *knp);

// ---------- Main Function ----------
int main() {

    double x = 10.0;

    /* ===============================
       1) bessik01a
       =============================== */
    double i0,i1,k0,k1,i0p,i1p,k0p,k1p;

    bessik01a(x,
              &i0,&i1,
              &k0,&k1,
              &i0p,&i1p,
              &k0p,&k1p);

    printf("bessik01a (x = %.2f)\n", x);
    printf("I0=%.10e  I1=%.10e\n", i0, i1);
    printf("K0=%.10e  K1=%.10e\n", k0, k1);
    printf("I0'=%g  I1'=%g\n\n", i0p, i1p);

    bessik01b(x,
              &i0,&i1,
              &k0,&k1,
              &i0p,&i1p,
              &k0p,&k1p);

    printf("bessik01b (x = %.2f)\n", x);
    printf("I0=%.10e  I1=%.10e\n", i0, i1);
    printf("K0=%.10e  K1=%.10e\n\n", k0, k1);

    int n = 5;
    int nm;
    double in[6], kn[6], inp[6], knp[6];

    bessikna(n, x, &nm, in, kn, inp, knp);

    printf("bessikna (n = %d, x = %.2f)\n", n, x);
    for (int k = 0; k <= nm; k++) {
        printf("k=%d  I=%.6e  K=%.6e\n", k, in[k], kn[k]);
    }
    printf("\n");

    /* ===============================
       4) bessiknb
       =============================== */
    bessiknb(n, x, &nm, in, kn, inp, knp);

    printf("bessiknb (n = %d, x = %.2f)\n", n, x);
    for (int k = 0; k <= nm; k++) {
        printf("k=%d  I=%.6e  K=%.6e\n", k, in[k], kn[k]);
    }

    
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
        fprintf(fp, "\nUSING K_PERP = %.12f. ITERATION %d\n\n", k_perp[i], i+1);
        compute_matrix_determinant(n_max, k_perp[i], fp);
        double k = k_perp[i]; // ο κυματικός αριθμός που βρήκες
        printf("Found root k_perp = %.12f\n", k_perp[i]);

        double omega = v * k;
        double f = omega / (2.0 * M_PI);
        printf("Frequency f = %.12f Hz\n", f);
        double lambda_m = 2.0 * M_PI / k;  
        printf("Wavelength λ = %.12f m\n", lambda_m);
    }

    if (N_roots == 0) {
        fprintf(stderr, "No roots found! Exiting.\n");
        return 1;
    }

    return 0;
}
