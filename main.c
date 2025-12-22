#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "head.h"
#include "globals.h"

int bessjyna(int n, double x, int *nm, double *jn, double *yn, double *jnp, double *ynp);
int bessjynb(int n, double x, int *nm, double *jn, double *yn, double *jnp, double *ynp);
int bessikna(int n, double x, int *nm, double *in, double *kn, double *inp, double *knp);
int bessiknb(int n, double x, int *nm, double *in, double *kn, double *inp, double *knp);

/* -------------------------
   Struct to hold value and derivatives
   ------------------------- */
typedef struct {
    double value;  // function value
    double dF;     // derivative w.r.t rho
    double Fp;     // derivative normalized
} FvResult;

FvResult F_v(int v, double complex k_perp, double rho) {
    int nm;
    FvResult res;

    if (cimag(k_perp) == 0.0) {  // real -> J_v
        double jn[20], yn[20], jnp[20], ynp[20];
        bessjyna(v, creal(k_perp) * rho, &nm, jn, yn, jnp, ynp);
        res.value = jn[v];
        res.dF = creal(k_perp) * jnp[v];
        res.Fp = jnp[v];
    } else {  // imaginary -> I_v
        double in[20], kn[20], inp[20], knp[20];
        double s = cimag(k_perp);  // positive imaginary part
        bessikna(v, s * rho, &nm, in, kn, inp, knp);
        res.value = in[v];
        res.dF = s * inp[v];
        res.Fp = inp[v];
    }

    return res;
}

FvResult G_v(int v, double complex k_perp, double rho) {
    int nm;
    FvResult res;

    if (cimag(k_perp) == 0.0) {  // real -> Y_v
        double jn[20], yn[20], jnp[20], ynp[20];
        bessjyna(v, creal(k_perp) * rho, &nm, jn, yn, jnp, ynp);
        res.value = yn[v];
        res.dF = creal(k_perp) * ynp[v];
        res.Fp = ynp[v];
    } else {  // imaginary -> K_v
        double in[20], kn[20], inp[20], knp[20];
        double s = cimag(k_perp);
        bessikna(v, s * rho, &nm, in, kn, inp, knp);
        res.value = kn[v];
        res.dF = s * knp[v];
        res.Fp = knp[v];
    }

    return res;
}

double complex I_alpha_beta(double alpha, double beta, double phi1, double phi2) {
    double complex result = 0.0 + 0.0*I;  // use 'result', not 'I'

    if (alpha + beta != 0.0) {
        result += (cexp(-I*beta*phi1) * 
                  (cexp(I*(alpha+beta)*phi2) - cexp(I*(alpha+beta)*phi1))) / 
                  (2.0*I*(alpha+beta));
    } else {
        result += (cexp(-I*beta*phi1) * (phi2 - phi1)) / 2.0;
    }

    if (alpha - beta != 0.0) {
        result += (cexp(I*beta*phi1) * 
                  (cexp(I*(alpha-beta)*phi2) - cexp(I*(alpha-beta)*phi1))) / 
                  (2.0*I*(alpha-beta));
    } else {
        result += (cexp(I*beta*phi1) * (phi2 - phi1)) / 2.0;
    }

    return result;
}

int delta(int l, int m) {
    return (l == m) ? 1 : 0;
}

double k_l(int l, double phi_c) {
    return l * M_PI / phi_c;
}

double k_n(int m, int n, double N) {
    return m + n * N;
}

double O_kl(int kl, double complex k_perp, double rho) {
    FvResult F_rho = F_v(kl, k_perp, rho);
    FvResult G_rho = G_v(kl, k_perp, rho);
    FvResult F_D   = F_v(kl, k_perp, D);
    FvResult G_D   = G_v(kl, k_perp, D);

    return (F_rho.value * G_D.Fp - F_D.Fp * G_rho.value) / G_D.Fp;
}

// Compute derivative dO_kl/drho using precomputed Fp and Gp
double O_klp(int kl, double complex k_perp, double rho) {
    FvResult F_rho = F_v(kl, k_perp, rho);
    FvResult G_rho = G_v(kl, k_perp, rho);
    FvResult F_D   = F_v(kl, k_perp, D);
    FvResult G_D   = G_v(kl, k_perp, D);

    return F_rho.Fp - (F_D.Fp * G_rho.Fp / G_D.Fp);
}

double Znt(double phi_c, double complex k_perp, double Fp) {
    return phi_c * cabs(k_perp) * Fp;
}

double Z_nqt(double n, double q, double complex k_perp, double a, double phi_c, double phi_i, int lmax) {
    double S = 0.0;
    double m = 1.0;
    double N = 5.0;

    double kn = k_n(m, n, N);
    double kq = k_n(m, q, N);

    for (int l = 0; l <= lmax; l++) {
        double kl = k_l(l, phi_c);

        double O_kl_a = O_kl(l, k_perp, a);
        double O_kl_a_p = O_klp(l, k_perp, a);

        double complex I1 = I_alpha_beta(kq, kl, phi_i, phi_i + phi_c);
        double complex I2 = I_alpha_beta(-kn, kl, phi_i, phi_i + phi_c);

        double delta_l0 = delta(l,0);

        S += (O_kl_a_p / O_kl_a) * (I1 * O_kl_a_p * I2) / (1.0 + delta_l0);
    }

    FvResult Fq = F_v(q, k_perp, a); // F_{k_q}(k_perp*a)
    double Znq = 2.0 * Fq.Fp * S / phi_c;

    return Znq;
}

// ---------- Main Function ----------
int main() {

    double complex k_perp = 3.0 + 0.0*I;  // purely real
    const double a = 1.0;    // example radius
    const double D = 1.5;   

    int l = 3;
    double phi_c = 1.0;
    double phi_i = 0.0;
    int q = 1;
    int n = 1;
    int lmax = 10; // π.χ. sum μέχρι 10
    int N_MAX = 5; // μέγεθος πίνακα

    double matrix[N_MAX][N_MAX];

    for (int n = 0; n < N_MAX; n++) {
        for (int q = 0; q < N_MAX; q++) {
            double Znq = Z_nqt(n, q, k_perp, a, phi_c, phi_i, lmax);
            double Zn = Z_nqt(n, n, k_perp, a, phi_c, phi_i, lmax);
            matrix[n][q] = Znq - delta(n, q) * Zn;
        }
    }

    // Print matrix
    printf("Matrix (Z_nq - delta_nq * Z_n):\n");
    for (int n = 0; n < N_MAX; n++) {
        for (int q = 0; q < N_MAX; q++) {
            printf("%10.6f ", matrix[n][q]);
        }
        printf("\n");
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
