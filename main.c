#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "head.h"
#include "globals.h"

#include <gsl/gsl_sf_bessel.h>

double Jnu(double n, double x) {
    if (fabs(x) < 1e-15) return (n == 0.0) ? 1.0 : 0.0;
    return gsl_sf_bessel_Jnu(n, x);
}

double Ynu(double n, double x) {
    if (x <= 0.0) {
        //fprintf(stderr, "Warning: Ynu not defined for x <= 0, returning NAN\n");
        return NAN;
    }
    return gsl_sf_bessel_Ynu(n, x);
}

double Inu(double n, double x) {
    if (fabs(x) < 1e-15) return (n == 0.0) ? 1.0 : 0.0;
    return gsl_sf_bessel_Inu(n, fabs(x)); // Iν(-x) = Iν(x)
}

double Knu(double n, double x) {
    if (x <= 0.0) {
        //fprintf(stderr, "Warning: Knu not defined for x <= 0, returning NAN\n");
        return NAN;
    }
    return gsl_sf_bessel_Knu(n, x);
}

/* ----------- Safe derivatives ----------- */

double Jnp(double n, double x) {
    if (fabs(x) < 1e-15) return 0.0; // derivative at 0
    return 0.5 * (Jnu(n-1, x) - Jnu(n+1, x));
}

double Ynp(double n, double x) {
    if (x <= 0.0) {
        //fprintf(stderr, "Warning: Ynp not defined for x <= 0, returning NAN\n");
        return NAN;
    }
    return 0.5 * (Ynu(n-1, x) - Ynu(n+1, x));
}

double Inp(double n, double x) {
    if (fabs(x) < 1e-15) return 0.0; // derivative at 0
    return 0.5 * (Inu(n-1, x) + Inu(n+1, x));
}

double Knp(double n, double x) {
    if (x <= 0.0) {
        //fprintf(stderr, "Warning: Knp not defined for x <= 0, returning NAN\n");
        return NAN;
    }
    return -0.5 * (Knu(n-1, x) + Knu(n+1, x));
}

/* -------------------------
   Struct to hold value and derivatives
   ------------------------- */
typedef struct {
    double value;  // function value
    double dF;     // derivative w.r.t rho
    double Fp;     // derivative normalized
} FvResult;

FvResult F_v(int v, double complex k_perp, double rho) {
    FvResult res;

    if (cimag(k_perp) == 0.0) {  // real -> J_v
        double k_real = creal(k_perp);
        res.value = Jnu(v, k_real * rho);
        res.dF    = k_real * Jnp(v, k_real * rho);
        res.Fp    = Jnp(v, k_real * rho);
    } else {  // imaginary -> I_v
        double s = cimag(k_perp);
        res.value = Inu(v, s * rho);
        res.dF    = s * Inp(v, s * rho);
        res.Fp    = Inp(v, s * rho);
    }

    return res;
}

FvResult G_v(int v, double complex k_perp, double rho) {
    FvResult res;

    if (cimag(k_perp) == 0.0) {  // real -> Y_v
        double k_real = creal(k_perp);
        res.value = Ynu(v, k_real * rho);
        res.dF    = k_real * Ynp(v, k_real * rho);
        res.Fp    = Ynp(v, k_real * rho);
    } else {  // imaginary -> K_v
        double s = cimag(k_perp);
        res.value = Knu(v, s * rho);
        res.dF    = s * Knp(v, s * rho);
        res.Fp    = Knp(v, s * rho);
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

double O_kl(double kl, double complex k_perp, double rho) {
    // Compute F and G at rho
    FvResult F_rho = F_v(kl, k_perp, rho);
    FvResult G_rho = G_v(kl, k_perp, rho);

    // Compute F and G at D (make sure D is defined)
    FvResult F_D   = F_v(kl, k_perp, 1.5);
    FvResult G_D   = G_v(kl, k_perp, 1.5);

    double result = (F_rho.value * G_D.Fp - F_D.Fp * G_rho.value) / G_D.Fp;

    return result;
}

// Compute derivative dO_kl/drho using precomputed Fp and Gp
double O_klp(double kl, double complex k_perp, double rho) {
    FvResult F_rho = F_v(kl, k_perp, rho);
    FvResult G_rho = G_v(kl, k_perp, rho);
    FvResult F_D   = F_v(kl, k_perp, 1.5);
    FvResult G_D   = G_v(kl, k_perp, 1.5);

    return F_rho.Fp - (F_D.Fp * G_rho.Fp / G_D.Fp);
}

double Znt(double phi_c, double complex k_perp, int n, double rho) {
    FvResult F = F_v(n, k_perp, rho);  // compute FvResult inside
    double Fp = F.Fp;                   // extract derivative

    double Zn = phi_c * cabs(k_perp) * Fp;

    return Zn;
}

double Z_nqt(double n, double q, double complex k_perp, double a, double phi_c, double phi_i, int lmax) {
    double S = 0.0;
    double m = 1.0;
    double N = 5.0;

    double kn = k_n(m, n, N);
    double kq = k_n(m, q, N);

    for (int l = 0; l <= lmax; l++) {
        double kl = k_l(l, phi_c);

        double O_kl_a = O_kl(kl, k_perp, a);
        double O_kl_a_p = O_klp(kl, k_perp, a);

        double complex I1 = I_alpha_beta(kq, kl, phi_i, phi_i + phi_c);
        double complex I2 = I_alpha_beta(-kn, kl, phi_i, phi_i + phi_c);

        double delta_l0 = delta(l,0);

        double complex term = (O_kl_a_p / O_kl_a) * (I1 * O_kl_a_p * I2) / (1.0 + delta_l0);
        S += creal(term);
    }

    FvResult Fq = F_v(q, k_perp, a); // F_{k_q}(k_perp*a)
    double Znq = 2.0 * Fq.Fp * S / phi_c;

    return Znq;
}

double det_matrix(int N, double mat[N][N]) {
    double A[N][N];
    int i, j, k;
    double det = 1.0;

    // Copy input matrix
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            A[i][j] = mat[i][j];

    // LU decomposition (no pivoting)
    for (k = 0; k < N; k++) {
        if (fabs(A[k][k]) < 1e-12) return 0.0; // nearly singular
        for (i = k + 1; i < N; i++) {
            double factor = A[i][k] / A[k][k];
            for (j = k; j < N; j++)
                A[i][j] -= factor * A[k][j];
        }
    }

    // determinant = product of diagonal
    for (i = 0; i < N; i++)
        det *= A[i][i];

    return det;
}

const double at = 1.0;       // example radius
const double phi_ct = 1.0;   // angle parameter
const double phi_it = 0.0;   // angle parameter
const int lmax = 10;        // maximum l for sums
const int N_MAX = 9;   // matrix dimension for real k

double compute_det_for_k_real(double k_real) {
    double complex k_perp = k_real + 0.0*I;

    double matrix[N_MAX][N_MAX];

    for (int n = 0; n < N_MAX; n++) {
        for (int q = 0; q < N_MAX; q++) {
            double Znq = Z_nqt(n, q, k_perp, at, phi_ct, phi_it, lmax);
            double Zn = Znt(phi_ct, k_perp, n, at);
            matrix[n][q] = Znq - delta(n, q) * Zn;
        }
    }

    return det_matrix(N_MAX, matrix);
}

double confirm_det_for_k(double complex k_perp) {

    double matrix[N_MAX][N_MAX];

    // Fill the matrix
    for (int n = 0; n < N_MAX; n++) {
        for (int q = 0; q < N_MAX; q++) {
            double Znq = Z_nqt(n, q, k_perp, at, phi_ct, phi_it, lmax);
            double Zn = Znt(phi_ct, k_perp, n, at);
            matrix[n][q] = Znq - delta(n, q) * Zn;
        }
    }

    // Optional: print the matrix
    /*printf("Matrix (Z_nq - delta_nq * Z_n):\n");
    for (int n = 0; n < N_MAX; n++) {
        for (int q = 0; q < N_MAX; q++) {
            printf("%10.6f ", matrix[n][q]);
        }
        printf("\n");
    }*/

    // Compute determinant
    double det = det_matrix(N_MAX, matrix);
    return det;
}

void riza(double (*func)(double),
          double down, double up,
          double step, double tol,
          int *nr, double rizes[]);

// ---------- Main Function ----------
int main() {

    int nr;
    double roots[100];

    // Example: scan for real k from 0.0 to 20.0
    riza(compute_det_for_k_real, 0.0, 100.0, 0.5, 1e-6, &nr, roots);

    printf("Found %d roots:\n", nr);

    for (int i = 0; i < nr; i++) {
        double k_real = roots[i];
        double complex k_perp = k_real + 0.0*I;

        printf("Matrix for root %d (k_real = %g):", i+1, k_real);

        // Use your confirm_det_for_k to build and print the matrix
        double det = confirm_det_for_k(k_perp);

        printf("Determinant at root %d: %g\n", i+1, det);
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
