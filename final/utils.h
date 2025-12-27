#ifndef UTILS_H
#define UTILS_H
#include <complex.h>
#include "globals.h"

double Jnu(double n, double x);
double Ynu(double n, double x);
double Inu(double n, double x);
double Knu(double n, double x);

/* Derivatives */
double Jnp(double n, double x);
double Ynp(double n, double x);
double Inp(double n, double x);
double Knp(double n, double x);

/* -------------------------
   Struct to hold value and derivatives
   ------------------------- */
typedef struct {
    double value;  // function value
    double dF;     // derivative w.r.t rho
    double Fp;     // derivative normalized
} FvResult;

/* ----------- Core radial functions ----------- */
FvResult F_v(int v, double complex k_perp, double rho);
FvResult G_v(int v, double complex k_perp, double rho);

/* ----------- Angular integrals ----------- */
double complex I_alpha_beta(double alpha, double beta,
                            double phi1, double phi2);

/* ----------- Utility math ----------- */
int delta(int l, int m);
double k_l(int l, double phi_c);
double k_n(int m, int n, double N);

/* ----------- Operators ----------- */
double O_kl(double kl, double complex k_perp, double rho);
double O_klp(double kl, double complex k_perp, double rho);

/* ----------- Z functions ----------- */
double Znu(double phi_c, double complex k_perp, int n, double rho);

double Znq(double n, double q, double complex k_perp,
           double a, double phi_c, double phi_i, int lmax);

#define NM (2*n_max + 1)

/* Matrix builders */
void build_Z_matrix(double matrix[NM][NM], double complex k_perp);

/* Determinant utilities */
double compute_det_matrix(int N, double mat[N][N]);
double compute_det_for_k(double k_perp);

/* Reduced system builder + solver */
double** create_mat_reduced(double complex k_perp);

#endif /* UTILS_H */
