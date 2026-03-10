#ifndef UTILS_H
#define UTILS_H
#include <complex.h>
#include "globals.h"

double Jnu(double n, double x);
double Ynu(double n, double x);
double Inu(double n, double x);
double Knu(double n, double x);

double Jnp(double n, double x);
double Ynp(double n, double x);
double Inp(double n, double x);
double Knp(double n, double x);

typedef struct {
    double value; 
    double dF;  
    double Fp;    
} FvResult;

FvResult F_v(int v, double k_perp, double rho);
FvResult G_v(int v, double k_perp, double rho);

double complex I_alpha_beta(double alpha, double beta,
                            double phi1, double phi2);

int delta(int l, int m);
double k_l(int l, double phi_c);
double k_n(int m, int n, double N);

double O_kl(double kl, double k_perp, double rho, double D);
double O_klp(double kl, double k_perp, double rho, double D);

double Zn(double phi_c, double k_perp, int n, double rho, double D);

double Znq(double n, double q, double k_perp,
           double a, double phi_c, double phi_i, int lmax, double D);

#define NM (2*n_max + 1)

void build_Z_matrix(double matrix[NM][NM], double k_perp, double D);

double compute_det_matrix(int N, double mat[N][N], double D);
double compute_det_for_k(double k_perp, double D);

double** create_mat_reduced(double k_perp, double D);

#endif 
