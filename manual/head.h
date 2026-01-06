#ifndef HEAD_H
#define HEAD_H

#include <complex.h>
#include <stdio.h>

#define EULER 0.57721566490153286060651209
#define PI 3.14159265358979323846
#define I_NPOINTS 10000
#define EPS_BESSEL 1e-8  
#define MAX_ITER 10000
#define TOL 1e-10

double J_n(double n, double x);
double Y_n(double n, double x);
double Jp_n(double n, double x);
double Yp_n(double n, double x);

double I_n(double n, double x);
double K_n(double n, double x);
double Ip_n(double n, double x);
double Kp_n(double n, double x);

double F_n(double n, double k_perp, double r);
double G_n(double n, double k_perp, double r);
double Fp_n(double n, double k_perp, double r); 
double Gp_n(double n, double k_perp, double r); 

double dF_n(double n, double k_perp, double r); 
double dG_n(double n, double k_perp, double r); 

double O_l(double l, double k_perp, double r);
double Op_l(double l, double k_perp, double r);
double dO_l(double l, double k_perp, double r);

double complex I_integral(double alpha, double beta, double phi1, double phi2);
double Z_n(double n, double m, double k_perp);
double complex S(double k_perp, double k_q, double k_n);
double complex Z_nq(double q, double n, double m, double k_perp);
double complex Z_diff(int q, int n, double m, double k_perp);

void getCofactor(double **mat, double **temp, int p, int q, int n);
double determinant(double **mat, int n);
double complex determinant_complex(double complex **mat, int N);
double complex **create_matrix(int N_mtx, double k_perp);

double f(double x, int n_max);
double brent(double a, double b, double tol, int n_max);
int find_roots(double xmin, double xmax, double dx, double tol, double tol_zero, double *roots, int max_roots, int n_max);

int delta(int n, int q);
double k_of_n(int n, int m);

void mat_vec_mul(int N, double **A, double *x, double *y);
double norm(int N, double *v);
void power_iteration(int N, double **A, double *eigvec, double *eigval);
void deflate_matrix(int N, double **A, double *eigvec, double eigval);
void find_all_eigenvalues(int N, double complex **A, FILE *fp);

void compute_matrix_determinant(int n_max, double k_perp, FILE *fp);

#endif
