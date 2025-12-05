#include "head.h"
#include "globals.h"
#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NPOINTS 10000  // number of points for numerical integration

int delta(int n, int q) {
    return (n == q) ? 1 : 0;
}

double k_of_n(int n, int m) {
    return m + n * N;
}

double complex I_integral(double alpha, double beta, double phi1, double phi2) {
    double complex sum = 0.0 + 0.0*I;
    double dphi = (phi2 - phi1) / NPOINTS;
    double phi;

    for (int i = 0; i <= NPOINTS; i++) {
        phi = phi1 + i * dphi;
        double weight = 1.0;
        if (i == 0 || i == NPOINTS) weight = 0.5;  // trapezoid endpoints
        sum += weight * cexp(I * alpha * phi) * cos(beta * (phi - phi1)) * dphi;
    }
    return sum;
}
double Z_n(double n, double m, double k_perp) {
    double k_n = k_of_n(n, m); 
    double phi_T = 2.0 * M_PI / (double)N;
    return phi_T * k_perp * Fp_n(k_n, k_perp, a);
}
double complex S(double k_perp, double k_q, double k_n) {
    double complex S = 0.0 + 0.0*I;

    for (int l = 0; l <= Lmax; l++) {
        double k_l = l * M_PI / phi_c; 

        double Op_l_a = Op_l(l, k_perp, a);        // O'_l(k_perp a)
        double O_l_a = O_l(l, k_perp, a);  // O_l(k_perp a)
        
        double complex I1 = I_integral(k_q, k_l, phi_i, phi_i + phi_c);
        double complex I2 = I_integral(-k_n, k_l, phi_i, phi_i + phi_c);

        S += (Op_l_a / O_l_a) * (I1 * Op_l_a * I2) / (1.0 + delta(l,0));
    }

    return S;
}
double complex Z_nq(double q, double n, double m, double k_perp) {
    double k_n = k_of_n(n, m); 
    double k_q = k_of_n(q, m); 

    double complex Sum = S(k_perp, k_q, k_n);
    double Fkq = F_n(k_q, k_perp, a);  // F_{k_q}(k_perp alpha)

    return (2.0 * Fkq * Sum) / phi_c;
}
double complex Z_diff(int q, int n, double m, double k_perp) {
    double complex Znq = Z_nq(q, n, m, k_perp);
    double complex Zn = Z_n(n, m, k_perp);  // Z_n
    return Znq - delta(n,q) * Zn;
}