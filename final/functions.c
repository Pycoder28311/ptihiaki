#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "globals.h"

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
    FvResult F_D   = F_v(kl, k_perp, D_radius);
    FvResult G_D   = G_v(kl, k_perp, D_radius);

    double result = (F_rho.value * G_D.Fp - F_D.Fp * G_rho.value) / G_D.Fp;

    return result;
}

// Compute derivative dO_kl/drho using precomputed Fp and Gp
double O_klp(double kl, double complex k_perp, double rho) {
    FvResult F_rho = F_v(kl, k_perp, rho);
    FvResult G_rho = G_v(kl, k_perp, rho);
    FvResult F_D   = F_v(kl, k_perp, D_radius);
    FvResult G_D   = G_v(kl, k_perp, D_radius);

    return F_rho.Fp - (F_D.Fp * G_rho.Fp / G_D.Fp);
}

double Znu(double phi_c, double complex k_perp, int n, double rho) {
    FvResult F = F_v(n, k_perp, rho);  // compute FvResult inside
    double Fp = F.Fp;                   // extract derivative

    double Zn = phi_c * cabs(k_perp) * Fp;

    return Zn;
}

double Znq(double n, double q, double complex k_perp, double a, double phi_c, double phi_i, int lmax) {
    double S = 0.0;

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