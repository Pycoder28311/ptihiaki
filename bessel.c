#include "head.h"
#include "globals.h"
#include <math.h>
#include <complex.h>
#include <stdio.h>

// Jn(x) using Taylor series
double J_n(double n, double x) {
    double sum = 0.0;
    double term;
    int k;
    for (k = 0; k < 500; k++) {
        term = pow(-1, k) * pow(x/2, 2*k + n) / (tgamma(k+1) * tgamma(n+k+1));
        sum += term;
        if (fabs(term) < 1e-15) break;
    }
    return sum;
}
double Y_n(double n, double x) {
    if (x == 0) return -INFINITY; // Y_n(0) -> -∞
    double Jn = J_n(n, x);
    double n_eps = n + EPS_BESSEL; // μικρή μετατόπιση για limit
    double Jn_eps = J_n(n_eps, x);
    double Jneg_eps = J_n(-n_eps, x);
    return (Jn_eps * cos(n_eps * PI) - Jneg_eps) / sin(n_eps * PI);
}
// Modified Bessel I_v(x) για δεκαδικό v
double I_n(double n, double x) {
    double sum = 0.0;
    double term;
    for (int k = 0; k < 50; k++) {
        term = pow(x/2, 2*k + n) / (tgamma(k+1) * tgamma(n+k+1));
        sum += term;
        if (fabs(term) < 1e-15) break;
    }
    return sum;
}
// Modified Bessel K_v(x) για δεκαδικό v
double K_n(double n, double x) {
    if (x == 0.0) return 1e300;  // K_n(0) -> πολύ μεγάλο
    double n_eps = n + EPS_BESSEL;
    double Ipos = I_n(n_eps, x);
    double Ineg = I_n(-n_eps, x);
    return (PI / 2.0) * (Ineg - Ipos) / sin(n_eps * PI);
}

// Fv using your J
double F_n(double n, double k_perp, double r) {
    if (k_perp * k_perp > 0.0) {
        return J_n(n, k_perp * r);  // Jv(k_perp * rho)
    } else {
        double s = sqrt(-k_perp*k_perp);  // s^2 = -k_perp^2
        return I_n(n, s * r);       // Iv(s * rho)
    }
}
// Gv using your Y
double G_n(double n, double k_perp, double r) {
    if (k_perp * k_perp > 0.0) {
        return Y_n(n, k_perp * r);  // Yv(k_perp * rho)
    } else {
        double s = sqrt(-k_perp*k_perp);  // s^2 = -k_perp^2
        return K_n(n, s * r);      // Kv(s * rho)
    }
}



// Derivative of J_v(x)
double Jp_n(double n, double x) {
    if (n == 0.0) return -J_n(1.0, x);   // J0' = -J1
    return 0.5 * (J_n(n-1.0, x) - J_n(n+1.0, x));
}
double Yp_n(double n, double x) {
    if (n == 0.0) return -Y_n(1.0, x);   // Y0' = -Y1
    return 0.5 * (Y_n(n-1.0, x) - Y_n(n+1.0, x));
}
// Derivative of I_v(x)
double Ip_n(double n, double x) {
    if (n == 0.0) return I_n(1.0, x);    // I0' = I1
    return 0.5 * (I_n(n-1.0, x) + I_n(n+1.0, x));
}
// Derivative of K_v(x)
double Kp_n(double n, double x) {
    if (n == 0.0) return -K_n(1.0, x);   // K0' = -K1
    return -0.5 * (K_n(n-1.0, x) + K_n(n+1.0, x));
}



// Derivative w.r.t argument (n double)
double Fp_n(double n, double k_perp, double r) {
    if (k_perp * k_perp > 0.0) {
        return Jp_n(n, k_perp * r);  // d/dρ Jv(k_perp * ρ)
    } else {
        double s = sqrt(-k_perp * k_perp);
        return Ip_n(n, s * r);       // d/dρ Iv(s * ρ)
    }
}
double Gp_n(double n, double k_perp, double r) {
    if (k_perp * k_perp > 0.0) {
        return Yp_n(n, k_perp * r);  // d/dρ Yv(k_perp * ρ)
    } else {
        double s = sqrt(-k_perp * k_perp);
        return Kp_n(n, s * r);       // d/dρ Kv(s * ρ)
    }
}
// Derivative w.r.t rho (chain rule)
double dF_n(double n, double k_perp, double r) {
    if (k_perp * k_perp > 0.0) {
        return k_perp * Jp_n(n, k_perp * r);  // d/dρ Jv(k_perp * ρ)
    } else {
        double s = sqrt(-k_perp * k_perp);
        return s * Ip_n(n, s * r);            // d/dρ Iv(s * ρ)
    }
}
double dG_n(double n, double k_perp, double r) {
    if (k_perp * k_perp > 0.0) {
        return k_perp * Yp_n(n, k_perp * r);  // d/dρ Yv(k_perp * ρ)
    } else {
        double s = sqrt(-k_perp * k_perp);
        return s * Kp_n(n, s * r);            // d/dρ Kv(s * ρ)
    }
}



double O_l(double l, double k_perp, double r) {
    double F_r = F_n(l, k_perp, r);
    double Gp_D  = Gp_n(l, k_perp, D);
    double Fp_D  = Fp_n(l, k_perp, D);
    double G_r = G_n(l, k_perp, r);

    return (F_r * Gp_D - Fp_D * G_r) / Gp_D;
}
double Op_l(double l, double k_perp, double r) {
    double Fp_r = Fp_n(l, k_perp, r);
    double Gp_D  = Gp_n(l, k_perp, D);
    double Fp_D  = Fp_n(l, k_perp, D);
    double Gp_r = Gp_n(l, k_perp, r);

    return (Fp_r * Gp_D - Fp_D * Gp_r) / Gp_D;
}
double dO_l(double l, double k_perp, double r) {
    double Op_r = Op_l(l, k_perp, r);

    return k_perp * Op_r; 
}