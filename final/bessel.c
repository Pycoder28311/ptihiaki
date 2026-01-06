#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_sf_bessel.h>

// Συναρτήσεις Bessel

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

// Παράγωγοι των συναρτήσεων Bessel

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