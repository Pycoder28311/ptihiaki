#include "head.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double complex **create_matrix(int N_mtx, double k_perp);
double complex determinant_complex(double complex **M_ptrs, int N_mtx);

double f(double k_perp, int n_max) {
    int N_mtx = 2 * n_max + 1;

    double complex **M = create_matrix(N_mtx, k_perp);
    if (!M) {
        fprintf(stderr, "Error: Could not allocate matrix.\n");
        return 0.0;  
    }

    double complex *M_ptrs[N_mtx];
    for (int i = 0; i < N_mtx; i++)
        M_ptrs[i] = M[i];

    double complex det = determinant_complex(M_ptrs, N_mtx);
    //printf("Debug: For k_perp = %.12f, Determinant = %.10lf + %.10lf i\n", k_perp, creal(det), cimag(det));

    for (int i = 0; i < N_mtx; i++)
        free(M[i]);
    free(M);

    // Επίστρεψε το πραγματικό μέρος της ρίζας
    return creal(det); 
}

// Για δοκιμή, δεν την χρησιμοποιώ
double f_Other(double k_perp, int n_max) {
    return k_perp - 1.0;  
}

double brent(double a, double b, double tol, int n_max) {
    double fa = f(a, n_max);
    double fb = f(b, n_max);

    if (fa * fb >= 0) return NAN;

    double c = a, fc = fa;
    double d = b - a, e = d;
    //printf("Debug: Starting Brent's method with a=%.12f (f=%.12e), b=%.12f (f=%.12e)\n", a, fa, b, fb);

    while (1) {
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        double m = 0.5 * (c - b);
        //printf("Debug: Current interval [%.20f, %.20f], m=%.20f, fb=%.12f %12f\n", b, c, m, f(b,n_max), f(c, n_max));
        double toler = 0.00001 * tol * fmax(1.0, fabs(b));
        //printf("Debug: Tolerance = %.20f\n", toler);

        if (fabs(m) <= toler || fb == 0.0) return b;

        double s, p, q;
        if (fabs(e) < toler || fabs(fa) <= fabs(fb)) {
            d = m; e = m;
        } else {
            s = fb / fa;
            if (a == c) {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                double r = fb / fc;
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0) q = -q;
            p = fabs(p);

            if ((2.0 * p) < fmin(3.0 * m * q - fabs(toler * q), fabs(e * q))) {
                e = d;
                d = p / q;
            } else {
                d = m;
                e = m;
            }
        }

        a = b;
        fa = fb;

        if (fabs(d) > toler)
            b += d;
        else
            b += (m > 0 ? toler : -toler);

        fb = f(b, n_max);

        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
    }
}

int find_roots(double xmin, double xmax, double dx, double tol, double tol_zero,
               double *roots, int max_roots, int n_max) {
    int found = 0;
    double x = xmin;
    double fx = f(x, n_max);

    while (x < xmax && found < max_roots) {
        double x2 = x + dx;
        if (x2 > xmax) x2 = xmax;
        double fx2 = f(x2, n_max);
        printf("Debug: Scanning interval [%.12f, %.12f] with f(x) = %.12e, f(x2) = %.12e %12f\n", x, x2, fx, fx2, tol_zero);

        if (fabs(fx) <= tol_zero) {
            if (found == 0 || fabs(roots[found-1] - x) > tol*0.001) {
                roots[found++] = x;
            }
        } else if (fx * fx2 < 0.0) {
            double r = brent(x, x2, tol, n_max);
            printf("Debug: Found root at x = %.30f with f(x) = %.12e\n", r, f(r, n_max));
            if (!isnan(r)) {
                int dup = 0;
                for (int k = 0; k < found; ++k) {
                    if (fabs(roots[k] - r) < tol * 10) { dup = 1; break; }
                }
                if (!dup && found < max_roots) roots[found++] = r;
            }
        }

        if (fabs(fx2) <= tol_zero) {
            if (found == 0 || fabs(roots[found-1] - x2) > tol*0.001) {
                roots[found++] = x2;
            }
        }

        x = x2;
        fx = fx2;
    }

    return found;
}