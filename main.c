#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "head.h"
#include "globals.h"

/* ---- f(x) = e^x - e, υπολογισμένο με σειρά Taylor ---- */
double f(double x) {
    return x - sqrt(2.0);
}

double brent(double a, double b, double tol) {
    double fa = f(a);
    double fb = f(b);

    if (fa * fb >= 0) {
        // Not bracketed
        return NAN;
    }

    double c = a, fc = fa;
    double d = b - a, e = d;

    while (1) {
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        double m = 0.5 * (c - b);
        double toler = 2.0 * tol * fmax(1.0, fabs(b));

        if (fabs(m) <= toler || fb == 0.0) return b;

        double s, p, q;
        if (fabs(e) < toler || fabs(fa) <= fabs(fb)) {
            // Bisection
            d = m;
            e = m;
        } else {
            s = fb / fa;
            if (a == c) {
                // Secant
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                // Inverse quadratic interpolation
                q = fa / fc;
                double r = fb / fc;
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0) q = -q;
            p = fabs(p);

            if ( (2.0 * p) < fmin(3.0 * m * q - fabs(toler * q), fabs(e * q)) ) {
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

        fb = f(b);

        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
    }
}

/* Find all roots in [xmin, xmax] by scanning and using brent on brackets */
int find_roots(double xmin, double xmax, double dx, double tol, double tol_zero,
               double *roots, int max_roots) {
    int found = 0;
    double x = xmin;
    double fx = f(x);

    while (x < xmax && found < max_roots) {
        double x2 = x + dx;
        if (x2 > xmax) x2 = xmax;
        double fx2 = f(x2);

        // Case 1: exact (or nearly exact) zero at x
        if (fabs(fx) <= tol_zero) {
            // avoid duplicates
            if (found == 0 || fabs(roots[found-1] - x) > tol*10) {
                roots[found++] = x;
            }
        }
        // Case 2: sign change -> bracketed root
        else if (fx * fx2 < 0.0) {
            double r = brent(x, x2, tol);
            if (!isnan(r)) {
                // avoid duplicates (close roots)
                int dup = 0;
                for (int k = 0; k < found; ++k) {
                    if (fabs(roots[k] - r) < tol * 10) { dup = 1; break; }
                }
                if (!dup && found < max_roots) roots[found++] = r;
            }
        }
        // Case 3: near-exact at x2
        if (fabs(fx2) <= tol_zero) {
            if (found == 0 || fabs(roots[found-1] - x2) > tol*10) {
                roots[found++] = x2;
            }
        }

        x = x2;
        fx = fx2;
    }

    return found;
}

int main() {

    double xmin = -1.0, xmax = 5.0;
    double dx = 0.25;         // scanning step (smaller -> more robust, slower)
    double tol = 1e-8;        // brent tolerance
    double tol_zero = 1e-12;  // consider f(x) ~ 0
    int max_roots = 100;
    double *roots = malloc(max_roots * sizeof(double));
    if (!roots) return 1;

    int n = find_roots(xmin, xmax, dx, tol, tol_zero, roots, max_roots);

    printf("Found %d root(s):\n", n);
    for (int i = 0; i < n; ++i) {
        printf("  root[%d] = %.12f  f(root)=%.12e\n", i, roots[i], f(roots[i]));
    }

    free(roots);




    int n_max = 2;
    int N_mtx = 2*n_max + 1; //Το μέγεθος του Πίνακα

    // --- CREATE MATRIX USING VARIABLES ---
    double complex **M = create_matrix(N_mtx, k_perp);

    // Μετατροπή σε pointer array για τη συνάρτηση
    double complex *M_ptrs[N_mtx];
    for (int i = 0; i < N_mtx; i++)
        M_ptrs[i] = M[i];

    double complex det = determinant_complex(M_ptrs, N_mtx);
    printf("Determinant = %.10lf + %.10lf i\n", creal(det), cimag(det));

    for (int i = 0; i < N_mtx; i++)
        free(M[i]);
    free(M);
    return 0;
}
