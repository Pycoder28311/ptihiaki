#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include "constants.h"

#include "Main.h"
#include "root.h"

extern double a, A,
              b,
              D,
              L,
              K, Kn, r,
              tol, err, UNITLEN, LD;
double nsp = 1.0; //10.0
double eps = 1e-12;

extern int Iresu, Jresu, M;
extern int Nmax, Mmax;
extern double e0, ei, m0, PI;
extern char workdir[300];
extern int maintype;

extern double complex w, w2, e;

void riza(double (*func)(double),
          double down, double up,
          double step1, double tol,
          int *nr, double rizes[])
{
    int nroot, nw, ier;
    int Nmaxt, Mmaxt;
    double ea, zx1, zx2;
    double xb1[500], xb2[500];
    double x1, x2, fa, fb, fc, root, step;

    Nmaxt = Nmax;
    Mmaxt = Mmax;

    *nr = 0;
    step = step1;

    if ((down + step > up) && (down < up))
        step = up - down;
    for (ea = down; (up - ea) / up > 1.e-10; ea += step) {

        if ((ea + step > up) && (ea < up))
            step = up - ea;

        zx1 = ea;
        zx2 = ea + step;
        nroot = 500;

        zbrak(func, zx1, zx2, nsp, xb1, xb2, &nroot);

        if (nroot != 0) {

            for (nw = 1; nw <= nroot; nw++) {

                x1 = xb1[nw];
                x2 = xb2[nw];
                root = zbrent(func, x1, x2, tol);

                if (root <= 0.0)
                    continue;

                fa = func(root - step / 20.0 / nsp);
                fb = func(root);
                fc = func(root + step / 20.0 / nsp);

                ier = 0;

                if ((fa * fb < 0.0) && (fabs(fb) < fabs(fc))) ier = 1;
                if ((fb * fc < 0.0) && (fabs(fb) < fabs(fa))) ier = 1;
                if (fb == 0.0) ier = 1;

                if (ier == 1) {
                    (*nr)++;
                    rizes[*nr] = root;
                }
            }
        }
    }
}

void zbrak(double (*fx)(double),
           double x1, double x2,
           int n,
           double xb1[], double xb2[],
           int *nb)
{
    int nbb, i;
    double x, fp, fc, dx;

    nbb = 0;
    dx = (x2 - x1) / n;     
    x = x1;
    fp = (*fx)(x);

    for (i = 1; i <= n; i++) {
        x += dx;
        fc = (*fx)(x);

        if (fabs(fc) < eps) fc = 0.0;
        if (fabs(fp) < eps) fp = 0.0;

        if (fp * fc < 0) {
            xb1[++nbb] = x - dx;
            xb2[nbb]   = x;

            if (*nb == nbb)
                return;
        }

        fp = fc;
    }

    *nb = nbb;
}

void nrerror(const char *msg) {
    exit(1);
}

static inline double sign(double a)
{
    return (a >= 0.0) ? 1.0 : -1.0;
}

double zbrent(double (*func)(double), double x1, double x2, double tol)
/*
Χρησιμοποιώντας τη μέθοδο Brent, βρίσκει τη ρίζα μιας συνάρτησης func,
η οποία είναι γνωστό ότι βρίσκεται μεταξύ των x1 και x2.
*/
{
    int iter;
    double aa = x1, b = x2, c = x2;
    double d = 0.0, e = 0.0;
    double fa = (*func)(aa);
    double fb = (*func)(b);
    double fc, p, q, r, s;
    double tol1, xm, min1, min2;

    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        //printf("Root must be bracketed in zbrent");
        return -11.0;
    }

    fc = fb;

    for (iter = 1; iter <= ITMAX; iter++) {

        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = aa;
            fc = fa;
            e = d = b - aa;
        }

        if (fabs(fc) < fabs(fb)) {
            aa = b;
            b = c;
            c = aa;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
        xm = 0.5 * (c - b);

        if (fabs(xm) <= tol1 || fb == 0.0)
            return b;

        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {

            s = fb / fa;

            if (aa == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - aa) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if (p > 0.0) q = -q;
            p = fabs(p);

            min1 = 3.0 * xm * q - fabs(tol1 * q);
            min2 = fabs(e * q);

            if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }

        aa = b;
        fa = fb;

        if (fabs(d) > tol1)
            b += d;
        else
            b += sign(tol1) * xm;

        fb = (*func)(b);
    }

    nrerror("Maximum number of iterations exceeded in zbrent");
    return -10.0;   
}

double rtflsp(double (*func)(double), double x1, double x2, double xacc)
/*
Χρησιμοποιώντας τη μέθοδο false position, βρίσκει τη ρίζα μιας συνάρτησης func,
η οποία είναι γνωστό ότι βρίσκεται μεταξύ των x1 και x2.
*/
{
    int j;
    double fl, fh, xl, xh, swap, dx, del, f, rtf;

    fl = (*func)(x1);
    fh = (*func)(x2);   

    if (fl * fh > 0.0)
        nrerror("Root must be bracketed in rtflsp");

    if (fl < 0.0) {
        xl = x1;
        xh = x2;
    } else {
        xl = x2;
        xh = x1;
        swap = fl;
        fl = fh;
        fh = swap;
    }

    dx = xh - xl;

    for (j = 1; j <= MAXIT; j++) {  
        rtf = xl + dx * fl / (fl - fh);
        f = (*func)(rtf);

        if (f < 0.0) {
            del = xl - rtf;
            xl = rtf;
            fl = f;
        } else {
            del = xh - rtf;
            xh = rtf;
            fh = f;
        }

        dx = xh - xl;

        if (fabs(del) < xacc || f == 0.0)
            return rtf;   
    }

    nrerror("Maximum number of iterations exceeded in rtflsp");
    return 0.0;   
}
