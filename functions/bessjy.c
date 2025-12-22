#include <math.h>
#include <stdlib.h>
#include "constants.h"

int bessjy01a(double x,
              double *j0, double *j1,
              double *y0, double *y1,
              double *j0p, double *j1p,
              double *y0p, double *y1p)
{
    double x2,r,ec,w0,w1,r0,r1,cs0,cs1;
    double cu,p0,q0,p1,q1,t1,t2;
    int k,kz;
    static double a[] = {
        -7.03125e-2,
         0.112152099609375,
        -0.5725014209747314,
         6.074042001273483,
        -1.100171402692467e2,
         3.038090510922384e3,
        -1.188384262567832e5,
         6.252951493434797e6,
        -4.259392165047669e8,
         3.646840080706556e10,
        -3.833534661393944e12,
         4.854014686852901e14,
        -7.286857349377656e16,
         1.279721941975975e19};
    static double b[] = {
         7.32421875e-2,
        -0.2271080017089844,
         1.727727502584457,
        -2.438052969955606e1,
         5.513358961220206e2,
        -1.825775547429318e4,
         8.328593040162893e5,
        -5.006958953198893e7,
         3.836255180230433e9,
        -3.649010818849833e11,
         4.218971570284096e13,
        -5.827244631566907e15,
         9.476288099260110e17,
        -1.792162323051699e20};
    static double a1[] = {
         0.1171875,
        -0.1441955566406250,
         0.6765925884246826,
        -6.883914268109947,
         1.215978918765359e2,
        -3.302272294480852e3,
         1.276412726461746e5,
        -6.656367718817688e6,
         4.502786003050393e8,
        -3.833857520742790e10,
         4.011838599133198e12,
        -5.060568503314727e14,
         7.572616461117958e16,
        -1.326257285320556e19};
    static double b1[] = {
        -0.1025390625,
         0.2775764465332031,
        -1.993531733751297,
         2.724882731126854e1,
        -6.038440767050702e2,
         1.971837591223663e4,
        -8.902978767070678e5,
         5.310411010968522e7,
        -4.043620325107754e9,
         3.827011346598605e11,
        -4.406481417852278e13,
         6.065091351222699e15,
        -9.833883876590679e17,
         1.855045211579828e20};

    if (x < 0.0) return 1;
    if (x == 0.0) {
        *j0 = 1.0;
        *j1 = 0.0;
        *y0 = -1e308;
        *y1 = -1e308;
        *j0p = 0.0;
        *j1p = 0.5;
        *y0p = 1e308;
        *y1p = 1e308;
        return 0;
    }

    x2 = x*x;

    if (x <= 12.0) {
        *j0 = 1.0;
        r = 1.0;
        for (k = 1; k <= 30; k++) {
            r *= -0.25 * x2 / (k*k);
            *j0 += r;
            if (fabs(r) < fabs(*j0)*1e-15) break;
        }

        *j1 = 1.0;
        r = 1.0;
        for (k = 1; k <= 30; k++) {
            r *= -0.25 * x2 / (k*(k+1));
            *j1 += r;
            if (fabs(r) < fabs(*j1)*1e-15) break;
        }
        *j1 *= 0.5*x;

        ec = log(0.5*x) + EL;

        cs0 = 0.0;
        w0 = 0.0;
        r0 = 1.0;
        for (k = 1; k <= 30; k++) {
            w0 += 1.0/k;
            r0 *= -0.25*x2/(k*k);
            r = r0 * w0;
            cs0 += r;
            if (fabs(r) < fabs(cs0)*1e-15) break;
        }
        *y0 = M_2_PI * (ec*(*j0) - cs0);

        cs1 = 1.0;
        w1 = 0.0;
        r1 = 1.0;
        for (k = 1; k <= 30; k++) {
            w1 += 1.0/k;
            r1 *= -0.25*x2/(k*(k+1));
            r = r1*(2.0*w1 + 1.0/(k+1));
            cs1 += r;
            if (fabs(r) < fabs(cs1)*1e-15) break;
        }
        *y1 = M_2_PI * (ec*(*j1) - 1.0/x - 0.25*x*cs1);
    } else {
        if (x >= 50.0) kz = 8;
        else if (x >= 35.0) kz = 10;
        else kz = 12;

        t1 = x - M_PI_4;
        p0 = 1.0;
        q0 = -0.125 / x;
        for (k = 0; k < kz; k++) {
            p0 += a[k] * pow(x, -2*k-2);
            q0 += b[k] * pow(x, -2*k-3);
        }
        cu = sqrt(M_2_PI / x);
        *j0 = cu * (p0*cos(t1) - q0*sin(t1));
        *y0 = cu * (p0*sin(t1) + q0*cos(t1));

        t2 = x - 0.75*M_PI;
        p1 = 1.0;
        q1 = 0.375 / x;
        for (k = 0; k < kz; k++) {
            p1 += a1[k]*pow(x,-2*k-2);
            q1 += b1[k]*pow(x,-2*k-3);
        }
        *j1 = cu*(p1*cos(t2) - q1*sin(t2));
        *y1 = cu*(p1*sin(t2) + q1*cos(t2));
    }

    *j0p = -(*j1);
    *j1p = *j0 - *j1 / x;
    *y0p = -(*y1);
    *y1p = *y0 - *y1 / x;

    return 0;
}

int bessjy01b(double x,
              double *j0, double *j1,
              double *y0, double *y1,
              double *j0p, double *j1p,
              double *y0p, double *y1p)
{
    double t, t2, dtmp, a0, p0, q0, p1, q1, ta0, ta1;

    if (x < 0.0) return 1;
    if (x == 0.0) {
        *j0 = 1.0;
        *j1 = 0.0;
        *y0 = -1e308;
        *y1 = -1e308;
        *j0p = 0.0;
        *j1p = 0.5;
        *y0p = 1e308;
        *y1p = 1e308;
        return 0;
    }

    if (x <= 4.0) {
        t = x/4.0;
        t2 = t*t;

        *j0 = ((((((-0.5014415e-3*t2 + 0.76771853e-2)*t2 -0.0709253492)*t2 +
                 0.4443584263)*t2 -1.7777560599)*t2 +3.9999973021)*t2 -
                 3.9999998721)*t2 + 1.0;

        *j1 = t*(((((((-0.1289769e-3*t2 + 0.22069155e-2)*t2 -0.0236616773)*t2 +
                    0.1777582922)*t2 -0.8888839649)*t2 +2.6666660544)*t2 -
                    3.999999971)*t2 + 1.9999999998);

        dtmp = (((((((-0.567433e-4*t2 +0.859977e-3)*t2 -0.94855882e-2)*t2 +
                   0.0772975809)*t2 -0.4261737419)*t2 +1.4216421221)*t2 -
                   2.3498519931)*t2 +1.0766115157)*t2 +0.3674669052;
        *y0 = M_2_PI*log(0.5*x)*(*j0) + dtmp;

        dtmp = (((((((0.6535773e-3*t2 -0.0108175626)*t2 +0.107657607)*t2 -
                   0.7268945577)*t2 +3.1261399273)*t2 -7.3980241381)*t2 +
                   6.8529236342)*t2 +0.3932562018)*t2 -0.6366197726;
        *y1 = M_2_PI*log(0.5*x)*(*j1) + dtmp/x;
    } else {
        t = 4.0/x;
        t2 = t*t;
        a0 = sqrt(M_2_PI/x);

        p0 = ((((-0.9285e-5*t2 +0.43506e-4)*t2 -0.122226e-3)*t2 +
               0.434725e-3)*t2 -0.4394275e-2)*t2 +0.999999997;
        q0 = t*(((((0.8099e-5*t2 -0.35614e-4)*t2 +0.85844e-4)*t2 -
                 0.218024e-3)*t2 +0.1144106e-2)*t2 -0.031249995);
        ta0 = x - M_PI_4;
        *j0 = a0*(p0*cos(ta0) - q0*sin(ta0));
        *y0 = a0*(p0*sin(ta0) + q0*cos(ta0));

        p1 = ((((0.10632e-4*t2 -0.50363e-4)*t2 +0.145575e-3)*t2 -
               0.559487e-3)*t2 +0.7323931e-2)*t2 +1.000000004;
        q1 = t*(((((-0.9173e-5*t2 +0.40658e-4)*t2 -0.99941e-4)*t2 +
                 0.266891e-3)*t2 -0.1601836e-2)*t2 +0.093749994);
        ta1 = x - 0.75*M_PI;
        *j1 = a0*(p1*cos(ta1) - q1*sin(ta1));
        *y1 = a0*(p1*sin(ta1) + q1*cos(ta1));
    }

    *j0p = -(*j1);
    *j1p = *j0 - *j1 / x;
    *y0p = -(*y1);
    *y1p = *y0 - *y1 / x;

    return 0;
}

int msta1(double x, int mp)
{
    double a0, f0, f1, f;
    int i, n0, n1, nn;

    a0 = fabs(x);
    n0 = (int)(1.1*a0) + 1;
    f0 = 0.5*log10(6.28*n0) - n0*log10(1.36*a0/n0) - mp;
    n1 = n0 + 5;
    f1 = 0.5*log10(6.28*n1) - n1*log10(1.36*a0/n1) - mp;

    for (i = 0; i < 20; i++) {
        nn = (int)(n1 - (n1 - n0)/(1.0 - f0/f1));
        f = 0.5*log10(6.28*nn) - nn*log10(1.36*a0/nn) - mp;
        if (fabs((double)(nn - n1)) < 1.0) break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }

    return nn;
}

int msta2(double x, int n, int mp)
{
    double a0, ejn, hmp, f0, f1, f, obj;
    int i, n0, n1, nn;

    a0 = fabs(x);
    hmp = 0.5*mp;
    ejn = 0.5*log10(6.28*n) - n*log10(1.36*a0/n);

    if (ejn <= hmp) {
        obj = mp;
        n0 = (int)(1.1*a0);
        if (n0 < 1) n0 = 1;
    } else {
        obj = hmp + ejn;
        n0 = n;
    }

    f0 = 0.5*log10(6.28*n0) - n0*log10(1.36*a0/n0) - obj;
    n1 = n0 + 5;
    f1 = 0.5*log10(6.28*n1) - n1*log10(1.36*a0/n1) - obj;

    for (i = 0; i < 20; i++) {
        nn = (int)(n1 - (n1 - n0)/(1.0 - f0/f1));
        f = 0.5*log10(6.28*nn) - nn*log10(1.36*a0/nn) - obj;
        if (fabs((double)(nn - n1)) < 1.0) break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }

    return nn + 10;
}

int bessjyna(int n, double x, int *nm,
             double *jn, double *yn,
             double *jnp, double *ynp)
{
    double bj0, bj1, f, f0, f1, f2, cs;
    int i, k, m, ecode;

    *nm = n;
    if ((x < 0.0) || (n < 0)) return 1;
    if (x < 1e-15) {
        for (i = 0; i <= n; i++) {
            jn[i] = 0.0;
            yn[i] = -1e308;
            jnp[i] = 0.0;
            ynp[i] = 1e308;
        }
        jn[0] = 1.0;
        jnp[1] = 0.5;
        return 0;
    }

    ecode = bessjy01a(x, &jn[0], &jn[1], &yn[0], &yn[1], &jnp[0], &jnp[1], &ynp[0], &ynp[1]);
    if (n < 2) return 0;

    bj0 = jn[0];
    bj1 = jn[1];

    if (n < (int)(0.9*x)) {
        for (k = 2; k <= n; k++) {
            jn[k] = 2.0*(k-1.0)*bj1/x - bj0;
            bj0 = bj1;
            bj1 = jn[k];
        }
    } else {
        m = msta1(x, 200);
        if (m < n) *nm = m;
        else m = msta2(x, n, 15);

        f2 = 0.0;
        f1 = 1.0e-100;
        for (k = m; k >= 0; k--) {
            f = 2.0*(k+1.0)/x*f1 - f2;
            if (k <= *nm) jn[k] = f;
            f2 = f1;
            f1 = f;
        }

        if (fabs(bj0) > fabs(bj1)) cs = bj0/f;
        else cs = bj1/f2;

        for (k = 0; k <= *nm; k++) {
            jn[k] *= cs;
        }
    }

    for (k = 2; k <= *nm; k++) {
        jnp[k] = jn[k-1] - k*jn[k]/x;
    }

    f0 = yn[0];
    f1 = yn[1];
    for (k = 2; k <= *nm; k++) {
        f = 2.0*(k-1.0)*f1/x - f0;
        yn[k] = f;
        f0 = f1;
        f1 = f;
    }

    for (k = 2; k <= *nm; k++) {
        ynp[k] = yn[k-1] - k*yn[k]/x;
    }

    return 0;
}

int bessjynb(int n, double x, int *nm,
             double *jn, double *yn,
             double *jnp, double *ynp)
{
    double t1, t2, f, f1, f2, bj0, bj1, bjk, by0, by1, cu, s0, su, sv;
    double ec, bs, byk, p0, p1, q0, q1;
    static double a[] = { -0.703125e-1, 0.1121520996093750, -0.5725014209747314, 6.074042001273483 };
    static double b[] = { 0.732421875e-1, -0.2271080017089844, 1.727727502584457, -2.438052969955606e1 };
    static double a1[] = { 0.1171875, -0.1441955566406250, 0.6765925884246826, -6.883914268109947 };
    static double b1[] = { -0.1025390625, 0.2775764465332031, -1.993531733751297, 2.724882731126854e1 };

    int i, k, m;

    *nm = n;
    if ((x < 0.0) || (n < 0)) return 1;
    if (x < 1e-15) {
        for (i = 0; i <= n; i++) {
            jn[i] = 0.0;
            yn[i] = -1e308;
            jnp[i] = 0.0;
            ynp[i] = 1e308;
        }
        jn[0] = 1.0;
        jnp[1] = 0.5;
        return 0;
    }

    if (x <= 300.0 || n > (int)(0.9*x)) {
        if (n == 0) *nm = 1;
        m = msta1(x, 200);
        if (m < *nm) *nm = m;
        else m = msta2(x, *nm, 15);

        bs = 0.0;
        su = 0.0;
        sv = 0.0;
        f2 = 0.0;
        f1 = 1.0e-100;

        for (k = m; k >= 0; k--) {
            f = 2.0*(k+1.0)/x*f1 - f2;
            if (k <= *nm) jn[k] = f;

            if ((k % 2 == 0) && (k != 0)) su += ((k & 2) - 1) * f / (double)k;
            else if (k > 1) sv += ((k & 2) - 1) * (double)k * f / (k*k - 1.0);

            f2 = f1;
            f1 = f;
        }

        s0 = bs + f;
        for (k = 0; k <= *nm; k++) {
            jn[k] /= s0;
        }

        ec = log(0.5*x) + 0.5772156649015329;
        by0 = M_2_PI * (ec*jn[0] - 4.0*su/s0);
        yn[0] = by0;
        by1 = M_2_PI * ((ec - 1.0)*jn[1] - jn[0]/x - 4.0*sv/s0);
        yn[1] = by1;
    } else {
        t1 = x - M_PI_4;
        p0 = 1.0;
        q0 = -0.125/x;
        for (k = 0; k < 4; k++) {
            p0 += a[k]*pow(x, -2*k-2);
            q0 += b[k]*pow(x, -2*k-3);
        }
        cu = sqrt(M_2_PI/x);
        bj0 = cu*(p0*cos(t1) - q0*sin(t1));
        by0 = cu*(p0*sin(t1) + q0*cos(t1));
        jn[0] = bj0;
        yn[0] = by0;

        t2 = x - 0.75*M_PI;
        p1 = 1.0;
        q1 = 0.375/x;
        for (k = 0; k < 4; k++) {
            p1 += a1[k]*pow(x, -2*k-2);
            q1 += b1[k]*pow(x, -2*k-3);
        }
        bj1 = cu*(p1*cos(t2) - q1*sin(t2));
        by1 = cu*(p1*sin(t2) + q1*cos(t2));
        jn[1] = bj1;
        yn[1] = by1;

        for (k = 2; k <= *nm; k++) {
            bjk = 2.0*(k-1.0)*bj1/x - bj0;
            jn[k] = bjk;
            bj0 = bj1;
            bj1 = bjk;
        }
    }

    jnp[0] = -jn[1];
    for (k = 1; k <= *nm; k++) {
        jnp[k] = jn[k-1] - k*jn[k]/x;
    }

    for (k = 2; k <= *nm; k++) {
        byk = 2.0*(k-1.0)*by1/x - by0;
        yn[k] = byk;
        by0 = by1;
        by1 = byk;
    }

    ynp[0] = -yn[1];
    for (k = 1; k <= *nm; k++) {
        ynp[k] = yn[k-1] - k*yn[k]/x;
    }

    return 0;
}