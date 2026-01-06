#include <math.h>
#include "constants.h"

int msta1(double x, int mp);
int msta2(double x, int n, int mp);

int bessik01a(double x,
              double *i0, double *i1,
              double *k0, double *k1,
              double *i0p, double *i1p,
              double *k0p, double *k1p)
{
    double r, x2, ca, cb, ct, ww, w0, xr, xr2;
    int k, kz;

    static double a[] = {
        0.125,
        7.03125e-2,
        7.32421875e-2,
        1.1215209960938e-1,
        2.2710800170898e-1,
        5.7250142097473e-1,
        1.7277275025845,
        6.0740420012735,
        2.4380529699556e1,
        1.1001714026925e2,
        5.5133589612202e2,
        3.0380905109224e3
    };

    static double b[] = {
        -0.375,
        -1.171875e-1,
        -1.025390625e-1,
        -1.4419555664063e-1,
        -2.7757644653320e-1,
        -6.7659258842468e-1,
        -1.9935317337513,
        -6.8839142681099,
        -2.7248827311269e1,
        -1.2159789187654e2,
        -6.0384407670507e2,
        -3.3022722944809e3
    };

    static double a1[] = {
        0.125,
        0.2109375,
        1.0986328125,
        1.1775970458984e1,
        2.1461706161499e2,
        5.9511522710323e3,
        2.3347645606175e5,
        1.2312234987631e7
    };

    if (x < 0.0) return 1;

    if (x == 0.0) {
        *i0  = 1.0;
        *i1  = 0.0;
        *k0  = 1e300;
        *k1  = 1e300;
        *i0p = 0.0;
        *i1p = 0.5;
        *k0p = -1e300;
        *k1p = -1e300;
        return 0;
    }

    x2 = x * x;

    if (x <= 18.0) {
        *i0 = 1.0;
        r = 1.0;
        for (k = 1; k <= 50; k++) {
            r *= 0.25 * x2 / (k * k);
            *i0 += r;
            if (fabs(r / (*i0)) < EPS) break;
        }

        *i1 = 1.0;
        r = 1.0;
        for (k = 1; k <= 50; k++) {
            r *= 0.25 * x2 / (k * (k + 1));
            *i1 += r;
            if (fabs(r / (*i1)) < EPS) break;
        }
        *i1 *= 0.5 * x;
    }
    else {
        if (x >= 50.0) kz = 7;
        else if (x >= 35.0) kz = 9;
        else kz = 12;

        ca = exp(x) / sqrt(2.0 * M_PI * x);
        xr = 1.0 / x;

        *i0 = 1.0;
        for (k = 0; k < kz; k++)
            *i0 += a[k] * pow(xr, k + 1);
        *i0 *= ca;

        *i1 = 1.0;
        for (k = 0; k < kz; k++)
            *i1 += b[k] * pow(xr, k + 1);
        *i1 *= ca;
    }

    if (x <= 9.0) {
        ct = -(log(0.5 * x) + EL);
        *k0 = 0.0;
        w0 = 0.0;
        r = 1.0;
        ww = 0.0;

        for (k = 1; k <= 50; k++) {
            w0 += 1.0 / k;
            r *= 0.25 * x2 / (k * k);
            *k0 += r * (w0 + ct);
            if (fabs((*k0 - ww) / (*k0)) < EPS) break;
            ww = *k0;
        }
        *k0 += ct;
    }
    else {
        cb = 0.5 / x;
        xr2 = 1.0 / x2;
        *k0 = 1.0;
        for (k = 0; k < 8; k++)
            *k0 += a1[k] * pow(xr2, k + 1);
        *k0 *= cb / (*i0);
    }

    *k1  = (1.0 / x - (*i1) * (*k0)) / (*i0);
    *i0p = *i1;
    *i1p = *i0 - (*i1) / x;
    *k0p = -(*k1);
    *k1p = -(*k0) - (*k1) / x;

    return 0;
}

int bessik01b(double x,
              double *i0, double *i1,
              double *k0, double *k1,
              double *i0p, double *i1p,
              double *k0p, double *k1p)
{
    double t, t2, dtmp, dtmp1;

    if (x < 0.0) return 1;

    if (x == 0.0) {
        *i0  = 1.0;
        *i1  = 0.0;
        *k0  = 1e300;
        *k1  = 1e300;
        *i0p = 0.0;
        *i1p = 0.5;
        *k0p = -1e300;
        *k1p = -1e300;
        return 0;
    }

    if (x < 3.75) {
        t = x / 3.75;
        t2 = t * t;
        *i0 = (((((0.0045813*t2 + 0.0360768)*t2 + 0.2659732)*t2 +
                1.2067492)*t2 + 3.0899424)*t2 + 3.5156229)*t2 + 1.0;
        *i1 = x*(((((0.00032411*t2+0.00301532)*t2+0.02658733*t2+
                0.15084934)*t2+0.51498869)*t2+0.87890594)*t2+0.5);
    }
    else {
        t = 3.75 / x;
        dtmp1 = exp(x) / sqrt(x);

        dtmp = (((((((0.00392377*t - 0.01647633)*t + 0.026355537)*t -
                  0.02057706)*t + 0.00916281)*t - 0.00157565)*t +
                  0.00225319)*t + 0.01328592)*t + 0.39894228;
        *i0 = dtmp * dtmp1;

        dtmp = (((((((-0.00420059*t + 0.01787654)*t - 0.02895312)*t +
                  0.02282967)*t - 0.01031555)*t + 0.00163801)*t -
                  0.00362018)*t - 0.03988024)*t + 0.39894228;
        *i1 = dtmp * dtmp1;
    }

    if (x < 2.0) {
        t = 0.5 * x;
        t2 = t * t;

        dtmp = (((((0.0000074*t2 + 0.0001075)*t2 + 0.00262698)*t2 +
                 0.0348859)*t2 + 0.23069756)*t2 + 0.4227842)*t2 - 0.57721566;
        *k0 = dtmp - (*i0) * log(t);

        dtmp = (((((-0.00004686*t2 - 0.00110404)*t2 - 0.01919402)*t2 -
                 0.18156897)*t2 - 0.67278578)*t2 + 0.15443144)*t2 + 1.0;
        *k1 = dtmp / x + (*i1) * log(t);
    }
    else {
        t = 2.0 / x;
        dtmp1 = exp(-x) / sqrt(x);

        dtmp = (((((0.00053208*t - 0.0025154)*t + 0.00587872)*t -
                 0.01062446)*t + 0.02189568)*t - 0.07832358)*t + 1.25331414;
        *k0 = dtmp * dtmp1;

        dtmp = (((((-0.00068245*t + 0.00325614)*t - 0.00780353)*t +
                 0.01504268)*t - 0.0365562)*t + 0.23498619)*t + 1.25331414;
        *k1 = dtmp * dtmp1;
    }

    *i0p = *i1;
    *i1p = *i0 - (*i1) / x;
    *k0p = -(*k1);
    *k1p = -(*k0) - (*k1) / x;

    return 0;
}

int bessikna(int n, double x, int *nm,
             double *in, double *kn,
             double *inp, double *knp)
{
    double bi0, bi1, bk0, bk1;
    double g, g0, g1, f, f0, f1;
    double h, h0, h1, s0;
    int k, m, ecode;

    if ((x < 0.0) || (n < 0)) return 1;

    if (x < EPS) {
        for (k = 0; k <= n; k++) {
            in[k]  = 0.0;
            kn[k]  = 1.0e300;
            inp[k] = 0.0;
            knp[k] = -1.0e300;
        }
        in[0]  = 1.0;
        inp[1] = 0.5;
        *nm = n;
        return 0;
    }

    *nm = n;

    ecode = bessik01a(x,
                      &in[0], &in[1],
                      &kn[0], &kn[1],
                      &inp[0], &inp[1],
                      &knp[0], &knp[1]);

    if (ecode != 0) return ecode;
    if (n < 2) return 0;

    bi0 = in[0];  bi1 = in[1];
    bk0 = kn[0];  bk1 = kn[1];

    if ((x > 40.0) && (n < (int)(0.25 * x))) {
        h0 = bi0;
        h1 = bi1;
        for (k = 2; k <= n; k++) {
            h = -2.0 * (k - 1.0) * h1 / x + h0;
            in[k] = h;
            h0 = h1;
            h1 = h;
        }
    }
    else {
        m = msta1(x, 200);
        if (m < n) *nm = m;
        else m = msta2(x, n, 15);

        f0 = 0.0;
        f1 = 1.0e-100;

        for (k = m; k >= 0; k--) {
            f = 2.0 * (k + 1.0) * f1 / x + f0;
            if (k <= *nm) in[k] = f;
            f0 = f1;
            f1 = f;
        }

        s0 = bi0 / f;
        for (k = 0; k <= *nm; k++)
            in[k] *= s0;
    }

    g0 = bk0;
    g1 = bk1;
    for (k = 2; k <= *nm; k++) {
        g = 2.0 * (k - 1.0) * g1 / x + g0;
        kn[k] = g;
        g0 = g1;
        g1 = g;
    }

    for (k = 2; k <= *nm; k++) {
        inp[k] = in[k - 1] - k * in[k] / x;
        knp[k] = -kn[k - 1] - k * kn[k] / x;
    }

    return 0;
}

int bessiknb(int n, double x, int *nm,
             double *in, double *kn,
             double *inp, double *knp)
{
    double s0, bs, f, f0, f1, sk0;
    double a0, bkl, vt, r, g, g0, g1;
    int k, kz, m, l;

    if ((x < 0.0) || (n < 0)) return 1;

    if (x < EPS) {
        for (k = 0; k <= n; k++) {
            in[k]  = 0.0;
            kn[k]  = 1.0e300;
            inp[k] = 0.0;
            knp[k] = -1.0e300;
        }
        in[0]  = 1.0;
        inp[1] = 0.5;
        *nm = n;
        return 0;
    }

    *nm = n;
    if (n == 0) *nm = 1;

    m = msta1(x, 200);
    if (m < *nm) *nm = m;
    else m = msta2(x, *nm, 15);

    bs = 0.0;
    sk0 = 0.0;
    f0 = 0.0;
    f1 = 1.0e-100;

    for (k = m; k >= 0; k--) {
        f = 2.0 * (k + 1.0) * f1 / x + f0;
        if (k <= *nm) in[k] = f;
        if ((k != 0) && (k == 2 * (k / 2))) {
            sk0 += 4.0 * f / k;
        }
        bs += 2.0 * f;
        f0 = f1;
        f1 = f;
    }

    s0 = exp(x) / (bs - f);
    for (k = 0; k <= *nm; k++)
        in[k] *= s0;

    if (x <= 8.0) {
        kn[0] = -(log(0.5 * x) + EL) * in[0] + s0 * sk0;
        kn[1] = (1.0 / x - in[1] * kn[0]) / in[0];
    }
    else {
        a0 = sqrt(M_PI_2 / x) * exp(-x);

        if (x >= 200.0) kz = 6;
        else if (x >= 80.0) kz = 8;
        else if (x >= 25.0) kz = 10;
        else kz = 16;

        for (l = 0; l < 2; l++) {
            bkl = 1.0;
            vt = 4.0 * l;
            r = 1.0;
            for (k = 1; k <= kz; k++) {
                r *= 0.125 * (vt - pow(2.0 * k - 1.0, 2)) / (k * x);
                bkl += r;
            }
            kn[l] = a0 * bkl;
        }
    }

    g0 = kn[0];
    g1 = kn[1];
    for (k = 2; k <= *nm; k++) {
        g = 2.0 * (k - 1.0) * g1 / x + g0;
        kn[k] = g;
        g0 = g1;
        g1 = g;
    }

    inp[0] = in[1];
    knp[0] = -kn[1];
    for (k = 1; k <= *nm; k++) {
        inp[k] = in[k - 1] - k * in[k] / x;
        knp[k] = -kn[k - 1] - k * kn[k] / x;
    }

    return 0;
}
