#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"

int cbessik01(double complex z,
              double complex *ci0, double complex *ci1,
              double complex *ck0, double complex *ck1,
              double complex *ci0p, double complex *ci1p,
              double complex *ck0p, double complex *ck1p)
{
    double complex z1, z2, zr, zr2, cr, ca, cb, cs, ct, cw=0.0+0.0*I;
    double a0, w0;
    int k, kz;

    static double a[] = {
        0.125, 0.0703125, 0.0732421875, 0.11215209960938,
        0.22710800170898, 0.57250142097473, 1.7277275025845,
        6.0740420012735, 24.380529699556, 110.01714026925,
        551.33589612202, 3038.0905109224
    };

    static double b[] = {
        -0.375, -0.1171875, -0.1025390625, -0.14419555664063,
        -0.27757644653320, -0.67659258842468, -1.9935317337513,
        -6.8839142681099, -27.248827311269, -121.59789187654,
        -603.84407670507, -3302.2722944809
    };

    static double a1[] = {
        0.125, 0.2109375, 1.0986328125, 11.775970458984,
        214.61706161499, 5951.1522710323, 233476.45606175,
        12312234.987631
    };

    a0 = cabs(z);
    z2 = z*z;
    z1 = z;

    if (a0 == 0.0) {
        *ci0 = CONE;
        *ci1 = CZERO;
        *ck0 = 1e308 + 0.0*I;
        *ck1 = 1e308 + 0.0*I;
        *ci0p = CZERO;
        *ci1p = 0.5 + 0.0*I;
        *ck0p = -1e308 + 0.0*I;
        *ck1p = -1e308 + 0.0*I;
        return 0;
    }

    if (creal(z) < 0.0) z1 = -z;

    if (a0 <= 18.0) {
        *ci0 = CONE;
        cr = CONE;
        for (k = 1; k <= 50; k++) {
            cr *= 0.25*z2/(double)(k*k);
            *ci0 += cr;
            if (cabs(cr / *ci0) < EPS) break;
        }

        *ci1 = CONE;
        cr = CONE;
        for (k = 1; k <= 50; k++) {
            cr *= 0.25*z2/(double)(k*(k+1));
            *ci1 += cr;
            if (cabs(cr / *ci1) < EPS) break;
        }
        *ci1 *= 0.5*z1;
    }
    else {
        if (a0 >= 50.0) kz = 7;
        else if (a0 >= 35.0) kz = 9;
        else kz = 12;

        ca = cexp(z1)/csqrt(2.0*M_PI*z1);
        *ci0 = CONE;
        zr = 1.0/z1;
        for (k = 0; k < kz; k++) {
            *ci0 += a[k]*cpow(zr,k+1.0);
        }
        *ci0 *= ca;

        *ci1 = CONE;
        for (k = 0; k < kz; k++) {
            *ci1 += b[k]*cpow(zr,k+1.0);
        }
        *ci1 *= ca;
    }

    if (a0 <= 9.0) {
        cs = CZERO;
        ct = -clog(0.5*z1)-EL;
        w0 = 0.0;
        cr = CONE;
        for (k = 1; k <= 50; k++) {
            w0 += 1.0/k;
            cr *= 0.25*z2/(double)(k*k);
            cs += cr*(w0+ct);
            if (cabs((cs-cw)/cs) < EPS) break;
            cw = cs;
        }
        *ck0 = ct + cs;
    }
    else {
        cb = 0.5/z1;
        zr2 = 1.0/z2;
        *ck0 = CONE;
        for (k = 0; k < 10; k++) {
            *ck0 += a1[k]*cpow(zr2,k+1.0);
        }
        *ck0 *= cb/(*ci0);
    }

    *ck1 = (1.0/z1 - (*ci1)*(*ck0))/(*ci0);

    if (creal(z) < 0.0) {
        if (cimag(z) < 0.0) {
            *ck0 += CII*M_PI*(*ci0);
            *ck1 = -(*ck1) + CII*M_PI*(*ci1);
        } else if (cimag(z) > 0.0) {
            *ck0 -= CII*M_PI*(*ci0);
            *ck1 = -(*ck1) - CII*M_PI*(*ci1);
        }
        *ci1 = -(*ci1);
    }

    *ci0p = *ci1;
    *ci1p = *ci0 - *ci1/z;
    *ck0p = -(*ck1);
    *ck1p = -(*ck0) - (*ck1)/z;

    return 0;
}

int msta1(double x, int mp); 
int msta2(double x, int n, int mp);

int cbessikna(int n, double complex z, int *nm,
              double complex *ci, double complex *ck,
              double complex *cip, double complex *ckp)
{
    double complex ci0, ci1, ck0, ck1, cf, cf1, cf2, cs, ckk;
    double a0;
    int k, m, ecode;

    a0 = cabs(z);
    *nm = n;

    if (a0 < 1.0e-100) {
        for (k = 0; k <= n; k++) {
            ci[k] = CZERO;
            ck[k] = -1e308 + 0.0*I;
            cip[k] = CZERO;
            ckp[k] = 1e308 + 0.0*I;
        }
        ci[0] = 1e308 + 0.0*I;
        cip[1] = 0.5 + 0.0*I;
        return 0;
    }

    ecode = cbessik01(z, &ci[0], &ci[1], &ck[0], &ck[1], &cip[0], &cip[1], &ckp[0], &ckp[1]);
    if (n < 2) return 0;

    ci0 = ci[0];
    ci1 = ci[1];
    ck0 = ck[0];
    ck1 = ck[1];

    m = msta1(a0, 200);
    if (m < n) *nm = m;
    else m = msta2(a0, n, 15);

    cf2 = CZERO;
    cf1 = 1.0e-100 + 0.0*I;

    for (k = m; k >= 0; k--) {
        cf = 2.0*(k+1)/z*cf1 + cf2;
        if (k <= *nm) ci[k] = cf;
        cf2 = cf1;
        cf1 = cf;
    }

    cs = ci0 / cf;
    for (k = 0; k <= *nm; k++) {
        ci[k] *= cs;
    }

    for (k = 2; k <= *nm; k++) {
        if (cabs(ci[k-1]) > cabs(ci[k-2])) {
            ckk = (1.0/z - ci[k]*ck[k-1]) / ci[k-1];
        } else {
            ckk = (ci[k]*ck[k-2] + 2.0*(k-1)/(z*z)) / ci[k-2];
        }
        ck[k] = ckk;
    }

    for (k = 2; k <= *nm; k++) {
        cip[k] = ci[k-1] - (double)k*ci[k]/z;
        ckp[k] = -ck[k-1] - (double)k*ck[k]/z;
    }

    return 0;
}

int cbessiknb(int n, double complex z, int *nm,
              double complex *ci, double complex *ck,
              double complex *cip, double complex *ckp)
{
    double complex z1, cbs, csk0, cf, cf0, cf1, ca0, cbkl;
    double complex cg, cg0, cg1, cs0, cr;
    double a0, vt, fac;
    int k, kz, l, m;

    a0 = cabs(z);
    *nm = n;

    if (a0 < 1.0e-100) {
        for (k = 0; k <= n; k++) {
            ci[k] = CZERO;
            ck[k] = 1e308 + 0.0*I;
            cip[k] = CZERO;
            ckp[k] = -1e308 + 0.0*I;
        }
        ci[0] = 1.0 + 0.0*I;
        cip[1] = 0.5 + 0.0*I;
        return 0;
    }

    z1 = z;
    if (creal(z) < 0.0) z1 = -z;

    if (n == 0) *nm = 1;

    m = msta1(a0, 200);
    if (m < *nm) *nm = m;
    else m = msta2(a0, *nm, 15);

    cbs = CZERO;
    csk0 = CZERO;
    cf0 = CZERO;
    cf1 = 1.0e-100 + 0.0*I;

    for (k = m; k >= 0; k--) {
        cf = 2.0*(k+1.0)/z1*cf1 + cf0;
        if (k <= *nm) ci[k] = cf;
        if ((k != 0) && (k == 2*(k>>1))) csk0 += 4.0*cf/(double)k;
        cbs += 2.0*cf;
        cf0 = cf1;
        cf1 = cf;
    }

    cs0 = cexp(z1)/(cbs - cf);
    for (k = 0; k <= *nm; k++) ci[k] *= cs0;

    if (a0 <= 9.0) {
        ck[0] = -(clog(0.5*z1) + EL)*ci[0] + cs0*csk0;
        ck[1] = (1.0/z1 - ci[1]*ck[0])/ci[0];
    } else {
        ca0 = csqrt(M_PI_2/z1)*cexp(-z1);
        if (a0 >= 200.0) kz = 6;
        else if (a0 >= 80.0) kz = 8;
        else if (a0 >= 25.0) kz = 10;
        else kz = 16;

        for (l = 0; l < 2; l++) {
            cbkl = CONE;
            vt = 4.0*l;
            cr = CONE;
            for (k = 1; k <= kz; k++) {
                cr *= 0.125*(vt - pow(2.0*k-1.0,2.0))/((double)k*z);
                cbkl += cr;
            }
            ck[l] = ca0*cbkl;
        }
    }

    cg0 = ck[0];
    cg1 = ck[1];
    for (k = 2; k <= *nm; k++) {
        cg = 2.0*(k-1.0)*cg1/z1 + cg0;
        ck[k] = cg;
        cg0 = cg1;
        cg1 = cg;
    }

    if (creal(z) < 0.0) {
        fac = 1.0;
        for (k = 0; k <= *nm; k++) {
            if (cimag(z) < 0.0) {
                ck[k] = fac*ck[k] + CII*M_PI*ci[k];
            } else {
                ck[k] = fac*ck[k] - CII*M_PI*ci[k];
            }
            ci[k] *= fac;
            fac = -fac;
        }
    }

    cip[0] = - ci[1];
    ck[0] = -ck[1];
    for (k = 1; k <= *nm; k++) {
        cip[k] = ci[k-1] - (double)k*ci[k]/z;
        ckp[k] = -ck[k-1] - (double)k*ck[k]/z;
    }

    return 0;
}

double complex besI(int n, double complex z) {
    double complex bsi[MAXN], bsk[MAXN], bsii[MAXN], bskk[MAXN];
    int nm;
    cbessikna(n, z, &nm, bsi, bsk, bsii, bskk);
    return bsi[n];
}

double complex besIa(int n, double complex z) {
    double complex sum = 1.0 + 0.0*I;
    double complex term = 1.0 + 0.0*I;
    int m = 4 * n * n;

    for(int i = 1; i <= 100; i++) {
        term = -term * ((double)m - pow(2.0*i-1.0, 2.0)) / (z * 8.0 * (double)i);
        sum += term;
        if(cabs(term/sum) < EPS) break;
    }

    return sum / csqrt(2.0 * M_PI * z);
}

double complex besII(int n, double complex z) {
    double complex bsi[MAXN], bsk[MAXN], bsii[MAXN], bskk[MAXN];
    int nm;
    cbessikna(n, z, &nm, bsi, bsk, bsii, bskk);
    return bsii[n];
}

double complex besIIa(int n, double complex z) {
    // Χρησιμοποιούμε την αναδιατυπωμένη σχέση: I_n' = I_{n+1} + n/z * I_n
    extern double complex besIa(int n, double complex z); // πρέπει να υπάρχει η besIa σε C
    return besIa(n+1, z) + (double)n / z * besIa(n, z);
}

double complex besK(int n, double complex z)
{
        double complex bsi[301],bsk[301],bsii[301],bskk[301];
        cbessikna(n,z,&n,bsi,bsk,bsii,bskk);
        return bsk[n];
}

double complex besKa(int n, double complex z)
{
        double complex sum=1.0 + 0.0*I;
        double complex term=1.0 + 0.0*I;
        
        int m=4*n*n;
        
        for(int i=1;i<=100;i++)
        {
                term = term * ((double)m - pow(2.0*(double)i - 1.0, 2.0)) / (z * 8.0 * (double)i);
                sum+=term;
                if(abs(term/sum)<1.e-16) break;
        }
        return sum*sqrt(M_PI/2.0/z);
}

double complex besKK(int n, double complex z) {
    double complex bsi[MAXN], bsk[MAXN], bsii[MAXN], bskk[MAXN];
    int nm;
    cbessikna(n, z, &nm, bsi, bsk, bsii, bskk);
    return bskk[n];
}

double complex besKKa(int n, double complex z) {
    // Χρησιμοποιούμε τη σχέση: K_n' = -K_{n+1} + n/z * K_n
    extern double complex besKa(int n, double complex z); // πρέπει να υπάρχει η besKa σε C
    return -besKa(n+1, z) + (double)n / z * besKa(n, z);
}

int besIK(int n, double complex z, double complex *bi, double complex *bk) {
    double complex bsi[MAXN], bsk[MAXN], bsii[MAXN], bskk[MAXN];
    int err;
    err = cbessikna(n, z, &n, bsi, bsk, bsii, bskk);
    *bi = bsi[n];
    *bk = bsk[n];
    return err;
}

int besIIKK(int n, double complex z, double complex *bii, double complex *bkk) {
    double complex bsi[MAXN], bsk[MAXN], bsii[MAXN], bskk[MAXN];
    int err;
    err = cbessikna(n, z, &n, bsi, bsk, bsii, bskk);
    *bii = bsii[n];
    *bkk = bskk[n];
    return err;
}