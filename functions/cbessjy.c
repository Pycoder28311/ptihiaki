#include <complex.h>
#include <math.h>

static const double complex cii = 0.0 + 1.0*I;
static const double complex cone = 1.0 + 0.0*I;
static const double complex czero = 0.0 + 0.0*I;
static const double eps = 1e-16;
static const double el = 0.5772156649015329; // Euler-Mascheroni

int cbessjy01(double complex z,
              double complex *cj0, double complex *cj1,
              double complex *cy0, double complex *cy1,
              double complex *cj0p, double complex *cj1p,
              double complex *cy0p, double complex *cy1p)
{
    double a0, w0, w1;
    int k, kz;
    double complex z1, z2, cr, cp, cs, cp0, cq0, cp1, cq1, ct1, ct2, cu;

    static double a[] = {
        -7.03125e-2, 0.112152099609375, -0.5725014209747314, 6.074042001273483,
        -1.100171402692467e2, 3.038090510922384e3, -1.188384262567832e5,
        6.252951493434797e6, -4.259392165047669e8, 3.646840080706556e10,
        -3.833534661393944e12, 4.854014686852901e14, -7.286857349377656e16,
        1.279721941975975e19
    };
    static double b[] = {
         7.32421875e-2, -0.2271080017089844, 1.727727502584457, -2.438052969955606e1,
         5.513358961220206e2, -1.825775547429318e4, 8.328593040162893e5,
        -5.006958953198893e7, 3.836255180230433e9, -3.649010818849833e11,
         4.218971570284096e13, -5.827244631566907e15, 9.476288099260110e17,
        -1.792162323051699e20
    };
    static double a1[] = {
         0.1171875, -0.1441955566406250, 0.6765925884246826, -6.883914268109947,
         1.215978918765359e2, -3.302272294480852e3, 1.276412726461746e5,
        -6.656367718817688e6, 4.502786003050393e8, -3.833857520742790e10,
         4.011838599133198e12, -5.060568503314727e14, 7.572616461117958e16,
        -1.326257285320556e19
    };
    static double b1[] = {
        -0.1025390625, 0.2775764465332031, -1.993531733751297, 2.724882731126854e1,
        -6.038440767050702e2, 1.971837591223663e4, -8.902978767070678e5,
        5.310411010968522e7, -4.043620325107754e9, 3.827011346598605e11,
        -4.406481417852278e13, 6.065091351222699e15, -9.833883876590679e17,
         1.855045211579828e20
    };

    a0 = cabs(z);
    z2 = z*z;
    z1 = z;

    if (a0 == 0.0) {
        *cj0 = cone;
        *cj1 = czero;
        *cy0 = -1e308 + 0.0*I;
        *cy1 = -1e308 + 0.0*I;
        *cj0p = czero;
        *cj1p = 0.5 + 0.0*I;
        *cy0p = 1e308 + 0.0*I;
        *cy1p = 1e308 + 0.0*I;
        return 0;
    }

    if (creal(z) < 0.0) z1 = -z;

    if (a0 <= 12.0) {
        *cj0 = cone;
        cr = cone;
        for (k=1; k<=40; k++) {
            cr *= -0.25*z2/(double)(k*k);
            *cj0 += cr;
            if (cabs(cr) < cabs(*cj0)*eps) break;
        }
        *cj1 = cone;
        cr = cone;
        for (k=1; k<=40; k++) {
            cr *= -0.25*z2/(k*(k+1.0));
            *cj1 += cr;
            if (cabs(cr) < cabs(*cj1)*eps) break;
        }
        *cj1 *= 0.5*z1;

        w0 = 0.0;
        cr = cone;
        cs = czero;
        for (k=1; k<=40; k++) {
            w0 += 1.0/k;
            cr *= -0.25*z2/(double)(k*k);
            cp = cr*w0;
            cs += cp;
            if (cabs(cp) < cabs(cs)*eps) break;
        }
        *cy0 = (2.0/M_PI)*((log(0.5*z1)+el)*(*cj0) - cs);

        w1 = 0.0;
        cr = cone;
        cs = cone;
        for (k=1; k<=40; k++) {
            w1 += 1.0/k;
            cr *= -0.25*z2/(k*(k+1.0));
            cp = cr*(2.0*w1 + 1.0/(k+1.0));
            cs += cp;
            if (cabs(cp) < cabs(cs)*eps) break;
        }
        *cy1 = (2.0/M_PI)*((log(0.5*z1)+el)*(*cj1) - 1.0/z1 - 0.25*z1*cs);
    } else {
        if (a0 >= 50.0) kz = 8;
        else if (a0 >= 35.0) kz = 10;
        else kz = 12;

        ct1 = z1 - M_PI_4;
        cp0 = cone;
        for (k=0; k<kz; k++) cp0 += a[k]*cpow(z1,-2.0*k-2.0);

        cq0 = -0.125/z1;
        for (k=0; k<kz; k++) cq0 += b[k]*cpow(z1,-2.0*k-3.0);

        cu = csqrt(2.0/(M_PI*z1));
        *cj0 = cu*(cp0*ccos(ct1) - cq0*csin(ct1));
        *cy0 = cu*(cp0*csin(ct1) + cq0*ccos(ct1));

        ct2 = z1 - 0.75*M_PI;
        cp1 = cone;
        for (k=0; k<kz; k++) cp1 += a1[k]*cpow(z1,-2.0*k-2.0);

        cq1 = 0.375/z1;
        for (k=0; k<kz; k++) cq1 += b1[k]*cpow(z1,-2.0*k-3.0);

        *cj1 = cu*(cp1*ccos(ct2) - cq1*csin(ct2));
        *cy1 = cu*(cp1*csin(ct2) + cq1*ccos(ct2));
    }

    if (creal(z) < 0.0) {
        if (cimag(z) < 0.0) {
            *cy0 -= 2.0*cii*(*cj0);
            *cy1 = -(*cy1 - 2.0*cii*(*cj1));
        } else if (cimag(z) > 0.0) {
            *cy0 += 2.0*cii*(*cj0);
            *cy1 = -(*cy1 + 2.0*cii*(*cj1));
        }
        *cj1 = -(*cj1);
    }

    *cj0p = -(*cj1);
    *cj1p = *cj0 - *cj1/z;
    *cy0p = -(*cy1);
    *cy1p = *cy0 - *cy1/z;

    return 0;
}

static const double complex cone = 1.0 + 0.0*I;

int msta1(double a0, int max) { /* Υλοποίηση κατάλληλη */ return max; }
int msta2(double a0, int n, int m) { /* Υλοποίηση κατάλληλη */ return n; }

int cbessjyna(int n, double complex z, int *nm,
              double complex *cj, double complex *cy,
              double complex *cjp, double complex *cyp)
{
    double a0, yak, ya1, ya0, wa;
    int m, k, lb, lb0;

    double complex cbj0, cbj1, cby0, cby1, cj0, cjk, cj1tmp;
    double complex cf, cf1, cf2, cs, cg0, cg1, cyk, cyl1, cyl2, cylk;
    double complex cp11, cp12, cp21, cp22;
    double complex ch0, ch1, ch2;

    if (n < 0) return 1;
    a0 = cabs(z);
    *nm = n;

    if (a0 < 1.0e-100) {
        for (k=0; k<=n; k++) {
            cj[k] = czero;
            cy[k] = -1e308 + 0.0*I;
            cjp[k] = czero;
            cyp[k] = 1e308 + 0.0*I;
        }
        cj[0] = cone;
        cjp[1] = 0.5 + 0.0*I;
        return 0;
    }

    cbessjy01(z, &cj[0], &cj[1], &cy[0], &cy[1], &cjp[0], &cjp[1], &cyp[0], &cyp[1]);

    cbj0 = cj[0];
    cbj1 = cj[1];
    cby0 = cy[0];
    cby1 = cy[1];

    if (n <= 1) return 0;

    if (n < (int)(0.25*a0)) {
        cj0 = cbj0;
        cj1tmp = cbj1;
        for (k=2; k<=n; k++) {
            cjk = 2.0*(k-1.0)*cj1tmp/z - cj0;
            cj[k] = cjk;
            cj0 = cj1tmp;
            cj1tmp = cjk;
        }
    } else {
        m = msta1(a0, 200);
        if (m < n) *nm = m;
        else m = msta2(a0, n, 15);

        cf2 = czero;
        cf1 = 1.0e-100 + 0.0*I;

        for (k=m; k>=0; k--) {
            cf = 2.0*(k+1.0)*cf1/z - cf2;
            if (k <= *nm) cj[k] = cf;
            cf2 = cf1;
            cf1 = cf;
        }

        if (cabs(cbj0) > cabs(cbj1)) cs = cbj0/cf;
        else cs = cbj1/cf2;

        for (k=0; k<=*nm; k++) cj[k] *= cs;
    }

    for (k=2; k<=*nm; k++) {
        cjp[k] = cj[k-1] - (double)k*cj[k]/z;
    }

    ya0 = cabs(cby0);
    lb = 0;
    cg0 = cby0;
    cg1 = cby1;

    for (k=2; k<=*nm; k++) {
        cyk = 2.0*(k-1.0)*cg1/z - cg0;
        yak = cabs(cyk);
        ya1 = cabs(cg0);
        if ((yak < ya0) && (yak < ya1)) lb = k;
        cy[k] = cyk;
        cg0 = cg1;
        cg1 = cyk;
    }

    lb0 = 0;
    if ((lb > 4) && (cimag(z) != 0.0)) {
        while (lb != lb0) {
            ch2 = cone;
            ch1 = czero;
            lb0 = lb;

            for (k=lb; k>=1; k--) {
                ch0 = 2.0*k*ch1/z - ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp12 = ch0;
            cp22 = ch2;

            ch2 = czero;
            ch1 = cone;

            for (k=lb; k>=1; k--) {
                ch0 = 2.0*k*ch1/z - ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp11 = ch0;
            cp21 = ch2;

            if (lb == *nm) cj[lb+1] = 2.0*lb*cj[lb]/z - cj[lb-1];

            if (cabs(cj[0]) > cabs(cj[1])) {
                cy[lb+1] = (cj[lb+1]*cby0 - 2.0*cp11/(M_PI*z))/cj[0];
                cy[lb]   = (cj[lb]*cby0 + 2.0*cp12/(M_PI*z))/cj[0];
            } else {
                cy[lb+1] = (cj[lb+1]*cby1 - 2.0*cp21/(M_PI*z))/cj[1];
                cy[lb]   = (cj[lb]*cby1 + 2.0*cp22/(M_PI*z))/cj[1];
            }

            cyl2 = cy[lb+1];
            cyl1 = cy[lb];

            for (k=lb-1; k>=0; k--) {
                cylk = 2.0*(k+1.0)*cyl1/z - cyl2;
                cy[k] = cylk;
                cyl2 = cyl1;
                cyl1 = cylk;
            }

            cyl1 = cy[lb];
            cyl2 = cy[lb+1];

            for (k=lb+1; k<n; k++) {
                cylk = 2.0*k*cyl2/z - cyl1;
                cy[k+1] = cylk;
                cyl1 = cyl2;
                cyl2 = cylk;
            }

            for (k=2; k<=*nm; k++) {
                wa = cabs(cy[k]);
                if (wa < cabs(cy[k-1])) lb = k;
            }
        }
    }

    for (k=2; k<=*nm; k++) {
        cyp[k] = cy[k-1] - (double)k*cy[k]/z;
    }

    return 0;
}

int cbessjynb(int n, double complex z, int *nm,
              double complex *cj, double complex *cy,
              double complex *cjp, double complex *cyp)
{
    int k, m;
    double a0, y0;
    double complex cf, cf0, cf1, cf2, cbs, csu, csv, cs0, ce;
    double complex ct1, cp0, cq0, cp1, cq1, cu, cbj0, cby0, cbj1, cby1;
    double complex cyy, cbjk, ct2;

    static double a[] = { -0.7031250000000000e-1, 0.1121520996093750, -0.5725014209747314, 6.074042001273483 };
    static double b[] = { 0.7324218750000000e-1, -0.2271080017089844, 1.727727502584457, -2.438052969955606e1 };
    static double a1[] = { 0.1171875, -0.1441955566406250, 0.6765925884246826, -6.883914268109947 };
    static double b1[] = { -0.1025390625, 0.2775764465332031, -1.993531733751297, 2.724882731126854e1 };

    y0 = cabs(cimag(z));
    a0 = cabs(z);
    *nm = n;

    if (a0 < 1.0e-100) {
        for (k = 0; k <= n; k++) {
            cj[k] = czero;
            cy[k] = -1e308 + 0.0*I;
            cjp[k] = czero;
            cyp[k] = 1e308 + 0.0*I;
        }
        cj[0] = cone;
        cjp[1] = 0.5 + 0.0*I;
        return 0;
    }

    if ((a0 <= 300.0) || (n > (int)(0.25*a0))) {
        if (n == 0) *nm = 1;
        m = msta1(a0, 200);
        if (m < *nm) *nm = m;
        else m = msta2(a0, *nm, 15);

        cbs = czero;
        csu = czero;
        csv = czero;
        cf2 = czero;
        cf1 = 1.0e-100 + 0.0*I;

        for (k = m; k >= 0; k--) {
            cf = 2.0*(k+1.0)*cf1/z - cf2;
            if (k <= *nm) cj[k] = cf;

            if (((k & 1) == 0) && (k != 0)) {
                if (y0 <= 1.0) {
                    cbs += 2.0*cf;
                } else {
                    cbs += (-1)*((k & 2)-1)*2.0*cf;
                }
                csu += (double)((-1)*((k & 2)-1))*cf/(double)k;
            } else if (k > 1) {
                csv += (double)((-1)*((k & 2)-1)*k)*cf/(double)(k*k-1.0);
            }

            cf2 = cf1;
            cf1 = cf;
        }

        if (y0 <= 1.0) cs0 = cbs + cf;
        else cs0 = (cbs + cf)/ccos(z);

        for (k = 0; k <= *nm; k++) cj[k] /= cs0;

        ce = clog(0.5*z) + el;
        cy[0] = M_2_PI*(ce*cj[0] - 4.0*csu/cs0);
        cy[1] = M_2_PI*(-cj[0]/z + (ce-1.0)*cj[1] - 4.0*csv/cs0);
    } else {
        ct1 = z - M_PI_4;
        cp0 = cone;
        for (k = 0; k < 4; k++) cp0 += a[k]*cpow(z, -2.0*k-2.0);

        cq0 = -0.125/z;
        for (k = 0; k < 4; k++) cq0 += b[k]*cpow(z, -2.0*k-3.0);

        cu = csqrt(M_2_PI/z);
        cbj0 = cu*(cp0*ccos(ct1) - cq0*csin(ct1));
        cby0 = cu*(cp0*csin(ct1) + cq0*ccos(ct1));
        cj[0] = cbj0;
        cy[0] = cby0;

        ct2 = z - 0.75*M_PI;
        cp1 = cone;
        for (k = 0; k < 4; k++) cp1 += a1[k]*cpow(z, -2.0*k-2.0);

        cq1 = 0.375/z;
        for (k = 0; k < 4; k++) cq1 += b1[k]*cpow(z, -2.0*k-3.0);

        cbj1 = cu*(cp1*ccos(ct2) - cq1*csin(ct2));
        cby1 = cu*(cp1*csin(ct2) + cq1*csin(ct2));
        cj[1] = cbj1;
        cy[1] = cby1;

        for (k = 2; k <= n; k++) {
            cbjk = 2.0*(k-1.0)*cbj1/z - cbj0;
            cj[k] = cbjk;
            cbj0 = cbj1;
            cbj1 = cbjk;
        }
    }

    cjp[0] = -cj[1];
    for (k = 1; k <= *nm; k++) {
        cjp[k] = cj[k-1] - (double)k*cj[k]/z;
    }

    if (cabs(cj[0]) > 1.0)
        cy[1] = (cj[1]*cy[0] - 2.0/(M_PI*z))/cj[0];

    for (k = 2; k <= *nm; k++) {
        if (cabs(cj[k-1]) >= cabs(cj[k-2]))
            cyy = (cj[k]*cy[k-1] - 2.0/(M_PI*z))/cj[k-1];
        else
            cyy = (cj[k]*cy[k-2] - 4.0*(k-1.0)/(M_PI*z*z))/cj[k-2];
        cy[k] = cyy;
    }

    cyp[0] = -cy[1];
    for (k = 1; k <= *nm; k++) {
        cyp[k] = cy[k-1] - (double)k*cy[k]/z;
    }

    return 0;
}

double complex besJ(int n, double complex z)
{
    double complex bsj[301], bsy[301], bsjj[301], bsyy[301];
    int nm;
    cbessjyna(n, z, &nm, bsj, bsy, bsjj, bsyy);
    return bsj[n];
}

double complex besJJ(int n, double complex z)
{
    double complex bsj[301], bsy[301], bsjj[301], bsyy[301];
    int nm;
    cbessjyna(n, z, &nm, bsj, bsy, bsjj, bsyy);
    return bsjj[n];
}

double complex besY(int n, double complex z)
{
    double complex bsj[301], bsy[301], bsjj[301], bsyy[301];
    int nm;
    cbessjyna(n, z, &nm, bsj, bsy, bsjj, bsyy);
    return bsy[n];
}

double complex besYY(int n, double complex z)
{
    double complex bsj[301], bsy[301], bsjj[301], bsyy[301];
    int nm;
    cbessjyna(n, z, &nm, bsj, bsy, bsjj, bsyy);
    return bsyy[n];
}

int besJY(int n, double complex z, double complex *bj, double complex *by)
{
    double complex bsj[301], bsy[301], bsjj[301], bsyy[301];
    int nm;
    int err = cbessjyna(n, z, &nm, bsj, bsy, bsjj, bsyy);
    *bj = bsj[n];
    *by = bsy[n];
    return err;
}

int besJJYY(int n, double complex z, double complex *bjj, double complex *byy)
{
    double complex bsj[301], bsy[301], bsjj[301], bsyy[301];
    int nm;
    int err = cbessjyna(n, z, &nm, bsj, bsy, bsjj, bsyy);
    *bjj = bsjj[n];
    *byy = bsyy[n];
    return err;
}