#ifndef BESSH
#define BESSH

#include <math.h>
#include <complex.h>

/* Real Bessel-related routines */
int msta1(double x, int mp);
int msta2(double x, int n, int mp);

int bessjy01a(double x, double *j0, double *j1, double *y0, double *y1,
              double *j0p, double *j1p, double *y0p, double *y1p);

int bessjy01b(double x, double *j0, double *j1, double *y0, double *y1,
              double *j0p, double *j1p, double *y0p, double *y1p);

int bessjyna(int n, double x, int *nm, double *jn, double *yn,
             double *jnp, double *ynp);

int bessjynb(int n, double x, int *nm, double *jn, double *yn,
             double *jnp, double *ynp);

int bessik01a(double x, double *i0, double *i1, double *k0, double *k1,
              double *i0p, double *i1p, double *k0p, double *k1p);

int bessik01b(double x, double *i0, double *i1, double *k0, double *k1,
              double *i0p, double *i1p, double *k0p, double *k1p);

int bessikna(int n, double x, int *nm, double *in, double *kn,
             double *inp, double *knp);

int bessiknb(int n, double x, int *nm, double *in, double *kn,
             double *inp, double *knp);

/* Complex Bessel routines */
int cbessjy01(double complex z,
              double complex *cj0, double complex *cj1,
              double complex *cy0, double complex *cy1,
              double complex *cj0p, double complex *cj1p,
              double complex *cy0p, double complex *cy1p);

int cbessjyna(int n, double complex z, int *nm,
              double complex *cj, double complex *cy,
              double complex *cjp, double complex *cyp);

int cbessjynb(int n, double complex z, int *nm,
              double complex *cj, double complex *cy,
              double complex *cjp, double complex *cyp);

int cbessik01(double complex z,
              double complex *ci0, double complex *ci1,
              double complex *ck0, double complex *ck1,
              double complex *ci0p, double complex *ci1p,
              double complex *ck0p, double complex *ck1p);

int cbessikna(int n, double complex z, int *nm,
              double complex *ci, double complex *ck,
              double complex *cip, double complex *ckp);

int cbessiknb(int n, double complex z, int *nm,
              double complex *ci, double complex *ck,
              double complex *cip, double complex *ckp);

/* Convenience wrapper functions */
double complex besJ(int n, double complex z);
double complex besJJ(int n, double complex z);
double complex besY(int n, double complex z);
double complex besYY(int n, double complex z);
double complex besI(int n, double complex z);
double complex besK(int n, double complex z);
double complex besIa(int n, double complex z);
double complex besKa(int n, double complex z);
double complex besII(int n, double complex z);
double complex besIIa(int n, double complex z);
double complex besKK(int n, double complex z);
double complex besKKa(int n, double complex z);

int besJY(int n, double complex z,
          double complex *bj, double complex *by);

int besJJYY(int n, double complex z,
            double complex *bjj, double complex *byy);

int besIK(int n, double complex z,
          double complex *bi, double complex *bk);

int besIIKK(int n, double complex z,
            double complex *bii, double complex *bkk);

#endif /* BESSH */
