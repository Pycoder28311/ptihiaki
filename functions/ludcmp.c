#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "root.h"


void ludcmp(double complex **a, int n, int *indx, double *d)
/*
Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of
a rowwise permutation of itself.
*/
{
    int i, imax, j, k;
    double big, temp;
    double complex dum, sum;
    double *vv;  /* scaling of each row */

    vv = (double *)malloc(n * sizeof(double));
    if (!vv) {
        fprintf(stderr, "Memory allocation failed in ludcmp\n");
        exit(1);
    }

    *d = 1.0;       /* No row interchanges yet */

    /* Loop over rows to get the implicit scaling information */
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++) {
            temp = cabs(a[i][j]);
            if (temp > big) big = temp;
        }
        if (big == 0.0)
            nrerror("Singular matrix in routine ludcmp");
        vv[i] = 1.0 / big;
    }

    /* Loop over columns of Crout's method */
    for (j = 1; j <= n; j++) {

        for (i = 1; i < j; i++) {
            sum = a[i][j];
            for (k = 1; k < i; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }

        big = 0.0;
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;

            temp = vv[i] * cabs(sum);
            if (temp >= big) {
                big = temp;
                imax = i;
            }
        }

        /* Pivoting */
        if (j != imax) {
            for (k = 1; k <= n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }

        indx[j] = imax;

        if (cabs(a[j][j]) == 0.0)
            a[j][j] = TINY + 0.0 * I;

        if (j != n) {
            dum = 1.0 / a[j][j];
            for (i = j + 1; i <= n; i++)
                a[i][j] *= dum;
        }
    }

    free(vv);
}

void lubksb(double complex **a, int n, int *indx, double complex b[])
/*
Solves the set of n linear equations A·X = B.
a[1..n][1..n] is the LU decomposition of A (from ludcmp).
indx[1..n] is the permutation vector.
b[1..n] is input as B and overwritten by the solution X.
*/
{
    int i, ii = 0, ip, j;
    double complex sum;

    /* Forward substitution */
    for (i = 1; i <= n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];

        if (ii) {
            for (j = ii; j <= i - 1; j++)
                sum -= a[i][j] * b[j];
        } else if (cabs(sum) != 0.0) {
            ii = i;
        }

        b[i] = sum;
    }

    /* Back substitution */
    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i + 1; j <= n; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

void Sbrlse(double complex **z1a, double complex *b, int n, double complex x[])
/* This function solves a linear system of equations */
{
    double complex a[205][205], bb[205], Csum;
    int i, i1, j, k, r, t, s, m, r1, n1, n0;
    double aa;

    /* Copy input matrix */
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            a[i][j] = z1a[i][j];
        }
    }

    n1 = n + 1;

    /* Append RHS vector */
    for (i = 1; i <= n; i++) {
        a[i][n1] = b[i];
    }

    n0 = n - 1;

    /* Forward elimination with partial pivoting */
    for (r = 1; r <= n0; r++) {
        aa = 0.0;
        k = r;

        do {
            if (aa < cabs(a[k][r])) {
                aa = cabs(a[k][r]);
                t = k;
            }
            if (k >= n) break;
            k++;
        } while (k - 1 < n);

        /* Row swap */
        for (s = r; s <= n1; s++) {
            bb[s] = a[r][s];
            a[r][s] = a[t][s];
            a[t][s] = bb[s];
        }

        r1 = r + 1;

        for (i = r1; i <= n1; i++) {
            for (m = r1; m <= n1; m++) {
                a[i][m] = a[i][m] - a[i][r] / a[r][r] * a[r][m];
            }
        }
    }

    /* Back substitution */
    x[n] = a[n][n1] / a[n][n];

    for (j = 1; j <= n0; j++) {
        i = n - j;
        Csum = 0.0 + 0.0 * I;
        i1 = i + 1;

        for (m = i1; m <= n; m++) {
            Csum += a[i][m] * x[m];
        }

        x[i] = (a[i][n1] - Csum) / a[i][i];
    }
}
