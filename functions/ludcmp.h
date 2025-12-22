#include <complex.h>

/* LU decomposition */
void ludcmp(double complex **aa, int n, int *indx, double *d);

/* LU back substitution */
void lubksb(double complex **a, int n, int *indx, double complex b[]);

/* Linear system solver */
void Sbrlse(double complex **z1a, double complex *b, int n, double complex x[]);
