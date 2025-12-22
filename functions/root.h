#ifndef ROOT_H
#define ROOT_H

/* Root search and root finding routines */

void nrerror(const char *msg);

void riza(double (*func)(double),
          double down, double up,
          double step, double tol,
          int *nr, double rizes[]);

void zbrak(double (*fx)(double),
           double x1, double x2,
           int n,
           double xb1[], double xb2[],
           int *nb);

double zbrent(double (*func)(double),
              double x1, double x2,
              double tol);

double rtflsp(double (*func)(double),
              double x1, double x2,
              double xacc);

#endif /* ROOT_H */
