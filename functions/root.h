#ifndef ROOT_H
#define ROOT_H

void nrerror(const char *msg);

void riza(double (*func)(double, double),
          double down, double up,
          double step, double tol,
          int *nr, double rizes[], double D);

void zbrak(double (*fx)(double, double),
           double x1, double x2,
           int n,
           double xb1[], double xb2[],
           int *nb, double D);

double zbrent(double (*func)(double, double),
              double x1, double x2,
              double tol, double D);

double rtflsp(double (*func)(double, double),
              double x1, double x2,
              double xacc, double D);

#endif 
