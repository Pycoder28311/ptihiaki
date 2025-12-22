#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>  /* for fabs, etc. */

/* Machine precision and tolerances */
#define EPS 1e-15

// constants.h
#ifndef EL
#define EL 0.5772156649015328606
#endif

/* Pi-related constants (define only if not already in math.h) */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_2_PI
#define M_2_PI 0.63661977236758134308   /* 2/pi */
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923  /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4 0.78539816339744830962  /* pi/4 */
#endif

/* Tiny number to avoid division by zero */
#define TINY 1.0e-20

#define CONE 1.0 + 0.0*I
#define CZERO 0.0 + 0.0*I
#define CII 1.0*I

#define MAXN 301
#define ITMAX 100        /* Maximum allowed number of iterations */
#define MAXIT 40        /* Maximum allowed number of iterations */

#endif /* CONSTANTS_H */
