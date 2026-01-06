#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>  

#define EPS 1e-15

#ifndef EL
#define EL 0.5772156649015328606
#endif

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

#define TINY 1.0e-20

#define CONE 1.0 + 0.0*I
#define CZERO 0.0 + 0.0*I
#define CII 1.0*I

#define MAXN 301
#define ITMAX 100    
#define MAXIT 40       

#endif 
