#include "head.h"
#include "globals.h"
#include <stdio.h>

// Kronecker delta function
int delta(int n, int q) {
    return (n == q) ? 1 : 0;
}

double k_of_n(int n, int m) {
    return m + n * N;
}
