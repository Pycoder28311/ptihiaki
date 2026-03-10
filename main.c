#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
//#include "head.h"

#include "final/globals.h"
#include "final/utils.h"
#include "functions/root.h"

#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

#define NM (2*n_max + 1) //Μέγεθος Πίνακα

double bessel(double x, int n) {
    return (gsl_sf_bessel_Jn(n-1, x) - gsl_sf_bessel_Jn(n+1, x)) / 2.0;
}

double my_function(double x, double D) {
    return compute_det_for_k(x, D);
}

double derivative(double (*f)(double, double), double x, double h, double D) {
    return (f(x + h, D) - f(x - h, D)) / (2.0 * h);
}

double my_function_derivative(double x, double D) {
    double h = 1e-6; // μικρό βήμα για τη διαφορική προσέγγιση
    return derivative(compute_det_for_k, x, h, D);
}

int cmp_double(const void *a, const void *b) {
    double da = *(double*)a;
    double db = *(double*)b;
    if (da < db) return -1;
    else if (da > db) return 1;
    else return 0;
}

#define MAX_ROOTS 8
#define TOL 1e-4
double roots[MAX_ROOTS];
int valid_nr = 0;
int nr = 0;

int compute_all_roots(double *all_roots, int max_size, double D) {
    int nr_deriv;
    double roots_deriv[100];

    // Ρίζες της παραγώγου
    riza(my_function_derivative, down, up, step, tol_riza, &nr_deriv, roots_deriv, D);

    // Φιλτράρισμα για f(x) ~ 0
    int valid_nr_deriv = 0;
    double roots_deriv_valid[100];
    for(int i = 0; i < nr_deriv; i++) {
        double fx = compute_det_for_k(roots_deriv[i], D);
        if(fabs(fx) < TOL && fabs(roots_deriv[i]) > 1e-15) {
            roots_deriv_valid[valid_nr_deriv++] = roots_deriv[i];
        }
    }

    // Ρίζες της αρχικής συνάρτησης
    //int nr;
    double roots_initial[100];
    riza(compute_det_for_k, down, up, step, tol_riza, &nr, roots_initial, D);

    int valid_nr_initial = 0;
    double roots_initial_valid[100];
    for(int i = 0; i < nr; i++) {
        if(fabs(roots_initial[i]) > 1e-15) {
            roots_initial_valid[valid_nr_initial++] = roots_initial[i];
        }
    }

    // Ενώνουμε τις δύο λίστες και αποθηκεύουμε global
    valid_nr = 0;
    for(int i = 0; i < valid_nr_deriv; i++) {
        roots[valid_nr++] = roots_deriv_valid[i];
    }
    for(int i = 0; i < valid_nr_initial; i++) {
        if(valid_nr >= MAX_ROOTS) break;
        roots[valid_nr++] = roots_initial_valid[i];
    }

    // Αν θέλεις ταξινόμηση
    for(int i = 0; i < valid_nr - 1; i++) {
        for(int j = i+1; j < valid_nr; j++) {
            if(roots[i] > roots[j]) {
                double tmp = roots[i];
                roots[i] = roots[j];
                roots[j] = tmp;
            }
        }
    }

    // Αν θέλεις να έχεις και το all_roots array ξεχωριστά
    for(int i = 0; i < valid_nr && i < max_size; i++) {
        all_roots[i] = roots[i];
    }

    return valid_nr;
}

int main() {
    double all_roots[MAX_ROOTS];
    double D_step = 0.0005; // Βήμα για D_radius
    double D_init = 0.01; // Αρχική τιμή για D_radius
    double D_final = 0.0105; // Τελική τιμή για D_radius

    double nu = 100.0;
    double xp = 1.0;
    double x = 3.0;

    // Υπολογισμός Yν(x) και Yν+1(x)
    double Ynu_1 = gsl_sf_bessel_Ynu(nu - 1, xp);
    double Ynu0 = gsl_sf_bessel_Ynu(nu, x);
    double Ynu1 = gsl_sf_bessel_Ynu(nu + 1, xp);
    // Παράγωγος με τύπο αναδρομής
    double Ynu_prime = 0.5 * (Ynu_1 - Ynu1);

    double Ynu0b = gsl_sf_bessel_Ynu(nu, xp);
    
    double Ynu_1c = gsl_sf_bessel_Ynu(nu - 1, x);
    double Ynu1c = gsl_sf_bessel_Ynu(nu + 1, x);
    double Ynu_primec = 0.5 * (Ynu_1c - Ynu1c);

    printf("Yν(%g, %g) = %g\n", nu, x, Ynu0/ Ynu_prime);
    printf("Yν'(%g, %g) = %g\n", nu, xp, Ynu_prime/Ynu0b);
    printf("Yνb(%g, %g) = %g\n", nu, xp, Ynu_primec/ Ynu_prime);
    double ratio1 = - (xp/nu) * pow(xp / x, nu); // Y/Y' dif x same n
    double ratio2 = - nu/xp; // Y/Y' same x same n
    double ratio3 = pow(xp/x, nu+1); // Y'/Y' dif x same n
    printf("Ratio n*(x/xp)^n = %g %g %g\n", ratio1, ratio2, ratio3);

    for (double D_radius = D_init; D_radius <= D_final; D_radius += D_step) {

        printf("\n===== D_radius = %lf =====\n", D_radius);

        compute_all_roots(all_roots, MAX_ROOTS, D_radius);

        printf("\nFound %d roots for k_perp\n\n", valid_nr);

        for (int i = 0; i < valid_nr; i++) { 
            double k_perp = roots[i]; 

            printf("Root %d: %lf", i+1, k_perp);
            printf(" frequency: %g", 3e8 * k_perp / (2 * M_PI));

            double matrix[NM][NM];

            build_Z_matrix(matrix, k_perp, D_radius);

            double det = compute_det_matrix(NM, matrix, D_radius);

            printf(" | Determinant: %g\n", det); 
        }

    }

    /*for (double i = 0.0; i < 10.0; i += 0.1) { 
        double k_perp = i; 
        printf("Testing k_perp = %g", i);

        double matrix[NM][NM];

        build_Z_matrix(matrix, k_perp);

        // Εκτύπωση του πίνακα για έλεγχο (προαιρετικό)
        printf("Matrix (Z_nq - delta_nq * Z_n):\n");
        for (int n = 0; n < N_MAX; n++) {
            for (int q = 0; q < N_MAX; q++) {
                printf("%10.6f ", matrix[n][q]);
            }
            printf("\n");
        }

        // Υπολογισμός και εκτύπωση του ορίζουσας
        double det = compute_det_matrix(NM, matrix);

        printf(" | Determinant: %g\n", det); 
    }*/

    // Παράδειγμα δημιουργίας πίνακα για την πρώτη ρίζα
    /*if (nr > 0) {
        double k_perp = roots[0];
        printf("\nMatrix for first root (k_perp = %g):\n", k_perp);

        double matrix[NM][NM];
        build_Z_matrix(matrix, k_perp, D_radius);

        for (int n = 0; n < NM; n++) {
            for (int q = 0; q < NM; q++) {
                printf("%10.6f ", matrix[n][q]);
            }
            printf("\n");
        }

        create_mat_reduced(k_perp, D_radius);
    }*/

    // Παλιός κώδικας με χρήση του head.h
    /*
        int max_roots = 10;
        int n_max = 1;

        double xmin = 0.0, xmax = 10.0; 
        double dx = 1.0;         
        double tol = 1e-10;      // Ακρίβεια εύρεσης ρίζας
        double tol_zero = 1e-12;  // Ακρίβεια συνάρτησης f(x) ~ 0
        double k_perp[max_roots]; 

        FILE *fp = fopen("output.txt", "w");
        if (!fp) {
            fprintf(stderr, "Error: Could not open file for writing.\n");
            return 1;
        }
        
        double N_roots = find_roots(xmin, xmax, dx, tol, tol_zero, k_perp, max_roots, n_max);

        for (int i = 0; i < N_roots; i++) {
            compute_matrix_determinant(n_max, k_perp[i], fp);
            printf("Found root k_perp = %.12f\n", k_perp[i]);
        }

        if (N_roots == 0) {
            fprintf(stderr, "No roots found! Exiting.\n");
            return 1;
        }
    */

    return 0;
}
