#include "globals.h"

const double a_radius = 1.0;       // Μικρή ακτίνα
const double D_radius = 1.5;       // Μεγάλη ακτίνα
const double N = 10.0;             // Αριθμός ραβδώσεων
const double phi_c = 0.314;          // Γωνία αυλάκωσης
const double phi_i = 0.0;          // Γωνία εισόδου
const int lmax = 5;               // Μέγιστο l για το S
const double m = 1.0;              // Σταθερά m ρυθμού
 
const int n_max = 2;               // Διαστάσεις Πίνακα
const double down = 0.0;           // Κάτω όριο για riza
const double up = 20.0;            // Άνω όριο για riza
const double step = 0.01;           // Βήμα για riza
const double tol_riza = 0.00001;      // Ακρίβεια για riza