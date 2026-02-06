#include "globals.h"

const double a_radius = 1.0;       // Μικρή ακτίνα 1.0e-2
const double D_radius = 1.5;       // Μεγάλη ακτίνα 1.00001e-2
const double N = 10.0;             // Αριθμός ραβδώσεων
const double phi_c = 0.314;          // Γωνία αυλάκωσης
const double phi_i = 0.0;          // Γωνία εισόδου
const int lmax = 5;               // Μέγιστο l για το S 2
const double m = 1.0;              // Σταθερά m ρυθμού
 
const int n_max = 2;               // Διαστάσεις Πίνακα 1
const double down = 0.0;           // Κάτω όριο για riza 0.0
const double up = 20.0;            // Άνω όριο για riza 3000.0
const double step = 0.01;           // Βήμα για riza 20.0
const double tol_riza = 0.00001;      // Ακρίβεια για riza