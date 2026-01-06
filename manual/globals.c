#include "globals.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int N = 5; // Αριθμός ραβδώσεων
double phi_i = 0.0; // Αρχική γωνία ράβδωσης
double phi_c = 0.1; // Η γωνία αυλάκωσης
//Οι γωνίες είναι σε ακτίνια

double D = 10.0; // Μεγάλη ακτίνα
double a = 9.0; // Μικρή ακτίνα

int Lmax = 5; // Το μέγιστο όριο σειράς S

double v = 3.0e8; // παράδειγμα: ταχύτητα κύματος (m/s)
