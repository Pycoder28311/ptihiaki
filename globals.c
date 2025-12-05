#include "globals.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int N = 10; // Αριθμός ραβδώσεων
double phi_i = 0.0; // Αρχική γωνία ράβδωσης
double phi_c = 0.1; // Η γωνία αυλάκωσης
//Οι γωνίες είναι σε ακτίνια
double D = 10.0; // Μεγάλη ακτίνα
double a = 9.0; // Μικρή ακτίνα

int Lmax = 3; // Το μέγιστο όριο σειράς S
