//momenta.cpp
#include <cmath>
#include "mom1D.h"

mom1D::mom1D(int max_, double sin_cutoff_):
max(max_), sin_cutoff(sin_cutoff_) {
	for (int i=0; i<3; i++) P[i] = 0;
}

int mom1D::index() {
	return P[2]*max*max + P[1]*max + P[0];
}

double mom1D::mod() {
	int a, b, c;

	if (P[2] >= max/2) a = max-P[2];
	else a = P[2];

        if (P[1] >= max/2) b = max-P[1];
        else b = P[1];

        if (P[0] >= max/2) c = max-P[0];
        else c = P[0]; 

	return M_PI*sqrt(a*a + b*b + c*c)/(1.0*max);
}

//Momenta:: Momenta(int MAX, double sinZcut, double sinYcut, double sinXcut) {
//	
//}
