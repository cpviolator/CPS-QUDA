//mom1D.h
#ifndef __momenta_h_
#define __momenta_h_

//class Momenta;

class mom1D {
//friend class Momenta;
public:
	mom1D(int max, double sin_cutoff_);
//private:
	int P[3];
	int index();
	double mod();
	const int max;
	const double sin_cutoff;
};

//class Momenta {
//public:
//	Momenta(int MAX, double sinZcut, double sinYcut, double sinXcut);
//	mom1D x;
//	mom1D y;
//	mom1D z;
//
//}

#endif
