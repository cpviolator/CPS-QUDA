//mom3D.h
#ifndef __momenta_h_
#define __momenta_h_

//class Momenta;

class mom3D {
//friend class Momenta;
public:
	mom3D(int max, double sin_cutoff_);
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
//	mom3D x;
//	mom3D y;
//	mom3D z;
//
//}

#endif
