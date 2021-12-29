//A spline interpolation object.  Pass support points, creates an interpolation
//Afterwards, can be passed x-values and will return interpolate y-values

#ifndef INTERP
#define INTERP

#include <math.h>

class Interpolator {
	public:
	virtual double interp(double) = 0;
	virtual double operator()(double) = 0;
	virtual double deriv(double) = 0;
	virtual ~Interpolator(){};
};


class Splinor : public Interpolator {
	public:
	Splinor(const double *const x, const double *const y, const int L);
	~Splinor();
	double interp(double);
	double operator()(double);
	double deriv(double);
	
	private:
	int len;
	int xlast, incr;
	double xa, xb; //boundaries
	void getSplineCoefficients(const double *const, const double *const);
	double *S;
	double *xarr;
	double *yarr;
};

#endif