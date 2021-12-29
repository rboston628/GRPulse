// *************************************************************************************
//					GENERAL RELATIVISTIC STAR ABSTRACT CLASS
// GRStar.h
//		This class is another abstract class
//		It extends the definition of a Star to Einstein's relativity
// *************************************************************************************

#ifndef GRSTARH
#define GRSTARH

#include "Star.h"

class GRStar : public Star {

public:
	//the elements of the metric and their derivatives
	//assumes a spherical metric of the form
	//	 ds^2 = -gtt dt^2 + 2 gtr dtdr + grr dr^2 + gθθ dΩ^2
	virtual double grr(int )=0, dgrrdr(int )=0;
	virtual double gtt(int )=0, dgttdr(int )=0;
	virtual double gtr(int )=0, dgtrdr(int )=0;
	virtual double gθθ(int )=0, dgθθdr(int )=0;
	//the special functions from Schwarzschild metric
	virtual double Nu(int) =0, dNudr(int) =0, d2Nudr2(int) =0;
	virtual double Lambda(int) =0, dLambdadr(int) =0, d2Lambdadr2(int) =0;;
	//define these in case of typing mixups
	double grt(int X){ return gtr(X);};
	double dgrtdr(int X) {return dgtrdr(X);};
	
	
	virtual void     getNuCenter(double*, int&) =0;
	virtual void getLambdaCenter(double*, int&) =0;
	virtual void     getNuSurface(double*, int&) =0;
	virtual void getLambdaSurface(double*, int&) =0;
	
	//the redshift
	virtual double getRedshift() =0;
};

#endif