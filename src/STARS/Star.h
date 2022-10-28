//**************************************************************************************
//							STAR ABSTRACT CLASS
// Star.h
//		A virtual class star for later models to inherit
//		Serves as the basis in code for other stellar models
//		Can be used in code as a generic star, without regard to particular physics
//**************************************************************************************

#ifndef STARH
#define STARH

#include "../constants.h"

class Star {
public:	
	//this is an identifier for the model of star
	//must be chosen to make UNIX accessible filename -- no spaces!
	char name[40];	
	virtual void graph_title(char*) =0;
	//this is the index where shooting methods for modes should join
	int indexFit;

	// boundary conditions
	//must return stellar structure expanded in powers of x=r/R near center
	virtual void getAstarCenter(double *, int&, double g=0) =0;
	virtual void getUCenter( double*, int&) =0;
	virtual void getVgCenter(double*, int&, double g=0) =0;
	virtual void getC1Center(double*, int&) =0;
	//must return stellar structure expanded in powers of t=1-r/R near surface
	virtual void getAstarSurface(double *, int&, double g=0) =0;
	virtual void getUSurface( double*, int&) =0;
	virtual void getVgSurface(double*, int&, double g=0) =0;
	virtual void getC1Surface(double*, int&) =0;

	//distance from center of star
	virtual double rad(int) =0;
	//density, pressure, potential, and their derivatives
	virtual double rho(int) =0, drhodr(int) =0;
	virtual double   P(int) =0,   dPdr(int) =0;
	virtual double Phi(int) =0, dPhidr(int) =0;
	//interior mass to r
	virtual double  mr(int) =0;
	//Schwarzschild discriminant
	virtual double Schwarzschild_A(int, double g=0.0) =0;
	virtual double getAstar(int, double g=0.0) = 0;
	virtual double getU( int) = 0;
	virtual double getVg(int, double g=0.0) = 0;
	virtual double getC( int) = 0;
	virtual double Gamma1(int) =0;
	virtual double sound_speed2(int, double g=0.0) =0;
	
	//destructor
	virtual ~Star(){};
	//return length of star
	virtual int length() =0;
	//dimensionfull quantities of interest
	virtual double Radius() =0;
	virtual double Mass() =0;
	virtual double Gee() =0;
	virtual double light_speed2() =0;
	//print relevant values of the star in .txt and gnuplot
	virtual void writeStar(char *c = NULL);
	virtual void printStar(char *c = NULL);
	virtual void printBV(  char *c = NULL, double const gam1=0.0);
	virtual void printCoefficients(char *c = NULL, double const gam1=0.0);
	virtual double SSR();
	
	//allow Modes to access private members of Star
	friend class ModeDriver;
};

#endif