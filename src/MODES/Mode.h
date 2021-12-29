//**************************************************************************************
//							The MODE Object
//	Mode.h
//		Equation agnostic 
//			-- all information specific to physics supplied by the driver
// 		Capable of nonradial, Cowling, 1PN modes by supplying different drivers
//**************************************************************************************

#ifndef MODEH
#define MODEH

#include "../constants.h"
#include "../STARS/Star.h"
#include "../MODES/ModeDriver.h"

template <size_t numvar> class Mode : public ModeBase {
public:
	static const unsigned int  num_var = numvar;
	//retrieve quantum numbers of pulsation
	int  modeOrder(){return k;};
	void modeNumbers(int& K, int& L, int& M){K=k; L=l; M=m;};
	//constructors: differ only in initial choice of frequency
	//initialize from quantum numbers and a background star
	Mode<numvar>(int k, int l, int m, ModeDriver*);
	//initialize from frequency guess and a background star
	Mode<numvar>(double omega2, int l, int m, ModeDriver*);
	//initialize from brackets on frequency and a background star
	Mode<numvar>(double omegMin, double omegMax, int l, int m, ModeDriver*);
	//deconstructor
	virtual ~Mode();
	
	//ways to access the frequency
	double getOmega2();
	double getFreq();
	double getPeriod();
	double SSR();
	double tidal_overlap();

	//methods to access mode functions
	double getRad(int x){ return rad[x];};
	double getY(int i, int x){ return pow(rad[x],l-2)*y[i][x];};
	double getYtilde(int i, int x){ return y[i][x];};
	
protected:
	//the equilibrium star
	Star *const star;
	//the mode driver specifying equations and boundaries
	ModeDriver *const driver;
	void basic_setup();	//handles the normal setup common to all constructors
	int k,l,m;			//the mode numbers
	double cee2, Gee;	//for specifying units
	double Gamma1;		//adiabatic exponent for perturbations	
	int len, len_star;	//length of mode arrays, and length of background model arrays
	double omega2;		//frequency squared
	double *rad;		//the radial coordinate, x = r/R
	//perturbation quantities
	double **y;			//the eigenfunctions y_1, y_2, ..., y_N
	double yCenter[numvar], ySurface[numvar];//central and surface values of y_1,...y_N
	
	bool converged;		//a flag for indicating if calculation converged -- no longer needed
	void convergeBisect(double);		//find frequency using a bisection search
	void convergeNewton(double, int);	//find frequency using a Newton method
	void linearMatch(double w2, double y0[numvar], double ys[numvar]);
	int verifyMode();
	void converge();
	//
	void RK4out(int, double, double y0[numvar]);
	void RK4in( int, double, double ys[numvar]);
	double RK4center(double, double y0[numvar], double ys[numvar]);
	
	double boundaryMatrix[numvar][numvar];	//the BCs specified in matrix form
	int    indexOrder[numvar];				//the order in which to multipy factors to match
	
	int xfit;			//the coordinate where inner/outer solutions are fit
	double sig2omeg;	//converstion from dimensionful frequency to dimensionless frequency
	double R;			//the radius of the star
public:
	//file output methods to write and view plots of mode
	void printMode(char *c = NULL);
	void writeMode(char *c = NULL);

};


//template classes in C++ cannot be split into multiple files, so we must include Mode.pp
#include "Mode.cpp"


#endif