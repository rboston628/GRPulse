//***************************************************************************************
//				POST-NEWTONIAN NONRADIAL PULSATION nlm MODE DRIVER
// PNNonradialMode.h
// 		Solves the 1PN LAWE in the 1PN Djiembowski variables
//		Modes are 8th order, reduced by harmonic coordinate condition
//***************************************************************************************

#ifndef PNFULLMODEDRIVERH
#define PNFULLMODEDRIVERH

#include "ModeDriver.h"
#include "../STARS/PNStar.h"

class PNNonradialModeDriver : public ModeDriver {
public:
	static const unsigned int num_var=8U;
	//constructor
	//initialize from a background star and adiabatic index
	PNNonradialModeDriver(PNStar*, double);
	//destructor
	~PNNonradialModeDriver();
			
	int length();
	double Gamma1();
	double rad(int);
	void getCoeff(double *CC, const int, const int, const double, const int);
	
private:
	PNStar *const star;
	int len;
	int len_star;
	double adiabatic_index;
	enum VarNames {y1=0, y2, y3, y4, z1, z2, z3, z5};
	
	double Gee, cee2;
	double zsurf;
	
	//additional perturbation quantities
	double *r, *A, *U, *C, *V, *beta, *Phi, *xstar;
	void initializeArrays();
	
	//central values needed for BC order 2
	double cc[2], Vgc2, Ac2, Uc2, betac0, phic0;
	//surface values needed for BC order 2
	double surfCoeff[3][num_var][num_var];//order 0, 1, 2 in t
	double surfCoeffN[num_var][num_var];//for polytropes of small non-integer 
	double n_surface;
	double cs[3], betas0, phis0;
	static const int BC_C = 2;
	static const int BC_S = 2;
	void setupBoundaries();
	void updateSurface(double (&A)[3][num_var][num_var], double omeg2, int l);
	int CentralBC(double **y, double *y0, double w2, int l, int m=0);
	int SurfaceBC(double **y, double *y0, double w2, int l, int m=0);
	
public:
	void getBoundaryMatrix(int, double *, double*, double**, int*);
	void varnames(std::string *names){
		names[y1] = "y1"; names[y2]="y2"; names[y3]="y3"; names[y4]="y4";
		names[z1] = "z1"; names[z2]="z2"; names[z3]="z3"; names[z5]="z5";
	}
	
	double SSR(double, int, ModeBase*);
	double DPfunc(int X, double freq, ModeBase*);
	double tidal_overlap(ModeBase*);
	double innerproduct(ModeBase*, ModeBase*);
};

#endif