//**************************************************************************************
//						COWLING PULSATION nlm MODE DRIVER
// Cowling ModeDriver.h
//		Solves the Newtonian LAWE in Dziembowski variables in the Cowling approximation
//		Does not consider perturbtations to gravitational field, making 2nd-order modes
//***************************************************************************************

#ifndef COWLINGDRIVERH
#define COWLINGDRIVERH

#include "ModeDriver.h"

class CowlingModeDriver : public ModeDriver {
public:
	static const unsigned int num_var=2U;
	
	//constructor
	//initialize from a background star and an adiabatic index
	CowlingModeDriver(Star*, double);
	//destructor
	virtual ~CowlingModeDriver();
	
	int length();
	double Gamma1();//the adiabatic index, if set; 0 if same as star
	double rad(int);//radial array x = r/R
	//return the coefficient matrix A from equation x*dy/dx = A*y
	void getCoeff(double *CC, const int, const int, const double, const int);

private:
	int len;		//number of grid points for mode
	int len_star;	//number of grid points in star
	double adiabatic_index;	//adiabatic index; set to 0 to use star's Gamma1
	
	//perturbation quantities
	double *r, *A, *U, *C, *V;
	void initializeArrays();
	
	static const int BC_C = 4;	//the desired order near the center
	static const int BC_S = 4;	//the desired order near the surface
	//expansion coefficients near surface;
	double As[BC_S+1], Vs[BC_S+1], cs[BC_S+1], Us[BC_S+1], cProds[BC_S+1], k_surface;
	//expansion coefficients near center;
	double Ac[BC_C/2+1], Vc[BC_C/2+1], cc[BC_C/2+1], Uc[BC_C/2+1], cProdc[BC_C/2+1];
	void setupBoundaries();
	int CentralBC(double **y, double *yo, double s2, int l, int m=0);
	int SurfaceBC(double **y, double *yo, double s2, int l, int m=0);

public:	
	void getBoundaryMatrix(int, double *, double*, double**, int*);
	void varnames(std::string *names){
		names[0] = "y1"; names[1]="y2";
	}
	
	//for the Mode passed, calculate sum-square residual
	double SSR(double, int l, ModeBase*);
	double tidal_overlap(ModeBase*);
	double innerproduct(ModeBase*, ModeBase*);
};

#endif