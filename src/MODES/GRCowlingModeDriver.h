//**************************************************************************************
//						GR COWLING PULSATION nlm MODE DRIVER
// GRCowlingModeDriver.h
//		Based on the equation of Yoshida & Lee (2002), using Dziembowski-like variables,
//		See their equations (42), (43), and definitions in (37) to (41), (44)
//		The metric component has been corrected, their 2nu -> nu, their 2lambda -> lambda
//			total metric: ds^2 = -e^nu dt^2 + e^lambda dr^2 + r^2 dOmega^2
//		In their (43), the "U" should be U_1
//		Only works for L>=2
//  	This method is equivalent to that of McDemrott, Van Horn and Scholl
//***************************************************************************************

#ifndef GRCOWLINGDRIVERH
#define GRCOWLINGDRIVERH



#include "ModeDriver.h"
#include "../STARS/GRStar.h"

class GRCowlingModeDriver : public ModeDriver {
public:
	//the order of equation set is second
	static const unsigned int num_var=2U;
	
	//constructor
	//initialize from a background star and an adiabatic index
	GRCowlingModeDriver(GRStar*, double);
	//destructor
	virtual ~GRCowlingModeDriver();

	int length();
	double Gamma1();//the adiabatic index, if set; 0 if same as star
	double rad(int);//radial array x = r/R
	//return the coefficient matrix A from equation x*dy/dx = A*y
	void getCoeff(double *CC, const int, const int, const double, const int);

private:
	GRStar *const star;
	int len;		//number of grid points for mode
	int len_star;	//number of grid points in star
	double adiabatic_index;	//adiabatic index; set to 0 to use star's Gamma1
	enum VarNames {y1=0, y2};
	
	double Gee, cee2;
	double zsurf;
	
	//perturbation quantities that appear in Dziembowski equations
	double *r, *A, *U1, *C, *V, *U2;
	void initializeArrays();
	
	static const int BC_C = 4;	//the desired order near the center
	static const int BC_S = 4;	//the desired order near the surface
	//expansion coefficients near surface;
	double nus[BC_S+2], ls[BC_S+2];
	double As[BC_S+1], Vgs[BC_S+1], cs[BC_S+1], U1s[BC_S+1], U2s[BC_S+1], k_surface;
	//expansion coefficients near center;
	double nuc[BC_C/2+3], lc[BC_C/2+1];
	double Ac[BC_C/2+1], Vgc[BC_C/2+1], cc[BC_C/2+1], U1c[BC_C/2+1], U2c[BC_C/2+1];
	double Asn, Vgn, U2sn, U1sn, c1sn;
	void setupBoundaries();
	int CentralBC(double **y, double *yo, double s2, int l, int m=0);
	int SurfaceBC(double **y, double *yo, double s2, int l, int m=0);

public:	
	void getBoundaryMatrix(int, double *, double*, double**, int*);
	void varnames(std::string *names){
		names[y1] = "y1"; names[y2]="y2";
	}
	
	//for the Mode passed, calculate sum-square residual
	double SSR(double, int l, ModeBase*);
	double tidal_overlap(ModeBase*);
	double innerproduct(ModeBase*, ModeBase*);
};

#endif