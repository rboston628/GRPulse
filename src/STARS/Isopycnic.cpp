//**************************************************************************************
//							ISOPYCNIC STAR
// Isopycnic.h
//		A simple stellar object in Newtonian physics of uniform density
//				rho = constant
//		Thi exploits the known analytic solution to Lane-Emden for n=0
//				y = 1 - x^2/6
//		Used as stellsr background for mode testing, to minimize background errors
//		THIS IS NOT PART OF GRPULSE -- ONLY FOR CERTAIN TESTS WITH NONRADIAL MODES
//**************************************************************************************


#ifndef ISOPYCNICCLASS
#define ISOPYCNICCLASS
#include "Isopycnic.h"

//initalize polytrope from length
Isopycnic::Isopycnic(int L)
	: len(L)
{
	GG=1.0;
	Gamma = 5./3.;

	sprintf(name, "isopycnic");
	
	//we find an appropriate grid spacing for array holding star data
	//an initial guess
	double dx = sqrt(6.0)/double(len-1), yS=1.0, ddx = dx;
	x = new double[len];
	y = new double[len];
	z = new double[len];
	double dxmax=1.0, dxmin=0, ySmax=-1.0, ySmin=1.0;
	
	populateValues(len, dx);
	
	//set initial density, pressure
	// specifying K, rho0 not necessary, just rescales the solutions
	rho0 = 1.0;
	P0 = 1.0;
	Rn = sqrt( 1.0/(4.*m_pi) );
	
	//now set physical properties of the polytrope
	mass = new double[len];
	mass[0] = 0.0;
	indexFit = len/2;
	double XS = 0.5*x[len-1];
	for(int X=1; X<len; X++){
		//as we scan through x,y,z, set matching point where y[X] = 0.5
		if(y[X-1]>0.5 & y[X+1]<=0.5) indexFit = X;
		mass[X] = -4.*m_pi*x[X]*x[X]*z[X];
	}
	mass[len-1] = 8.0*m_pi*sqrt(6.0);
	indexFit /= 2;
}

Isopycnic::~Isopycnic(){
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] mass;
}

int Isopycnic::populateValues(const int Len, double dx){
	//set our initial conditions
	x[0] = 0.0; y[0] = 1.0; z[0] = 0.0;
	for(int X = 0; X<Len; X++){
		x[X] = double(X)*dx;
		y[X] = 1.0 - x[X]*x[X]/6.0;
		z[X] = -x[X]/3.0;
	}
	x[len-1] =  sqrt(6.0);
	y[len-1] = 0.0;
	z[len-1] = -sqrt(6.0)/3.0;
	return 0;
}

//Here we define functions to access radius, pressure, etc.
double Isopycnic::rad(int X){
	return (Rn*x[X]);
}
double Isopycnic::rho(int X){
	return rho0;
}
double Isopycnic::drhodr(int X){
	return 0.0;
}
inline double Isopycnic::P(int X){
	return P0*y[X];
}
double Isopycnic::dPdr(int X){
	return P0*(z[X]/Rn);
}
double Isopycnic::Phi(int X){
	return -P0/rho0 *(y[X]-1.0);
}
 double Isopycnic::dPhidr(int X){
	return -P0/rho0*z[X]/Rn;
}
double Isopycnic::mr(int X){
	return pow(Rn,3)*rho0*mass[X];
}

double Isopycnic::Schwarzschild_A(int X, double GamPert){
	if(GamPert == 0.0) return -0.6*z[X]/(y[X]*Rn);
	else               return -z[X]/(y[X]*Rn)/GamPert;
}

double Isopycnic::getAstar(int X, double GamPert){
	if(GamPert == 0.0) return 0.6*x[X]*z[X]/(y[X]);
	else               return x[X]*z[X]/(y[X])/GamPert;
}

double Isopycnic::getU(int X){
	return 3.0;
}
double Isopycnic::getVg(int X, double GamPert){
	if(GamPert == 0.0) return -0.6*x[X]*z[X]/(y[X]);
	else               return -x[X]*z[X]/(y[X])/GamPert;
}
double Isopycnic::getC(int X){
	return 1.0;
}
double Isopycnic::Gamma1(int X){
	return 5./3.;
}

double Isopycnic::sound_speed2(int X, double GamPert){
	if(GamPert == 0.0) return Gamma  *P0*y[X]/rho0;
	else               return GamPert*P0*y[X]/rho0;
}


double Isopycnic::Radius(){return rad(len-1);}	//total radius
double Isopycnic::Mass(){  return mr(len-1);}//total mass


// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
//  Note that because of extremely simple structure, coefficiens available to arbitrary precision
void Isopycnic::setupCenter(){}

void Isopycnic::getAstarCenter(double *Ac, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	for(int k=0; k<=maxPow/2; k++){
		Ac[k] = -2./Gam1;
	}
}

void Isopycnic::getVgCenter(double *Vc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	for(int k=0; k<=maxPow/2; k++){
		Vc[k] = 2./Gam1;
	}
}

void Isopycnic::getUCenter(double *Uc, int& maxPow){
	if(maxPow>=0) Uc[0] = 3.0;
	for(int k=1; k<= maxPow/2; k++){
		Uc[k] = 0.0;
	}
}

void Isopycnic::getC1Center(double *cc, int& maxPow){
	if(maxPow>=0) cc[0] = 1.0;
	for(int k=1; k<= maxPow/2; k++){
		cc[k] = 0.0;
	}
}


// **************************  SURFACE BOUNDARY  ***************************************
// the following provide coefficients for surface expansions of A*, Vg, U, c1 in terms of t=1-r/R
//	Note that A*, Vg require a power -1
//	Note that up to order N, we only require:
//		A*, Vg, c1 up to order N-1
//		U          up to order N
//  Note that because of extremely simple structure, coefficiens available to arbitrary precision
void Isopycnic::setupSurface(){}

void Isopycnic::getAstarSurface(double *As, int& maxPow, double g){
	//we make use of the fact  that A* and Vg are simply related in polytropes
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	int O=1;
	if(maxPow>= -1) As[O-1] = -1./Gam1;
	if(maxPow>=  0) As[O  ] = 1.5/Gam1;
	for(int k=1; k<= maxPow; k++){
		As[O+k] = -pow(2,k+1)/Gam1;
	}
}

void Isopycnic::getVgSurface(double *Vs, int& maxPow, double g){
// coefficients of Vg must extend up to maxPow-1
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	int O=1;
	if(maxPow>= -1) Vs[O-1] = 1./Gam1;
	if(maxPow>=  0) Vs[O  ] = -1.5/Gam1;
	for(int k=1; k<= maxPow; k++){
		Vs[O+k] = pow(2,k+1)/Gam1;
	}
}

void Isopycnic::getUSurface(double *Us, int& maxPow){
// coefficients of U must extend up to order maxPow
	if(maxPow>=0) Us[0]  = 3.0;
	for(int k=1; k<= maxPow; k++){
		Us[k] = 0.0;
	}
}

void Isopycnic::getC1Surface(double *cs, int& maxPow){
// coefficients of c1 are only needed up to order maxPow-1
	if(maxPow>=0) cs[0]  = 1.0;
	for(int k=1; k<= maxPow; k++){
		cs[k] = 0.0;
	}
}

#endif