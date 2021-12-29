//Class to calculate read in MESA data for pulsation calculations d

#ifndef MESACLASSH
#define MESACLASSH

#include "Star.h"

class MESA : public Star {
public:	
	MESA(const char*, int);	//constructor
	virtual ~MESA(); //destructor - clears all space in memory
	int length(){return len;}
	double Mass();
	double Radius();
	double Gee();
	double light_speed2();
	//return a scale factor to put frequencies in seconds
	//double getF0(){ return F0;}
	void graph_title(char* buff){
		sprintf(buff, "MESA model w/ Mass=%lg", Mtot);
	}
	
	double rad(int);
	double rho(int), drhodr(int);
	double   P(int),   dPdr(int);
	double Phi(int), dPhidr(int);
	double mr(int);
	
	double Gamma1(int);
	double Schwarzschild_A(int, double GamPert=0.0);
	double getAstar(int, double GamPert=0.0);
	double getVg(int, double GamPert=0.0);
	double getU(int);
	double getC(int);
	double sound_speed2(int, double GamPert=0);
	
private:
	void printSection(int, int);
	void subgridCubicSpline(const int, const int, const int*);
	void spline(
		const int, const int, const int *const,
		const double *const,  const double *const,
		double*, double*
	);
	void getSplineCoefficients(const int, double*, const double *const, const double *const);
	
	//void printBruntVaisala(const char*);
		
	int Ntot, len, subgrid;  //Ntot = grid size from MESA
	double G, c2;
	//in units of g, cm, erg/s, respectively
	double Mtot, Rtot, Ltot;
	double Mstar, Rstar, Lstar;
	double Dscale, Pscale, Gscale;
	//arrays for density, pressure, mass, and gravitational field
	double *radi;	
	Splinor *dens, *pres, *mass, *grav, *Gam1, *BVfq;
	Splinor *aSpline, *vSpline, *uSpline, *cSpline;
	
	
	//for BCs of modes
	double nc;// effective polytropic index at center
	double ac[8]; //coefficients of theta near center
	double dc[3],pc[3];
	double A0[3],V0[3],c0[3],U0[3];
	//surface
	double ds[5], ps[5], ts[5], dels;
	double A1[5], V1[5], c1[5], U1[5];
	
	//double as[6]; //coefficients of theta near surface
	//double Vs[5];
	//double cs[4];
	//double Us[4];
	void setupCenter();
	void setupSurface();
	
public:
	void getAstarCenter(double *, int&, double g=0);
	void getUCenter(double*, int&);
	void getVgCenter(double*, int&, double g=0);
	void getC1Center(double*, int&);
	void getAstarSurface(double *, int&, double g=0);
	void getUSurface(double*, int&);
	void getVgSurface(double*, int&, double g=0);
	void getC1Surface(double*, int&);
	void writeStar(char *c=NULL);
	double SSR();
private:
	void printBV(char *c);
	void printCoefficients(char *c);
};

#endif