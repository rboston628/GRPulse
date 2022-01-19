// *************************************************************************************
//					GENERAL RELATIVISTIC POLYTROPIC STAR
// GRPolytrope.h
//		A simple stellar object in GR obeying everywhere a polytropic equation of state
//				P~rho^Gamma
//		Here rho = total mass-energy density, following Tooper, 1964
//		Solves Relativivistic Lane-Emden equation, in eqn 2.25, 2.26 in (Tooper 1964)
//		y=theta, x=xi, v=v (a unitless mass) in Tooper's notation
// *************************************************************************************

#ifndef GRPOLYTROPEH
#define GRPOLYTROPEH

#include "GRStar.h"

class GRPolytrope : public GRStar {
public:

	void graph_title(char* buff){
		sprintf(buff, "GR polytrope n=%0.2f, sigma=%0.2f", n, sigma);
	}
	
	GRPolytrope(double, double, int);	//constructor
	GRPolytrope(double, double, int, const double);	//constructor
	virtual ~GRPolytrope(); //destructor - clears all space in memory
	int length(){return len;}
	//these three functions specify units
	double Radius();
	double Mass();
	double Gee(){return 1.0;}
	double light_speed2(){return 1.0;}
	
	double getX(int X){return Y[X][x];}
	double getY(int X){return Y[X][y];}
	double getV(int X){return Y[X][v];}
	//this is proper distance from c enter
	double getXproper(int X);
	double getSigma(){return sigma;};
	double getRedshift(){return zsurf;};
	
	double rad(int);
	//Note: here rho is the total mass-energy density, different from classical case
	double rho(int), drhodr(int);
	double   P(int),   dPdr(int);
	double Phi(int), dPhidr(int);
	double mr(int);
	//GR-specific functions
	double gtt(int), dgttdr(int);
	double grr(int), dgrrdr(int);
	double gtr(int), dgtrdr(int);
	double gθθ(int), dgθθdr(int);
	double Nu(int);
	double dNudr(int);
	double d2Nudr2(int);
	double Lambda(int);
	double dLambdadr(int);
	double d2Lambdadr2(int);
	
	double Schwarzschild_A(int, double GamPert=0.0);
	double getAstar(int, double GamPert=0.0);
	double getVg(int, double GamPert=0.0);
	double getU(int);
	double getC(int);
	double Gamma1(int);
	double sound_speed2(int, double GamPert=0.0);
		
private:
	double n;		//polytropic index
	double Gamma;	//polytropic exponent
	double sigma;	//relativistic paramter as per Toopeer (1964)
	double zsurf;	//surface redshift
	double rho0;	//central density
	double P0;		//central pressure
	double Rn;		//radius scale factor
	int len;
	//lane-emden solution functions
	enum VarNames {x=0, y, v, nu, la, numvar};
	double **Y;
	double *base; //= pow(y, n-1), used to speed of calculations of quantities
	double *dydx; //need to store derivative separately
	
	//this deriv comes from eqn 2.25 in Tooper, eqn 2.8 in Bludman
	double deriv(double sigma, double yy[numvar]);
	double setSigma(double ysurface[numvar]);
	double setZsurf(double ysurface[numvar]);
	void   centerInit(double s, double z, double ycenter[numvar]);
	void   RK4step(double dx, double s, double yin[numvar], double yout[numvar]);
	
	
	//integrate Lane-Emden using RK4 method
	double RK4integrate(double s, double& dx);
	//to find the correct step size
	double RK4integrate(const int, double);
	//when we know the correct step-size
	int RK4integrate(const int, double, int);
	
	//methods for handling the BCs
	void setupCenter( );
	double thc[4+1], nuc[4+3], lc[4+1];
	void setupSurface();
	double ths[4+1], nus[4+3], ls[4+1];
public:
	//must return stellar structure expanded in powers of x=r/R near center
	void getAstarCenter(double *, int&, double g=0);
	void getUCenter( double*, int&);
	void getVgCenter(double*, int&, double g=0);
	void getC1Center(double*, int&);
	void getNuCenter(double*, int&);
	void getLambdaCenter(double*, int&);
	//must return stellar structure expanded in powers of t=1-r/R near surface
	void getAstarSurface(double *, int&, double g=0);
	void getUSurface( double*, int&);
	void getVgSurface(double*, int&, double g=0);
	void getC1Surface(double*, int&);
	void getNuSurface(double*, int&);
	void getLambdaSurface(double*, int&);
	void writeStar(char *c=NULL);
	double SSR();
	private:
//	void printCoefficients(char *c);
};

#endif