//**************************************************************************************
//						POST-NEWTONIAN POLYTROPIC STAR
// PNPolytrope.h
//		A simple stellar object in 1PN physics obeying everywhere a polytropic
//			equation of state of the form 
//				P~rho^Gamma
//			where rho = mass-energy density
// 		compare to Tooper 1964
//		derivations in 1pn4_polytrope.nb
//		Based on my own work with perturbed 1PN equations of stellar structure
// 		We use units where G=c=1
//***************************************************************************************


#ifndef PNPOLYTROPEH
#define PNPOLYTROPEH

#include "PNStar.h"

class PNPolytrope : public PNStar {
public:
	
	void graph_title(char* buff){
		sprintf(buff, "1PN polytrope n=%0.2f, sigma=%0.2le", n, sigma);
	}
	
	PNPolytrope(double, double, double, int);	//constructor, M, R, n and length
	PNPolytrope(double, double, int);	//constructor
	PNPolytrope(double, double, int, const double);	//constructor
	virtual ~PNPolytrope(); //destructor - clears all space in memory
	int length(){return len;}
	//these three functions specify units
	double Radius();
	double Mass();
	virtual double Gee(){return GG;};
	virtual double light_speed2(){return cc;};
	
	//these return value of indicated variable
	double getX(int X){return Y[X][x];}
	double getY(int X){return Y[X][y];}
	double getYderiv(int X){return Y[X][z];}
	//this is proper distance from center
	double getXproper(int X);
	double getSigma(){return sigma;};
	double getRedshift(){return zsurf;};
	
	double rad(int);
	//Note: here rho is the total mass-energy density, different from classical case
	double rho(int), drhodr(int);
	double   P(int),   dPdr(int);
	double Phi(int), dPhidr(int);
	double mr(int);
	//the 1PN functions
	double Psi(int), dPsidr(int);
	double Wx(int ), dWxdr(int );
	double Wy(int ), dWydr(int );
	double Wz(int ), dWzdr(int );
	
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
	double zsurf;	//the surface redshift -- not to be confused with z1
	//units and dimensionalty
	double GG;
	double cc;
	double rho0;	//central density
	double P0;		//central pressure
	double Rn;		//radius scale factor
	int len;
	//lane-emden solution functions
	enum VarNames {x=0, y,z,u,v,phi, numvar};
	double **Y;
//	double *x; //normalized radius
//	double *y; //lane-emden solution (theta θ)
//	double *z; //derivative, z=dy/dx
//	double *u; //unitless g, u = dφ/dx
//	double *v; //unitless gPN, v = dψ/dx
//	double *φ; //unitless newtonian potential, Phi = (n+1)P0/rho0 φ
	double *base; //= pow(y, n-1), used to speed of calculations of quantities
	double *dydx; //need to store derivative separately
	double φc0;// value of φ at center -- has to be derived from surface values
	
	inline double deriv1PN(double, double yy[numvar]);
	double setSigma(double ysurface[numvar]);
	std::complex<double> setZsurf(double ysurface[numvar]);
	void   centerInit(double ycenter[numvar]);
	void   RK4step(  double dx, double s, double yin[numvar], double yout[numvar]);
	
	//integrate Lane-Emden using RK4 method
	std::complex<double> RK4integrate(double, double& dx);
	//to find the correct step size
	double RK4integrate(const int, double);
	//when we know the correct step-size
	int RK4integrate(const int, double, int);
	
	//methods for handling the BCs
	double θc[3], φc[3], ψc[3]; //coefficients of theta near center
	double θs[4], φs[4], ψs[4]; //coefficients of theta near surface
	double igs[4];				//coefficients of 1/theta ner surface
	void setupCenter( );
	void setupSurface();	
public:
	//must return stellar structure expanded in powers of x=r/R near center
	void getAstarCenter(double *, int&, double g=0);
	void getUCenter(double*, int&);
	void getVgCenter(double*, int&, double g=0);
	void getC1Center(double*, int&);
	void getBetaCenter(double*, int&, double g=0);
	void getPhiCenter(double*, int&);
	//must return stellar structure expanded in powers of t=1-r/R near surface
	void getAstarSurface(double *, int&, double g=0);
	void getUSurface(double*, int&);
	void getVgSurface(double*, int&, double g=0);
	void getC1Surface(double*, int&);
	void getBetaSurface(double*, int&, double g=0);
	void getPhiSurface(double*, int&);
	void writeStar(char *c=NULL);
	double SSR();
private:
	void printCoefficients(char *c);
};

#endif