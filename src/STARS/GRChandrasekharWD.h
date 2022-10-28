//**************************************************************************************
//							CHANDRASEKHAR WHITE DWARF
// ChandrasekharWD++.cpp
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//		This model assumes T=0, and ignores Coulombic and other effects
//		The surface is not treated in any special way
//		Updated to include a chemical profile, indicated by mu_electron
// **************************************************************************************/

#ifndef GRChandrasekharWDH
#define GRChandrasekharWDH

#include "GRStar.h"
#include "../../lib/chandra.h"


class GRChandrasekharWD : public GRStar {
public:

	void graph_title(char* buff){
		sprintf(buff, "GR Chandrasekhar WD with y_0=%1.2f, sigma=%0.2f", Y0, sigma);
	}
	
	//the constructors
	void basic_setup();
	void init_arrays();
	GRChandrasekharWD(double, int,               double MU0=1);
	GRChandrasekharWD(double, int,               double mu,  double AN, double BN);
	GRChandrasekharWD(double, int, const double, double MU0=1);
	virtual ~GRChandrasekharWD();   //destructor
	int length(){return len;}
	//these three functions specify units
	double Radius();	//total radius
	double Mass();//total mass
	double Gee(); //{return G_CGS;};
	//in Newtonian, light speed is infinity...
	double light_speed2();//{return C_CGS*C_CGS;};
	
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
	double Y0;		// central value of y, y^2=1+x^2
	double X0, X02, Y02;
	double A0;		//pressure scale
	double B0;		//density scale
	double Rn;		//radius scale
	//
	double zsurf;	//surface redshift
	double sigma;	//relativistic parameter as per Tooper & Chandrasekhar (1969)
	
	int len;
	double dx;
	//lane-emden solution functions
	enum VarNames {xi=0, y, v, nu, la, numrk, x, f, g, numvar};
	double **Y;
//	double *xi;	//normalized radius
//	double *x;  //the relativity factor x = pF/mc
//	double *y;	//Chandrasekhar's y, y^2=1+x^2
//	double *z;	//derivative (dy/dxi)  note: dx/dxi = (dy/dxi)/x
//	double *mass;
//	double *f;  //factor relating x to pressure
//	double *g;  //factor relating x to energy
		
	//parameters of the chemical profile
	double mu0;
	//integrate using basic RK4
	double RK4integrate(double s, double& dx);
	double RK4integrate(const int, double);
	int RK4integrate(const int, double, int);
	
	//the T=0 Fermi function
	//double factor_f(double x);
	
	//methods for handling the BCs
	double yc[4], xc[3], fc[2];	//series coefficients of y,x,f near center
	void setupCenter();		//prepare values near center
	void setupSurface();	//prepare values near surface
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
	//
	void writeStar(char *c=NULL);
	double SSR();
};

#endif