//**************************************************************************************
//						POST-NEWTONIAN CHANDRASEKHAR WHITE DWARF
// PNChandrasekharWD.h
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//  	Builds off of work in Tooper (1964) to implement CHWD in GR
//		This model assumes T=0, and ignores Coulombic and other effects
//		The surface is not treated in any special way
//		Assumes mu_e = 1
//
//		Utilizes my own work with perturbed 1PN equations of stellar structure
// 		We use units where G=c=1
//***************************************************************************************


#ifndef PNChandrasekharWDH
#define PNChandrasekharWDH

#include "PNStar.h"

class PNChandrasekharWD : public PNStar {
public:
	
	void graph_title(char* buff){
		sprintf(buff, "Chandrasekhar WD with y_0=%1.2f", Y0);
	}
	
	PNChandrasekharWD(double, double, int);	//constructor
	PNChandrasekharWD(double, double, int, const double);	//constructor
	virtual ~PNChandrasekharWD(); //destructor - clears all space in memory
	int length(){return len;}
	//these three functions specify units
	double Radius();
	double Mass();
	virtual double Gee(){return G_CGS;};
	virtual double light_speed2(){return C_CGS*C_CGS;};
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
	double getVp(int);
	double getEta(int);
		
	void writeStar(char *c=NULL);
	
	double SSR();
	
private:
	double Y0;		// central value of y, y^2=1+x^2
	double X0, X02, Y02;
	double sigma;	//relativistic paramter as per Toopeer (1964)
	double zsurf;	//the surface redshift -- not to be confused with z1
	double A0;		//central pressure
	double B0;		//central density
	double Rn;		//radius scale factor
	int len;
	double dxi;
	

	//lane-emden solution functions
	double *xi;	//normalized radius
	double *x;  //the relativity factor x = pF/mc
	double *dx; //derivative dx/dxi
	double *φ;  //unitless newtonian potential, Phi = (n+1)P0/rho0 p
	double *u;  //unitless dPhi/dr, u = dφ/dxi
	double *v;  //unitless dPsi/dr, v = dψ/dxi
	double *y;	//Chandrasekhar's y
	double *f;
	double *h;
	double φc0;// value of φ at center -- has to be derived from surface values
	
	//integrate Lane-Emden using RK4 method
	void   RK4integrate(double z[2], double f[2], double&, int b=1);
	double RK4integrate(const int, double);
	//a grid-multiplying RK4 method
	int RK4integrate(const int, double, int);
	
	//the T=0 Fermi function
	double factor_f(double x);
	double factor_g(double x);
	double factor_h(double x);
	
	//methods for handling the BCs
	double xc[2], yc[2], φc[3], ψc[3], fc[2], hc[2]; //coefficients near center
	double xs[4], ys[4], φs[4], ψs[4];               //coefficients near surface
	void setupCenter( );
	void setupSurface();	
public:
	//must return stellar structure expanded in powers of x=r/R near center
	void getAstarCenter(double*, int&, double g=0);
	void getUCenter(    double*, int&);
	void getVgCenter(   double*, int&, double g=0);
	void getC1Center(   double*, int&);
	void getBetaCenter( double*, int&, double g=0);
	void getPhiCenter(  double*, int&);
	//must return stellar structure expanded in powers of t=1-r/R near surface
	void getAstarSurface(double*, int&, double g=0);
	void getUSurface(   double*, int&);
	void getVgSurface(  double*, int&, double g=0);
	void getC1Surface(  double*, int&);
	void getBetaSurface(double*, int&, double g=0);
	void getPhiSurface( double*, int&);
};

#endif