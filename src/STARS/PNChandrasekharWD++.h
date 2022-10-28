//**************************************************************************************
//						POST-NEWTONIAN CHANDRASEKHAR WHITE DWARF
// PNChandrasekharWD.h
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//  	Builds off of work in Tooper (1964) to implement CHWD in GR
//		This model assumes T=0, and ignores Coulombic and other effects
//		The surface is not treated in any special way
//		Updated to include a chemical profile, indicated by mu_electron
//
//		Based on my own work with perturbed 1PN equations of stellar structure
// 		We use units where G=c=1
//***************************************************************************************


#ifndef PNChandrasekharWDH
#define PNChandrasekharWDH

#include "PNStar.h"
#include "../../lib/chandra.h"

class PNChandrasekharWD : public PNStar {
public:
	
	void graph_title(char* buff){
		sprintf(buff, "1PN Chandrasekhar WD with y_0=%1.2f", Y0);
	}
	
	//the constructors
	void basic_setup();
	void init_arrays();
	PNChandrasekharWD(double, int,               double MU0, double K, double AC, double AS); //constructor
	PNChandrasekharWD(double, int, const double, double MU0, double K, double AC, double AS); //constructor
	virtual ~PNChandrasekharWD(); //destructor - clears all space in memory
	int length(){return len;}
	//these three functions specify units
	double Radius();
	double Mass();
	virtual double Gee(){return GG;};
	virtual double light_speed2(){return cee2;};
	double getRedshift(){return zsurf;};
	
	//these return value of indicated variable -- used in testing
	double getXi(int X){return Y[X][xi];}
	double getX(int X){return  Y[X][x];}
	double getY(int X){return  Y[X][y];}
		
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
	void printDegeneracy(char *c);
	void printChemicalGradient(char *c);
	
	double SSR();
	
private:
	int len;
	double dxi;
	double Y0;		// central value of y, y^2=1+x^2
	double X0, X02, Y02;
	double sigma;	//relativistic paramter as per Tooper (1964), set by physical constants
	double zsurf;	//the surface redshift -- cannot be freely set once Y0 is set
	// units and dimensionality
	double GG;
	double cee2;
	double A0;		//central pressure
	double B0;		//central density
	double Rn;		//radius scale factor

	

	//solution functions
	enum VarName {
		xi=0,	//normalized radius
		x, 		//the relativity factor x = pF/mc
		u,		//unitless dPhi/dr, u = dφ/dxi
		v,		//unitless dPsi/dr, v = dψ/dxi
		phi,	//unitless newtonian potential, Phi = 8A/B* phi
		y,		//Chandrasekhar's y
		f,		//Factor f, for P = A f(x)
		h,		//Factor h = (mu*x^3 + sigma*g(x), rho = B h(x)
		numvar
	};
	double **Y;
	double *mass;
	double *dxdxi;	//derivative dx/dxi
	double *dfdx;   //derivative df/dx = 8x^4/y
	double *dhdx;   //derivative dh/dx = 3x^2(mu + 8*sigma*(y-1))
	double phic0;// value of φ at center -- has to be derived from surface values
		
	//chemical gradient
	double mu0, k, acore, aswap;
	double* mue;   //mean atomic mass per electron
	double* dmue;  //derivative of above
	void chemical_gradient(const double, const double, double&, double&);
	
	double setZsurf(double ysurface[numvar]);
	
	enum X0Y0setState {setX0, setY0};
	void setX0Y0(double xy, X0Y0setState);
	void centerInit(double[numvar]);
	void RK4step(double, double[numvar], double[numvar]);
	//integrate Lane-Emden using RK4 method
	double RK4integrate(double, double&, bool f=true);
	double RK4integrate(const int, double);
	//a grid-multiplying RK4 method
	int RK4integrate(const int, double, int);
	
	//methods for handling the BCs
	//double xc[2], yc[2], φc[3], ψc[3], fc[2], hc[2]; //coefficients near center
	double φc[3], ψc[3], hc[2];
	double xs[4], ys[4], φs[4], ψs[4];               //coefficients near surface
	double yc[4], xc[3], fc[2];	//series coefficients of y,x,f near center
	
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
