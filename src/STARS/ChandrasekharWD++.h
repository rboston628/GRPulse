//**************************************************************************************
//							CHANDRASEKHAR WHITE DWARF
// ChandrasekharWD++.h
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//		This model assumes T=0, and ignores Coulombic and other effects
//		The surface is not treated in any special way
//		Updated to include a chemical profile, indicated by mu_electron
// **************************************************************************************/

#ifndef ChandrasekharWDH
#define ChandrasekharWDH

#include "Star.h"
#include "../../lib/chandra.h"

class ChandrasekharWD : public Star {
public:

	void graph_title(char* buff){
		sprintf(buff, "Chandrasekhar WD with y_0=%1.2f", Y0);
	}
	
	//the constructors
	void basic_setup();
	void init_arrays();
	ChandrasekharWD(double, int,               double MU0, double K, double AC, double AS);
	ChandrasekharWD(double, int,               double mu,  double AN, double BN);
	ChandrasekharWD(double, int, const double, double MU0, double K, double AC, double AS);
	virtual ~ChandrasekharWD();   //destructor
	int length(){return len;}
	//these three functions specify units
	double Radius();	//total radius
	double Mass();//total mass
	double Gee();//
	//in Newtonian, light speed is infinity...
	double light_speed2();
	
	//these return value of indicated variable -- used in testing
	double getXi(int X){return Y[X][xi];}
	double getX(int X){return Y[X][x];}
	double getY(int X){return Y[X][y];}
	double getYderiv(int X){return Y[X][z];}
	
	double rad(int);
	double rho(int), drhodr(int);
	double   P(int),   dPdr(int);
	double Phi(int), dPhidr(int);
	double mr(int);
	
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
	int len;
	double dxi;
	
	//solution functions
	enum VarName {
		xi=0,	//normalized radius
		x, 		//the relativity factor x = pF/mc
		y,		//Chandrasekhar's y
		z,		//derivative dy/dxi
		f,		//factor f, for P = A f(x)
		numvar
	};
	double **Y;
	//for the mass
	double *mass;
		
	//parameters of the chemical profile
	double mu0, k, acore, aswap;
	void chemical_gradient(const double, const double, double&, double&);
	
	
	//to handle chemical profile
	double* mue;   //mean atomic mass per electron
	double* dmue;  //derivative of above
	//integrate using basic RK4
	void centerInit(double ycenter[numvar]);
	void RK4step(double dx, double yin[numvar], double yout[numvar]);
	double RK4integrate(const int, double);
	int RK4integrate(const int, double, int);
	
	//methods for handling the BCs
	double yc[4], xc[3], fc[2];	//series coefficients of y,x,f near center
	void setupCenter();		//prepare values near center
	void setupSurface();	//prepare values near surface
	
public:
	//methods to find central, surfae power series expansions of key variables in pulsation
	void getAstarCenter(double *, int&, double g=0);
	void getUCenter(double*, int&);
	void getVgCenter(double*, int&, double g=0);
	void getC1Center(double*, int&);
	void getAstarSurface(double *, int&, double g=0);
	void getUSurface(double*, int&);
	void getVgSurface(double*, int&, double g=0);
	void getC1Surface(double*, int&);

	//a particular output generation for this model of white dwarf
	void writeStar(char *c=NULL);
private:
	void printChemicalGradient(char *c);
};

#endif