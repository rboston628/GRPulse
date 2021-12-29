//**************************************************************************************
//							CHANDRASEKHAR WHITE DWARF
// ChandrasekharWD.h
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//		This model assumes T=0, and ignores Coulombic and other effects
//		Future step will be to adjust the composition, include finite T effects
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
	ChandrasekharWD(double, double, int);
	ChandrasekharWD(double, int, const double);
	virtual ~ChandrasekharWD(); //destructor
	int length(){return len;}
	//these three functions specify units
	double Radius();	//total radius
	double Mass();//total mass
	double Gee();//
	//in Newtonian, light speed is infinity...
	double light_speed2();
	
	//these return value of indicated variable -- used in testing
	double getXi(int X){return xi[X];}
	double getX(int X){return x[X];}
	double getY(int X){return y[X];}
	double getYderiv(int X){return z[X];}
	
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
	//the values of K, rho0, P0 relate to the choice of units and scale factors
	//going to set all to 1
	//to rescale, basically "finish" the part of Rn we left off, sqrt[Pc/(Grho^2)]
	// radius scales like Rn
	// mass scales like rho0*Rn3
	double A0;		//central pressure
	double B0;		//central density
	double Rn;		//radius scale factor
	int len;
	//lane-emden solution functions
	double *xi;	//normalized radius
	double *x;  //the relativity factor x = pF/mc
	double *y;	//Chandrasekhar's y, y^2=1+x^2
	double *z;	//derivative (dy/dxi)  note: dx/dxi = (dy/dxi)/x
	double *mass;
	double *f;
		
	//integrate using basic RK4
	double RK4integrate(const int, double);
	int RK4integrate(const int, double, int);
	
	//methods for handling the BCs
	double yc[4], xc[2], fc[2];	//series coefficients of y,x,f near center
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
	void printCoefficients(char *c);
};

#endif