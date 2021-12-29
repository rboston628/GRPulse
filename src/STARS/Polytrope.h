//**************************************************************************************
//							POLYTROPIC STAR
// Polytrope.h
//		A simple stellar object in Newtonian physics obeying everywhere a polytropic
//			equation of state of the form 
//				P~rho^Gamma
//		Solves Lane-Emden equation, y=theta, x=xi, z=dy/dxi in more usual notation
//		See Hansen & Kawaler Chapter 7 for further information
//**************************************************************************************

#ifndef POLYTROPEH
#define POLYTROPEH

#include "Star.h"

class Polytrope : public Star {
public:

	void graph_title(char* buff){
		sprintf(buff, "polytrope n=%1.2f", n);
	}
	//char* name(){return polyname;};
	
	Polytrope(double, double, double, int);	//constructor, M, R, n and length
	Polytrope(double, int);					//constructor, n and length
	Polytrope(double, int, const double);	//constructor, n and length, dx
	virtual ~Polytrope();  //destructor
	int length(){return len;}
	//these three functions specify units
	double Radius();	//total radius
	double Mass();		//total mass
	double Gee(){return GG;};
	//in Newtonian, light speed is infinity... just use the max value to represent this
	virtual double light_speed2(){return std::numeric_limits<double>::max();};
	
	//these return value of indicated variable -- used in testing
	double getX(int X){     return Y[X][x];}
	double getY(int X){     return Y[X][y];}
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
	int len;
	double n;		//polytropic index
	double Gamma;	//polytropic exponent
	//the values of K, rho0, P0 relate to the choice of units and scale factors
	//going to set all to 1
	//to rescale, basically "finish" the part of Rn we left off, sqrt[Pc/(Grho^2)]
	// radius scales like finished Rn
	// mass scales like rho0*Rn3
	double rho0;	//central density
	double P0;		//central pressure
	double Rn;		//radius scale factor
	//lane-emden solution functions
	enum VarName {x=0, y, z, numvar};
	double **Y;   //at each point X, variables x, y, z
	double *mass;
	double *base; //= pow(y,n-1), avoids repeated calls to pow(y,n)
	double GG;
	
	void centerInit(double ycenter[numvar]);
	void RK4step(double dx, double yin[numvar], double yout[numvar]);
	
	//integrate Lane-Emden using basic RK4
	double RK4integrate(const int, double);
	//a grid-multiplying RK4 method
	int RK4integrate(const int, double, int);
	//find the location of edge of star for last step of solution
	void findEdge(int);
	double integrateEdge(int Nedge, double dx);
	
	//methods for handling the BCs
	double ac[4], as[6];	//expansion coefficients of y near center, surface
	void setupCenter();		//prepare values of ac[]
	void setupSurface();	//prepare values of as[]
	char polyname[40];		//a name of the polytrope... not yet implemented properly
public:
	//methods to find central, surfae power series expansions of key variables in pulsation
	void getAstarCenter(double *, int&, double g=0);
	void getUCenter(double*, int&);
	void getVgCenter(double*, int&, double g=0);
	void getC1Center(double*, int&);
	void getAstarSurface(double*, int&, double g=0);
	void getUSurface(double*, int&);
	void getVgSurface(double*, int&, double g=0);
	void getC1Surface(double*, int&);
};

#endif