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

#ifndef ISOPYCNICH
#define ISOPYCNICH

#include "Star.h"

class Isopycnic : public Star {
public:
	void graph_title(char* buff){
		sprintf(buff, "uniform density");
	}
	
	Isopycnic(int);			//constructor, n and length
	virtual ~Isopycnic();	//destructor
	int length(){return len;}
	//these three functions specify units
	double Radius();		//total radius
	double Mass();			//total mass
	double Gee(){return GG;};
	//in Newtonian, light speed is infinity... just use the max value to represent this
	virtual double light_speed2(){return std::numeric_limits<double>::max();};
			
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
	double Gamma;	//polytropic exponent
	double rho0;	//central density
	double P0;		//central pressure
	double Rn;		//radius scale factor
	//lane-emden solution functions
	double *x;	//normalized radius (xi)
	double *y;	//lane-emden solution (theta)
	double *z;	//derivative (dtheta/dxi)
	double *mass;
	double GG;
	
	//integrate Lane-Emden using basic RK4
	int populateValues(const int, double);
	
	//methods for handling the BCs
	void setupCenter();		//for conformity
	void setupSurface();	//for conformity
	char isoname[40];		//a name of the star
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