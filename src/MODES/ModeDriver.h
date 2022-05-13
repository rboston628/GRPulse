#ifndef MODEDRIVERH
#define MODEDRIVERH

#include "../constants.h"
#include "../STARS/Star.h"


//**************************************************************************************
//							MODE ABSTRACT BASE CLASS
//		This is an abstract base class for all Mode objects, allowing polymorphism
//		Can be used in code as a generic mode, without regard to particular physics
//		The generic Mode<N> defined in Mode.h cannot always be used polymorphically
//		When a mode is needed and exact physics unknown, use ModeBase as type
//**************************************************************************************
class ModeBase {
public:
	//file output methods to write and view plots of mode
	virtual void printMode(char *c = NULL) = 0;
	virtual void writeMode(char *c = NULL) = 0;
		
	virtual ~ModeBase(){};
	
	virtual int modeOrder() = 0;
	virtual void modeNumbers(int&, int&, int&) =0;
	virtual double getOmega2() = 0;
	virtual double SSR() = 0;
	virtual double tidal_overlap() =0;
	virtual double getFreq() =0;
	virtual double getPeriod()=0;
	
	virtual double getRad(int x) =0;
	virtual double getY(int i, int x) =0;
	virtual double getYtilde(int i, int x) =0;
};

//**************************************************************************************
//							MODE DRIVER ABSTRACT CLASS
// ModeDriver.h
//		An abstract class for all sets of mode equations to inherit
//		A ModeDriver object will provide the physical equations
//			and he boundary conditions obeyed by the individual modes
//		A single ModeDriver is capable of serving many Modes
//		Serves as the basis in code for other models of pulsations
//**************************************************************************************
class ModeDriver {
public:
	const size_t num_var;
	//constructor
	ModeDriver(int nv, Star *s) : num_var(nv), star(s)  {};
	virtual ~ModeDriver(){};
	virtual int length() =0;
	virtual double Gamma1() =0;
	virtual double rad(int) =0;

	virtual int CentralBC(double **y, double *yo, double s2, int l, int m=0) =0;
	virtual int SurfaceBC(double **y, double *ys, double s2, int l, int m=0) =0;
	virtual void getCoeff(double *CC, const int, const int, const double, const int) =0;
	virtual void setupBoundaries() =0;
	
	virtual double SSR(double, int, ModeBase*) =0;
	virtual double tidal_overlap(ModeBase*) =0;
	virtual double innerproduct(ModeBase*, ModeBase*) =0;
		
	//the following two methods are added to make Mode agnostic
	//provide the order of BCs to check in forming Wronskian
	virtual void getBoundaryMatrix(int, double *, double*, double **, int*) =0;
	virtual void varnames(std::string*) =0;	//names of variables to print out

protected:
	Star *const star;
	int central_bc_order;
	int surface_bc_order;
	
	friend class ModeBase;
	template <size_t N> friend class Mode;
	friend class Star;
};

#endif