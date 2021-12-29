// *************************************************************************************
//					POST-NEWTONIAN STAR ABSTRACT CLASS
// PNStar.h
//		This class is another abstract class
//		It extends the definition of a Star to Post-Newtonian physics
// *************************************************************************************

#ifndef PNSTARH
#define PNSTARH

#include "Star.h"

class PNStar : public Star {

public:
	//the potentials defined only for the post-newtonian case
	virtual double Psi(int )=0, dPsidr(int )=0;
	virtual double  Wx(int )=0,  dWxdr(int )=0;
	virtual double  Wy(int )=0,  dWydr(int )=0;
	virtual double  Wz(int )=0,  dWzdr(int )=0;
	
	//must return stellar structure expanded in powers of x=r/R near center
	virtual void getBetaCenter( double*, int&, double g=0) =0;
	virtual void getPhiCenter(  double*, int&) =0;
	//must return stellar structure expanded in powers of t=1-r/R near surface
	virtual void getBetaSurface(double*, int&, double g=0) =0;
	virtual void getPhiSurface( double*, int&) =0;
	
	//the redshift
	virtual double getRedshift() =0;
};

#endif