#ifndef SPLINORCLASS
#define SPLINORCLASS

#include "Splinor.h"

Splinor::Splinor(const double *const x, const double *const y, const int L){
	len = L;
	xlast = 0;
	incr = +1;
	xa = x[0];
	xb = x[L-1];
	
	S = new double[len];
	xarr = new double[len];
	yarr = new double[len];
	
	
	//the coefficient and interpolatio methods assume that x is strictly increasing
	//if x is decreasing, then read it in backwards.
	if(xa<xb){
		for(int i=0; i<len; i++){
			xarr[i] = x[i];
			yarr[i] = y[i];	
		}
	}
	else if (xa>xb){
		double xt = xb;
		xb = xa;
		xa = xt;
		for(int i=0; i<len; i++){
			xarr[i] = x[len-1-i];
			yarr[i] = y[len-1-i];
		}
	
	}	
	getSplineCoefficients(xarr, yarr);
}

Splinor::~Splinor(){
	delete[] S;
	delete[] xarr;
	delete[] yarr;
}

double Splinor::operator()(double xpos){
	return interp(xpos);
}

double Splinor::interp(double xpos){
	int xind = xlast;	//begin at position we last checked
	if(xarr[xind] > xpos) xind =0;
	bool wrap = false;	//flag, if we have already wrapped around array
	if      (xpos == xb) return yarr[len-1];	//if at the end, give the end
	else if (xpos == xa) return yarr[0];		//if at beginning, give beginning
	else if (xpos  < xa) xind = 0;
	else if (xpos  > xb) xind = len-2;
	//if(xpos > xb || xpos < xa) return 0.0;//outside range, return 0
	else {
		while(xarr[xind] <= xpos){	//scan through range until we pass it	
			xind += incr;
			if (xind >len-1 | xind < 0) {
				//if we reach end without finding it,
				//wrap back around one time
				if(!wrap) xind = (xind >=0 ? 0 : len-2);
				//if we've already wrapped, quit with error
				else if (wrap){
					return nan("array wrap");
				}
				//set flag to note we've wrapped through array once
				wrap = true;
			}
		}
		//we passed the point we need, so rewind once 
		xind -= incr;
	}
	xlast = xind; //remember position for next time we look - normally fit in order
	
	double a, b, c, d, h0, h1;
	double t, dx;
	h1 = xarr[xind+1]-xarr[xind  ];
	a = (S[xind+1]-S[xind])/(6.*h1);
	b = S[xind]/2;
	c = (yarr[xind+1]-yarr[xind])/h1 - h1*(2.*S[xind]+S[xind+1])/6.;
	d = yarr[xind];
	t = xpos - xarr[xind];
	return t*( t*(a*t+b)+c ) + d;
}

double Splinor::deriv(double xpos){
	int xind = xlast;	//begin at position we last checked
	if(xarr[xind] > xpos) xind = 0;
	bool wrap = false;	//flag, if we have already wrapped around array
	if      (xpos  <= xa) xind = 0;
	else if (xpos  >= xb) xind = len-2;
	else {
		//while(xpos >= xarr[xind]){	//scan through range until we pass it	
		while(xarr[xind] <= xpos){
			xind++;
			if (xind >len-1) {
				//if we reach end without finding it,
				//wrap back around one time
				if(!wrap) xind = 0;
				//if we've already wrapped, quit with error
				else if (wrap){
					return nan("");
				}
				//set flag to note we've wrapped through array once
				wrap = true;
			}
		}
		//we passed the point we need, so rewind once 
		xind--;
	}
	xlast = xind; //remember position for next time we look - normally fit in order
		
	double a, b, c, d, h0, h1;
	double t, dx;
	h1 = xarr[xind+1]-xarr[xind  ];
	a = (S[xind+1]-S[xind])/(6.*h1);
	b = S[xind]/2;
	c = (yarr[xind+1]-yarr[xind])/h1 - h1*(2.*S[xind]+S[xind+1])/6.;
	t = xpos - xarr[xind];
	return t*(3.*a*t+2.*b)+c;//formula for derivative of splinor
}

void Splinor::getSplineCoefficients(
	const double *const x,
	const double *const y
){
	int N = len-1;
	//set up spline equations
	S[0] = S[N] = 0.0;
	double *a = new double[len];
	double *b = new double[len];
	double *c = new double[len];
	double *d = new double[len];
	double *gamma = new double[len];
	double beta;
	//top and bottom rows of coefficient arrays are empty
	a[0] = b[0] = c[0] = d[0] = 0.0;
	a[N] = b[N] = c[N] = d[N] = 0.0;
	//fill in coefficient arrays for solution
	a[1] = 0.0;
	b[1] = 2.*(x[2]-x[0]);
	c[1] = x[2]-x[1];
	d[1] = 6.0*( (y[2]-y[1])/(x[2]-x[1])-(y[1]-y[0])/(x[1]-x[0]) );
	for(int k=2; k<N; k++){
		a[k] = x[k]-x[k-1];
		b[k] = 2.*( x[k+1]-x[k-1] );
		c[k] = x[k+1]-x[k];
		d[k] = 6.*( (y[k+1]-y[k])/(x[k+1]-x[k])-(y[k]-y[k-1])/(x[k]-x[k-1]) );
	}
	a[N-1] = x[N-1]-x[N-2];
	b[N-1] = 2.*( x[N]-x[N-2] );
	c[N-1] = 0.0;
	d[N-1] = 6.*( (y[N]-y[N-1])/(x[N]-x[N-1])-(y[N-1]-y[N-2])/(x[N-1]-x[N-2]) );
		
	//backsubstitute
	gamma[1] = 0.0;
	beta = b[1];
	S[1] = d[1]/beta;
	for(int k=2; k<N; k++){
		gamma[k] = c[k-1]/beta;
		beta = b[k] - a[k]*gamma[k];
		S[k] =(d[k] - a[k]*S[k-1])/beta;
	}
	for(int k=(N-2); k>=1; k--){
		S[k] -= gamma[k+1]*S[k+1];
	}
	delete[] gamma;
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d;
}


#endif