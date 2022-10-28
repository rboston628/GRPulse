#ifndef ROOTFINDINGMETHODS_C
#define ROOTFINDINGMETHODS_C

#include "rootfind.h"
#include <complex>

double pseudo_unif(){
	static int a=53, b=122, r=17737, dart = 39;//for pseudorandom positioning
	dart = (a*dart+b)%r; //generates a psuedo-random integer in (0,r)
	return (double(dart)/double(r));
}

// *************************************************************************************
//					BISECTION METHODS
//  These methods are for implementing bisection searches in one parameter
// *************************************************************************************

// this will find brackets using a simple expansion/movement search
// one bracket is moved to one side until the function changes sign, to find brackets
void bisection_find_brackets_move(
	std::function<double(double)> func,		//the function to find zero of
	double x0,								//an initial guess for the zero
	double& xmin,							//the lower bracket -- will be returned
	double& xmax							//the upper bracket -- will be returned
){
	double x=x0, y = func(x), ymin, ymax;
	double dx = x0;
	
	//find brackets on x that bound a zero in y
	
	if(y > 0){
		xmin = x; ymin = y;
		while(y > 0 && !isnan(y)){
			x += dx;
			y = func(x);
			if(isnan(y)){
				x = x0; dx *= 0.1; y = 1.0;
			}
		}
		xmax = x; ymax = y;
	}
	else if (y < 0){
		xmax = x; ymax = y;
		while(y < 0){
			x *= 0.5;
			y = func(x);
		}
		xmin = x; ymin = y;
	}

	//swap if they are backwards
	if(xmin > xmax){
		double temp = xmax;
		xmax = xmin; xmin = temp;
		temp = ymax;
		ymax = ymin; ymin = temp;
	}
}

// this will find brackets using newton's method to look for the next zero, then a bit beyond it
// why use one zero-finding method to prepare another zero-finding method? 
//    because Newton's method can fail, but bisection searches can go to nearly arbitrary accuracy
void bisection_find_brackets_newton(
	std::function<double(double)> func,		//the function to find zero of
	double x,								//an initial guess for the zero
	double& xmin,							//the lower bracket -- will be returned
	double& xmax							//the upper bracket -- will be returned
){
	double y1=func(x), y2=y1;	
	double x1=x, x2=x;
	double xdx, ydx;
	double ymax, ymin;
	//while the two ys are on same side of axis, keep reposition until zero is bound
	//we use Newton's method to the nearest zero
	while(y1*y2 >= 0.0){
		//compute numerical derivative
		xdx = 1.01*x2;
		ydx = func(xdx);
		xdx = x2 - y2*(xdx-x2)/(ydx-y2);
		//Newton's method can fail at this if it only approaches the zero from one direction
		//In such a case, slightly broaden the bracket to get to other side of zero
		if( xdx == x2 ) {
			if(xdx>x1) xdx *= 1.01;
			if(xdx<x1) xdx *= 0.99;
		}
		//limit amount of change allowed in single step
		if(xdx > 1.1*x2) xdx = 1.1*x2;	//if we increased, don't increase too much
		if(xdx < 0.9*x2) xdx = 0.9*x2;	//if we decreased, don't decrease too much
		ydx = y2;
		y2 = func(xdx);
		x2 = xdx;
	}
	//sometimes the brackets will be on either end of hyperbolic divergence
	//check that solution changes continuously between the two brackets
	double x3 = 0.5*(x1+x2), y3 = func(x3);
	double scale = 1.01;
	while( (y3*y1>=0.0 & fabs(y3)>fabs(y1)) | (y3*y2>=0.0 & fabs(y3)>fabs(y2))){
		//this can sometimes be solved by taking broader steps in secant method
		scale *= 1.1;
		x2=x1;
		y2=y1;
		while(y1*y2 >= 0.0){
			//compute numerical derivative
			xdx = scale*x2;
			ydx = func(xdx);
			xdx = x2 - y2*(xdx-x2)/(ydx-y2);
			//Newton's method can fail at this if it only approaches the zero from one direction
			//In such a case, slightly broaden the bracket to get to other side of zero
			if( xdx == x2 ) {
				if(xdx>x1) xdx *= 1.01;
				if(xdx<x1) xdx *= 0.99;
			}
			ydx = y2;
			y2 = func(xdx);
			x2 = xdx;
		}
		x3 = 0.5*(x1+x2);
		y3 = func(x3);
	}
	
	if(x2 > x1){
		xmax = x2; ymax = y2;
		xmin = x1; ymin = y1; 
	}
	else {
		xmax = x1; ymax = y1;
		xmin = x2; ymin = y2;
	}
	//swap if they are backwards
	if(xmin > xmax){
		double temp = xmax;
		xmax = xmin; xmin = temp;
		temp = ymax;
		ymax = ymin; ymin = temp;
	}


}

//given brackets bounding a single zero, find the zero
// the brackets xmin, xmax MUST bound a single zero
double bisection_search(
	std::function<double(double)> func,		//the function to find zero of
	double &x,								//the location of zero -- will be returned
	double xmin,							//the lower bracket
	double xmax								//the upper bracket
){
	//now use bisection to find dx so that y=0.0
	x = 0.5*(xmin+xmax);
	double y = func(x), y2 = y;
	double ymin = func(xmin), ymax=func(xmax);
	double xold = x;
	while( fabs(y)>0.0 || isnan(y) ){	
		//printf("BISECT [%le %le], %le\n", xmin, xmax, x);	
		if( (y*ymax>0.0) ){
			xmax = x;
			ymax = y;
		}
		else if( (y*ymin>0.0) ){
			xmin = x;
			ymin = y;
		}
		x = 0.5*(xmin+xmax);
		y = func(x);
		
		//it can happen that the search becomes stuck
		// in this case, picking pseudorandom location within brackets can sometimes help
		if(y2==y){//if the problem becomes stuck in a loop
			x = xmin + pseudo_unif()*fabs(xmax-xmin);
			y = func(x);
		}
		y2 = y;		
		//if the brackets are not moving, stop the search
		if(xold == fabs(xmax-xmin)) break;
		xold = fabs(xmax-xmin);
	}
	x = 0.5*(xmin+xmax);
	return func(x);
}



// *************************************************************************************
//					NEWTON METHODS
//  These methods are for gradient-descent methods
// *************************************************************************************

template<>
double gen_abs<double> (double x){
	return fabs(x);
}

template<>
std::complex<double> gen_abs(std::complex<double> x){
	return std::abs(x);
}

#endif