//**************************************************************************************
//						POST-NEWTONIAN POLYTROPIC STAR
// PNPolytrope.cpp
//		A simple stellar object in 1PN physics obeying everywhere a polytropic
//			equation of state of the form 
//				P~rho^Gamma
//			where rho = total mass-energy density
// 		compare to Tooper 1964
//		derivations in 1pn4_polytrope.nb
//		Based on my own work with perturbed 1PN equations of stellar structure
// 		We use units where G=c=1
//**************************************************************************************

#ifndef PNPOLYTROPECLASS
#define PNPOLYTROPECLASS

#include "PNPolytrope.h"

//initalize polytrope from index and length, fit to mass M, radius R
PNPolytrope::PNPolytrope(double BigM, double BigR, double n, int L)
	: n(n), len(L), zsurf(G_CGS*BigM/C_CGS/C_CGS/BigR), Gamma(1.0+1.0/n),
	GG(G_CGS), cc(C_CGS*C_CGS)
{
	sigma = (1.-pow(1.+zsurf,-2))/(2.*(n+1));
	//if index is out of range, fail
	if(n >= 5 || n < 0){
		n = nan(""); len = 0;
		printf("Invalid polytropic index.  Polytrope not initialized.\n");
		exit(EXIT_FAILURE);
	}
	if(zsurf >= 0.5){
		printf("\tYour specified star would collapse to a black hole.\n");
		printf("\tSwitching to a smaller value.\n");
		zsurf = 0.49;
		sigma = zsurf/(1.+n);
	}
	if(sigma> n/(n+1.) &&n!=0){
		printf("\tOverlarge value for sigma, resorting to max value.\n");
		sigma = n/(n+1.);
	}
	if(n==0) Gamma = 0.0;
	//name this polytrope for files
	sprintf(name, "PNPolytrope%1.1f_%1.2le", n, sigma);
	
	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with y[len-1]=0.0
	//we will find the proper dx with a bisection search
	
	//an initial guess -- based on 0pn constant-volume model
	double dx = sqrt(6.0)/len, yS=1.0, ddx = dx;
	φc0 = -1.0;
	Y = new double*[len];
	for(int i=0; i<len; i++) 
		Y[i] = new double[numvar];

	double dxmax=1.0, dxmin=0, ySmax=-1.0, ySmin=1.0;
	
	
	//we first find reasonable values of the parameters φc0 and sigma
	// this step is necessary for large values of z
	int tracker = 0;
	double s1=sigma, ds =0.01*sigma, dds=sigma;
	std::complex<double> target(zsurf, zsurf);
	std::complex<double> f1 = RK4integrate(s1, dx), f2 = RK4integrate(s1+ds,dx);
	printf("target = \t%0.16le\n", target.real());
	printf("result = \t%0.16le\n", f1.real());
	printf("       = \t%0.16le\n", f1.imag());
	while(fabs(dds) > 1.e-16 && tracker < 100){
		dds = ds*(target-f1).real()/(f2-f1).real();
		s1 = s1 + dds;
		ds = dds;
		f1 = RK4integrate(s1,   dx);
		f2 = RK4integrate(s1+ds,dx);
		printf("result = \t%0.16le\n", f1.real());
		printf("       = \t%0.16le\n", f1.imag());
	}
	sigma = s1;
	f1 = RK4integrate(sigma, dx);
	tracker=0;
	while(abs(target-f1) > 1.e-16 && tracker<100){
		f1 = RK4integrate(sigma, dx);
		tracker++;
	}
	printf("target = \t%0.16le\n", target.real());
	printf("result = \t%0.16le\n", f1.real());
	printf("       = \t%0.16le\n", f1.imag());
	
	//now find dx so that surface occurs in len steps
	yS = RK4integrate(len, dx);
	//find brackets on dx that bound a zero in yS
	if(yS > 0){
		dxmin = dx; ySmin = yS;
		while(yS > 0 && !isnan(yS)){
			dx += ddx;
			yS = RK4integrate(len, dx);
			if(isnan(yS)){
				if(n!= int(n)) yS = -1.0;
				else {
					dx = 1.0/len; ddx *= 0.1; yS = 1.0;
				}
			}
		}
		dxmax = dx; ySmax = yS;
	}
	else if (yS < 0){
		dxmax = dx; ySmax = yS;
		while(yS < 0){
			dx = 0.99*dx;
			yS = RK4integrate(len,dx);
		}
		dxmin = dx; ySmin = yS;
	}
	dx = 0.5*(dxmin+dxmax);
	yS = RK4integrate(len, dx);
	
	//some possible errors that might occur
	if(dx<0) {printf("somehow dx is negative...\n"); dx=-dx;}
	if(ySmin*ySmax > 0.0) {printf("big problem, chief\n"); exit(EXIT_FAILURE);}
	
	//use bisection to find dx so that yS=0.0
	double dxold = 1.0;
	while( fabs(dxmin-dxmax)>0.0 || isnan(yS) ){
		dx = 0.5*(dxmin+dxmax);
		yS = RK4integrate(len, dx);
		
		if(isnan(yS)){
			yS = -1.0;
		}
		if( (yS*ySmax>0.0) ){
			dxmax = dx;
			ySmax = yS;
		}
		else if( (yS*ySmin>0.0) ){
			dxmin = dx;
			ySmin = yS;
		}
		//if the brackets are not moving, stop the search
		if(dxold == fabs(dxmin-dxmax)) break;
		dxold = fabs(dxmin-dxmax);
	}
	
	sigma = setSigma(Y[len-1]);
	
	//check that sigma, phi0 lead to correct surface values
	f1 = setZsurf(Y[len-1]);
	printf("target = \t%0.16le\n", target.real());
	printf("result = \t%0.16le\n", f1.real());
	printf("       = \t%0.16le\n", f1.imag());
	printf("sigma=%le\n", sigma);
	RK4integrate(len, dx,1);
	f1 = setZsurf(Y[len-1]);
	//check that end is at surface
	printf("sigma=%le\n", sigma);
	printf("y1 = %le\n", Y[len-1][y]);
	
	//continue adjusting the value of sigma until the surface redshift is correct
	tracker = 0;
	while(fabs( zsurf-f1.imag() ) > 1e-16  & tracker<100){
		φc0  -= (Y[len-1][u]*Y[len-1][x] + Y[len-1][phi]);
		sigma = setSigma(Y[len-1]);
		RK4integrate(len, dx,1);
		f1 = setZsurf(Y[len-1]);
		tracker++;
	}
	RK4integrate(len, dx,1);
	f1 = setZsurf(Y[len-1]);
	printf("target = \t%0.16le\n", target.real());
	printf("result = \t%0.16le\n", f1.real());
	printf("       = \t%0.16le\n", f1.imag());
	printf("sigma=%le\n", sigma);
	printf("y1 = %le\n", Y[len-1][y]);
	
	//now set physical properties of the polytrope
	rho0 = BigM/(4.*m_pi)*pow(BigR,-3)*(-Y[len-1][x]/Y[len-1][z]);
	P0   = GG*pow(BigM,2)*pow(BigR,-4)/(4.*m_pi)*pow(Y[len-1][z],-2)/(n+1.);
	Rn   = BigR/Y[len-1][x];
	
	//calculate derivative and other useful thngs
	dydx = new double[len];
	base = new double[len];
	dydx[0] = 0.0;
	base[0] = 1.0;
	//
	indexFit = len/2;
	for(int X=1; X<len; X++){
		dydx[X] = -((1.+sigma*Y[X][y])*Y[X][u] + sigma*Y[X][v]);
		//as we scan through x,y,z, set matching point as where y[X] = 0.5
		if(Y[X-1][y]>0.5 & Y[X][y]<=0.5) indexFit = X;
		//the array base saves some calclations
		base[X] = pow(Y[X][y],n-1.0);
	}
	base[len-1] = 0.0;
	indexFit /= 2;
	setupCenter();
	setupSurface();
}


//initalize polytrope from index, redshift, and length
PNPolytrope::PNPolytrope(double n, double zsurf, int Len)
	: n(n), zsurf(zsurf), sigma( (1.-pow(1.+zsurf,-2))/(2.*(n+1)) ), len(Len), Gamma(1.0+1.0/n),
	GG(1.), cc(1.)
{
	//if index is out of range, fail
	if(n >= 5 || n < 0){
		n = nan(""); len = 0;
		printf("Invalid polytropic index.  Polytrope not initialized.\n");
		exit(EXIT_FAILURE);
	}
	if(zsurf >= 0.5){
		printf("\tYour specified star would collapse to a black hole.\n");
		printf("\tSwitching to a smaller value.\n");
		zsurf = 0.49;
		sigma = zsurf/(1.+n);
	}
	if(sigma> n/(n+1.) &&n!=0){
		printf("\tOverlarge value for sigma, resorting to max value.\n");
		sigma = n/(n+1.);
	}
	if(n==0) Gamma = 0.0;
	//name this polytrope for files
	sprintf(name, "PNPolytrope%1.1f_%1.2le", n, sigma);

	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with y[len-1]=0.0
	//we will find the proper dx with a bisection search

	//an initial guess -- based on 0pn constant-volume model
	double dx = sqrt(6.0)/len, yS=1.0, ddx = dx;
	φc0 = -1.0;//2.*zsurf/(n+1.);
	Y = new double*[len];
	for(int i=0; i<len; i++) 
		Y[i] = new double[numvar];
	double dxmax=1.0, dxmin=0, ySmax=-1.0, ySmin=1.0;
	
	//we first find reasonable values of the parameters φc0 and sigma
	// this step is necessary for large values of z
	int tracker = 0;
	std::complex<double> target(zsurf, zsurf);
	std::complex<double> f1 = RK4integrate(sigma, dx);
	printf("target = \t%0.16le\n", target.real());
	printf("result = \t%0.16le\n", f1.real());
	printf("       = \t%0.16le\n", f1.imag());
	while(abs(target-f1) > 0.0 && tracker<100){
		f1 = RK4integrate(sigma, dx);
		tracker++;
	}
	printf("target = \t%0.16le\n", target.real());
	printf("result = \t%0.16le\n", f1.real());
	printf("       = \t%0.16le\n", f1.imag());
	
	//now find dx so that surface occurs in len steps
	yS = RK4integrate(len, dx);
	//find brackets on dx that bound a zero in yS
	if(yS > 0){
		dxmin = dx; ySmin = yS;
		while(yS > 0 || isnan(yS)){
			dx += ddx;
			yS = RK4integrate(len, dx);
			if(isnan(yS)){
				if(n!= int(n)) yS = -1.0;
				else {
					dx = 1.0/len; ddx *= 0.1; yS = 1.0;
				}
			}
		}
		dxmax = dx; ySmax = yS;
	}
	else if (yS < 0){
		dxmax = dx; ySmax = yS;
		while(yS < 0){
			dx = 0.99*dx;
			yS = RK4integrate(len,dx);
		}
		dxmin = dx; ySmin = yS;
	}
	dx = 0.5*(dxmin+dxmax);
	yS = RK4integrate(len, dx);
	
	if(dx<0) {printf("somehow dx is negative...\n"); dx=-dx;}
	if(ySmin*ySmax > 0.0) {printf("big problem, chief\n"); exit(EXIT_FAILURE);}

	//use bisection to find dx so that yS=0.0
	double dxold = 1.0;
	while( fabs(dxmin-dxmax)>0.0 || isnan(yS) ){
		dx = 0.5*(dxmin+dxmax);
		yS = RK4integrate(len, dx);
		
		if(isnan(yS)){
			yS = -1.0;
		}
		if( (yS*ySmax>0.0) ){
			dxmax = dx;
			ySmax = yS;
		}
		else if( (yS*ySmin>0.0) ){
			dxmin = dx;
			ySmin = yS;
		}
		//if the brackets are not moving, stop the search
		if(dxold == fabs(dxmin-dxmax)) break;
		dxold = fabs(dxmin-dxmax);
	}
	//check that sigma, phi0 lead to correct surface values
	f1 = setZsurf(Y[len-1]);
	printf("target = \t%0.16le\n", target.real());
	printf("result = \t%0.16le\n", f1.real());
	printf("       = \t%0.16le\n", f1.imag());
	printf("sigma=%le\n", sigma);
	RK4integrate(len, dx,1);
	//check that end is at surface
	printf("sigma=%le\n", sigma);
	printf("y1 = %le\n", Y[len-1][y]);
	
	//continue adjusting the value of sigma until the surface redshift is correct
	tracker = 0;
	while(fabs(zsurf-setZsurf(Y[len-1]).real()) > 1e-16  & tracker<100){
		φc0  -= (Y[len-1][u]*Y[len-1][x] + Y[len-1][phi]);
		sigma = setSigma(Y[len-1]);
		RK4integrate(len, dx,1);
		tracker++;
	}
	RK4integrate(len, dx,1);
	f1 = setZsurf(Y[len-1]);
	printf("target = \t%0.16le\n", target.real());
	printf("result = \t%0.16le\n", f1.real());
	printf("       = \t%0.16le\n", f1.imag());
	printf("sigma=%le\n", sigma);
	printf("y1 = %le\n", Y[len-1][y]);
	//now set physical properties of the polytrope
	//set initial density, pressure
	//if we use units c2=1, then P0 = sigma
	rho0 = 1.0;
	P0 = sigma;
	if(sigma==0) P0=1.0;
	Rn = sqrt( (n+1.0)*P0/(4.*m_pi*rho0*rho0) ); // 1/A as defined by Tooper (1964)
	printf("Rn=%le\n", Rn);
	printf("kappa=%le\n", P0*pow(rho0,-2));
	
	//calculate derivative and other useful thngs
	dydx = new double[len];
	base = new double[len];
	dydx[0] = 0.0;
	base[0] = 1.0;
	//
	indexFit = len/2;
	for(int X=1; X<len; X++){
		dydx[X] = -((1.+sigma*Y[X][y])*Y[X][u] + sigma*Y[X][v]);
		//as we scan through x,y,z, set matching point as where y[X] = 0.5
		if(Y[X-1][y]>0.5 & Y[X][y]<=0.5) indexFit = X;
		//the array base saves some calclations
		base[X] = pow(Y[X][y],n-1.0);
	}
	base[len-1] = 0.0;
	indexFit /= 2;
	setupCenter();
	setupSurface();
}

//initalize polytrope from index, sigma, length, and dx
PNPolytrope::PNPolytrope(double n, double sigma, int Len, const double dx)
	: n(n), sigma(sigma), zsurf(1./sqrt(1.-2.*(n+1.)*sigma) - 1.),len(Len), 
	Gamma(1.0+1.0/n), GG(1.), cc(1.)
{
	zsurf = sigma;
	//if index is out of range, fail
	if( n > 5.0 || n < 0.0 ){
		n = nan(""); len = 0;
		printf("\nInvalid polytropic index.  Polytrope not initialized.\n");
		exit(EXIT_FAILURE);
	}
	if(n==0) Gamma = 0.0;

	//name this polytrope for files
	sprintf(name, "PNPolytrope%1.1f_%1.2le", n, sigma);
	//initialize arrays
	Y = new double*[len];
	for(int i=0; i<len; i++) 
		Y[i] = new double[numvar];
	
	φc0 = -1.0;
	RK4integrate(len, dx, 1);
	//continue adjusting the value of sigma until the surface redshift is correct
	int tracker = 0;
	while(fabs(zsurf-setZsurf(Y[len-1]).real()) > 1e-16  & tracker<100){
		φc0  -= (Y[len-1][u]*Y[len-1][x] + Y[len-1][phi]);
		sigma = setSigma(Y[len-1]);
		RK4integrate(len, dx,1);
		tracker++;
	}

	//set initial density, pressure
	rho0 = 1.0;
	P0 = sigma;
	if(sigma==0) P0=1.0;
	Rn = sqrt( (n+1.0)/(4.*m_pi) ); // 1/A as defined by Tooper (1964)
	
	//calculate derivative and other useful thngs
	dydx = new double[len];
	base = new double[len];
	dydx[0] = 0.0;
	base[0] = 1.0;
	for(int X=1; X<len; X++){
		dydx[X] = -((1.+sigma*Y[X][y])*Y[X][u] + sigma*Y[X][v]);
		//the array base saves some calclations
		base[X] = pow(Y[X][y], n-1.0);
		
	}
	indexFit = 512*round(double(len)/1024.0);
	indexFit /= 2;
	printf("  indexFit  = %d\n", indexFit);
	printf("r[indexFit] = %0.32le\t", rad(indexFit));
	printf("y[indexFit] = %0.32le\n", Y[indexFit][y]);
	printf("zsurf = %0.32le \t sigma = %0.32le\n", zsurf, sigma);
	setupCenter();
	setupSurface();
}


PNPolytrope::~PNPolytrope(){
	for(int i=0; i<len; i++)
		delete[] Y[i];
	delete[] base;
	delete[] dydx;
}


//  dz/dx = d^2(y)/dx^2
inline double PNPolytrope::deriv1PN(double s, double yy[numvar]){
	std::complex<double> YN(yy[y],0.0);
	YN = pow(YN,n);
	return -YN.real()*( 1. - (n+n+2.)*s*yy[phi] + s*4.*yy[y]  ) 
				- s*yy[u]*yy[z] - 2.*yy[z]/yy[x]; 
}

double PNPolytrope::setSigma(double ysurface[numvar]){
	return fabs( (1.-pow(1.+zsurf,-2))/(2.*(n+1.)*ysurface[u]*ysurface[x]) );
}
std::complex<double> PNPolytrope::setZsurf(double ysurface[numvar]){
	return std::complex<double> (
		1./sqrt(1.+2.*(n+1.)*sigma*ysurface[phi]          )-1.,
		1./sqrt(1.-2.*(n+1.)*sigma*ysurface[u]*ysurface[x])-1.
	);//*/
}

void PNPolytrope::centerInit(double ycenter[numvar]){
	ycenter[x] = 0.0;
	ycenter[y] = 1.0;
	ycenter[phi] = φc0;
	ycenter[z] = 0.0;	//dtheta/dx
	ycenter[u] = 0.0;   //dphi/dx
	ycenter[v] = 0.0;   //dpsi/dx	
}


void PNPolytrope::RK4step(double dx, double s, double yin[numvar], double yout[numvar]){
	double YC[numvar] = {yin[x], yin[y], yin[z], yin[u], yin[v], yin[phi]};
	double K[numvar][4];
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	std::complex<double> YCN(0.,0.);
	//calculate all Ks
	for(int a = 0; a<4; a++){
		//now from these, calculate next shift
		YCN = YC[y];
		YCN = std::pow(YCN,n);
		//K = dy = dx*(dy/dx)
		K[y][a] = -dx*((1.+sigma*YC[y])*YC[u] + s*YC[v]); //dtheta/dxi
		K[z][a] = dx*deriv1PN(s, YC);      //dz/dxi
		K[u][a] = dx*(YCN.real() - 2.*YC[u]/YC[x]);       //du/dxi
		K[v][a] = dx*(YCN.real()*( 3.*YC[y] - 2.*(n+1.)*YC[phi] ) - 2.*YC[v]/YC[x] ); //dv/dxi
		K[phi][a] = dx*YC[u];                          //dphi/dx
		//these boundary values at bottom of 1pnpolytrope1.nb
		if(YC[x]==0) {
			K[z][a] =-dx/3.0*(1.0 + 4.0*sigma - 2.*(n+1.)*sigma*φc0);
			K[u][a] = dx/3.0; K[v][a]=dx*(1.-2./3.*(n+1.)*φc0); 
			K[y][a] =K[phi][a]=0.0;}
		//calculate intermediate positions using previous shift vectors
		YC[x] = yin[x] + B[a]*dx;   //radius xi
		for(int b=1; b<numvar; b++)
			YC[b] = yin[b] + B[a]*K[b][a];
	}
	yout[x] = yin[x] + dx;
	for(int b=1; b<numvar; b++)
		yout[b] = yin[b] + K[b][0]/6.0 + K[b][1]/3.0 + K[b][2]/3.0 + K[b][3]/6.0;
}


//integrate the polytrope using RK4 until the end is found
std::complex<double> PNPolytrope::RK4integrate(double s, double& dx){
	double y1[numvar], y2[numvar];
	
	//set the initial conditions
	centerInit(y2);
	
	int X=0;
	while(y2[y] > 0.0 && !isnan(y2[y])){
		for(int b=0; b<numvar; b++) y1[b] =  y2[b];
		RK4step(dx, s, y1, y2);
		X++;
	}
	
	if(s==sigma) sigma = setSigma(y1);
	φc0 -= (y1[u]*y1[x] + y1[phi]);
	dx = y2[x]/double(len+1);
	return setZsurf(y1);
}

//integrate the polytrope up to Len using RK4
double PNPolytrope::RK4integrate(const int Len, double dx){
	//set the initial conditions
	centerInit(Y[0]);
	
	//integrate with RK4
	int X=0;
	for(X = 0; X<Len-1; X++){
		RK4step(dx, sigma, Y[X], Y[X+1]);
		//sometimes, for noninteger n, if array is too large values y<0 lead to nans
		//if these occur, it is safe to terminate the integration
		if(Y[X+1][y] <0.0){
			φc0 -= (Y[X][u]*Y[X][x] + Y[X][phi]);
			sigma = setSigma(Y[X]);
			return Y[X+1][y];
		}
	}
	//adjust sigma and φc0 based on surface values
	φc0 -= (Y[X][u]*Y[X][x] + Y[X][phi]);
	sigma = setSigma(Y[X]);
	return Y[X][y];
}

int PNPolytrope::RK4integrate(const int Len, double dx, int grid){
	grid=1;	
	//set the initial conditions
	centerInit(Y[0]);

	//integrate with RK4
	for(int X = 0; X<Len-1; X++){
		RK4step(dx, sigma, Y[X], Y[X+1]);
	}
	return Len;
}

//Here we define functions to access radius, pressure, etc.
double PNPolytrope::rad(int X){
	return Rn*Y[X][x];
}
double PNPolytrope::rho(int X){
	if(n==0) return rho0;
	else     return rho0*base[X]*Y[X][y];
}
double PNPolytrope::drhodr(int X){
	if(n==0) return 0.0;
	else     return n*rho0*base[X] * (dydx[X]/Rn);
}
double PNPolytrope::P(int X){
	if(n==0) return P0*Y[X][y];
	else     return P0*base[X]*Y[X][y]*Y[X][y];
}
double PNPolytrope::dPdr(int X){
	if(n==0) return P0 * dydx[X]/Rn;
	else     return (n+1.)*P0*base[X]*Y[X][y] * (dydx[X]/Rn);
}
double PNPolytrope::Phi(int X){
	return (n+1.)*P0/rho0 * Y[X][phi];
}
double PNPolytrope::dPhidr(int X){
	return (n+1.)*P0/rho0 * Y[X][u]/Rn;
}
double PNPolytrope::mr(int X){
	//m = 4 pi Int[r^2 rho,r] 
	//  = 4 pi Rn^3 rho0 Int[x^2 theta^n,x] 
	//  = 4 pi Rn^3 rho0 x^2 dphi/dx
	return 4.*m_pi*rho0*pow(Rn,3)*Y[X][u]*pow(Y[X][x],2);
}
double PNPolytrope::Psi(int){
	return 0.0; //TO BE CONTINUED...
}
double PNPolytrope::dPsidr(int X){
	return (n+1.)*P0*P0/rho0/rho0 * Y[X][v]/Rn;
}
//the gravitomagnetic potential vansihes for non-rotating star
double PNPolytrope::Wx(int X){return 0.0;}
double PNPolytrope::dWxdr(int ){return 0.0;}
double PNPolytrope::Wy(int ){return 0.0;}
double PNPolytrope::dWydr(int ){return 0.0;}
double PNPolytrope::Wz(int ){return 0.0;}
double PNPolytrope::dWzdr(int ){return 0.0;}

double PNPolytrope::sound_speed2(int X, double GamPert){
	// vs2 = Gamma1 * P/(rho+P/c^2) 
	//     = Gamma1* P0 y^n+1 /(rho0 y^n + P0 y^n+1/c2)
	//     = Gamma1*P0*y/(rho0 + P0 y/c^2)
	if(GamPert==0.0) return Gamma  *P0/rho0*Y[X][y]/(1.0+sigma*Y[X][y]);
	else             return GamPert*P0/rho0*Y[X][y]/(1.0+sigma*Y[X][y]);
}

//Schwarzschild discriminant as for GR case
double PNPolytrope::Schwarzschild_A(int X, double GamPert){
	if(GamPert==0.0) return dydx[X]/Y[X][y]/Rn* (n/(1.+sigma*Y[X][y])-(n+1.)/Gamma);
	else             return dydx[X]/Y[X][y]/Rn* (n/(1.+sigma*Y[X][y])-(n+1.)/GamPert);	
}

double PNPolytrope::getAstar(int X, double GamPert){
	if(GamPert==0.0) return Y[X][x]*dydx[X]/Y[X][y]* ((n+1.)/Gamma   - n/(1.+sigma*Y[X][y]));
	else             return Y[X][x]*dydx[X]/Y[X][y]* ((n+1.)/GamPert - n/(1.+sigma*Y[X][y]));	
}

double PNPolytrope::getU(int X){
	// U = 4piG*r*(rho + 3P - 2 rho*Phi)/g
	if(X==0) return 3.0;
	return Y[X][x]*pow(Y[X][y],n)*(1.+sigma*(3.*Y[X][y]-2.*(n+1.)*Y[X][phi]))/(Y[X][u]+sigma*Y[X][v]);
}

double PNPolytrope::getVg(int X, double GamPert){
	if(GamPert==0.0) return -(n+1.)*Y[X][x]*dydx[X]/Y[X][y]/Gamma;
	else			 return -(n+1.)*Y[X][x]*dydx[X]/Y[X][y]/GamPert;
}

double PNPolytrope::getC(int X){
	if(X==0) return 0.5*Y[len-1][u]*Y[len-1][x]/(φc[1]+sigma*ψc[1]);
	return (Y[len-1][u]/Y[len-1][x])*Y[X][x]/(Y[X][u]+sigma*Y[X][v]);
}

double PNPolytrope::Gamma1(int X){
	return Gamma;
}

//integrate eq 3.8 of Tooper (1964) using trapezoidal rule
double PNPolytrope::getXproper(int X){
	double xbar = 0.0;
	double twosign = 2.0*sigma*(n+1.);
	double dx, int1 = 1.0, int2 = 1.0; //in lim x->0, z/x->0
	for(int i=0; i<X-1; i++){
		dx = Y[i+1][x] - Y[i][x];
		int1 = int2;
		int2 = 1.0/sqrt(1.0 + twosign*Y[i+1][u]*Y[i+1][x]);
		xbar += 0.5*dx*(int2+int1);
	}
	return xbar;
}

double PNPolytrope::Radius(){return rad(len-1);}	//total radius
double PNPolytrope::Mass(){return mr(len-1);}		//total mass



// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
void PNPolytrope::setupCenter(){
	//the central coefficients -- expanded in terms of xi
	double x1 = Y[len-1][x];
	double x12 = x1*x1;

	//the central coefficients
	φc[0] = φc0;
	φc[1] = x12/6.;
	φc[2] =-x12*x12*n/120.
			+ sigma*x12*x12*n/60.*(-2.+φc0*(n+1.));
	θc[0] = 1.0;
	θc[1] =-x12/6.
			+ sigma*x12/3.*((n+1.)*φc0-2.);
	θc[2] = x12*x12*n/120.
			+ sigma*x12*x12/180.*(10.+n*(15.-6.*φc0)-6.*n*n*φc0);
	ψc[0] = 0.0;
	ψc[1] = x12/6.*(3.-2*(n+1.)*φc0);
	ψc[2] = x12*x12/120.*(n+1.)*(2.*n*φc0-5.) 
			- sigma*x12*x12/60.*(n+1.)*(6-(3.+7.*n)*φc0+2.*n*(n+1.)*φc0*φc0);
}
void PNPolytrope::getAstarCenter(double *AC, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	//depending on power requested, return appropriate number of terms
	if(maxPow>=0) AC[0] = 0.0;
	if(maxPow>=2) AC[1] = 2.*θc[1]*( (n+1.)/Gam1 - n/(1.+sigma) );
	//notice that θc[1] is already carrying a factor x12
	//if more  terms than this requested, cap number of terms
	if(maxPow> 2) maxPow = 2;
}
void PNPolytrope::getVgCenter(   double *Vc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] = -2.*θc[1]*(n+1.)/Gam1; 
	if(maxPow> 2) maxPow = 2; 
}
void PNPolytrope::getUCenter(    double *Uc, int& maxPow){
	double x12 = Y[len-1][x]*Y[len-1][x];
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] =-x12*n/5.
						 + sigma*x12/5.*(2.*n*(n+1.)*φc0 - 6.*n - 5.);
	if(maxPow> 2) maxPow = 2;
}
void PNPolytrope::getC1Center(   double *cc, int& maxPow){
	double x1 = Y[len-1][x];
	double u1 = Y[len-1][u];
	double u1x1 = x1*u1;
	if(maxPow>=0) cc[0] = 0.5*u1x1/(φc[1]+sigma*ψc[1]);
	if(maxPow>=2) cc[1] =    -u1x1*(φc[2]+sigma*ψc[2])*pow(φc[1]+sigma*ψc[1],-2);
	if(maxPow> 2) maxPow = 2;
}

void PNPolytrope::getBetaCenter( double *bc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double Phi1 = (n+1.)*P0/rho0 * Y[len-1][phi];
	if(maxPow>=0) bc[0] = P0/rho0/Phi1*θc[0]/(1.+sigma);
	if(maxPow>=2) bc[1] = P0/rho0/Phi1*θc[1]*pow(1.+sigma,-2);
	if(maxPow> 2) maxPow = 2; 
}
void PNPolytrope::getPhiCenter(  double *pc, int& maxPow){
	if(maxPow>=0) pc[0] = φc0  /φs[0];
	if(maxPow>=2) pc[1] = φc[1]/φs[0];
	if(maxPow> 2) maxPow = 2;
}


// **************************  SURFACE BOUNDARY  ***************************************
// the following provide coefficients for surface expansions of A*, Vg, U, c1 in terms of t=1-r/R
//	Note that A*, Vg require a power -1
//	Note that up to order N, we only require:
//		A*, Vg, c1 up to order N-1
//		U          up to order N
//	Here maxPow represents the maixmum order of expansions of y1,..,y4 in LAWE
//	If maxPow = 0, we need terms -1
//	If maxPow = 1, we need terms -1, 0
//	If maxPow = 2, we need terms -1, 0, 1
void PNPolytrope::setupSurface(){
	double x1  =  Y[len-1][x];
	double θs1 = -Y[len-1][z]*x1;
	double φs1 = -Y[len-1][u]*x1;
	double ψs1 = -Y[len-1][v]*x1;
	double φs0 = ( pow(1.+zsurf,-2) - 1. )/(2.*(n+1.)*sigma);
	
	//first, we must find the coefficients in expansion theta = a1*t+a2*t^2 + ...
	θs[0] = 0.0;
	ψs[0] = 0.0;
	φs[0] = φs0;
	for(int i=1; i<4;i++) {
		θs[i] = θs1*(1. - double(i-1)/2.*sigma*φs1);
		φs[i] = φs1;
		ψs[i] = ψs1;
	}
	if(n==0){
		θs[2] += x1*x1*( -1./2. + sigma*φs0);
		θs[3] += x1*x1*( -1./3. + sigma*(2./3.*φs0 - 5./6.*θs[1] + φs[1]/2.) );
		φs[2] += x1*x1/2.;
		φs[3] += x1*x1/3.;
		ψs[2] -= x1*x1*φs0;
		ψs[3] += 0.5*x1*x1*θs[1] - 2./3.*x1*x1*φs0 - x1*x1*φs[1]/3.;
	}
	else if(n==1){
		θs[3] += θs[1]*x1*x1*( -1./6. + 2./3.*sigma*φs0);
		φs[3] += θs[1]*x1*x1/6.;
		ψs[3] -= 2./3.*x1*x1*θs[1]*φs0;
	}
	double gs[4];
	gs[0] = -   (φs[1]+sigma*ψs[1]);
	gs[1] = -2.*(φs[2]+sigma*ψs[2]);
	gs[2] = -3.*(φs[3]+sigma*ψs[3]);
	
	igs[0] = 1./gs[0];
	igs[1] = -gs[1]/(gs[0]*gs[0]);
	igs[2] = (gs[1]*gs[1]-gs[0]*gs[2])/(gs[0]*gs[0]*gs[0]);
}
void PNPolytrope::getAstarSurface(double *As, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma : g);
	double AV = (n*Gam1-n-1.)/Gam1;
	int O=1;
	if(maxPow>=O-1) As[O-1] = AV;
	if(maxPow>=O+0) As[O  ] = AV*(θs[2]/θs[1] -1.) - n*sigma*θs[1];
	if(maxPow>=O+1) As[O+1] = AV*(
					( θs[1]*(2.*θs[3]-θs[2]) - θs[2]*θs[2] )/(θs[1]*θs[1])
					+ n*sigma*(θs[1]-2.*θs[2]) );
	if(maxPow >O+1) maxPow = O+1;
}
void PNPolytrope::getVgSurface(   double *Vs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma : g);
	int O=1;
	if(maxPow>=O-1) Vs[O-1] = (n+1.)/Gam1;
	if(maxPow>=O+0) Vs[O  ] = (n+1.)/Gam1*( θs[2] - θs[1] )/θs[1];
	if(maxPow>=O+1) Vs[O+1] = (n+1.)/Gam1*( 
				 θs[1]*(2.*θs[3]-θs[2]) - θs[2]*θs[2]
			)/(θs[1]*θs[1]);
	if(maxPow >O+1) maxPow = O+1;
}
void PNPolytrope::getUSurface(    double *Us, int& maxPow){
	double x12 = Y[len-1][x]*Y[len-1][x];
	for(int a=0; a<maxPow; a++) Us[a] = 0.0;
	
	if(n==0){
		if(maxPow>=0) Us[0] = x12*(1.-2.*sigma*φs[0])*igs[0];
		if(maxPow>=1) Us[1] = x12*(1.-2.*sigma*φs[0])*(igs[1]-igs[0])
							 + sigma*x12*(3.*θs[1]-2.*φs[1])*igs[0];
		if(maxPow>=2) Us[2] = x12*(1.-2.*sigma*φs[0])*(igs[2]-igs[1])
							 + sigma*x12*(3.*θs[1]-2.*φs[1])*(igs[1]-igs[0])
							 + sigma*x12*(-2.5*x12)*igs[0];
	}
	else if(n==1){
		if(maxPow>=0) Us[0] = 0.0;
		if(maxPow>=1) Us[1] = x12*θs[1]*(1.-4.*sigma*φs[0])*igs[0];
		if(maxPow>=2) Us[2] = x12*θs[1]*(1.-4.*sigma*φs[0])*(igs[1]-igs[0])
							 +x12*θs[2]*(1.-4.*sigma*φs[0])*igs[0]
							 +sigma*x12*θs[1]*(3.*θs[1]-4.*φs[1])*igs[0];
	}
	else if(n==2){
		if(maxPow>=0) Us[0] = 0.0;
		if(maxPow>=1) Us[1] = 0.0;
		if(maxPow>=2) Us[2] = x12*igs[0]*θs[1]*θs[1]*(1.-6.*sigma*φs[0]);
	}
	if(maxPow > 2) maxPow = 2;
}
void PNPolytrope::getC1Surface(   double *cs, int& maxPow){
	double x1u1 = Y[len-1][x]*Y[len-1][u];
	if(maxPow>=0) cs[0]  = x1u1*(igs[0]);//1.0
	if(maxPow>=1) cs[1]  = x1u1*(igs[0]-igs[1]);//-3.0;
	if(maxPow>=2) cs[2]  = x1u1*(igs[1]-igs[2]);//3.0;
	if(maxPow> 2) maxPow = 2;
}

void PNPolytrope::getBetaSurface( double *bs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma : g);
	double Phi1 = (n+1.)*P0/rho0 * Y[len-1][phi];
	if(maxPow>=0) bs[0]  = 0.0;
	if(maxPow>=1) bs[1]  = Gam1*θs[1]/Phi1;
	if(maxPow>=2) bs[2]  = Gam1*θs[2]/Phi1 - sigma*Gam1*θs[1]*θs[1]/Phi1;
	//the pattern would continue, if needed...
	if(maxPow> 2) maxPow = 2;
}
void PNPolytrope::getPhiSurface(  double *ps, int& maxPow){
	if(maxPow>=0) ps[0]  = 1.;
	if(maxPow>=1) ps[1]  = φs[1]/φs[0];
	if(maxPow>=2) ps[2]  = φs[2]/φs[0];
	//pattern would continue for higher terms...
	if(maxPow> 2) maxPow = 2;
}

//method to print pertinent values of star to .txt, and plot them in gnuplot
void PNPolytrope::writeStar(char *c){
	//create names for files to be opened
	char pathname[256];
	char rootname[256];
	char txtname[256];
	char outname[256];
	if(c==NULL)	sprintf(pathname, "./out/%s", name);
	else{
		sprintf(pathname, "./%s/star/", c);
	}
	sprintf(txtname, "%s/%s.txt", pathname, name);
	sprintf(outname, "%s/%s.png", pathname, name);

	FILE *fp;
	if(!(fp = fopen(txtname, "w")) ){
		char command[256];
		sprintf(command, "mkdir -p %s", pathname);
		system(command);
		fp = fopen(txtname, "w");
	}
	
	
	//print results to text file
	// radius rho pressure gravity
	double irc=1./rho(0), ipc=1./P(0), R=Radius(), ig=1./dPhidr(length()-1), iPN=1./dPsidr(length()-1);
	double iq = 1./(dPhidr(length()-1)+dPsidr(length()-1)), MT = Mass();
	for(int X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			rad(X)/R, 
			rho(X)*irc, -drhodr(X)*irc*R, 
			P(X)*ipc, -dPdr(X)*ipc*R, 
			mr(X)/MT, 
			(dPhidr(X)+dPsidr(X))*iq,
			Phi(X));
	}
	fclose(fp);
	//plot file in png in gnuplot, and open png
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	char title[256]; graph_title(title);
	fprintf(gnuplot, "set title 'Profile for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'rho/rho_c, P/P_c, m/M, g/g_S'\n");
	fprintf(gnuplot, "plot ");
	fprintf(gnuplot, " '%s' u 1:2 w l t 'rho'", txtname);
	fprintf(gnuplot, ", '%s' u 1:4 w l t 'P'", txtname);
	fprintf(gnuplot, ", '%s' u 1:6 w l t 'm'", txtname);
	fprintf(gnuplot, ", '%s' u 1:7 w l t 'q'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s/%s.png'\n", pathname, "phi");
	fprintf(gnuplot, "set title 'Newtonian Potential for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'Phi'\n");
	fprintf(gnuplot, "plot ");
	fprintf(gnuplot, " '%s' u 1:8 w l t 'Phi'", txtname);
	fprintf(gnuplot, "\n");
	
	//print the pulsation coeffcients frequency
	sprintf(txtname, "%s/coefficients.txt", pathname);
	sprintf(outname, "%s/coefficients.png", pathname);
	fp  = fopen(txtname, "w");
	int maxpow=2;
	double A,U,V,C, read[2], x0 = Y[1][x]/Y[len-1][x];
	getAstarCenter(read, maxpow, 5./3.);
	A = read[0] + read[1]*x0*x0;
	getVgCenter(read, maxpow, 5./3.);
	V = read[0] + read[1]*x0*x0;
	getC1Center(read, maxpow);
	C = read[0] + read[1]*x0*x0;
	getUCenter(read, maxpow);
	U = read[0] + read[1]*x0*x0;
	fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
		x0, A, U, V, C);
	for(int X=1; X< length()-1; X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			Y[X][x]/Y[len-1][x],
			getAstar(X, 5./3.),
			getU(X),
			getVg(X, 5./3.),
			getC(X)
		);
	}
	double reads[3], t1 = 1.-Y[len-2][x]/Y[len-1][x];
	getAstarSurface(reads, maxpow, 5./3.);
	A = reads[0]/t1 + reads[1] + reads[2]*t1;
	getVgSurface(reads, maxpow, 5./3.);
	V = reads[0]/t1 + reads[1] + reads[2]*t1;
	getC1Surface(reads, maxpow);
	C = reads[0] + reads[1]*t1 + reads[2]*t1*t1;
	getUSurface(reads, maxpow);
	U = reads[0] + reads[1]*t1 + reads[2]*t1*t1;
	fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
		1.0, A, U, V, C);
	fclose(fp);	
	//plot file in png in gnuplot
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 800,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Pulsation Coefficients for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'A*, U, V_g, c_1'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'A*'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'U'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:4 w l t 'V_g'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:5 w l t 'c_1'", txtname);
	fprintf(gnuplot, "\n");	
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	
	printCoefficients(pathname);
}

void PNPolytrope::printCoefficients(char *c){
	char pathname[256];
	char rootname[256];
	char txtname[256];
	char outname[256];

	char command[300];
	sprintf(pathname, "%s/wave_coefficient", c);
	sprintf(command, "mkdir %s", pathname);
	system(command);
	
	char title[256]; graph_title(title);
	
	//print the coefficients of the center and surface, for series analysis
	sprintf(txtname, "%s/center.txt", pathname);
	FILE *fp = fopen(txtname, "w");
	double gam1 = 5./3.;
	int bc = 2;
	double A0[2], V0[2], c0[2], U0[2];
	getAstarCenter(A0, bc, gam1);
	getVgCenter(V0, bc, gam1);
	getUCenter(U0, bc);
	getC1Center(c0, bc);
	fprintf(fp, "A*  :\t%0.16le\t%0.16le\n", A0[0],A0[1]);
	fprintf(fp, "U   :\t%0.16le\t%0.16le\n", U0[0],U0[1]);
	fprintf(fp, "Vg  :\t%0.16le\t%0.16le\n", V0[0],V0[1]);
	fprintf(fp, "c1  :\t%0.16le\t%0.16le\n", c0[0],c0[1]);
	fclose(fp);
	sprintf(txtname, "%s/surface.txt", pathname);
	fp = fopen(txtname, "w");
	double A1[3], V1[3], c1[3], U1[3];
	getAstarSurface(A1, bc, gam1);
	getVgSurface(V1, bc, gam1);
	getUSurface(U1, bc);
	getC1Surface(c1, bc);
	fprintf(fp, "A*  :\t%0.16le\t%0.16le\t%0.16le\n", A1[0],A1[1],A1[2]);
	fprintf(fp, "U   :\t%0.16le\t%0.16le\t%0.16le\n", U1[0],U1[1],U1[2]);
	fprintf(fp, "Vg  :\t%0.16le\t%0.16le\t%0.16le\n", V1[0],V1[1],V1[2]);
	fprintf(fp, "c1  :\t%0.16le\t%0.16le\t%0.16le\n", c1[0],c1[1],c1[2]);
	fclose(fp);
	
	//print fits to those coefficients at center and surface
	int NC=15, NS=15;
	sprintf(txtname, "%s/centerfit.txt", pathname);
	fp = fopen(txtname, "w");
	double x2 = 0.0, xi1 = Y[len-1][x];
	fprintf(fp, "x           \tA*          \tA*_fit      \tU           \tU_fit       \tVg          \tVg_fit      \tc1          \tc1_fit\n");
	for(int X=0; X<NC; X++){
		x2 = Y[X][x]/xi1;
		fprintf(fp, "%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\n",
			x2,
			getAstar(X, gam1),
			A0[0] + A0[1]*x2*x2,
			getU(X),
			U0[0] + U0[1]*x2*x2,
			getVg(X, gam1),
			V0[0] + V0[1]*x2*x2,
			getC(X),
			c0[0] + c0[1]*x2*x2
		);
	}
	fclose(fp);
	sprintf(txtname, "%s/surfacefit.txt", pathname);
	fp = fopen(txtname, "w");
	double t = 1.;
	fprintf(fp, "x           \tA*          \tA*_fit      \tU           \tU_fit       \tVg          \tVg_fit      \tc1          \tc1_fit\n");
	for(int X=len-1; X>=len-NS-1; X--){
		t = 1. - Y[X][x]/xi1;
		fprintf(fp, "%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\n",
			Y[X][x]/xi1,
			getAstar(X, gam1),
			A1[0]/t + A1[1] + A1[2]*t,
			getU(X),
			U1[0] + U1[1]*t + U1[2]*t*t,
			getVg(X, gam1),
			V1[0]/t + V1[1] + V1[2]*t,
			getC(X),
			c1[0] + c1[1]*t + c1[2]*t*t
		);
	}
	fclose(fp);
	
	//print the pulsation coeffcients frequency
	sprintf(txtname, "%s/coefficients.txt", pathname);
	sprintf(outname, "%s/coefficients.png", pathname);
	fp  = fopen(txtname, "w");
	for(int X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			Y[X][x]/xi1,
			getAstar(X),
			getU(X),
			getVg(X),
			getC(X)
		);
	}
	fclose(fp);	
	//plot file in png in gnuplot;
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Pulsation Coefficients for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'A*, U, V_g, c_1'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 100\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'A*'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'U'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:4 w l t 'V_g'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:5 w l t 'c_1'", txtname);
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	//fits
	gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	sprintf(txtname, "%s/centerfit.txt", pathname);
	sprintf(outname, "%s/centerfit.png", pathname);
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Central Fitting by Power Series for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'difference'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 10\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xrange [0:%le]\n", Y[NC][x]/xi1);
	fprintf(gnuplot, "plot '%s' u 1:(abs($2-$3)/$2) w lp t 'A*'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($4-$5)/$4) w lp t 'U'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($6-$7)/$6) w lp t 'Vg'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($8-$9)/$8) w lp t 'c1'", txtname);
	fprintf(gnuplot, "\n");
	sprintf(txtname, "%s/surfacefit.txt", pathname);
	sprintf(outname, "%s/surfacefit.png", pathname);
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Surface Fitting by Power Series for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'difference'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 100\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xrange [%le:1]\n", Y[len-NS-1][x]/xi1);
	fprintf(gnuplot, "plot '%s' u 1:(abs($2-$3)/abs($2)) w lp t 'A*'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($4-$5)/abs($4)) w lp t 'U'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($6-$7)/abs($6)) w lp t 'Vg'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($8-$9)/abs($8)) w lp t 'c1'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
}

double PNPolytrope::SSR(){
	double checkEuler=0.0, checkPoiss=0.0, checkPsi=0.0;
		
	//sum up errors in equations
	double r = 0.0;
	double R = Radius();
	double d2Phi=0.0;
	double d2Psi=0.0;
	double e1, e2, e3, n1, n2, n3;
	char txtname[100], outname[100];
	for(int X=3; X<len-4; X++){
		r = rad(X);
		//Poisson equation for Phi
		double  h=rad(X+1)-rad(X);
		double  b3=dPhidr(X-3),
				b2=dPhidr(X-2),
				b1=dPhidr(X-1),
				a1=dPhidr(X+1),
				a2=dPhidr(X+2),
				a3=dPhidr(X+3);
		d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h);
		e1 = fabs( 4.*m_pi*rho(X)*r - 2.*dPhidr(X) - d2Phi*r );
		n1 = fabs( 4.*m_pi*rho(X)*r ) + fabs( 2.*dPhidr(X) ) + fabs(d2Phi*r);
		//Euler equation
		e2 = fabs( dPdr(X)
				+ dPhidr(X)*(rho(X)+P(X))
				+ rho(X)*dPsidr(X) );
		n2 = fabs( dPdr(X) )
				+ fabs(dPhidr(X)*(rho(X)+P(X)))
				+ fabs(rho(X)*dPsidr(X));
		// Equation for Psi
		b3=dPsidr(X-3);
		b2=dPsidr(X-2);
		b1=dPsidr(X-1);
		a1=dPsidr(X+1);
		a2=dPsidr(X+2);
		a3=dPsidr(X+3);
		d2Psi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h);
		e3 = fabs( 12.*m_pi*P(X) - 8.*m_pi*rho(X)*Phi(X) - 2.*dPsidr(X)/r - d2Psi ); 
		n3 = fabs( 12.*m_pi*P(X) )
				+ fabs( 8.*m_pi*rho(X)*Phi(X) )
				+ fabs(2.*dPsidr(X)/r) + fabs(d2Psi);
		e1 /= n1;
		e2 /= n2;
		e3 /= n3;
		checkPoiss += e1*e1;
		checkEuler += e2*e2;
		checkPsi   += e3*e3;	
	}
	return sqrt((checkPoiss + checkEuler + checkPsi)/double(3*len));
}

#endif
