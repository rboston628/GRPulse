//**************************************************************************************
//							CHANDRASEKHAR WHITE DWARF
// ChandrasekharWD++.h
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//		This model assumes T=0, and ignores Coulombic and other effects
//		The surface is not treated in any special way
//		Updated to include composition gradient, indicated by mu_e
// **************************************************************************************/

#ifndef GRChandrasekharWDCLASS
#define GRChandrasekharWDCLASS


#include "GRChandrasekharWD.h"

void GRChandrasekharWD::basic_setup(){
	//exclude unphysical values of Y0
	if(Y0 < 1.) Y0 = 1.;
	Y02 = Y0*Y0;
	X0  = sqrt(Y02-1.);
	X02 = Y02-1.;
	//name this model for files
	sprintf(name, "ChandrasekharWD.%1.1f", Y0);
	Rn = sqrt(2.*A0/(m_pi*Gee()))/B0;
	
	Y = new double*[len];
	for(int i=0; i<len; i++)
		Y[i] = new double[numvar];
}

void GRChandrasekharWD::init_arrays(){
	mass = new double[len];
	dydx = new double[len];
	mass[0] = 0.0;
	dydx[0] = 0.0;
	for(int X=1; X<len; X++){
		dydx[X] = deriv(sigma, Y[X]);
		mass[X] = -4.*m_pi*Y[X][v];
	}
}

//initalize white dwarf from central value of y and length
GRChandrasekharWD::GRChandrasekharWD(double Y0, int L, double mu0)
	: Y0(Y0), len(L), mu0(mu0), A0(Chandrasekhar::A0), B0(Chandrasekhar::B0*mu0)
{
	basic_setup();
	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with x[len-1]=0.0
	//we will find the proper dx with a bisection search

	//an initial guess -- based on constant volume star	
	dx = sqrt(6.0)/(len-1);
	double yS=1.0, ddx = dx;
	double dxmax=1.0, dxmin=0, ySmax=-1.0, ySmin=1.0;
	//this is used to control bisection search to guarantee it terminates
	double dxold = 1.0;
		
	//we first find reasonable values of the parameters sigma
	int tracker = 0;
	double s1 = sigma, ds = 0.01*sigma, dds=sigma;
	double target = zsurf, f1 = RK4integrate(s1, dx), f2 = RK4integrate(s1+ds, dx);
	printf("target = \t%0.16le\n", target);
	printf("result = \t%0.16le\n", f1);
	while(fabs(dds) > 1.e-16 && tracker<100){
		dds = ds*(target-f1)/(f2-f1);
		s1 = s1 + dds;
		ds = dds;
		f1 = RK4integrate(s1,   dx);
		f2 = RK4integrate(s1+ds,dx);
		printf("result = \t%0.26le\t %0.26le\t %0.26le\n", f1, s1, dds);
		tracker++;
	}
	sigma = s1;
	f1 = RK4integrate(sigma, dx);
	tracker=0;
	while(fabs(f1-target) > 1.e-16 && tracker<100){
		f1 = RK4integrate(sigma,dx);
		printf("result = \t%0.26le\t %0.26le\n", f1, sigma);
		tracker++;
	}
	printf("done\n");
	printf("target = \t%0.16le\n", target);
	printf("result = \t%0.16le\n", f1);	
	
		
	yS = RK4integrate(len, dx);
	//find brackets on dx that bound a zero in yS
	if(yS > 0){
		dxmin = dx; ySmin = yS;
		while(yS > 0 && !isnan(yS)){
			dx += ddx;
			yS = RK4integrate(len, dx);
			if(isnan(yS)){
				dx = 1.0/len; ddx *= 0.1; yS = 1.0;
			}
		}
		dxmax = dx; ySmax = yS;
	}
	else if (yS < 0){
		dxmax = dx; ySmax = yS;
		while(yS < 0){
			dx *= 0.5;
			yS = RK4integrate(len,dx);
		}
		dxmin = dx; ySmin = yS;
	}
	dx = 0.5*(dxmin+dxmax);
	yS = RK4integrate(len, dx);
	
	if(dx<0) {printf("somehow dx is negative...\n"); dx=-dx;}
	if(ySmin*ySmax > 0.0) {printf("big problem, chief\n"); exit(EXIT_FAILURE);}
	
	//now use bisection to find dx so that yS=0.0
	while( fabs(yS)>0.0 || isnan(yS) ){
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
	RK4integrate(len, dx, 1);
	
	printf("sigma = \t%le\n", sigma);
	printf("zsurf = \t%le %le\n", zsurf, setZsurf(Y[len-1]));
	zsurf = setZsurf(Y[len-1]);

	//now set physical properties of the white dwarf
	init_arrays();
	indexFit = len/2;
	for(int X=1; X<len; X++){
		//scan through x, set matching point where x[X] = 0.5
		if(x[X-1]>0.5 & x[X+1]<=0.5) indexFit = X;
	}
	indexFit /= 2;
	setupCenter();
	setupSurface();
	printf("%0.8lf\t%le\t%le\n", 1./Y02, mr(len-1)/MSOLAR, rad(len-1));
}


//initalize white dwarf from central value of y and length, with specific step size
//this should only be used for testing scaling issues
GRChandrasekharWD::GRChandrasekharWD(double Y0, int L, const double dx, double mu0)
	: Y0(Y0), len(L), dx(dx), mu0(mu0), A0(Chandrasekhar::A0), B0(Chandrasekhar::B0*mu0)
{	
	basic_setup();	
	//we know dx, so no need for bisection search
	RK4integrate(len, dx, 1);
	//continue integrating until the surface redshift is self-consistent
	int tracker = 0;
	while(fabs(zsurf-setZsurf(Y[len-1])) > 1e-16  & tracker<100){
		zsurf = setZsurf(Y[len-1]);
		RK4integrate(len, dx,1);
		tracker++;
	}

	//now set physical properties of the white dwarf
	init_arrays();
	indexFit = 512*round(double(len)/1024.0);
	indexFit /= 2;
	printf("  indexFit  = %d\n", indexFit);
	printf("r[indexFit] = %0.32le\n", rad(indexFit));
	setupCenter();
	setupSurface();
	printf("%0.2lf\t%le\t%le\n", 1./Y02, mr(len-1)/MSOLAR, rad(len-1));
}


//initialize a sphere of degenerate fermions using 
//  mu_e = mu
//  A0 = AN
//  B0 = BN
//  Use cases are:
//    1. testing results against tables with old values of A,B
//    2. implementing non-electron degenerate gasesous spheres (neutron stars)
//  This will use a constant value of mu_e throughout
//  When using a neutron star, mu_e will range from 2 (C12) to 1 (pure-neutron matter
GRChandrasekharWD::GRChandrasekharWD(double Y0, int L, double mu, double AN, double BN)
	: Y0(Y0), len(L), mu0(mu), A0(AN), B0(BN*mu)
{
	basic_setup();
	printf("A0: %le\n", A0);
	printf("B0: %le\n", B0);
	printf("Rn: %le\n", Rn);

	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with x[len-1]=0.0
	//we will find the proper dx with a bisection search

	//an initial guess -- based on constant volume star
	dx = sqrt(6.0)/(len-1);
	double yS=1.0, ddx = dx;
	double dxmax=1.0, dxmin=0, ySmax=-1.0, ySmin=1.0;
	//this is used to control bisection search to guarantee it terminates
	double dxold = 1.0;
	
	//we first find reasonable values of the parameters sigma
	int tracker = 0;
	double s1 = sigma, ds = 0.01*sigma, dds=sigma;
	double target = zsurf, f1 = RK4integrate(s1, dx), f2 = RK4integrate(s1+ds, dx);
	printf("target = \t%0.16le\n", target);
	printf("result = \t%0.16le\n", f1);
	while(fabs(dds) > 1.e-16 && tracker<100){
		dds = ds*(target-f1)/(f2-f1);
		s1 = s1 + dds;
		ds = dds;
		f1 = RK4integrate(s1,   dx);
		f2 = RK4integrate(s1+ds,dx);
		printf("result = \t%0.26le\t %0.26le\t %0.26le\n", f1, s1, dds);
		tracker++;
	}
	sigma = s1;
	f1 = RK4integrate(sigma, dx);
	tracker=0;
	while(fabs(f1-target) > 1.e-16 && tracker<100){
		f1 = RK4integrate(sigma,dx);
		printf("result = \t%0.26le\t %0.26le\n", f1, sigma);
		tracker++;
	}
	printf("done\n");
	printf("target = \t%0.16le\n", target);
	printf("result = \t%0.16le\n", f1);
	
	yS = RK4integrate(len, dx);	
	//find brackets on dx that bound a zero in yS
	if(yS > 0){
		dxmin = dx; ySmin = yS;
		while(yS > 0 && !isnan(yS)){
			dx += ddx;
			yS = RK4integrate(len, dx);
			if(isnan(yS)){
				dx = 1.0/len; ddx *= 0.1; yS = 1.0;
			}
		}
		dxmax = dx; ySmax = yS;
	}
	else if (yS < 0){
		dxmax = dx; ySmax = yS;
		while(yS < 0){
			dx *= 0.5;
			yS = RK4integrate(len,dx);
		}
		dxmin = dx; ySmin = yS;
	}
	dx = 0.5*(dxmin+dxmax);
	yS = RK4integrate(len, dx);
	
	//some possible errors that might occur
	if(dx<0)              {printf("somehow dx is negative...\n"); dx=-dx;}
	if(ySmin*ySmax > 0.0) {printf("big problem, chief\n"); exit(EXIT_FAILURE);}
	
	//now use bisection to find dx so that yS=0.0
	while( fabs(yS)>0.0 || isnan(yS) ){
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
	RK4integrate(len, dx, 1);
	printf("sigma = \t%le\n", sigma);
	printf("zsurf = \t%le %le\n", zsurf, setZsurf(Y[len-1]));
	zsurf = setZsurf(Y[len-1]);
		
	init_arrays();
	indexFit = len/2;
	for(int X=1; X<len; X++){
		//scan through x, set matching point where x[X] = 0.5
		if(Y[X-1][xi]>0.5 & Y[X+1][xi]<=0.5) indexFit = X;
	}
	indexFit /= 2;
	setupCenter();
	setupSurface();
	printf("%0.8lf\t%le\t%le\n", 1./Y02, mr(len-1)/MSOLAR, exp(-Y[len-1][la]) );
}

GRChandrasekharWD::~GRChandrasekharWD(){
	for(int i=0; i<len; i++)
		delete[] Y[i];
	delete[] dydx;
	delete[] mass;
}

//equation for dy/dx
//see Tooper & Chandrasekhar (1969) eq 7 -- solve for dy/dxi
double GRChandrasekharWD::deriv(double s, double yy[numvar]){
	if(yy[xi]==0.0) return 0.0;
	return -( 1.+s*(yy[y]-1.) )*std::exp(yy[la])*( yy[v] + 0.125*s*pow(yy[xi],3)*yy[f] )*pow(yy[xi],-2);
}

//see Tooper 1964, eq 3.7
double GRChandrasekharWD::setSigma(double ysurface[numvar]){
	// exp(la) = (1+zsurf)^2 = 1/[1 - 2GM/c^2R]
	// exp(la) = (1+zsurf)^2 = 1/[1 - 2(n+1)*sigma*v1/x1]
	// 1 - 2(n+1)sigma*v1/x1 = exp(-la) = (1+zsurf)^{-2}
	// sigma   = [1-exp(-la)] * x1/(2(n+1)v1) = [1 - (1+zsurf)^{-2}]*x1/[2(n+1)v1]
	return fabs( (1.-std::exp(-ysurface[la]))*ysurface[xi]/(2.*(n+1.)*ysurface[v]) );
	//return fabs( (1.-pow(1.+zsurf,-2))*ysurface[x]/(2.*(n+1.)*ysurface[v]) );
}

double GRChandrasekharWD::setZsurf(double ysurface[numvar]){
	// (1 + z)^2 = exp(la) = 1/[1 - 2Gm/c^2r]
	// zsurf = [1 - 2GM/c^2R]^(-1/2) - 1 
	// zsurf = exp(la/2) - 1
	return std::exp(0.5*ysurface[la]) - 1.;
}

double GRChandrasekharWD::centerInit(double s, double z, double ycenter[numvar]){
	ycenter[xi] = 0.0;
	ycenter[y]  = Y0;
	ycenter[v]  = 0.0;
	// Tooper 1964, eq 2.19 with 1-2GM/c^2R = (1+z)^{-2}
	ycenter[nu] = -2.*(n+1.)*std::log(1.+s) - 2.*std::log(1.+z);
	// Tooper 1964, eq 2.9, e^(-lambda) = 1 --> lambda = 0.0
	ycenter[la] = 0.0;
	// these variables can be determined analytically from y
	ycenter[x] = X0;
	ycenter[f] = Chandrasekhar::factor_f(X0);
	ycenter[g] = Chandrasekhar::factor_g(X0);
}

double GRChandrasekharWD::RK4step(double dx, double s, double yin[numvar], double yout[numvar]){
	double YC[numvar]={yin[xi],yin[y],yin[v],yin[nu],yin[la],0,yin[x],yin[f],yin[g]};
	double K[numrk][4];
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	std::complex<double> YCN(0.,0.);
	double twosig, expL;
	//calculate all Ks
	for(int a = 0; a<4; a++){
		//now from these, calculate next shift
		YCN = YC[y];
		YCN = std::pow(YCN,n);
		expL = std::exp(YC[la]);
		twosig = 2.*s*(n+1.)*YCN.real()*YC[xi];
		//K = dy = dx*(dy/dx)
		//Tooper & Chandrasekhar 1969 eq 7
		K[y][a]  = dx*deriv(s, YC);
		//Tooper & Chandrasekhar 1969 eq 8
		K[v][a]  = dx*( YC[xi]*YC[xi]*pow(YC[x],3) +  0.125*s*YC[g] ); 
		//Tooper 1964, eq 2.4
		K[la][a] = dx*( expL*(twosig - 1./YC[xi]) + 1./YC[xi] );
		twosig *= s*YC[y];
		//Tooper 1964, eq 2.3
		K[nu][a] = dx*( expL*(twosig + 1./YC[xi]) - 1./YC[xi] );
		//at center, must use analytic expansion
		if(YC[xi]==0) K[la][a] = K[nu][a] = K[y][a] = 0.0;
		//calculate intermediate positions using previous shift vectors
		YC[xi] = yin[xi] + B[a]*dx;   //radius xi
		for(int b=1; b<numrk; b++)
			YC[b] = yin[b] + B[a]*K[b][a];
		YC[x] = sqrt(YC[y]*YC[y] - 1.);
		YC[f] = Chandrasekhar::factor_f(YC[x]);
		YC[g] = Chandrasekhar::factor_g(YC[x]);
	}
	yout[xi] = yin[xi] + dx;
	for(int b=1; b<numrk; b++)
		yout[b] = yin[b] + K[b][0]/6.0 + K[b][1]/3.0 + K[b][2]/3.0 + K[b][3]/6.0;
	yout[x] = sqrt(yout[y]*yout[y] - 1.);
	yout[f] = Chandrasekhar::factor_f(yout[x]);
	yout[g] = Chandrasekhar::factor_g(yout[x]);
}


//integrate the WD using RK4 until the end is found
double GRChandrasekharWD::RK4integrate(double s, double& dx){
	double y1[numvar], y2[numvar];

	//set our initial conditions
	centerInit(s, zsurf, y2);
	
	//integrate with RK4
	int X=0;
	while(y2[y] > 0.0 && !isnan(y2[y]) ){
		for(int b=0; b<numvar; b++) y1[b] = y2[b];
		RK4step(dx, s, y1, y2);
		X++;
		//a possible error occurs when theta never reaches 0
		if(X>2*len) return setZsurf(y1);
	}
	//correct the value of sigma using the new surface values
	//if(s==sigma) sigma = setSigma(y1);
	//a possible error occurs if the second step is zero, making dx=0
	if(y1[x] != 0.0) dx = y1[x]/(len+1);
	//return guess for z in terms of field
	return setZsurf(y1);
}

//integrate the polytrope up to Len using RK4
double GRChandrasekharWD::RK4integrate(const int Len, double dx){

	//set our initial conditions
	centerInit(sigma, zsurf, Y[0]);
	
	//integrate with RK4
	for(int X = 0; X<Len-1; X++){
		RK4step(dx, sigma, Y[X], Y[X+1]);
		//sometimes, for noninteger n, if array is too large values y<0 lead to nans
		//if these occur, it is safe to terminate the integration
		if(Y[X+1][y] <0.0 || isnan(Y[X+1][y])){
			return Y[X+1][y]; //used in bisection search, dx too large should give negative
		}
	}
	//adjust value of sigma here
	return Y[Len-1][y];
}

int    GRChandrasekharWD::RK4integrate(const int Len, double dx, int grid){
	if(grid<1) grid=1;

	//set our initial conditions
	centerInit(sigma, zsurf, Y[0]);

	//integrate with RK4
	for(int X = 0; X<Len-1; X++){
		RK4step(dx, sigma, Y[X], Y[X+1]);
	}
	return Len;
}

//Here we define functions to access radius, pressure, etc.
double GRChandrasekharWD::rad(int X){
	return Rn*Y[X][xi];
}
double GRChandrasekharWD::rho(int X){
	return B0*pow(Y[X][x],3);
}
double GRChandrasekharWD::drhodr(int X){
	//note: dx/dxi = z*y/x
	return 3.*B0*Y[X][x]*Y[X][y]*dydx[X]/Rn;
}
double GRChandrasekharWD::P(int X){
	return A0*Y[X][f];
}
double GRChandrasekharWD::dPdr(int X){
	//see eqn (9) in Chandrasekhar 1935
	//note: dx/dxi = z*y/x
	return 8.*A0*pow(Y[X][x],3)*dydx[X]/Rn;
}
//*** THESE NEED SPECIAL CONSIDERATION ***
double GRPolytrope::Phi(int X){
	//return the Newtonian potential Gm/r
	return 8.*A0/B0*pow(Rn,2)*Y[X][v]/Y[X][x];
}
double GRPolytrope::dPhidr(int X){
	//return the Newtonian gravity Gm/r^2
	return 8.*A0/B0*Rn*Y[X][v]/Y[X][x]/Y[X][x];
}

double GRChandrasekharWD::mr(int X){
	return B0*pow(Rn,3)*mass[X];
}

//METRIC FIELDS
double GRPolytrope::Nu(int X){
	return Y[X][nu];
}
double GRPolytrope::dNudr(int X){
	//equation 2.19.5 in Tooper
	return -2.0*sigma*(n+1.)/(1.+sigma*Y[X][y])*dydx[X]/Rn;
}
double GRPolytrope::d2Nudr2(int X){
	//take derivative of equation 2.19.5 in Tooper
	double twonpsig = 2.*(n+1.)*sigma;
	double syp1 = sigma*Y[X][y] + 1.;
	return twonpsig*(sigma*dydx[X]*dydx[X]/syp1/syp1 
			- deriv(sigma, Y[X])/(1.+sigma*Y[X][y]))/Rn/Rn;
}

double GRPolytrope::Lambda(int X){
	return Y[X][la];
}
double GRPolytrope::dLambdadr(int X){
	double twonpsig = 2.*(n+1.)*sigma;
	return twonpsig*(pow(Y[X][x],3)*pow(Y[X][y],n)-Y[X][v])/(pow(Y[X][x],2)-twonpsig*Y[X][v]*Y[X][x])/Rn;
}
double GRPolytrope::d2Lambdadr2(int X){
	double twonpsig = 2.*(n+1.)*sigma;
	return 0.0;
}
//metric components
double GRPolytrope::gtt(int X){
	return -exp(Nu(X));
}
double GRPolytrope::dgttdr(int X){
	return dNudr(X)*gtt(X);
}
double GRPolytrope::gtr(int X){
	return 0.0;
}
double GRPolytrope::dgtrdr(int X){
	return 0.0;
}
double GRPolytrope::grr(int X){
	return exp(Lambda(X));
}
double GRPolytrope::dgrrdr(int X){
	return dLambdadr(X)*grr(X);
}
double GRPolytrope::gθθ(int X){
	return pow(rad(X),2);
}
double GRPolytrope::dgθθdr(int X){
	return 2.*rad(X);
}

double GRChandrasekharWD::Schwarzschild_A(int X, double GamPert){
	if(GamPert==0.0) return  0.0;
	else        	 return -dydx[X]*(8.*pow(Y[X][x],3)/Y[X][f]/GamPert - 3.*Y[X][y]/pow(Y[X][x],2))/Rn;
}

double GRChandrasekharWD::getAstar(int X, double GamPert){
	if(GamPert==0.0) return  0.0;
	else        	 return  dydx[X]*(8.*pow(Y[X][x],3)/Y[X][f]/GamPert - 3.*Y[X][y]/pow(Y[X][x],2))*Y[X][xi];
}

double GRChandrasekharWD::getU(int X){
	if(X==0) return 3.;
	return pow(Y[X][xi],3)*pow(Y[X][y],n)/Y[X][v];
}

double GRChandrasekharWD::getVg(int X, double GamPert){
	if(GamPert==0.0) return -8.*pow(Y[X][x],3)*Y[X][xi]*dydx[X]/Y[X][f]/Gamma1(X);
	else			 return -8.*pow(Y[X][x],3)*Y[X][xi]*dydx[X]/Y[X][f]/GamPert;
}

double GRChandrasekharWD::getC(int X){
	if(X==0) return z[len-1]*xi[len-1]/2./yc[1];//-3.*z[len-1]/xi[len-1]/X02/X0;
	return (z[len-1]/z[X])*(xi[X]/xi[len-1])*(mue[X]/mue[len-1]);
}

double GRChandrasekharWD::Gamma1(int X){
	//see Tooper & Chandrasekhar (1969) eqn 15
	return 8.*pow(x[X],5)/(3.*y[X]*f[X]);
}

double GRChandrasekharWD::sound_speed2(int X, double GamPert){
	if(GamPert == 0.0) return Gamma1(X)*A0/B0*f[X]/pow(x[X],3)/mue[X];
	else               return GamPert  *A0/B0*f[X]/pow(x[X],3)/mue[X];
}


double GRChandrasekharWD::Radius(){return Rn*xi[len-1];}	//total radius
double GRChandrasekharWD::Mass(){return mr(len-1);}//total mass
double GRChandrasekharWD::Gee(){return G_CGS;};
//in Newtonian, light speed is infinity...
double GRChandrasekharWD::light_speed2(){return C_CGS*C_CGS};



// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
void GRChandrasekharWD::setupCenter(){	
	//the central coefficients -- expanded in terms of r/R
	double z1 = Y[len-1][z];
	double x1 = Y[len-1][xi];
	//
	yc[0]=	Y0;
	yc[1]=	-X02*X0/6.*pow(x1,2);
	yc[2]=	X02*X02*Y0/40.*pow(x1,4);
	yc[3]=	X0*(-14.*X02*X02-19.*X02*X02*X02)/5040.*pow(x1,6);
	//
	xc[0]=	X0;
	xc[1]=	Y0*yc[1]/X0;
	xc[2]=  yc[2]*pow(yc[0],3)/pow(xc[0],3) - yc[0]*yc[2]/xc[0] - 0.5*yc[1]*yc[1]/xc[0];
	//
	fc[0]=	Y0*X0*(2.*X0-3.) + 3.*asinh(X0);
	fc[1]=	2.*xc[1]*(3.*xc[0]*xc[0]*(xc[0]-1.)+2.*xc[0]-1.);
}

void GRChandrasekharWD::getAstarCenter(double *Ac, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	//depending on power requested, return appropriate number of terms
	if(maxPow>=0) Ac[0] = 0.0;
	if(maxPow>=2) Ac[1] = 16.*X0*X02*yc[1]/Gam1/fc[0] - 6.*x[1]/X0;
	if(maxPow>=4) Ac[2] = 
		-16.*fc[1]*X02*X0*yc[1]/Gam1/fc[0]/fc[0]
			+ 48.*X02*xc[1]*yc[1]/Gam1/fc[0]
			- 12.*xc[2]/X02
			+ 6.*xc[1]*xc[1]/X02
			+ 32.*X02*X0*yc[2]/Gam1/fc[0];		
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}

void GRChandrasekharWD::getVgCenter(double *Vc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);		
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] = -16.*yc[1]*X02*X0/Gam1/fc[0];
	if(maxPow>=4) Vc[2] = 
		-16.*(2.*fc[0]*X02*X0*yc[2]
			+3.*fc[0]*X02*xc[1]*yc[1]
			-fc[1]*X02*X0*yc[1]
		)/Gam1/fc[0]/fc[0];
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4; 
}

void GRChandrasekharWD::getUCenter(double *Uc, int& maxPow){
	double x1 = xi[len-1];
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] = 9.0*xc[1]/X0 + 0.9*pow(x1,2)*Y0*X0;
	if(maxPow>=4) Uc[2] = 
		9.0*xc[1]*xc[1]/X0/X0
		+9.*xc[2]/X0+2.7*x1*x1*Y0*xc[1]
		+3./25.*X02*pow(x1*x1,2)
		+93./1400.*pow(x1,4)*X02*X02;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}

void GRChandrasekharWD::getC1Center(double *cc, int& maxPow){
	double x1 = xi[len-1];
	double z1 = z[len-1];
	if(maxPow>=0) cc[0] = z1*x1/2./yc[1];//-3.*z1/x1/X02/X0;//
	if(maxPow>=2) cc[1] = -z1*x1*yc[2]/yc[1]/yc[1];//-9.*z1*x1*Y0/10./X02;//
	if(maxPow>=4) cc[2] = -z1*x1*(1.5*yc[1]*yc[3]-2.*yc[2]*yc[2])/pow(yc[1],3);//-3.*z1*pow(x1,3)*(31.*Y02+25.)/1400./X0;//
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
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
void GRChandrasekharWD::setupSurface(){}

void GRChandrasekharWD::getAstarSurface(double *As, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double x1 = xi[len-1];
	double a1 = -dydx[len-1]*x1;
	int O=1;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=-1) As[O-1] = 1.5;
	if(maxPow>= 0) As[O  ] = 0.75*a1;
	if(maxPow>= 1) As[O+1] = 0.75*a1 + 8.*a1*a1/Gam1 - 3.*a1*a1/8.;
	if(maxPow>= 2) As[O+2] = a1*3./16.*pow(a1-2.,2)+4./3.*a1*a1*(12.+5.*a1)/Gam1;
	if(maxPow>= 3) As[O+3] = -3.*a1*pow(a1-2,3)/32.+(724.*a1*a1/45. + 20.*a1 + 24.)*a1*a1/Gam1;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 3) maxPow = O+3;
}

void GRChandrasekharWD::getVgSurface(double *Vs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double x1 = xi[len-1];
	double a1 =-x1*z[len-1];
	int O=1;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=-1) Vs[O-1] = 0.0;
	if(maxPow>= 0) Vs[O  ] = 0.0;
	if(maxPow>= 1) Vs[O+1] = -8.*a1*a1/Gam1;
	if(maxPow>= 2) Vs[O+2] = -4./3.*a1*a1*(12.+5.*a1)/Gam1;
	if(maxPow>= 3) Vs[O+3] = -(724.*a1*a1/45. + 20.*a1 + 24.)*a1*a1/Gam1;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 3) maxPow = O+3;
}

void GRChandrasekharWD::getUSurface(double *Us, int& maxPow){
	//these happen to all be zero at surface
	for(int j=0; j<=maxPow; j++){
		Us[j] = 0.0;
	}
}

void GRChandrasekharWD::getC1Surface(double *cs, int& maxPow){
	if(maxPow >= 0) cs[0] =  1.;
	if(maxPow >= 1) cs[1] = -3.;
	if(maxPow >= 2) cs[2] =  3.;
	if(maxPow >= 3) cs[3] = -1.;
	if(maxPow >= 4) cs[4] =  0.; //does not actually appear in equations
	//if more  terms than this requested, cap number of terms
	if(maxPow  > 4) maxPow = 4;
}


void GRChandrasekharWD::writeStar(char *c){
	//create names for files to be opened
	char filename[256];
	char rootname[256];
	char txtname[256];
	char outname[256];
	if(c==NULL)	sprintf(filename, "./out/%s", name);
	else{
		sprintf(filename, "./%s/star", c);
	}
	sprintf(txtname, "%s/%s.txt", filename, name);
	sprintf(outname, "%s/%s.png", filename, name);

	FILE *fp;
	if(!(fp = fopen(txtname, "w")) ){
		if(c==NULL){
			system("mkdir ./out");
			char command[256]; sprintf(command, "mkdir ./out/%s", name);
			system(command);
			fp = fopen(txtname, "w");
		}
		else {
			char command[256]; sprintf(command, "mkdir ./%s", c);
			system(command);
			sprintf(command, "mkdir %s", filename);
			system(command);
			if(!(fp = fopen(filename, "w"))) printf("big trouble, boss\n");		
		}
	}
	//print results to text file
	// radius rho pressure gravity
	double irc=1./rho(0), ipc=1./P(0), R=Radius(), ig=1./dPhidr(length()-1);
	for(int X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t",
			rad(X)/R, rho(X)*irc, -drhodr(X)*irc*R,
			P(X)*ipc, -dPdr(X)*ipc*R,
			mr(X)/Mass(), dPhidr(X)*ig);
		fprintf(fp, "%le\t%le\t%le", x[X], y[X], f[X]);
		fprintf(fp, "\n");
		//fflush(fp);
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
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'rho'", txtname);
	//fprintf(gnuplot, ", '%s' u 1:3 w l t '|drho/dr|'", txtname);
	fprintf(gnuplot, ", '%s' u 1:4 w l t 'P'", txtname);
	//fprintf(gnuplot, ", '%s' u 1:5 w l t '|dP/dr|'", txtname);
	fprintf(gnuplot, ", '%s' u 1:6 w l t 'm'", txtname);
	fprintf(gnuplot, ", '%s' u 1:7 w l t 'g'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	
	//print the pulsation coeffcients frequency
	sprintf(txtname, "%s/coefficients.txt", filename);
	sprintf(outname, "%s/coefficients.png", filename);
	fp  = fopen(txtname, "w");
	for(int X=1; X< length()-1; X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			xi[X],
			-xi[X]*Rn*Schwarzschild_A(X, 0.),
			getU(X),
			getVg(X, 0.),
			getC(X)
		);
	}
	fclose(fp);	
	//plot file in png in gnuplot
	gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 800,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Pulsation Coefficients for %s'\n", title);
	//fprintf(gnuplot, "set xlabel 'log_{10} r/R'\n");
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
	
	
	//print the Brunt-Vaisala frequency
	sprintf(txtname, "%s/BruntVaisala.txt", filename);
	sprintf(outname, "%s/BruntVaisala.png", filename);
	fp  = fopen(txtname, "w");
	double N2 = -1.0;
	for(int X=1; X< length()-1; X++){
		//N^2 = -g*A
		N2 = -dPhidr(X)*Schwarzschild_A(X,0.);
		//if(N2<0.0) N2 = 0.0;	//show a drop for negative N^2
		//else N2 = log10(N2);
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\n",
			//1./(1.-mr(X)/Mass()),
			P(X),
			//xi[X]/xi[len-1],
			N2,
			2.*sound_speed2(X,0.)*pow(Rn*xi[X],-2));
	}
	fclose(fp);	
	//plot file in png in gnuplot
	gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 800,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Brunt-Vaisala for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'log_{10} r/R'\n");
	fprintf(gnuplot, "set logscale x 10\n");
	fprintf(gnuplot, "set format x '10^{%%L}'\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xlabel 'log_{10} P\n");
	fprintf(gnuplot, "set ylabel 'log_{10} N^2 & log_{10} L_1^2 (Hz^2)\n");\
	fprintf(gnuplot, "set yrange [1e-6:1e2]\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'N^2'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'L_1^2'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	
	//print the chemical gradient
	sprintf(txtname, "%s/chemical.txt", filename);
	sprintf(outname, "%s/chemical.png", filename);
	fp  = fopen(txtname, "w");
	double H,He;
	for(int X=1; X< length()-1; X++){
		He = 2.*(mue[X]-1.)/mue[X];
		H = 1.-He;
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\n",
			1./(1.-mr(X)/Mass()),
			H,
			He
		);
	}
	fclose(fp);	
	//plot file in png in gnuplot
	gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 800,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Chemical composition for %s'\n", title);
	fprintf(gnuplot, "set logscale x 10\n");
	fprintf(gnuplot, "set format x '10^{%%L}'\n");
	fprintf(gnuplot, "set xlabel '-log_{10} (1-m/M)\n");
	fprintf(gnuplot, "set ylabel 'abundance\n");
	fprintf(gnuplot, "set yrange [-0.01: 1.01]\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'X (H1)'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'Y (He4)'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);	
}

#endif