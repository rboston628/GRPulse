// *************************************************************************************
//					GENERAL RELATIVISTIC POLYTROPIC STAR
// GRPolytrope.cpp
//		A simple stellar object in GR obeying everywhere a polytropic equation of state
//				P~rho^Gamma
//		Here rho = total mass-energy density, following Tooper, 1964
//		Using metric form
//			ds^2 = -exp(nu) c^2dt^2 + exp(lambda) dr^2 + r^2 dOmega^2
//		Solves Relativivistic Lane-Emden equation, in eqn 2.25, 2.26 in 
//				Tooper (1964), ApJ 140, 434T
//		y=theta, x=xi, v=v (a unitless mass) in Tooper's notation
// *************************************************************************************

#ifndef GRPOLYTROPECLASS
#define GRPOLYTROPECLASS

#include "GRPolytrope.h"

//initalize polytrope from index and length
GRPolytrope::GRPolytrope(double n, double zsurf, int L)
	: n(n), zsurf(zsurf), len(L), Gamma(1.0+1.0/n), 
			sigma( (1.-pow(1.+zsurf,-2))/(2.*(n+1)) )
{
	printf("Making GR polytrope N=%1.1f, z=%le\n", n,zsurf);
	//if index is out of range, fail
	if(n >= 5 || n < 0){
		n = nan(""); len = 0;
		printf("Invalid polytropic index.  Polytrope not initialized.\n");
		exit(EXIT_FAILURE);
	}
	if(zsurf >= 0.5){
		printf("\tYour specified star would collapse to a black hole.\n");
		zsurf = 0.499;
	}
	if(sigma> n/(n+1.) && n!=0){
		printf("Overlarge value for sigma, resorting to max value.\n");
		sigma = n/(n+1.);
	}
	if(n==0) Gamma = 0.0;
	//name this polytrope for files
	sprintf(name, "GRpolytrope%1.1f_%1.1f", n, sigma);
	
	Y = new double*[len];
	for(int i=0;i<len;i++)
		Y[i] = new double[numvar];	
	
	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with y[len-1]=0.0
	//we will find the proper dx with a bisection search
	
	//an initial guess based on uniform density star -- slightly larger
	double dx = 2.*sqrt(6.0)/(len-1);
	//run one test-run to adjust dx
	RK4integrate(sigma, dx);
	
	//use a 2D Newton's method to find sigma, dx simultaneously
	size_t const np=2;
	double z1[np] = {sigma, dx};
	double dz[np] = {0.1*sigma, dx};
	double target[np] = {zsurf, 0.0};
	// require the zsurf to match at surface, and theta to vanish there
	std::function<void(double[np],double[np])> match_zsurf = [this](double f[np], double x[np]){
		this->RK4integrate(x[0], x[1]);
		f[0] = std::exp(0.5*this->Y[len-1][la]) - 1.;
		f[1] = this->Y[len-1][y];
	};
	// both sigma and dx are limited to be positive
	std::function<bool(double[np])> sign_limit = [](double z[np])->bool{
		return (z[0]>0.0)&&(z[1]>0.0);
	};
	//now run the newtonian searc
	newton_search<np>(match_zsurf, target, z1, dz, 1.e-13, 1000, sign_limit);
	sigma = z1[0];
	dx = z1[1];	
	
	RK4integrate(len, dx, 1);
	printf("sigma = \t%le\n", sigma);
	printf("zsurf = \t%le %le\n", zsurf, setZsurf(Y[len-1]));
	printf("edge  = \t%le\n", Y[len-1][y]);
		
	//set initial density, pressure
	//if we use units c2=1, then P0 = sigma
	rho0 = 1.0;
	P0 = sigma;
	if(sigma==0) P0 = 1.0;
	// 1/A as defined by Tooper (1964)
	Rn = sqrt( (n+1.)*P0/(4.*m_pi*rho0*rho0) );
	
	//now set physical properties of the polytrope
	dydx = new double[len];
	base = new double[len];
	dydx[0] = 0.0;
	base[0] = 1.0;
	//
	indexFit = len/2;
	for(int X=1; X<len-1; X++){
		dydx[X] = deriv(sigma, Y[X]);
		//as we scan through x,y,z, set matching point where y[X] = 0.5
		if(Y[X-1][y]>0.5 & Y[X+1][y]<=0.5) indexFit = X;
		//the array base saves some calclations
		base[X] = pow(Y[X][y], n-1.0);	
	}
	dydx[len-1] = deriv(sigma, Y[len-1]);
	base[len-1] = 0.0;
	indexFit /= 2;
	setupCenter();
	setupSurface();
}

//initalize polytrope from index, length, and dx
GRPolytrope::GRPolytrope(double n, double sigma, int L, const double dx)
	: n(n), sigma(sigma), len(L), Gamma(1.0+1.0/n),
			zsurf( 1./sqrt(1.-2.*(n+1.)*sigma) - 1. )
{
	//if index is out of range, fail
	if( n >= 5.0 || n < 0.0 ){
		n = nan(""); len = 0;
		printf("\nInvalid polytropic index.  Polytrope not initialized.\n");
		exit(EXIT_FAILURE);
	}
	if(sigma> n/(n+1.) && n!=0){
		printf("Overlarge value for sigma, resorting to max value.\n");
		sigma = n/(n+1.);
	}
	//name this polytrope for files
	sprintf(name, "GRpolytrope%1.1f_%1.1f", n, sigma);

	//initialize arrays
	Y = new double*[len];
	for(int i=0;i<len;i++) Y[i] = new double[numvar];
	
	RK4integrate(len, dx, 1);
	//continue integrating until the surface redshift is self-consistent
	int tracker = 0;
	while(fabs(zsurf-setZsurf(Y[len-1])) > 1e-16  & tracker<100){
		zsurf = setZsurf(Y[len-1]);
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
		dydx[X] = deriv(sigma, Y[X]);
		//the array base saves some calclations
		base[X] = pow(Y[X][y], n-1.0);
	}
	indexFit = 512*round(double(len)/1024.0);
	indexFit /= 2;
	printf("  indexFit  = %d\n", indexFit);
	printf("r[indexFit] = %0.32le\n", rad(indexFit));
	printf("zsurf = %le\t sigma=%le\n", zsurf,sigma);
	
	setupCenter();
	setupSurface();	
}


GRPolytrope::~GRPolytrope(){
	for(int i=0; i<len; i++)
		delete[] Y[i];
	delete[] base;
	delete[] dydx;
}


double GRPolytrope::deriv(double s, double yy[numvar]){
	//equation for dy/dx
	//see Tooper 1964, eq 2.25
	if(yy[x]==0.0) return 0.0;
	std::complex<double> yn = pow(yy[y],n+1);
	return -( 1.+s*yy[y] )*std::exp(yy[la])*( yy[v] + s*pow(yy[x],3)*yn.real() )*pow(yy[x],-2);
}

//see Tooper 1964, eq 3.7
double GRPolytrope::setSigma(double ysurface[numvar]){
	// exp(la) = (1+zsurf)^2 = 1/[1 - 2GM/c^2R]
	// exp(la) = (1+zsurf)^2 = 1/[1 - 2(n+1)*sigma*v1/x1]
	// 1 - 2(n+1)sigma*v1/x1 = exp(-la) = (1+zsurf)^{-2}
	// sigma   = [1-exp(-la)] * x1/(2(n+1)v1) = [1 - (1+zsurf)^{-2}]*x1/[2(n+1)v1]
	return fabs( (1.-std::exp(-ysurface[la]))*ysurface[x]/(2.*(n+1.)*ysurface[v]) );
}

double GRPolytrope::setZsurf(double ysurface[numvar]){
	// (1 + z)^2 = exp(la) = 1/[1 - 2Gm/c^2r]
	// zsurf = [1 - 2GM/c^2R]^(-1/2) - 1 
	// zsurf = exp(la/2) - 1
	return std::exp(0.5*ysurface[la]) - 1.;
}

void GRPolytrope::centerInit(double s, double z, double ycenter[numvar]){
	ycenter[x] = 0.0;
	ycenter[y] = 1.0;
	ycenter[v] = 0.0;
	// Tooper 1964, eq 2.19 with 1-2GM/c^2R = (1+z)^{-2}
	ycenter[nu] = -2.*(n+1.)*std::log(1.+s) - 2.*std::log(1.+z);
	// Tooper 1964, eq 2.9, e^(-lambda) = 1 --> lambda = 0.0
	ycenter[la] = 0.0; 
}

void GRPolytrope::RK4step(double dx, double s, double yin[numvar], double yout[numvar]){
	double YC[numvar]={yin[x],yin[y],yin[v],yin[nu],yin[la]};
	double K[numvar][4];
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	std::complex<double> YCN(0.,0.);
	double twosig, expL;
	//calculate all Ks
	for(int a = 0; a<4; a++){
		//now from these, calculate next shift
		YCN = YC[y];
		YCN = std::pow(YCN,n);
		expL = std::exp(YC[la]);
		twosig = 2.*s*(n+1.)*YCN.real()*YC[x];
		//K = dy = dx*(dy/dx)
		K[v][a]  = dx*( YC[x]*YC[x]*YCN.real() ); //Tooper 1964, eq 2.27
		K[la][a] = dx*( expL*(twosig - 1./YC[x]) + 1./YC[x] ); //Tooper 1964, eq 2.4
		twosig *= s*YC[y];
		K[nu][a] = dx*( expL*(twosig + 1./YC[x]) - 1./YC[x] ); //Tooper 1964, eq 2.3
		K[y][a]  = -(1.+s*YC[y])/s/(n+1.)/2.*K[nu][a]; //Tooper 1964, eq. 2.14
		if(YC[x]==0) K[la][a] = K[nu][a] = K[y][a] = 0.0; //at center, must use analytic expansion
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
double GRPolytrope::RK4integrate(double s, double& dx){
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
		if(X>2+len) break; 
	}
	//a possible error occurs if the second step is zero, making dx=0
	if(y1[x] != 0.0) dx = y1[x]/(len+1);
	//return guess for z in terms of field
	for(int b=0; b<numvar; b++){
		Y[len-1][b] = y1[b];
	}
	return setZsurf(y1);
}

//integrate the polytrope up to Len using RK4
double GRPolytrope::RK4integrate(const int Len, double dx){

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

int GRPolytrope::RK4integrate(const int Len, double dx, int grid){
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
double GRPolytrope::rad(int X){
	return (Rn*Y[X][x]);
}
double GRPolytrope::rho(int X){
	if(n==0) return rho0;
	else     return rho0*Y[X][y]*base[X];//pow(y[X], n);
}
double GRPolytrope::drhodr(int X){
	if(n==0) return 0.0;
	else     return n*rho0*base[X]*(dydx[X]/Rn);//n*pow(y[X], n-1)*(z[X]/Rn);
}
double GRPolytrope::P(int X){
	if(n==0) return P0*Y[X][y];
	else     return P0*base[X]*Y[X][y]*Y[X][y];//pow(y[X], n+1);
}
double GRPolytrope::dPdr(int X){
	if(n==0) return P0*dydx[X]/Rn;
	else     return (n+1)*P0*base[X]*Y[X][y]*(dydx[X]/Rn);//(n+1)*pow(y[X], n)*(z[X]/Rn);
}
inline double GRPolytrope::mr(int X){
	return 4.*m_pi*rho0*pow(Rn,3)*Y[X][v];
}

//*** THESE NEED SPECIAL CONSIDERATION ***
double GRPolytrope::Phi(int X){
	//return the Newtonian potential Gm/r
	return 4.*m_pi*rho0*pow(Rn,2)*Y[X][v]/Y[X][x];
}
double GRPolytrope::dPhidr(int X){
	//return the Newtonian gravity Gm/r^2
	return 4.*m_pi*rho0*Rn*Y[X][v]/Y[X][x]/Y[X][x];
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
	return twonpsig*(pow(Y[X][x],3)*pow(Y[X][y],n)-Y[X][v])/pow(Y[X][x],2)*std::exp(Y[X][la])/Rn;
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

double GRPolytrope::sound_speed2(int X, double GamPert){
	//see Thorne RSSD eq 2.10
	// vs2 = Gamma1 * P/(rho+P/c^2) = Gamma1* P0 y^n+1 /(rho0 y^n + P0 y^n+1/c2)
	// vs2 = Gamma1*P0*y/(rho0 + P0 y/c^2)
	if(GamPert==0.0) return Gamma  *P0/rho0*Y[X][y]/(1.0+sigma*Y[X][y]);
	else             return GamPert*P0/rho0*Y[X][y]/(1.0+sigma*Y[X][y]);
}

double GRPolytrope::Schwarzschild_A(int X, double GamPert){
	if(GamPert==0.0) return dydx[X]/Y[X][y]/Rn* (n/(1.+sigma*Y[X][y])-(n+1.)/Gamma);
	else             return dydx[X]/Y[X][y]/Rn* (n/(1.+sigma*Y[X][y])-(n+1.)/GamPert);
}

double GRPolytrope::getAstar(int X, double GamPert){
	if(GamPert==0.0) return Y[X][x]*dydx[X]/Y[X][y]* ((n+1.)/Gamma   - n/(1.+sigma*Y[X][y]));
	else             return Y[X][x]*dydx[X]/Y[X][y]* ((n+1.)/GamPert - n/(1.+sigma*Y[X][y]));
}

double GRPolytrope::getU(int X){
	if(X==0) return 3.;
	return pow(Y[X][x],3)*pow(Y[X][y],n)/Y[X][v];
}
double GRPolytrope::getVg(int X, double GamPert){
	if(GamPert==0.0) return -Y[X][x]*dydx[X]/Y[X][y]* (n+1.)/Gamma;
	else			 return -Y[X][x]*dydx[X]/Y[X][y]* (n+1.)/GamPert;
}
double GRPolytrope::getC(int X){
	double vxi3 = Y[len-1][v]*pow(Y[len-1][x],-3);
	if(X==0) return 3.*vxi3*std::exp(-Y[0][nu])/(1.+3.*sigma);
	return -(1. + sigma*Y[X][y])*(vxi3)*(Y[X][x]/dydx[X])*std::exp(-Y[X][nu]);
}

double GRPolytrope::Gamma1(int X){
	return Gamma;
}


//integrate eq 3.8 of Tooper (1964) using trapezoidal rule
double GRPolytrope::getXproper(int k){
	double xbar = 0.0;
	double twosign = 2.0*sigma*(n+1.);
	double dx, int1 = 1.0, int2 = 1.0; //in lim x->0, z/x->0
	for(int i=0; i<k-1; i++){
		dx = Y[i+1][x] - Y[i][x];
		int1 = int2;
		int2 = std::exp(0.5*Y[i+1][la]);
		xbar += 0.5*dx*(int2+int1);
	}
	return xbar;
}

double GRPolytrope::Radius(){return rad(len-1);} //total radius
double GRPolytrope::Mass(){return mr(len-1);}    //total mass

// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
void GRPolytrope::setupCenter(){
	double xi1 = Y[len-1][x];
	double xi2 = xi1*xi1;
	double s=sigma, s2=sigma*sigma;
	double sigma_factor = (1.+sigma)*(1.+3.*sigma);
	//expansion coefficients of theta
	thc[0] = 1.0;
	thc[2] = -xi2/6.0*sigma_factor;
	thc[4] = xi2*xi2/360.0*sigma_factor*(30.0*s2 + n*(3.-2.*s+15.*s2));
	
	//expansion coefficients of lambda
	lc[0] = 0.0;
	lc[2] = 2.*(n+1.)/3.*xi2*s;
	lc[4] = -(n+1.)*xi2*xi2/45.*s*(-10.*s + n*(3. + 2.*s + 9.*s2));
	
	//expansion coefficients of nu
	nuc[0] = Y[0][nu];
	nuc[2] = xi2*(n + 1.)/3.*s*(1. + 3.*s);
	nuc[4] = -xi2*xi2*(n + 1.)/180.*s*(1. + 3.*s)*(5.*s*(3.*s - 1.) + n*(3. - 2.*s + 15.*s2));
	nuc[6] = xi2*xi2*xi2*sigma*(1. + 3.*sigma)/22680.*(
		70.*s2*(4. - 3.*s + 9.*s2) + 2.*pow(n,3)*(12. - 45.*s + 65.*s2  - 123.*s*s2 + 315.*s2*s2)
	+   3.*n*n*(3. - 111.*s + 127.*s2 - 249.*s*s2 + 630.*s2*s2)
	+   3.*n*( -5. - 81.*s  + 177.*s2 - 237.*s*s2 + 630.*s2*s2)
	);

}
void GRPolytrope::getAstarCenter(double *AC, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	if(maxPow>=0) AC[0] = 0.0;
	if(maxPow>=2) AC[1] = 2.*( (n+1.)/Gam1 - n/(1.+sigma) )*thc[2];
	if(maxPow>=4) AC[2] = 2.*( (n+1.)/Gam1 - n/(1.+sigma) )*(4.*thc[4] - 2.*thc[2]*thc[2]) 
						+ 2.*n*sigma*thc[2]*thc[2]/(1.+sigma)/(1.+sigma); 
	if(maxPow> 4) maxPow = 4;
}
void GRPolytrope::getVgCenter(double *Vc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] = -2.*(n+1.)* thc[2]/Gam1;
	if(maxPow>=4) Vc[2] =  2.*(n+1.)*(thc[2]*thc[2]-2.*thc[4])/Gam1;
	if(maxPow> 4) maxPow = 4; 
}
void GRPolytrope::getC1Center(double *cc, int& maxPow){
	if(maxPow>=0) cc[0] = getC(0);
	if(maxPow>=2) cc[1] = -3.*Y[len-1][v]/Y[len-1][x]*(1.+sigma)*(5.*sigma + n*(5.*sigma - 1.))/(10.+30.*sigma)*std::exp(-Y[0][nu]);
	if(maxPow>=4) cc[2] = 0.0;
	if(maxPow> 4) maxPow = 4;
}
void GRPolytrope::getUCenter(double *Uc, int& maxPow){
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] = 4.*nuc[4]/nuc[2];
	if(maxPow>=4) Uc[2] = 4.*(3.*nuc[2]*nuc[6]-2.*nuc[4]*nuc[4])/(nuc[2]*nuc[2]);
	if(maxPow> 4) maxPow = 4;
}
void GRPolytrope::getNuCenter(double *nc, int& maxPow){
	if(maxPow>=0) nc[0] = nuc[0];
	if(maxPow>=2) nc[1] = nuc[2]; 
	if(maxPow>=2) nc[2] = nuc[4];
	if(maxPow>=4) nc[3] = nuc[6];
	if(maxPow> 4) maxPow = 4;
}
void GRPolytrope::getLambdaCenter(double *ll, int& maxPow){
	if(maxPow>=0) ll[0] = lc[0];
	if(maxPow>=2) ll[1] = lc[2];
	if(maxPow>=4) ll[2] = lc[4];
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
void GRPolytrope::setupSurface(){
	double xi1 = Y[len-1][x];
	double vs0 = Y[len-1][v];
	
	ths[0] = 0.0;
	ths[1] = -xi1*dydx[len-1];
	ths[2] = 0.0;
	ths[3] = 0.0;
	ths[4] = 0.0;
	
	nus[0] = log(xi1 - 2.*(n+1.)*vs0*sigma) - log(xi1); 
	nus[1] = -xi1*dNudr(len-1)*Rn;
	nus[2] = 0.0;
	nus[3] = 0.0;
	nus[4] = 0.0;
	
	ls[0] = log(xi1) -log(xi1 - 2.*(n+1.)*vs0*sigma); 
	ls[1] = -xi1*dLambdadr(len-1)*Rn;
	ls[2] = 0.0;
	ls[3] = 0.0;
	ls[4] = 0.0;

}
void GRPolytrope::getAstarSurface(double *As, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma : g);
	double AV = (n*Gam1-n-1.)/Gam1;
	int O=1;
	if(maxPow>=O-1) As[O-1] = AV;
	if(maxPow>=O+0) As[O  ] = AV*(ths[2]/ths[1] -1.) - n*sigma*ths[1];
	if(maxPow>=O+1) As[O+1] = AV*(
					( ths[1]*(2.*ths[3]-ths[2]) - ths[2]*ths[2] )/(ths[1]*ths[1])
					+ n*sigma*(ths[1]-2.*ths[2]) );
	if(maxPow >O+1) maxPow = O+1;
}
void GRPolytrope::getVgSurface(double *Vs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma : g);
	int O=1;
	if(maxPow>=O-1) Vs[O-1] = (n+1.)/Gam1;
	if(maxPow>=O+0) Vs[O  ] = (n+1.)/Gam1*( ths[2] - ths[1] )/ths[1];
	if(maxPow>=O+1) Vs[O+1] = (n+1.)/Gam1*( 
				 ths[1]*(2.*ths[3]-ths[2]) - ths[2]*ths[2]
			)/(ths[1]*ths[1]);
	if(maxPow >O+1) maxPow = O+1;
}
void GRPolytrope::getC1Surface(double *cs, int& maxPow){
	double x1 = Y[len-1][x];
	double v1 = Y[len-1][v];
	if(maxPow>=0) cs[0]  = v1/x1*std::exp(-nus[0])/ths[1];
	if(maxPow>=1) cs[1]  = v1/x1*std::exp(-nus[0])*(sigma*ths[1]*ths[1] - 2.*ths[2] - ths[1]*(1.+nus[1]))/ths[1]/ths[1];
	if(maxPow>=2) cs[2]  =  3.;
	if(maxPow> 2) maxPow = 2;
}
void GRPolytrope::getUSurface(double *Us, int& maxPow){
/*	double x1 = x[len-1];
	double θs1 = -z[len-1]*x1;
	double φs1 = -u[len-1]*x1;
	for(int a=0; a<maxPow; a++) Us[a] = 0.0;
	
	if(n==0){
		if(maxPow>=0) Us[0] = -x1*x1/φs1;
		if(maxPow>=1) Us[1] =  x1*x1*(3.*φs1 + x1*x1)/(φs1*φs1);
		if(maxPow>=2) Us[2] = -x1*x1*(3.*φs1*φs1 +4.*x1*x1*φs1 + x1*x1*x1*x1)/(φs1*φs1*φs1);
	}
	else if(n==1){
		if(maxPow>=0) Us[0] = 0.0;
		if(maxPow>=1) Us[1] = -x1*x1*θs1/φs1;
		if(maxPow>=2) Us[2] = x1*x1*(2.*θs1)/φs1 + 0.5*x1*x1*sigma*θs1;
	}
	else if(n==2){
		if(maxPow>=0) Us[0] = 0.0;
		if(maxPow>=1) Us[1] = 0.0;
		if(maxPow>=2) Us[2] = -x1*x1*θs1*θs1/φs1;
	}*/
	if(maxPow>=0) Us[0] = 0.0;
	if(maxPow>=1) Us[1] = 0.0;
	if(maxPow>=2) Us[2] = 0.0;
	if(maxPow >2) maxPow = 2;
	
}

void GRPolytrope::getNuSurface(double *ns, int& maxPow){
	if(maxPow>=0) ns[0] = nus[0];
	if(maxPow>=1) ns[1] = nus[1];
	if(maxPow>=2) ns[2] = nus[2];
	if(maxPow>=3) ns[3] = nus[3];
	if(maxPow>=4) ns[4] = nus[4];
	if(maxPow> 4) maxPow = 4;
}
void GRPolytrope::getLambdaSurface(double *ll, int& maxPow){
	if(maxPow>=0) ll[0] = ls[0];
	if(maxPow>=1) ll[1] = ls[1];
	if(maxPow>=2) ll[2] = ls[2];
	if(maxPow>=3) ll[3] = ls[3];
	if(maxPow>=4) ll[4] = ls[4];
	if(maxPow> 4) maxPow = 4;
}


double GRPolytrope::SSR(){
	double checkTOV1=0.0, checkTOV2=0.0, checkEuler=0.0;
		
	//sum up errors in equations
	double r = 0.0;
	double R = Radius();
	double expLambda;
	double e1, e2, e3, n1, n2, n3;
	char txtname[100], outname[100];
	for(int X=1; X<len-1; X++){
		r = rad(X);
		expLambda = exp(-Lambda(X));
		//TOV1 -- equation 2.3 of Tooper1
		e1 = fabs( expLambda*(r*dNudr(X)     + 1.0) - 1.0 - 8.0*m_pi*P(X)*pow(r,2) );
		n1 = fabs( expLambda*r*dNudr(X) )    + fabs(expLambda) + 1.0 + fabs(8.0*m_pi*P(X)*pow(r,2));
		//TOV2 -- equation 2.4 of Tooper1
		e2 = fabs( expLambda*(r*dLambdadr(X) - 1.0) + 1.0 - 8.0*m_pi*rho(X)*pow(r,2) );
		n2 = fabs( expLambda*r*dLambdadr(X)) + fabs(expLambda) + 1.0 + fabs(8.0*m_pi*rho(X)*pow(r,2));
		//Euler equation -- equation 2.6 of Tooper1
		e3 = fabs( dPdr(X) + 0.5*(rho(X)+P(X))*dNudr(X) );
		n3 = fabs(dPdr(X)) + fabs( 0.5*(rho(X)+P(X))*dNudr(X) );
		e1 /= n1;
		e2 /= n2;
		e3 /= n3;
		checkTOV1  += e1*e1;
		checkTOV2  += e2*e2;
		checkEuler += e3*e3;	
	}
	return sqrt((checkTOV1 + checkTOV2 + checkEuler)/double(3*len));
}

//method to print pertinent values of star to .txt, and plot them in gnuplot
void GRPolytrope::writeStar(char *c){
	//create names for files to be opened
	char pathname[256];
	char txtname[256];
	char outname[256];
	if(c==NULL)	sprintf(pathname, "./out/%s", name);
	else{
		sprintf(pathname, "./%s/star/", c);
	}
	char command[300];
	sprintf(command, "mkdir -p %s", pathname);
	
	
	sprintf(txtname, "%s/%s.txt", pathname, name);
	sprintf(outname, "%s/%s.png", pathname, name);

	FILE *fp = fopen(txtname, "w");

	//print results to text file
	// radius rho pressure gravity
	double irc=1./rho(0), ipc=1./P(0), R=Radius(), ig=1./dPhidr(length()-1), MT = Mass();
	for(int X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			rad(X)/R, 
			rho(X)*irc, -drhodr(X)*irc*R, 
			P(X)*ipc, -dPdr(X)*ipc*R, 
			mr(X)/MT, 
			Nu(X),
			Lambda(X));
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
	fprintf(gnuplot, "  '%s' u 1:2 w l t 'rho'", txtname);
	fprintf(gnuplot, ", '%s' u 1:4 w l t 'P'", txtname);
	fprintf(gnuplot, ", '%s' u 1:6 w l t 'm'", txtname);
	fprintf(gnuplot, ", '%s' u 1:7 w l t 'nu'", txtname);
	fprintf(gnuplot, ", '%s' u 1:8 w l t 'lambda'", txtname);
	fprintf(gnuplot, "\n");
	
	printBV(pathname, 5./3.);
	printCoefficients(pathname, 5./3.);
}

#endif