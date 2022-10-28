//**************************************************************************************
//							CHANDRASEKHAR WHITE DWARF
// ChandrasekharWD++.cpp
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//		This model assumes T=0, and ignores Coulombic and other effects
//		The surface is not treated in any special way
//		Updated to include composition gradient, indicated by mu_e
// **************************************************************************************/

#ifndef ChandrasekharWDCLASS
#define ChandrasekharWDCLASS

#include "ChandrasekharWD++.h"

/*void ChandrasekharWD::chemical_gradient(
	const double xi, 	// the radial position
	const double dx, 	// the radial step size -- used to scale total radius
	double& mu, 		// return for the chemical potential
	double& dmu			// return for the derivative of the chemical potential
){
	//if the relative core is set to 1, the chem gradient is constant
	if(acore==1.){
		mu = mu0;
		dmu = 0.0;
		return;
	}
	//otherwise calculate the chem gradient as a sigmoidal curve from mu0 to 1
	double xcore = dx*(acore*len);	// position of the 
	double xswap = dx*(aswap*len);	// position of sigmoid midpoint
	double EXP = exp(k*(xi-xswap));	// the expoential
	// within the core, constant distribution
	if(xi<xcore){
		mu = mu0;
		dmu= 0.0;
	}
	// outside the core, sigmoidal distribution
	else {
		mu = (mu0-1.) + 1.0/(1.+EXP);
		dmu= -k*EXP*pow(1.+EXP,-2);
	}
}*/

void ChandrasekharWD::chemical_gradient(
	const double x, 	// the degeneracy factor
	const double dydxi, // the derivative dy/dxi, needed to find dmue/dxi
	double& mu, 		// return for the chemical potential mue
	double& dmu			// return for the derivative of the chemical potential, dmue/dxi
){
	//if the relative core is set to 1, the chem gradient is constant
	if(mu0==1.0){
		mu =  1.0;
		dmu = 0.0;
		return;
	}
	//otherwise calculate the chem gradient as a sigmoidal curve from mu0 to 1
	k = 2.0;//
	double F0 = Chandrasekhar::factor_f(X0);
	double F  = Chandrasekhar::factor_f(x);
	double z = -log(F/F0);
	double zeec = -log(1e-5);
	double EXP = exp(k*(z-zeec));	// the expoential
	mu = (mu0-1.) + 1.0/(1.+EXP);
	// IF USING LOG X
	//   dmudx = k/x * EXP/(1+EXP)^2
	//   dmudxi = dmudx * sqrt(1+x^2)/x * dydxi
	//dmu = k/x * sqrt(1.+x*x)/x * EXP*pow(1.+EXP,-2) * dydxi;
	// IF USING X
	//   dmudx = -k * EXP/(1+EXP)^2
	//   dmudxi = dmudx * sqrt(1+x^2)/x * dydxi
	//dmu = -k * sqrt(1.+x*x)/x * EXP*pow(1.+EXP,-2) * dydxi;
	// IF USING F
	//   dmudx = - k f'(x) * EXP/(1+EXP)^2 = - 8kx^4/y * EXP/(1+EXP)^2
	//   dmudxi = dmudx * sqrt(1+x^2)/x * dydxi = -8kx^3 * EXP/(1+EXP)^2
	//dmu  = -8.*k*pow(x,3) * EXP*pow(1.+EXP,-2) * dydxi;
	// IF USING LOG F
	//   dmudx = k f'(x)/f(x) * EXP/(1+EXP)^2  = 8k x^4/y/f * EXP/(1+EXP)^2
	//   dmudxi = dmudx * y/x * dydxi  = 8k x^3/f * dydxi * EXP/(1+EXP)^2
	dmu = 8.*k*pow(x,3)/F * EXP*pow(1.+EXP,-2) * dydxi;
	//printf("%le %le %le %le %le\n", x, Chandrasekhar::factor_f(x), z, EXP, mu); 	
}

void ChandrasekharWD::basic_setup(){
	//exclude unphysical values of Y0
	if(Y0 < 1.) Y0 = 1.;
	Y02 = Y0*Y0;
	X02 = Y02-1.;
	X0  = sqrt(X02);
	//name this model for files
	sprintf(name, "ChandrasekharWD.%1.3f", Y0);
	//set initial density, pressure
	//A0 = 6.01e22;	//this value taken from Chandraeskhar's book
	//B0 = 9.82e5;	//this value taken from Chandrasekhar's book
	Rn = sqrt(2.*A0/(m_pi*Gee()))/B0;
	printf("sigma=%le\n", A0/(B0*C_CGS*C_CGS));
	
	Y = new double*[len];
	for(int i=0; i<len; i++)
		Y[i] = new double[numvar];
}

void ChandrasekharWD::init_arrays(){
	mass = new double[len];
	mue  = new double[len];
	dmue = new double[len];
	for(int X=0; X<len; X++){
		Y[X][x] = sqrt(std::complex<double>(Y[X][y]*Y[X][y]-1.0,0.)).real();
		Y[X][f] = Chandrasekhar::factor_f(Y[X][x]);
		chemical_gradient(Y[X][x], Y[X][z], mue[X], dmue[X]);
		mass[X] = -4.*m_pi*Y[X][xi]*Y[X][xi]*Y[X][z]/mue[X];
	}
}

//initalize white dwarf from central value of y and length
ChandrasekharWD::ChandrasekharWD(
	double Y0, int L,
	double mu0, double k, double acore, double aswap
)
	: Y0(Y0), len(L),
	A0(Chandrasekhar::A0), B0(Chandrasekhar::B0),
	mu0(mu0), k(k), acore(acore), aswap(aswap)
{
	basic_setup();
	
	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with y[len-1]=0.0
	//we will find the proper dx with a bisection search

	//an initial guess -- based on constant volume star	
	dxi = 2.*sqrt(6.0)/(len-1);
	double dxmax=1.0, dxmin=0.0;
	
	//the bisection search
	// the function to search
	std::function<double(double)> find_surface = [this](double x)->double{return RK4integrate(len, x);};
	// use a Newton's mthod to find brackets around initial guess
	bisection_find_brackets_move(find_surface, dxi, dxmin, dxmax);
	// now contract the brackets using bisection
	bisection_search(find_surface, dxi, dxmin, dxmax);
	RK4integrate(len, dxi, 1);
	printf("edge  =\t%le\n", Y[len-1][y]-1.0);

	//now set physical properties of the white dwarf
	init_arrays();
	indexFit = len/2;
	for(int X=1; X<len-1; X++){
		//scan through x, set matching point where Y[X][x] = 0.5
		if(Y[X-1][x]>0.5 & Y[X+1][x]<=0.5) indexFit = X;
	}
	indexFit /= 2;
	setupCenter();
	setupSurface();
	printf("%0.2lf\t%le\t%le\n", 1./Y02, Mass()/MSOLAR, Radius());
}


//initalize white dwarf from central value of y and length, with specific step size
//this should only be used for testing scaling issues
ChandrasekharWD::ChandrasekharWD(
	double Y0, int L, const double dxi, 
	double mu0, double k, double acore, double aswap
)
	: Y0(Y0), len(L), dxi(dxi), 
	A0(Chandrasekhar::A0), B0(Chandrasekhar::B0),
	mu0(mu0), k(k), acore(acore), aswap(aswap)
{	
	basic_setup();	
	//we know dx, so no need for bisection search
	RK4integrate(len, dxi, 1);

	//now set physical properties of the white dwarf
	init_arrays();
	indexFit = 512*round(double(len)/1024.0);
	indexFit /= 2;
	printf("  indexFit  = %d\n", indexFit);
	printf("r[indexFit] = %0.32le\n", rad(indexFit));
	setupCenter();
	setupSurface();
	printf("%0.2lf\t%le\t%le\n", 1./Y02, Mass()/MSOLAR, Radius());
}


//initialize a sphere of degenerate fermions using 
//  mu_e = mu = constant
//  A0 = AN
//  B0 = BN
//  Use cases are:
//    1. testing results against tables with old values of A,B
//    2. implementing non-electron degenerate gasesous spheres (neutron stars)
//  This will use a constant value of mu_e throughout
ChandrasekharWD::ChandrasekharWD(
	double Y0, int L, 
	double mu, double AN, double BN
)
	: Y0(Y0), len(L), 
	A0(AN), B0(BN),
	mu0(mu), k(0), acore(1.0), aswap(0)
{
	basic_setup();
	printf("A0: %le\n", A0);
	printf("B0: %le\n", B0);
	printf("Rn: %le\n", Rn);

	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with x[len-1]=0.0
	//we will find the proper dx with a bisection search

	//an initial guess -- based on constant volume star
	dxi = sqrt(6.0)/(len-1);
	double dxmax=1.0, dxmin=0;
		
	std::function<double(double)> find_surface = [this](double x)->double{return RK4integrate(len, x);};
	bisection_find_brackets_newton(find_surface, dxi, dxmin, dxmax);
	bisection_search(find_surface, dxi, dxmin, dxmax);
	RK4integrate(len, dxi, 1);
		
	init_arrays();
	indexFit = len/2;
	for(int X=1; X<len; X++){
		//scan through x, set matching point where Y[X][x] = 0.5
		if(Y[X-1][x]>0.5 & Y[X+1][x]<=0.5) indexFit = X;
	}
	indexFit /= 2;
	setupCenter();
	setupSurface();
	printf("%0.8lf\t%le\t%le\t%le\n", 1./Y02, mr(len-1)/MSOLAR, rad(len-1)/1.e5, 1.-2.*G_CGS*mr(len-1)/(C_CGS*C_CGS*rad(len-1))  );
}

ChandrasekharWD::~ChandrasekharWD(){
	for(int i=0; i>len; i++)
		delete[] Y[i];
	delete[] mass;
	delete[] mue;
	delete[] dmue;
}

void ChandrasekharWD::centerInit(double ycenter[numvar]){
	ycenter[xi] = 0.0;
	ycenter[x] = X0;
	ycenter[y] = Y0;
	ycenter[z] = 0.0;
	ycenter[f] = Chandrasekhar::factor_f(X0);
}

void ChandrasekharWD::RK4step(double dx, double yin[numvar], double yout[numvar]){
	double YC[numvar] = {0.};
	double K[4][numvar] = {0.};
	static const double B[4] = {0.5,0.5,1.0,0.0};
	
	std::complex<double> YCN(0.,0.);
	double MUC = 1.0, DMUC = 0.0;
	
	for(int b=xi; b<=z; b++)
		YC[b] = yin[b];
	YC[x] = sqrt(yin[y]*yin[y]-1.);
	chemical_gradient(YC[x], YC[z], MUC, DMUC);
	
	for(int a=0; a<4; a++){
		// find the shift vector
		YCN = YC[y]*YC[y]-1.0;
		YCN = pow(YCN, 1.5);
		K[a][xi] = dx;
		K[a][y] =  dx*YC[z];
		K[a][z] = -dx*( YCN.real()*MUC*MUC + 2.0*YC[z]/YC[xi] - YC[z]*DMUC/MUC);
		// possible divide-by-zero errors at pure middle
		if(YC[xi]==0) 
			K[a][z] = -dx/3.0*YCN.real()*MUC*MUC;
		// prepare for next step
		for(int b=xi; b<=z; b++)
			YC[b] = yin[b] + B[a]*K[a][b];
		if(YC[y]>=1.0) YC[x] = sqrt(YC[y]*YC[y]-1.0);
		chemical_gradient(YC[x], YC[z], MUC, DMUC);
	}
	yout[xi] = yin[xi] + dx;
	for(int b=y; b<=z; b++)
		yout[b] = yin[b] + K[0][b]/6.0 + K[1][b]/3.0 + K[2][b]/3.0 + K[3][b]/6.0;
	yout[x] = sqrt(yout[y]*yout[y]-1.);
}


//integrate the polytrope up to Len using RK4
double ChandrasekharWD::RK4integrate(const int Len, double dx){
	//set our initial conditions
	centerInit(Y[0]);
	
	int X=0;
	for(X = 0; X<Len-1; X++){
		RK4step(dx, Y[X], Y[X+1]);
		if(Y[X+1][y]<1.0) {return Y[X+1][y]-1.;}
	}
	return Y[Len-1][y]-1.;
}

int    ChandrasekharWD::RK4integrate(const int Len, double dx, int grid){
	grid=1;
	
	//set our initial conditions
	centerInit(Y[0]);
	
	//integrate with RK4
	int X=0;
	for(X = 0; X<Len-1; X++){
		RK4step(dx, Y[X], Y[X+1]);
	}
	return Len;
}

//Here we define functions to access radius, pressure, etc.
double ChandrasekharWD::rad(int X){
	return Rn*Y[X][xi];
}
double ChandrasekharWD::rho(int X){
	return B0*mue[X]*pow(Y[X][x],3);
}
double ChandrasekharWD::drhodr(int X){
	return 3.*B0*mue[X]*Y[X][x]*Y[X][y]*Y[X][z]/Rn + B0*pow(Y[X][x],3)*dmue[X]/Rn; //note: dx/dxi = z*y/x
}
double ChandrasekharWD::P(int X){
	return A0*Y[X][f];
}
double ChandrasekharWD::dPdr(int X){
	return 8.*A0*pow(Y[X][x],3)*Y[X][z]/Rn; //see eqn (9) in Chandrasekhar 1935
}
double ChandrasekharWD::Phi(int X){
	//zeroed to join exterior solution at surface, where Phi->0 at infty
	return -8.*A0/B0 * (Y[X][y] - 1. + Y[len-1][z]*Y[len-1][xi])/mue[X];
}
double ChandrasekharWD::dPhidr(int X){
	return -8.*A0/B0 * Y[X][z]/Rn/mue[X];
}
double ChandrasekharWD::mr(int X){
	return B0*pow(Rn,3)*mass[X];
}

double ChandrasekharWD::Schwarzschild_A(int X, double GamPert){
	if(GamPert==0.0) return  dmue[X]/mue[X]/Rn;
	else        	 return  dmue[X]/mue[X]/Rn 
								- Y[X][z]*(8.*pow(Y[X][x],3)/Y[X][f]/GamPert-3.*Y[X][y]/pow(Y[X][x],2))/Rn;
}

double ChandrasekharWD::getAstar(int X, double GamPert){
	if(GamPert==0.0) return -dmue[X]/mue[X]*Y[X][xi];
	else        	 return -dmue[X]/mue[X]*Y[X][xi] 
								+ Y[X][z]*(8.*pow(Y[X][x],3)/Y[X][f]/GamPert-3.*Y[X][y]/pow(Y[X][x],2))*Y[X][xi];
}

double ChandrasekharWD::getU(int X){
	if(X==0) return 3.0;
	return -mue[X]*mue[X]*Y[X][xi]*pow(Y[X][x],3)/Y[X][z];
}

double ChandrasekharWD::getVg(int X, double GamPert){
	if(GamPert==0.0) return -8.*pow(Y[X][x],3)*Y[X][xi]*Y[X][z]/Y[X][f]/Gamma1(X);
	else			 return -8.*pow(Y[X][x],3)*Y[X][xi]*Y[X][z]/Y[X][f]/GamPert;
}

double ChandrasekharWD::getC(int X){
	if(X==0) return Y[len-1][z]*Y[len-1][xi]/2./yc[1]*mue[0]/mue[len-1];//-3.*Y[len-1][z]/Y[len-1][xi]/X02/X0;
	return (Y[len-1][z]/Y[X][z])*(Y[X][xi]/Y[len-1][xi])*(mue[X]/mue[len-1]);
}

double ChandrasekharWD::Gamma1(int X){
	if(Y[X][x]==0.) return 0.0;
	return 8./3.*pow(Y[X][x],5)/(Y[X][y]*Y[X][f]);
}

double ChandrasekharWD::sound_speed2(int X, double GamPert){
	if(GamPert == 0.0) return Gamma1(X)*A0/B0*Y[X][f]/pow(Y[X][x],3)/mue[X];
	else               return GamPert  *A0/B0*Y[X][f]/pow(Y[X][x],3)/mue[X];
}


double ChandrasekharWD::Radius(){return Rn*Y[len-1][xi];}	//total radius
double ChandrasekharWD::Mass(){return mr(len-1);}//total mass
double ChandrasekharWD::Gee(){return G_CGS;};
//in Newtonian, light speed is infinity...
double ChandrasekharWD::light_speed2(){return 1.0;};



// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
void ChandrasekharWD::setupCenter(){	
	//the central coefficients -- expanded in terms of r/R
	double z1 = Y[len-1][z];
	double x1 = Y[len-1][xi];
	double mu2 = mue[0]*mue[0];
	//
	yc[0]=  Y0;
	yc[1]= -X02*X0/6.*pow(x1,2)*mu2;
	yc[2]=  X02*X02*Y0/40.*pow(x1,4)*mu2*mu2;
	yc[3]=  X0*(-14.*X02*X02-19.*X02*X02*X02)/5040.*pow(x1,6)*pow(mu2,3);
	//
	xc[0]=  X0;
	xc[1]= -X02*Y0/6.*pow(x1,2)*mu2;
	xc[2]=  X02*X0*(4.+9.*X02)/360.*pow(x1,4);
	//
	fc[0]=  Chandrasekhar::factor_f(X0);
	fc[1]= -X02*X0*(2. - 3.*X0 + 3.*X02)/3.*pow(x1,2);
}

void ChandrasekharWD::getAstarCenter(double *Ac, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double xi2 = Y[len-1][xi]*Y[len-1][xi];
	double ac1 = (8.*X0*X02/Gam1/fc[0] - 3.*Y0/X02);
	//depending on power requested, return appropriate number of terms
	if(maxPow>=0) Ac[0] = 0.0;
	if(maxPow>=2) Ac[1] = 2.*yc[1]*ac1;
	if(maxPow>=4) Ac[2] = 4.*yc[2]*ac1
			- yc[1]*16.*X02*X0*fc[1]/Gam1/fc[0]/fc[0]
			+ yc[1]*48.*X02*xc[1]/Gam1/fc[0]
			+ yc[1]*12.*xc[1]*Y0/pow(X0,3)
			- 6.*yc[1]*yc[1]/X02;		
	//if more  terms than this requested, cap number of terms
	for(int i=6; i<=maxPow; i+=2) Ac[i] = 0.0;
}

void ChandrasekharWD::getVgCenter(double *Vc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double x1 = Y[len-1][xi];
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] = -16.*X02*X0*yc[1]/Gam1/fc[0];
	if(maxPow>=4) Vc[2] =-8.*   pow(x1,2)*pow(X0,6)*(5.*fc[1]+4.*pow(x1,2)*fc[0]*X0*Y0)/(15.*Gam1*fc[0]*fc[0]);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4; 
}

void ChandrasekharWD::getUCenter(double *Uc, int& maxPow){
	double x1 = Y[len-1][xi];
	double mu2 = mue[0]*mue[0];
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] =-0.6*mu2*x1*x1*X0*Y0;
	if(maxPow>=4) Uc[2] = 
		9.*xc[2]/X0 - pow(x1*x1*mu2,2)*X02*(0.08+0.134*X02);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}

void ChandrasekharWD::getC1Center(double *cc, int& maxPow){
	double x1 = Y[len-1][xi];
	double z1 = Y[len-1][z];
	double murat = mue[0]/mue[len-1];
	if(maxPow>=0) cc[0] = z1*x1/2./yc[1]*murat;
	if(maxPow>=2) cc[1] = -z1*x1*yc[2]/yc[1]/yc[1]*murat;
	if(maxPow>=4) cc[2] = -z1*x1*(1.5*yc[1]*yc[3]-2.*yc[2]*yc[2])/pow(yc[1],3)*murat;
	//if more terms than this requested, cap number of terms
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
void ChandrasekharWD::setupSurface(){}

void ChandrasekharWD::getAstarSurface(double *As, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-2) : g);
	double x1 =  Y[len-1][xi];
	double a1 = -Y[len-1][z]*x1;
	int O=1;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=-1) As[O-1] = 0.0;//1.5;
	if(maxPow>= 0) As[O  ] = 0.75*a1;
	if(maxPow>= 1) As[O+1] = 0.75*a1 - (8./Gam1 + 0.375)*a1*a1;
	if(maxPow>= 2) As[O+2] = 0.0;//a1*3./16.*pow(a1-2.,2)+4./3.*a1*a1*(12.+5.*a1)/Gam1;
	if(maxPow>= 3) As[O+3] = 0.0;//-3.*a1*pow(a1-2,3)/32.+(724.*a1*a1/45. + 20.*a1 + 24.)*a1*a1/Gam1;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 3) maxPow = O+3;
}

void ChandrasekharWD::getVgSurface(double *Vs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-2) : g);
	double x1 =  Y[len-1][xi];
	double a1 = -Y[len-1][z]*x1;
	int O=1;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=-1) Vs[O-1] = 1.5;
	if(maxPow>= 0) Vs[O  ] = 0.75*a1;
	if(maxPow>= 1) Vs[O+1] = 0.75*a1 + 8.*a1*a1/Gam1 - 3.*a1*a1/8.;
	if(maxPow>= 2) Vs[O+2] = a1*3./16.*pow(a1-2.,2)+4./3.*a1*a1*(12.+5.*a1)/Gam1;
	if(maxPow>= 3) Vs[O+3] = -3.*a1*pow(a1-2,3)/32.+(724.*a1*a1/45. + 20.*a1 + 24.)*a1*a1/Gam1;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 3) maxPow = O+3;
}

void ChandrasekharWD::getUSurface(double *Us, int& maxPow){
	//these happen to all be zero at surface
	for(int j=0; j<=maxPow; j++){
		Us[j] = 0.0;
	}
}

void ChandrasekharWD::getC1Surface(double *cs, int& maxPow){
	if(maxPow >= 0) cs[0] =  1.;
	if(maxPow >= 1) cs[1] = -3.;
	if(maxPow >= 2) cs[2] =  3.;
	if(maxPow >= 3) cs[3] = -1.;
	if(maxPow >= 4) cs[4] =  0.; //does not actually appear in equations
	//if more  terms than this requested, cap number of terms
	if(maxPow  > 4) maxPow = 4;
}


void ChandrasekharWD::writeStar(char *c){
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
	double irc=1./rho(0), ipc=1./P(0), R=Radius(), ig=1./dPhidr(length()-1);
	for(int X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t",
			rad(X)/R, rho(X)*irc, -drhodr(X)*irc*R,
			P(X)*ipc, -dPdr(X)*ipc*R,
			1.-mr(X)/Mass(), dPhidr(X)*ig);
		fprintf(fp, "%le\t%le\t%le", Y[X][x], Y[X][y], Y[X][f]);
		fprintf(fp, "\n");
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
	fprintf(gnuplot, "set xrange [0.9:1]\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'rho'", txtname);
	fprintf(gnuplot, ", '%s' u 1:4 w l t 'P'", txtname);
	fprintf(gnuplot, ", '%s' u 1:6 w l t '1-m'", txtname);
	fprintf(gnuplot, ", '%s' u 1:7 w l t 'g'", txtname);
	fprintf(gnuplot, "\n");
	
	sprintf(outname, "%s/degenerate.png", pathname);
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Degenerate functions for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'xi'\n");
	fprintf(gnuplot, "set ylabel 'x, y, f'\n");
 	fprintf(gnuplot, "set xrange [0.9:1]\n");
 	fprintf(gnuplot, "set logscale y 10\n");
 	fprintf(gnuplot, "plot '%s' u 1:8 w l t 'x'", txtname);
 	fprintf(gnuplot, ",    '%s' u 1:9 w l t 'y'", txtname);
 	fprintf(gnuplot, ",    '%s' u 1:10 w l t 'f'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	
	
	printBV(pathname);
	printCoefficients(pathname);
	printChemicalGradient(pathname);
}
	
void ChandrasekharWD::printChemicalGradient(char *pathname){
	char txtname[256];
	char outname[256];
	char title[256]; graph_title(title);
	
	//print the Brunt-Vaisala frequency
	//print the chemical gradient
	sprintf(txtname, "%s/chemical.txt", pathname);
	sprintf(outname, "%s/chemical.png", pathname);
	FILE *fp  = fopen(txtname, "w");
	double H,He;
	for(int X=1; X< length()-1; X++){
		He = 2.*(mue[X]-1.)/mue[X];
		H = 1.-He;
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\n",
			(1.-mass[X]/mass[len-1]),
			H,
			He
		);
	}
	fclose(fp);
	
	//plot chem gradient
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 800,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Chemical composition for %s'\n", title);
	fprintf(gnuplot, "set logscale x 10\n");
	fprintf(gnuplot, "set format x '10^{%%L}'\n");
	fprintf(gnuplot, "set xlabel 'log_{10} (1-m/M)\n");
	fprintf(gnuplot, "set ylabel 'abundance\n");
	fprintf(gnuplot, "set yrange [-0.01: 1.01]\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'X (H1)'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'Y/Z (He4/C12/O16)'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);	
}

#endif