//**************************************************************************************
//							CHANDRASEKHAR WHITE DWARF
// ChandrasekharWD.cpp
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//		This model assumes T=0, and ignores Coulombic and other effects
//		Future step will be to adjust the composition, include finite T effects
// **************************************************************************************/

#ifndef ChandrasekharWDCLASS
#define ChandrasekharWDCLASS


#include "ChandrasekharWD.h"

//initalize white dwarf from central value of y and length
ChandrasekharWD::ChandrasekharWD(double Y0, double mu_electron, int L)
	: Y0(Y0), len(L)
{
	//exclude unphysical values of Y0
	if(Y0 < 1.) Y0 = 1.;
	Y02 = Y0*Y0;
	X02 = Y02-1.;
	X0  = sqrt(X02);
	//name this model for files
	sprintf(name, "ChandrasekharWD.%1.1f", Y0);

	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with y[len-1]=0.0
	//we will find the proper dx with a bisection search

	//an initial guess -- based on constant volume star
	double dx = sqrt(6.0)/(len-1), yS=1.0, ddx = dx;
	xi= new double[len]; //normalized radius
	x = new double[len]; //the relativity factor x = pF/mc
	y = new double[len]; //Chandrasekhar's y, y^2=1+x^2
	z = new double[len]; //derivative (dy/dxi)  note: dx/dxi = (dy/dxi)/x
	double dxmax=1.0, dxmin=0, ySmax=-1.0, ySmin=1.0;
	//this is used to control bisection search to guarantee it terminates
	double dxold = 1.0;
		
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
	RK4integrate(len, dx, 1);

	//now set physical properties of the white dwarf
	
	//set initial density, pressure
	//the values of A0, B0, as defined in lib/chandra.h
	//  these have changed greatly since Chandrasekhar's time, 
	//  so to compare to the table in his 1939 book, change to use his values
	A0 = Chandrasekhar::A0;
	B0 = Chandrasekhar::B0*mu_electron;
	Rn = sqrt(2.*A0/(m_pi*Gee()))/B0;
	
	mass = new double[len];
	f    = new double[len];
	mass[0] = 0.0;
	x[0] = X0;
	f[0] = Chandrasekhar::factor_f(x[0]);
	for(int X=1; X<len; X++){
		x[X] = sqrt(y[X]*y[X]-1.);
		if(y[X]<1.) x[X] = 0.0;	
		mass[X] = -4.*m_pi*xi[X]*xi[X]*z[X];
		f[X] = Chandrasekhar::factor_f(x[X]);
	}
	indexFit = len/2;
	for(int X=1; X<len; X++){
		//scan through x, set matching point where x[X] = 0.5
		if(x[X-1]>0.5 & x[X+1]<=0.5) indexFit = X;
	}
	indexFit /= 2;
	setupCenter();
	setupSurface();
	printf("%0.2lf\t%le\t%le\n", 1./Y02, Mass()/MSOLAR, Radius());
}


//initalize white dwarf from central value of y and length, with specific step size
ChandrasekharWD::ChandrasekharWD(double Y0, int L, const double dx)
	: Y0(Y0), len(L)
{
	//exclude unphysical values of Y0
	if(Y0 < 1.) Y0 = 1.;
	Y02 = Y0*Y0;
	X0  = sqrt(Y02-1.);
	X02 = Y02-1.;
	//name this model for files
	sprintf(name, "ChandrasekharWD.%1.1f", Y0);

	//reserve room for arrays
	xi= new double[len]; //normalized radius
	x = new double[len]; //the relativity factor x = pF/mc
	y = new double[len]; //Chandrasekhar's y, y^2=1+x^2
	z = new double[len]; //derivative (dy/dxi)  note: dx/dxi = (dy/dxi)/x
		
	//we know dx, so no need for bisection search
	RK4integrate(len, dx, 1);

	//now set physical properties of the white dwarf
	A0 = 1.0; //set to 1 for better numerical accuracy
	B0 = 1.0; //set to 1 for better numerical accuracy
	Rn = sqrt(2.*A0/(m_pi*Gee()))/B0;	//the radial scale factor
	
	mass = new double[len];
	f    = new double[len];
	mass[0] = 0.0;
	x[0] = X0;
	f[0] = Chandrasekhar::factor_f(x[0]);
	for(int X=1; X<len; X++){
		x[X] = sqrt(y[X]*y[X]-1.);
		if(y[X]<1.) x[X] = 0.0;	
		mass[X] = -4.*m_pi*xi[X]*xi[X]*z[X];
		f[X] = Chandrasekhar::factor_f(x[X]);
	}
	indexFit = 512*round(double(len)/1024.0);
	indexFit /= 2;
	printf("  indexFit  = %d\n", indexFit);
	printf("r[indexFit] = %0.32le\n", rad(indexFit));
	setupCenter();
	setupSurface();
	printf("%0.2lf\t%le\t%le\n", 1./Y02, Mass()/MSOLAR, Radius());
}


ChandrasekharWD::~ChandrasekharWD(){
	delete[] xi;
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] mass;
	delete[] f;
}

//integrate the polytrope up to Len using RK4
double ChandrasekharWD::RK4integrate(const int Len, double dx){
	//the YC, ZC, XC are corrected values of y, z, x to be used in equations
	double YC,ZC,XC;
	//the arrays K and L contain the corrections to y and z to be applied
	// when generating the next step of the correction
	double K[4], L[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	
	std::complex<double> YCN(0.0,0.0);//helps avoid NaNs near edge
	
	//set our initial conditions
	xi[0] = 0.0; y[0] = Y0; z[0] = 0.0;
	for(int X = 0; X<Len-1; X++){
		XC = xi[X]; YC = y[X]; ZC = z[X];
		for(int a = 0; a<4; a++){
			//now from these, calculate next shift from derivatives
			YCN = YC*YC - 1.;
			YCN = pow(YCN, 1.5);
			// K = dy = dxi*(dy/dxi) = dxi*z
			K[a] =  dx*ZC;
			// L = dz = dxi*(dy^2/dxi^2) = dxi*[ -(y^2-1)^3/2 - 2(dy/dxi)/xi ]
			L[a] = -dx*( YCN.real() + 2.0*ZC/XC );
			if(XC==0) L[0] = -dx/3.0*YCN.real(); //must use analytic expression at very center
			//calculate "corrected" positions using previous shift vectors
			XC = (double(X)+B[a])*dx;//x[X] + B[a]*dx;
			YC = y[X] + B[a]*K[a];
			ZC = z[X] + B[a]*L[a];
		}		
		xi[X+1] = double(X+1)*dx;
		y[X+1] = y[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		z[X+1] = z[X] + L[0]/6.0 + L[1]/3.0 + L[2]/3.0 + L[3]/6.0;
		//sometimes, for noninteger n, if array is too large values y<0 lead to nans
		//if these occur, it is safe to terminate the integration
		if(y[X+1]<1.0) {return y[X+1]-1.;}
	}
	return y[Len-1]-1.;
}

int ChandrasekharWD::RK4integrate(const int Len, double dx, int grid){
	grid=1;
	//the YC, ZC, XC are mid-point values of y, z, x to be used in equations
	double YC,ZC,XC;
	//the arrays K and L contain the corrections to y and z
	double K[4], L[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};

	std::complex<double> YCN(0.0,0.0);
	
	//set our initial conditions
	x[0] = 0.0; y[0] = Y0; z[0] = 0.0;
	for(int X = 0; X<Len-1; X++){
		XC = double(X)*dx; YC = y[X]; ZC = z[X];
		for(int a = 0; a<4; a++){
			//to handle NaNs at surface, use real part of complex root
			YCN = YC*YC - 1.;
			YCN = pow(YCN, 1.5);
			K[a] = dx*ZC;
			L[a] = -dx*( YCN.real() + 2.0*ZC/XC );
			if(X==0) L[0] = -dx/3.0*YCN.real();
			//calculate "corrected" positions using previous shift vectors
			XC = (double(X)+B[a])*dx;
			YC = y[X] + B[a]*K[a];
			ZC = z[X] + B[a]*L[a];
		}
		xi[X+1] = double(X+1)*dx;
		y[X+1] = y[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		z[X+1] = z[X] + L[0]/6.0 + L[1]/3.0 + L[2]/3.0 + L[3]/6.0;
	}
	return Len;
}

//Here we define functions to access radius, pressure, etc.
double ChandrasekharWD::rad(int X){
	return Rn*xi[X];
}
double ChandrasekharWD::rho(int X){
	return B0*pow(x[X],3);
}
double ChandrasekharWD::drhodr(int X){
	return 3.*B0*x[X]*y[X]*z[X]/Rn; //note: dx/dxi = z*y/x
}
double ChandrasekharWD::P(int X){
	return A0*f[X];
}
double ChandrasekharWD::dPdr(int X){
	return 8.*A0*pow(x[X],3)*z[X]/Rn; //see eqn (9) in Chandrasekhar 1935
}
double ChandrasekharWD::Phi(int X){
	//zeroed to join exterior solution at surface, where Phi->0 at infty
	return -8.*A0/B0 * (y[X] - 1. + z[len-1]*xi[len-1]);
}
double ChandrasekharWD::dPhidr(int X){
	return -8.*A0/B0 * z[X]/Rn;
}
double ChandrasekharWD::mr(int X){
	return B0*pow(Rn,3)*mass[X];
}

double ChandrasekharWD::Schwarzschild_A(int X, double GamPert){
	if(GamPert==0.0) return 0.0;
	else        	 return -z[X]*(8.*pow(x[X],3)/f[X]/GamPert-3.*y[X]/pow(x[X],2))/Rn;
}

double ChandrasekharWD::getAstar(int X, double GamPert){
	if(GamPert==0.0) return 0.0;
	else        	 return xi[X]*z[X]*(8.*pow(x[X],3)/f[X]/GamPert-3.*y[X]/pow(x[X],2));
}

double ChandrasekharWD::getU(int X){
	if(X==0) return 3.0;
	return -xi[X]*pow(x[X],3)/z[X];
}

double ChandrasekharWD::getVg(int X, double GamPert){
	if(GamPert==0.0) return -3.*y[X]*xi[X]*z[X]/pow(x[X],2);
	else			 return -8.*pow(x[X],3)*xi[X]*z[X]/f[X]/GamPert;
}

double ChandrasekharWD::getC(int X){
	if(X==0) return z[len-1]*xi[len-1]/2./yc[1];
	return (z[len-1]/z[X]) * (xi[X]/xi[len-1]);
}

double ChandrasekharWD::Gamma1(int X){
	return 8./3.*pow(x[X],5)/(y[X]*f[X]);
}

double ChandrasekharWD::sound_speed2(int X, double GamPert){
	if(GamPert == 0.0) return 8./3.  *A0/B0*x[X]*x[X]/y[X];
	else               return GamPert*A0/B0*f[X]/pow(x[X],3);
}


double ChandrasekharWD::Radius(){return rad(len-1);}	//total radius
double ChandrasekharWD::Mass(){return mr(len-1);}		//total mass
// the units of this model are set by our use of A0,B0 in CGS
double ChandrasekharWD::Gee(){return G_CGS;};
//in Newtonian, light speed is infinity...
double ChandrasekharWD::light_speed2(){return C_CGS*C_CGS;};



// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
void ChandrasekharWD::setupCenter(){	
	//the central coefficients -- expanded in terms of r/R
	double z1 =  z[len-1];
	double x1 = xi[len-1];
	//
	Y02 = Y0*Y0;
	X0  = sqrt(Y02-1.);
	X02 = Y02-1.;
	//
	yc[0]=	Y0;
	yc[1]=	-X02*X0/6.*pow(x1,2);
	yc[2]=	X02*X02*Y0/40.*pow(x1,4);
	yc[3]=	X0*(5.-29.*Y02+43.*Y02*Y02-19.*Y02*Y02*Y02)/5040.*pow(x1,6);
	//
	xc[0]=	X0;
	xc[1]=	Y0*yc[1]/X0;
	//
	fc[0]=	Y0*X0*(2.*X0-3.) + asinh(X0);
	fc[1]=	2.*xc[1]*(3.*xc[0]*xc[0]*(xc[0]-1.)+2.*xc[0]-1.);
}

void ChandrasekharWD::getAstarCenter(double *Ac, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	//depending on power requested, return appropriate number of terms
	if(maxPow>=0) Ac[0] = 0.0;
	if(maxPow>=2) Ac[0] = 16.*X0*X02*yc[1]/Gam1/fc[0] - 6.*Y0*yc[1]/X02;
	if(maxPow>=4) Ac[0] = 
		-16.*fc[1]*X02*X0*yc[1]/Gam1/fc[0]/fc[0]
			+ 48.*X02*xc[1]*yc[1]/Gam1/fc[0]
			+ 12.*xc[1]*Y0*yc[1]/X02/X0
			- 6.*yc[1]*yc[1]/X02
			+ 32.*X02*X0*yc[2]/Gam1/fc[0] - 12.*Y0*yc[2]/X02;		
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}

void ChandrasekharWD::getVgCenter(double *Vc, int& maxPow, double g){
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

void ChandrasekharWD::getUCenter(double *Uc, int& maxPow){
	double x1 = xi[len-1];
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] = -0.6*Y0*X0*pow(x1,2);
	if(maxPow>=4) Uc[2] = (25.-57.*Y02+32.*Y02*Y02)/350.*pow(x1,4);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}

void ChandrasekharWD::getC1Center(double *cc, int& maxPow){
	double x1 = xi[len-1];
	double z1 = z[len-1];
	if(maxPow>=0) cc[0] = z1*x1/2./yc[1];
	if(maxPow>=2) cc[1] = -z1*x1*yc[2]/yc[1]/yc[1];
	if(maxPow>=4) cc[2] = -z1*x1*(1.5*yc[1]*yc[3]-2.*yc[2]*yc[2])/pow(yc[1],3);
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
void ChandrasekharWD::setupSurface(){}

void ChandrasekharWD::getAstarSurface(double *As, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-1) : g);
	double x1 = xi[len-1];
	double a1 = -z[len-1]*x1;
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

void ChandrasekharWD::getVgSurface(double *Vs, int& maxPow, double g){
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
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			rad(X)/R, rho(X)*irc, -drhodr(X)*irc*R,
			P(X)*ipc, -dPdr(X)*ipc*R,
			mr(X)/Mass(), dPhidr(X)*ig);
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
	fprintf(gnuplot, ", '%s' u 1:4 w l t 'P'", txtname);
	fprintf(gnuplot, ", '%s' u 1:6 w l t 'm'", txtname);
	fprintf(gnuplot, ", '%s' u 1:7 w l t 'g'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	
	sprintf(txtname, "%s/degenerate.txt", pathname);
	sprintf(outname, "%s/degenerate.png", pathname);
	fp = fopen(txtname, "w");
	for(int X=0; X< length(); X++){
		fprintf(fp, "%le\t%le\t%le\t%le\n", xi[X], x[X], y[X], f[X]);
	}
	fclose(fp);	
	gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Degenerate functions for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'xi'\n");
	fprintf(gnuplot, "set ylabel 'x, y, f'\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'rho'", txtname);
 	fprintf(gnuplot, "set xrange [0:%le]\n", xi[len-1]);
 	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'x'", txtname);
 	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'y'", txtname);
 	fprintf(gnuplot, ",    '%s' u 1:4 w l t 'f'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	
}

void ChandrasekharWD::printCoefficients(char *c){
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
	int bc = 4;
	double A0[bc/2+1], V0[bc/2+1], c0[bc/2+1], U0[bc/2+1];
	getAstarCenter(A0, bc, gam1);
	getVgCenter(V0, bc, gam1);
	getUCenter(U0, bc);
	getC1Center(c0, bc);
	fprintf(fp, "A*  :\t%0.16le\t%0.16le\t%0.16le\n", A0[0],A0[1],A0[2]);
	fprintf(fp, "U   :\t%0.16le\t%0.16le\t%0.16le\n", U0[0],U0[1],U0[2]);
	fprintf(fp, "Vg  :\t%0.16le\t%0.16le\t%0.16le\n", V0[0],V0[1],V0[2]);
	fprintf(fp, "c1  :\t%0.16le\t%0.16le\t%0.16le\n", c0[0],c0[1],c0[2]);
	fclose(fp);
	sprintf(txtname, "%s/surface.txt", pathname);
	fp = fopen(txtname, "w");
	double A1[bc+1], V1[bc+1], c1[bc+1], U1[bc+1];
	getAstarSurface(A1, bc, gam1);
	getVgSurface(V1, bc, gam1);
	getUSurface(U1, bc);
	getC1Surface(c1, bc);
	fprintf(fp, "A*  :\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", A1[0],A1[1],A1[2],A1[3],A1[4]);
	fprintf(fp, "U   :\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", U1[0],U1[1],U1[2],U1[3],U1[4]);
	fprintf(fp, "Vg  :\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", V1[0],V1[1],V1[2],V1[3],V1[4]);
	fprintf(fp, "c1  :\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", c1[0],c1[1],c1[2],c1[3],c1[4]);
	fclose(fp);
	
	//print fits to those coefficients at center and surface
	int NC=15, NS=15;
	sprintf(txtname, "%s/centerfit.txt", pathname);
	fp = fopen(txtname, "w");
	double x2 = 0.0, xi1 = xi[len-1];
	fprintf(fp, "x           \tA*          \tA*_fit      \tU           \tU_fit       \tVg          \tVg_fit      \tc1          \tc1_fit\n");
	for(int X=0; X<NC; X++){
		x2 = xi[X]/xi1;
		fprintf(fp, "%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\n",
			x2,
			getAstar(X, gam1),
			A0[0] + A0[1]*x2*x2 + A0[2]*x2*x2*x2*x2,
			getU(X),
			U0[0] + U0[1]*x2*x2 + U0[2]*x2*x2*x2*x2,
			getVg(X, gam1),
			V0[0] + V0[1]*x2*x2 + V0[2]*x2*x2*x2*x2,
			getC(X),
			c0[0] + c0[1]*x2*x2 + c0[2]*x2*x2*x2*x2
		);
	}
	fclose(fp);
	sprintf(txtname, "%s/surfacefit.txt", pathname);
	fp = fopen(txtname, "w");
	double t = 1.;
	fprintf(fp, "x           \tA*          \tA*_fit      \tU           \tU_fit       \tVg          \tVg_fit      \tc1          \tc1_fit\n");
	for(int X=len-1; X>=len-NS-1; X--){
		t = 1. - xi[X]/xi1;
		fprintf(fp, "%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\n",
			xi[X]/xi1,
			getAstar(X, gam1),
			A1[0]/t + A1[1] + A1[2]*t + A1[3]*t*t + A1[4]*t*t*t,
			getU(X),
			U1[0] + U1[1]*t + U1[2]*t*t + U1[3]*t*t*t + U1[4]*t*t*t*t,
			getVg(X, gam1),
			V1[0]/t + V1[1] + V1[2]*t + V1[3]*t*t + V1[4]*t*t*t,
			getC(X),
			c1[0] + c1[1]*t + c1[2]*t*t + c1[3]*t*t*t + c1[4]*t*t*t*t
		);
	}
	fclose(fp);
	
	//print the pulsation coeffcients frequency
	sprintf(txtname, "%s/coefficients.txt", pathname);
	sprintf(outname, "%s/coefficients.png", pathname);
	fp  = fopen(txtname, "w");
	for(int X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			xi[X]/xi1,
			getAstar(X),
			getU(X),
			getVg(X),
			getC(X)
		);
	}
	fclose(fp);	
	//plot file in png in gnuplot
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
	fprintf(gnuplot, "set xrange [0:%le]\n", xi[NC]/xi1);
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
	fprintf(gnuplot, "set xrange [%le:1]\n", xi[len-NS-1]/xi1);
	fprintf(gnuplot, "plot '%s' u 1:(abs($2-$3)/abs($2)) w lp t 'A*'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($4-$5)/abs($4)) w lp t 'U'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($6-$7)/abs($6)) w lp t 'Vg'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($8-$9)/abs($8)) w lp t 'c1'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
}

#endif