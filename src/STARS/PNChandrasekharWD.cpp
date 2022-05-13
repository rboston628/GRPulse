//**************************************************************************************
//						POST-NEWTONIAN CHANDRASEKHAR WHITE DWARF
// PNChandrasekharWD.h
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar (1935)
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//  	Builds off of work in Tooper (1964) to implement CHWD in GR
//		This model assumes T=0, and ignores Coulombic and other effects
//		The surface is not treated in any special way
//		Assumes mu_e = 1
//
//		Utilizes my own work with perturbed 1PN equations of stellar structure
// 		We use units where G=c=1
//**************************************************************************************

#ifndef PNChandrasekharWDCLASS
#define PNChandrasekharWDCLASS

#include "PNChandrasekharWD.h"

//initalize polytrope from index and length
PNChandrasekharWD::PNChandrasekharWD(double Y0, double zsurf, int Len)
	: Y0(Y0), zsurf(zsurf), len(Len)
{
	//exclude unphysical values of Y0
	if(Y0 < 1.) Y0 = 1.;
	Y02 = Y0*Y0;
	X0  = sqrt(Y02-1.);
	X02 = Y02-1.;
	
	//now set physical properties of the polytrope
	//set initial density, pressure
	A0 = 6.01e22;	//this value taken from Chandraeskhar's book
	B0 = 9.82e5;	//this value taken from Chandrasekhar's book
	Rn = sqrt(2.*A0/(m_pi*G_CGS))/B0;
	sigma = A0/(B0*C_CGS*C_CGS);
	printf("sigma=%le\n", sigma);
	
	//check surface redshift and sigma
	if(zsurf >= 0.5){
		printf("\tYour specified star would collapse to a black hole.\n");
		printf("\tSwitching to a smaller value.\n");
		zsurf = 0.49;
	}
	//name this star for files
	sprintf(name, "PNChandrasekharWD%1.3f_%1.2le", Y0, zsurf);
	
	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with y[len-1]=0.0
	//we will find the proper dx with a bisection search
	
	//an initial guess -- based on 0pn constant-volume model
	double dxi = sqrt(6.0)/len, yS=1.0, ddxi = dxi;
	φc0 = -zsurf/sigma/2.;
	xi= new double[len];
	x = new double[len];
	u = new double[len];
	v = new double[len];
	φ = new double[len];
	
	printf("dxi=%le\n", dxi);
	//we first find reasonable values of the parameters φc0
	// this step is necessary for large values of z
	int tracker = 0;
	int const np = 2;
	double dfdz[2][2];
	double target[2] = {zsurf, zsurf};
	double z1[2] = {X0,φc0};
	double z2[2] = {0.,0.0};
	double dz[2] = {0.,0.0};
	double f1[2], f2[2];
	double abs1, abs2;
	RK4integrate(z1, f1, dxi);
	printf("x0 = %le\n", X0);
	printf("x1 = %le\n", x[len-1]);
	RK4integrate(z1, f1, dxi);
	RK4integrate(z1, f1, dxi);
	RK4integrate(z1, f1, dxi);
	printf("x0 = %le\n", X0);
	printf("x1 = %le\n", x[len-1]);
	printf("target = \t%0.16le\t%0.16le\n", zsurf, zsurf);
	printf("result = \t%0.16le\t%0.16le\n", f1[0], f1[1]);
	abs1=abs2=0.0;
	for(int i=0; i<np; i++) abs1 += pow(target[i]-f1[i],2);
	for(int i=0; i<np; i++) abs2 += pow(target[i],2);
	while(sqrt(abs1)/zsurf> 1.e-10 & tracker<100){
		//calculate the matrix of derivatived (df/dz)
		for(int i=0; i<np; i++) z2[i] = z1[i];
		for(int i=0; i<np; i++){
			z2[i] = 1.01*z1[i];
			RK4integrate(z2, f2, dxi,0);
			for(int j=0;j<np;j++){
				dfdz[j][i] = (f2[j]-f1[j])/(z2[i]-z1[i]);
			}
			z2[i] = z1[i];
		}
		for(int i=0; i<np; i++) f2[i] = (target[i]-f1[i]);
		//solve the equation (T-f1) = (dfdx)*(dz) for dz
		invertMatrix(dfdz, f2, dz);
		for(int i=0; i<np; i++) z2[i] = z1[i] + dz[i];
		//if the increase is too large, throttle it
		abs1=abs2=0.0;
		for(int i=0; i<np; i++) abs1 += z1[i]*z1[i];
		for(int i=0; i<np; i++) abs2 += z2[i]*z2[i];
		if( abs2 < 0.9*abs1 | 1.1*abs1 < abs2){
			for(int i=0; i<np; i++) z2[i] = z1[i] + 0.01*dz[i];
		}
		//use updated variable, and try again
		for(int i=0; i<np; i++) z1[i] = z2[i];
		RK4integrate(z1, f1, dxi);
		//prepare for next step
		abs1=abs2=0.0;
		for(int i=0; i<np; i++) abs1 += pow(target[i]-f1[i],2);
		for(int i=0; i<np; i++) abs2 += pow(target[i],2);
	}
	printf("target = \t%0.16le\t%0.16le\n", zsurf, zsurf);
	printf("result = \t%0.16le\t%0.16le\n", f1[0], f1[1]);
	X0 = z1[0];
	φc0= z1[1];
	X02 = X0*X0;
	Y02 = 1.0 + X02;
	Y0 = sqrt(Y02);
	printf("PRE-BISECT\n");
	printf("\tdxi=%le\n", dxi);
	printf("\tx0 = %le\n", X0);
	printf("\tx1 = %le\n", x[len-1]);
	
	double dximax=1.0, dximin=0, ySmax=-1.0, ySmin=1.0;
	yS = RK4integrate(len, dxi);
	//find brackets on dx that bound a zero in yS
	if(yS > 0){
		dximin = dxi; ySmin = yS;
		while(yS > 0 || isnan(yS)){
			dxi+= ddxi;
			yS = RK4integrate(len, dxi);
			if(isnan(yS)){
					dxi = 1.0/len; ddxi *= 0.1; yS = 1.0;
			}
		}
		dximax = dxi; ySmax = yS;
	}
	else if (yS < 0){
		dximax = dxi; ySmax = yS;
		while(yS < 0){
			dxi = 0.99*dxi;
			yS = RK4integrate(len,dxi);
		}
		dximin = dxi; ySmin = yS;
	}
	dxi = 0.5*(dximin+dximax);
	yS = RK4integrate(len, dxi);
	
	if(dxi<0) {printf("somehow dx is negative...\n"); dxi=-dxi;}
	if(ySmin*ySmax > 0.0) {printf("big problem, chief\n"); exit(EXIT_FAILURE);}

	//now use bisection to find dx so that yS=0.0
	double dxiold = 1.0;
	while( fabs(dximin-dximax)>0.0 || isnan(yS) ){
		dxi= 0.5*(dximin+dximax);
		yS = RK4integrate(len, dxi);
		
		if(isnan(yS)){
			yS = -1.0;
		}
		if( (yS*ySmax>0.0) ){
			dximax = dxi;
			ySmax  = yS;
		}
		else if( (yS*ySmin>0.0) ){
			dximin = dxi;
			ySmin  = yS;
		}
		//if the brackets are not moving, stop the search
		if(dxiold == fabs(dximin-dximax)) break;
		dxiold = fabs(dximin-dximax);
	}

	printf("POST-BISECT\n");
	printf("\tdxi= %le\n", dxi);
	printf("\tx0 = %le\n", X0);
	printf("\tx1 = %le\n", x[len-1]);
	z1[0] = X0; z1[1] = φc0;
	f1[0] = -8.*sigma*φ[len-1];
	f1[1] =  8.*sigma*u[len-1]*xi[len-1];
	printf("target = \t%0.16le\t%0.16le\n", zsurf, zsurf);
	printf("result = \t%0.16le\t%0.16le\n", f1[0], f1[1]);


	//check that end is at surface
	ddxi = dxi;
	RK4integrate(len, ddxi, 1);
	printf("\tx1 = %le\n", x[len-1]);
	z1[0] = X0; z1[1] = φc0;
	RK4integrate(z1, f1, ddxi, 0);
	abs1=abs2 = 0.0;
	for(int i=0; i<np; i++) abs1 += pow(target[i]-f1[i],2);
	for(int i=0; i<np; i++) abs2 += pow(target[i],2);
	while(sqrt(abs1)/zsurf>1.e-10 & tracker<100){
		for(int i=0; i<np; i++) z2[i] = z1[i];
		for(int i=0; i<np; i++){
			z2[i] = 1.01*z1[i];
			RK4integrate(z2, f2, ddxi, 0);
			for(int j=0;j<np;j++){
				dfdz[j][i] = (f2[j]-f1[j])/(z2[i]-z1[i]);
			}
			z2[i] = z1[i];
		}
		for(int i=0; i<np; i++) f2[i] = (target[i]-f1[i]);
		invertMatrix(dfdz, f2, dz);
		for(int i=0; i<np; i++) z2[i] = z1[i] + dz[i];
		abs1=abs2 = 0.0;
		for(int i=0; i<np; i++) abs1 += z1[i]*z1[i];
		for(int i=0; i<np; i++) abs2 += z2[i]*z2[i];
		if( abs2 < 0.9*abs1 | 1.1*abs1 < abs2){
			for(int i=0; i<np; i++) z2[i] = z1[i] + 0.01*dz[i];
		}
		for(int i=0; i<np; i++) z1[i] = z2[i];
		RK4integrate(z1, f1, ddxi,1);
		abs1=abs2 = 0.0;
		for(int i=0; i<np; i++) abs1 += pow(target[i]-f1[i],2);
		for(int i=0; i<np; i++) abs2 += pow(target[i],2);
		//tracker++;
	}
	X0 = z1[0];
	φc0= z1[1];
	X02 = X0*X0;
	Y02 = 1.0 + X02;
	Y0 = sqrt(Y02);
	RK4integrate(len, dxi, 1);
	printf("\tdxi= %le\n", dxi);
	printf("\tx0 = %le\n", X0);
	printf("\tx1 = %le\n", x[len-1]);
	f1[0] = -8.*sigma*φ[len-1];
	f1[1] =  8.*sigma*u[len-1]*xi[len-1];
	printf("target = \t%0.16le\t%0.16le\n", zsurf, zsurf);
	printf("result = \t%0.16le\t%0.16le\n", f1[0], f1[1]);

	//calculate derivative and other useful thngs
	y = new double[len];
	f = new double[len];
	h = new double[len];
	dx= new double[len];
	//
	indexFit = len/2;
	for(int X=0; X<len; X++){
		 y[X] = sqrt(1.0 + x[X]*x[X]);
		 f[X] = factor_f(x[X]);
		 h[X] = factor_h(x[X]);
		dx[X] = -((h[X] + sigma*f[X])*u[X] + sigma*h[X]*v[X])*y[X]*pow(x[X],-4);
		//as we scan through x,y,z, set matching point as where x[X] = 0.5
		if(x[X-1]>0.5 & x[X+1]<=0.5) indexFit = X;
	}
	indexFit /= 2;
	setupCenter();
	setupSurface();
}

//initalize polytrope from index, length, and dx
PNChandrasekharWD::PNChandrasekharWD(double Y0, double zsurf, int Len, const double dxi)
	: Y0(Y0), zsurf(zsurf), len(Len), dxi(dxi)
{
	//exclude unphysical values of Y0
	if(Y0 < 1.) Y0 = 1.;
	Y02 = Y0*Y0;
	X0  = sqrt(Y02-1.);
	X02 = Y02-1.;
	//check surface redshift and sigma
	if(zsurf >= 0.5){
		printf("\tYour specified star would collapse to a black hole.\n");
		printf("\tSwitching to a smaller value.\n");
		zsurf = 0.49;
	}
	//name this star for files
	sprintf(name, "PNChandrasekharWD%1.3f_%1.2le", Y0, sigma);

	//now set physical properties of the polytrope
	//set initial density, pressure
	A0 = 6.01e22;	//this value taken from Chandraeskhar's book
	B0 = 9.82e5;	//this value taken from Chandrasekhar's book
	Rn = sqrt(2.*A0/(m_pi*Gee()))/B0;
	sigma = A0/(B0*C_CGS*C_CGS);


	//initialize arrays
	xi= new double[len];
	x = new double[len];
	u = new double[len];
	v = new double[len];
	φ = new double[len];
	
	φc0 = -1.0;
	RK4integrate(len, dxi, 1);
	//continue adjusting the value of sigma until the surface redshift is correct
	int tracker = 0;
	while(fabs(φ[len-1] + zsurf/8./sigma) > 1e-16  & tracker<100){
		φc0 -= (zsurf/8./sigma + φ[len-1]);
		RK4integrate(len, dxi,1);
		tracker++;
	}

	//calculate derivative and other useful thngs
	 y = new double[len];
	 f = new double[len];
	 h = new double[len];
	dx = new double[len];
	for(int X=0; X<len; X++){
		 y[X] = sqrt(1.0 + x[X]*x[X]);
		 f[X] = factor_f(x[X]);
		 h[X] = factor_h(x[X]);
		dx[X] =-( (h[X] + sigma*f[X])*u[X] + sigma*h[X]*v[X]  )*y[X]*pow(x[X],-4);		
	}
	indexFit = 512*round(double(len)/1024.0);
	indexFit /= 2;
	printf("  indexFit  = %d\n", indexFit);
	printf("r[indexFit] = %0.32le\n", rad(indexFit));
}


PNChandrasekharWD::~PNChandrasekharWD(){
	delete[] x;
	delete[] y;
	delete[] dx;
	delete[] u;
	delete[] v;
	delete[] φ;
}


//integrate the polytrope up to Len using RK4
void PNChandrasekharWD::RK4integrate(double Z[2], double F[2], double& dxi, int f){
	//the YC, ZC, XC are corrected values of y, z, x to be used in equations
	double xi1,  x1,u1,v1,p1,  y1,f1,h1;
	double xi2,  x2,u2,v2,p2,  y2,f2,h2;
	double xiC,  xC,uC,vC,pC,  yC,fC,hC;
	//the arrays K and L contain the corrections to y and z to be applied
	// when generating the next step of the correction
	double K[4], L[4], M[4], N[4], P[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	std::complex<double> YCN(0.0,0.0);	
	//set our initial conditions
	xi2 = 0.0; x2=Z[0]; u2=0.0; v2=0.0; p2=Z[1]; 
	y2 = sqrt(x2*x2+1.);
	f2 = factor_f(x2);
	h2 = factor_h(x2);
	//integrate with RK4
	int X = 0;
	while(x2 >= 0.0){
		xi1=xi2; x1=x2; u1=u2; v1=v2; p1=p2; y1=y2; f1=f2; h1=h2;
		xiC=xi1; xC=x1; uC=u1; vC=v1; pC=p1; yC=y1; fC=f1; hC=h1;
		//calculate all Ks
		for(int a = 0; a<4; a++){
			//now from these, calculate next shift
			K[a] =-dxi*(  (hC + sigma*fC)*uC + sigma*hC*vC  )*yC*pow(xC,-4);
			M[a] = dxi*(   hC             - 2.*uC/xiC );
			N[a] = dxi*(3.*fC - 16.*hC*pC - 2.*vC/xiC );
			P[a] = dxi*uC;
			//at the very center, M and N need special treatment
			if(xiC==0) { 
				M[0]= dxi/3.0*factor_h(xC);
				N[0]= dxi/3.0*(3.0*factor_f(xC)-16.*factor_h(xC)*pC); K[0]=P[0]=0.0;}
			//calculate intermediate positions using previous shift vectors
			xiC= xi1+ B[a]*dxi;
			xC = x1 + B[a]*K[a]; // x(xi)
			uC = u1 + B[a]*M[a]; // dphi/dxi
			vC = v1 + B[a]*N[a]; // dpsi/dxi
			pC = p1 + B[a]*P[a]; // phi(xi)
			yC = sqrt(1.0 + xC*xC);// y(x)
			fC = factor_f(xC);     // f(x)
			hC = factor_h(xC);     // h(x)
		}
		xi2 =xi1+ dxi;
		x2 = x1 + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		u2 = u1 + M[0]/6.0 + M[1]/3.0 + M[2]/3.0 + M[3]/6.0;
		v2 = v1 + N[0]/6.0 + N[1]/3.0 + N[2]/3.0 + N[3]/6.0;
		p2 = p1 + P[0]/6.0 + P[1]/3.0 + P[2]/3.0 + P[3]/6.0;
		y2 = sqrt(1.0+x2*x2);
		f2 = factor_f(x2);
		h2 = factor_h(x2);
		X++;
	}
	x[len-1] = x1; 
	if(f) dxi = xi2/(len-1);
	F[0] =-8.*sigma*p1;
	F[1] = 8.*sigma*u1*xi1;
}

//integrate the polytrope up to Len using RK4
double PNChandrasekharWD::RK4integrate(const int Len, double dxi){
	//the YC, ZC, XC are corrected values of y, z, x to be used in equations
	double xiC, XC,UC,VC,ΦC,   YC,FC,HC;
	//the arrays K and L contain the corrections to y and z to be applied
	// when generating the next step of the correction
	double K[4], M[4], N[4], P[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	//set our initial conditions
	xi[0] = 0.0; x[0]=X0; u[0]=0.0; v[0]=0.0; φ[0]=φc0;
	//integrate with RK4
	int X=0;
	for(X = 0; X<Len-1; X++){
		xiC= xi[X];
		XC = x[X]; UC = u[X]; VC = v[X]; ΦC=φ[X];
		YC = sqrt(1.0 + XC*XC);
		FC = factor_f(XC);
		HC = factor_h(XC);
		//calculate all Ks
		for(int a = 0; a<4; a++){
			//now from these, calculate next shift
			K[a] =-dxi*(  (HC + sigma*FC)*UC + sigma*HC*VC  )*YC*pow(XC,-4);
			M[a] = dxi*(   HC             - 2.*UC/xiC );
			N[a] = dxi*(3.*FC - 16.*HC*ΦC - 2.*VC/xiC );
			P[a] = dxi*UC;
			//at the very center, M and N need special treatment
			if(xiC==0) { 
				M[0]= dxi/3.0*factor_h(X0);
				N[0]= dxi/3.0*(3.0*factor_f(X0)-16.*factor_h(X0)*φc0); K[0]=P[0]=0.0;}
			//calculate intermediate positions using previous shift vectors
			xiC= (double(X) + B[a])*dxi;
			XC = x[X] + B[a]*K[a]; // x(xi)
			UC = u[X] + B[a]*M[a]; // dphi/dxi
			VC = v[X] + B[a]*N[a]; // dpsi/dxi
			ΦC = φ[X] + B[a]*P[a]; // phi(xi)
			YC = sqrt(1.0 + XC*XC);// y(x)
			FC = factor_f(XC);     // f(x)
			HC = factor_h(XC);     // h(x)
		}
		xi[X+1]= double(X+1)*dxi;
		x[X+1] = x[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		u[X+1] = u[X] + M[0]/6.0 + M[1]/3.0 + M[2]/3.0 + M[3]/6.0;
		v[X+1] = v[X] + N[0]/6.0 + N[1]/3.0 + N[2]/3.0 + N[3]/6.0;
		φ[X+1] = φ[X] + P[0]/6.0 + P[1]/3.0 + P[2]/3.0 + P[3]/6.0;
		//sometimes, for noninteger n, if array is too large values y<0 lead to nans
		//if these occur, it is safe to terminate the integration
		if(x[X+1] <0.0){
			return x[X+1];
		}
	}
	return x[X];
}

int    PNChandrasekharWD::RK4integrate(const int Len, double dxi, int grid){
	grid=1;
	//the YC, ZC, XC are mid-point values of y, z, x to be used in equations
	double xiC, XC,UC,VC,ΦC,   YC,FC,HC;
	//the arrays K and L contain the corrections to y and z to be applied
	// when generating the next step of the correction
	double K[4], M[4], N[4], P[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	//set our initial conditions
	xi[0] = 0.0; x[0]=X0; u[0]=0.0; v[0]=0.0; φ[0]=φc0;
	
	for(int X = 0; X<Len-1; X++){
		xiC= xi[X];
		XC = x[X]; UC = u[X]; VC = v[X]; ΦC=φ[X];
		YC = sqrt(1.0 + XC*XC);
		FC = factor_f(XC);
		HC = factor_h(XC);
		for(int a = 0; a<4; a++){
			//now from these, calculate next shift
			K[a] =-dxi*(  (HC + sigma*FC)*UC + sigma*HC*VC  )*YC*pow(XC,-4);
			M[a] = dxi*(   HC             - 2.*UC/xiC );
			N[a] = dxi*(3.*FC - 16.*HC*ΦC - 2.*VC/xiC );
			P[a] = dxi*UC;
			//at the very center, M and N need special treatment
			if(xiC==0) { 
				M[0]= dxi/3.0*factor_h(X0);
				N[0]= dxi/3.0*(3.0*factor_f(X0)-16.*factor_h(X0)*φc0); K[0]=P[0]=0.0;}
			//calculate intermediate positions using previous shift vectors
			xiC= xi[X]+ B[a]*dxi;
			XC = x[X] + B[a]*K[a]; // x(xi)
			UC = u[X] + B[a]*M[a]; // dphi/dxi
			VC = v[X] + B[a]*N[a]; // dpsi/dxi
			ΦC = φ[X] + B[a]*P[a]; // phi(xi)
			YC = sqrt(1.0 + XC*XC);// y(x)
			FC = factor_f(XC);     // f(x)
			HC = factor_h(XC);     // h(x)
		}
		xi[X+1]= xi[X]+ dxi;
		x[X+1] = x[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		u[X+1] = u[X] + M[0]/6.0 + M[1]/3.0 + M[2]/3.0 + M[3]/6.0;
		v[X+1] = v[X] + N[0]/6.0 + N[1]/3.0 + N[2]/3.0 + N[3]/6.0;
		φ[X+1] = φ[X] + P[0]/6.0 + P[1]/3.0 + P[2]/3.0 + P[3]/6.0;
	}
	return Len;
}

//Here we define functions to access radius, pressure, etc.
double PNChandrasekharWD::rad(int X){
	return Rn*xi[X];
}
double PNChandrasekharWD::rho(int X){
	return B0*h[X];
}
double PNChandrasekharWD::drhodr(int X){
	return 3.*B0*x[X]*x[X]*(1. + 8.*sigma*(y[X]-1))*(dx[X]/Rn);
}
double PNChandrasekharWD::P(int X){
	return A0*f[X];
}
double PNChandrasekharWD::dPdr(int X){
	return 8.*A0*pow(x[X],4)/y[X]*(dx[X]/Rn);
}
double PNChandrasekharWD::Phi(int X){
	return 8.*A0/B0*φ[X];
}
double PNChandrasekharWD::dPhidr(int X){
	return 8.*A0/B0*(u[X]/Rn);
}
double PNChandrasekharWD::mr(int X){
	//m = 4 pi Int[r^2 rho,r] 
	//  = 4 pi Rn^3 B Int[xi^2 h(x),xi] 
	//  = 4 pi Rn^3 B xi^2 dphi/dxi
	return 4.*m_pi*B0*pow(Rn,3)*u[X]*pow(xi[X],2);
}
double PNChandrasekharWD::Psi(int){
	return 0.0; //TO BE CONTINUED...
}
double PNChandrasekharWD::dPsidr(int X){
	return 8.*A0*A0/B0/B0*(v[X]/Rn);
}
//the gravitomagnetic potential vansihes for non-rotating star
double PNChandrasekharWD::Wx(int X){return 0.0;}
double PNChandrasekharWD::dWxdr(int ){return 0.0;}
double PNChandrasekharWD::Wy(int ){return 0.0;}
double PNChandrasekharWD::dWydr(int ){return 0.0;}
double PNChandrasekharWD::Wz(int ){return 0.0;}
double PNChandrasekharWD::dWzdr(int ){return 0.0;}

//functions for various thermodynamic quantities of degenerate electron stars
double PNChandrasekharWD::factor_f(double x){
	return x*(2.*x*x-3.)*sqrt(1.+x*x) + 3.*asinh(x);
}
double PNChandrasekharWD::factor_g(double x){
	return 8.0*pow(x,3)*(sqrt(1.+x*x)-1.) - factor_f(x);
}
double PNChandrasekharWD::factor_h(double x){
	return pow(x,3) + sigma*factor_g(x);
}

double PNChandrasekharWD::Gamma1(int X){
	return 8./3.*(x[X]*x[X]*h[X])/f[X]/y[X]/(1.+8.*sigma*(y[X]-1));   
}

double PNChandrasekharWD::sound_speed2(int X, double GamPert){
	// vs2 = Gamma1 * P/(rho+P/c^2) 
	//     = Gamma1* P0 y^n+1 /(rho0 y^n + P0 y^n+1/c2)
	//     = Gamma1*P0*y/(rho0 + P0 y/c^2)
	if(GamPert==0.0) return Gamma1(X)*A0/B0*f[X]/(h[X]+sigma*f[X]);
	else             return GamPert  *A0/B0*f[X]/(h[X]+sigma*f[X]);
}

//Schwarzschild discriminant as for GR case
double PNChandrasekharWD::Schwarzschild_A(int X, double GamPert){
	if(GamPert==0.0) return 0.0;
	else        	 return ( -8.*pow(x[X],4)/(f[X]*y[X]*GamPert)
							  +3.*pow(x[X],2)*( 1.+8.*sigma*(y[X]-1.) )/(h[X]+sigma*f[X]) )*(dx[X]/Rn);	
}

double PNChandrasekharWD::getAstar(int X, double GamPert){
	if(GamPert==0.0) return 0.0;
	else        	 return ( 8.*pow(x[X],4)/(f[X]*y[X]*GamPert)
							  -3.*pow(x[X],2)*( 1.+8.*sigma*(y[X]-1.) )/(h[X]+sigma*f[X]) )*(xi[X]*dx[X]);	
}

double PNChandrasekharWD::getU(int X){
	if(X==0) return 3.0;
	return xi[X]*(3.*sigma*f[X] + h[X]*(1. - 16.*sigma*φ[X]))/(u[X]+sigma*v[X]);
}

double PNChandrasekharWD::getVg(int X, double GamPert){
	if(GamPert==0.0) return -8.*xi[X]*pow(x[X],4)/( f[X]*y[X]*Gamma1(X) )*dx[X];
	else			 return -8.*xi[X]*pow(x[X],4)/( f[X]*y[X]*GamPert   )*dx[X];
}

double PNChandrasekharWD::getC(int X){
	if(X==0) return 0.5*u[len-1]*xi[len-1]/(φc[1]+sigma*ψc[1]);
	return u[len-1]/(u[X]+sigma*v[X])*xi[X]/xi[len-1];
}

double PNChandrasekharWD::Radius(){return rad(len-1);}	//total radius
double PNChandrasekharWD::Mass(){return mr(len-1);}//total mass

double PNChandrasekharWD::SSR(){
	double checkEuler=0.0, checkPoiss=0.0, checkPsi=0.0;
		
	//sum up errors in equations
	double r = 0.0;
	double R = Radius();
	double cee2 = light_speed2();
	double d2Phi=0.0;
	double d2Psi=0.0;
	double e1, e2, e3, n1, n2, n3;
	char txtname[100], outname[100];
	for(int X=6; X<len-7; X++){
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
		e1 = fabs( 4.*m_pi*G_CGS*rho(X)*r - 2.*dPhidr(X) - d2Phi*r );
		n1 = fabs( 4.*m_pi*G_CGS*rho(X)*r ) + fabs( 2.*dPhidr(X) ) + fabs(d2Phi*r);
		//Euler equation	
		e2 = fabs( dPdr(X)
				+ dPhidr(X)*(rho(X)+P(X)/cee2)
				+ rho(X)*dPsidr(X)/cee2 );
		n2 = fabs( dPdr(X) )
				+ fabs(dPhidr(X)*(rho(X)+P(X)/cee2))
				+ fabs(rho(X)*dPsidr(X)/cee2);
		// Equation for Psi
		b3=dPsidr(X-3);
		b2=dPsidr(X-2);
		b1=dPsidr(X-1);
		a1=dPsidr(X+1);
		a2=dPsidr(X+2);
		a3=dPsidr(X+3);
		d2Psi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h);
		e3 = fabs( 12.*m_pi*G_CGS*P(X) - 8.*m_pi*G_CGS*rho(X)*Phi(X) - 2.*dPsidr(X)/r - d2Psi ); 
		n3 = fabs( 12.*m_pi*G_CGS*P(X) )
				+ fabs( 8.*m_pi*G_CGS*rho(X)*Phi(X) )
				+ fabs( 2.*dPsidr(X)/r) + fabs(d2Psi);
		e1 /= n1;
		e2 /= n2;
		e3 /= n3;
		checkPoiss += e1*e1;
		checkEuler += e2*e2;
		checkPsi   += e3*e3;	
	}
	
	return sqrt((checkPoiss + checkEuler + checkPsi)/double(3*len));
}

// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
void PNChandrasekharWD::setupCenter(){
	//the central coefficients -- expanded in terms of xi
	double xi1  = xi[len-1];
	double xi12 = xi1*xi1;
	X02 = X0*X0;
	Y02 = 1.0 + X0*X0;
	Y0 = sqrt(1.0 + X0*X0);
	//first order
	xc[0] = X0;
	yc[0] = Y0;
	fc[0] = factor_f(X0);
	hc[0] = factor_h(X0);
	φc[0] = φc0;
	ψc[0] = 0.0;
	//second order
	xc[1] = xi12*(4.*hc[0]*yc[0]*sigma*(4.*hc[0]*φc[0]-fc[0]) - hc[0]*hc[0]*yc[0])/6.*pow(X0,-4);
	yc[1] = X0*xc[1]/Y0;
	fc[1] = 8.*(xc[1]/yc[0] - yc[0]*xc[1] + xc[0]*xc[0]*xc[1]*yc[0]);
	hc[1] = 3.*X02*xc[1] + 24.*sigma*X02*xc[1]*(Y0-1.);
	φc[1] = xi12/6.*hc[0];
	ψc[1] = xi12/6.*(3.*fc[0] - 16.*hc[0]*φc[0]);
	//fourth order
	φc[2] = xi12/20.*hc[1];
	ψc[2] = xi12/20.*(3.*fc[1]- 16.*hc[1]*φc[0] - 16.*hc[0]*φc[1]);
}
void PNChandrasekharWD::getAstarCenter(double *AC, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	//depending on power requested, return appropriate number of terms
	if(maxPow>=0) AC[0] = 0.0;
	if(maxPow>=2) AC[1] = 16.*X02*X02*xc[1]/(Gam1*fc[0]*Y0)
						     - 6.*X02*xc[1]*(1.+8.*sigma*(Y0-1))/(hc[0]+sigma*fc[0]);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 2) maxPow = 2;
}
void PNChandrasekharWD::getVgCenter(   double *Vc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double x12 = xi[len-1]*xi[len-1];
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] =-16.*X02*X02*xc[1]/(Gam1*fc[0]*Y0);
	if(maxPow> 2) maxPow = 2; 
}
void PNChandrasekharWD::getUCenter(    double *Uc, int& maxPow){
	double x12 = xi[len-1]*xi[len-1];
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] = 6.*hc[1]/(5.*hc[0])
		- 2.*sigma*(8.*pow(hc[0],3)*x12 - 9.*hc[0]*fc[1] + 9.*fc[0]*hc[1])/(5.*hc[0]*hc[0]);
	if(maxPow> 2) maxPow = 2;
}
void PNChandrasekharWD::getC1Center(   double *cc, int& maxPow){
	double u1x1 = u[len-1]*xi[len-1];
	if(maxPow>=0) cc[0] = 0.5*u1x1/(φc[1]+sigma*ψc[1]);
	if(maxPow>=2) cc[1] =    -u1x1*(φc[2]+sigma*ψc[2])*pow(φc[1]+sigma*ψc[1],-2);
	if(maxPow> 2) maxPow = 2;
}
void PNChandrasekharWD::getBetaCenter( double *bc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double phi1 = φ[len-1];
	if(maxPow>=0) bc[0] = Gam1* fc[0]/(hc[0]+sigma*fc[0])/8./phi1;
	if(maxPow>=2) bc[1] = Gam1*(fc[1]*hc[0]-fc[0]*hc[1])*pow(hc[0]+sigma*fc[0],-2)/8./phi1;
	//second order not appear in equations
	if(maxPow> 2) maxPow = 2; 
}
void PNChandrasekharWD::getPhiCenter(  double *pc, int& maxPow){
	if(maxPow>=0) pc[0] = φc[0]/φ[len-1];
	if(maxPow>=2) pc[1] = φc[1]/φ[len-1];//this order does not appear in equations
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
void PNChandrasekharWD::setupSurface(){
	double xi1 =  xi[len-1];
	double xs1 = -dx[len-1]*xi1;
	double φs1 = - u[len-1]*xi1;
	double ψs1 = - v[len-1]*xi1;
	//
	ψs[0] = 0.0;
	φs[0] = -zsurf/8./sigma;
	for(int i=1; i<4;i++) {
		φs[i] = φs1;
		ψs[i] = ψs1;
	}
	//
	xs[0] = 0.0;
	xs[1] = xs1;
	xs[2] = 2.*xs1 + pow(xs1,3)/(φs[1]+sigma*ψs[1]);
	xs[3] = 3.*xs1 + 3.*pow(xs1,5)/pow(φs[1]+sigma*ψs[1],2)
			+ pow(xs1,3)*(60.+5.*(φs[1]+sigma*ψs[1]) + 40.*sigma*φs[1])/10./(φs[1]+sigma*ψs[1]);
	//
	ys[0] = 1.0;
	ys[1] = 0.0;
	ys[2] = xs[1]*xs[1]/2.;
	ys[3] = xs[1]*xs[2];
}

void PNChandrasekharWD::getAstarSurface(double *As, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-1) : g);
	double AV = 3. - 5./Gam1;
	int O=1;
	if(maxPow>=-1) As[O-1] = AV;
	if(maxPow>= 0) As[O  ] = AV*(xs[2]-xs[1])/xs[1];
	if(maxPow>= 1) As[O+1] = AV*( -xs[1]*xs[1]/7. - pow(xs[2]/xs[1],2) - (xs[2]-2.*xs[3])/xs[1] )
							  + xs[1]*xs[1]*(3./7.);
	if(maxPow > 1) maxPow = O+1;
}
void PNChandrasekharWD::getVgSurface(double *Vs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-1) : g);
	double VG = 5./Gam1;
	int O=1;
	if(maxPow>=-1) Vs[O-1] = VG;
	if(maxPow>= 0) Vs[O  ] = VG*(xs[2]-xs[1])/xs[1];
	if(maxPow>= 1) Vs[O+1] =-5./Gam1*( 
			pow(xs[1],4) + 7.*xs[2]*xs[2] + 7.*xs[2]*xs[1] - 14.*xs[3]*xs[1]
		)/(7.*pow(xs[1],2));
	if(maxPow > 1) maxPow = O+1;
}
void PNChandrasekharWD::getUSurface(double *Us, int& maxPow){
	if(maxPow >=0) Us[0] = 0.0;
	if(maxPow >=1) Us[1] = 0.0;
	if(maxPow >=2) Us[2] = 0.0;
	if(maxPow > 2) maxPow = 2;
}
void PNChandrasekharWD::getC1Surface(double *cs, int& maxPow){
	double x1u1 = xi[len-1]*u[len-1];
	if(maxPow>=0) cs[0]  =  - x1u1/(φs[1]+sigma*ψs[1]);
	if(maxPow>=1) cs[1]  = 3.*x1u1/(φs[1]+sigma*ψs[1]);
	if(maxPow>=2) cs[2]  =-3.*x1u1/(φs[1]+sigma*ψs[1]);
	if(maxPow> 2) maxPow = 2;
}

void PNChandrasekharWD::getBetaSurface(double *bs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-1) : g);
	if(maxPow>=0) bs[0]  = 0.0;
	if(maxPow>=1) bs[1]  = 0.0;
	if(maxPow>=2) bs[2]  = 0.0;
	if(maxPow> 2) maxPow = 2;
}
void PNChandrasekharWD::getPhiSurface(double *ps, int& maxPow){
	double phi1 = -zsurf/8./sigma;
	if(maxPow>=0) ps[0]  = 1.;
	if(maxPow>=1) ps[1]  = φs[1]/φs[0];
	if(maxPow>=2) ps[2]  = φs[2]/φs[0];
	//the pattern would continue for higher terms...
	if(maxPow> 2) maxPow = 2;
}


//method to print pertinent values of star to .txt, and plot them in gnuplot
void PNChandrasekharWD::writeStar(char *c){
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
	double irc=1./rho(0), ipc=1./P(0), R=Radius(), ig=1./dPhidr(length()-1), iPN=1./dPsidr(length()-1);
	double iq = 1./(dPhidr(length()-1)+dPsidr(length()-1)), MT = Mass();
	for(int X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le",
			rad(X)/R, rho(X)*irc, -drhodr(X)*irc*R,
			P(X)*ipc, -dPdr(X)*ipc*R,
			mr(X)/MT, (dPhidr(X)+dPsidr(X))*iq);
		fprintf(fp, "\t%le\t%le\t%le", x[X], y[X], f[X]);
		fprintf(fp, "\t%le", φ[X]);
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
	fprintf(gnuplot, ", '%s' u 1:4 w l t 'P'", txtname);
	fprintf(gnuplot, ", '%s' u 1:6 w l t 'm'", txtname);
	fprintf(gnuplot, ", '%s' u 1:7 w l t 'g'", txtname);
	fprintf(gnuplot, "\n"); 
	
	//print fthe phi
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s/%s.png'\n", filename, "phi");
	fprintf(gnuplot, "set title 'Newtonian Potential for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'Phi'\n");
	fprintf(gnuplot, "plot ");
	fprintf(gnuplot, " '%s' u 1:11 w l t 'Phi'", txtname);
	fprintf(gnuplot, "\n");
	
	//print the pulsation coeffcients frequency
	sprintf(txtname, "%s/coefficients.txt", filename);
	sprintf(outname, "%s/coefficients.png", filename);
	fp  = fopen(txtname, "w");
	int maxpow=2;
	double A,U,V,C, read[2], x = xi[1]/xi[len-1];
	getAstarCenter(read, maxpow, 5./3.);
	A = read[0] + read[1]*x*x;
	getVgCenter(read, maxpow, 5./3.);
	V = read[0] + read[1]*x*x;
	getC1Center(read, maxpow);
	C = read[0] + read[1]*x*x;
	getUCenter(read, maxpow);
	U = read[0] + read[1]*x*x;
	fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
		xi[0], A, U, V, C);
	for(int X=1; X< length()-1; X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			xi[X],
			-xi[X]*Rn*Schwarzschild_A(X, 5./3.),
			getU(X),
			getVg(X, 5./3.),
			getC(X)
		);
	}
	double reads[3], t1 = 1.-xi[len-2]/xi[len-1];
	getAstarSurface(reads, maxpow, 5./3.);
	A = reads[0]/t1 + reads[1] + reads[2]*t1;
	getVgSurface(reads, maxpow, 5./3.);
	V = reads[0]/t1 + reads[1] + reads[2]*t1;
	getC1Surface(reads, maxpow);
	C = reads[0] + reads[1]*t1 + reads[2]*t1*t1;
	getUSurface(reads, maxpow);
	U = reads[0] + reads[1]*t1 + reads[2]*t1*t1;
	fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
		xi[len-1], A, U, V, C);
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
	
	//print the Brunt-Vaisala frequency
	sprintf(txtname, "%s/BruntVaisala.txt", filename);
	sprintf(outname, "%s/BruntVaisala.png", filename);
	fp  = fopen(txtname, "w");
	double N2 = -1.0;
	for(int X=1; X< length()-1; X++){
		//N^2 = -g*A
		N2 = -dPhidr(X)*Schwarzschild_A(X,5./3.);
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\n",
			P(X),
			N2,
			2.*sound_speed2(X,0.)*pow(Rn*xi[X],-2));
	}
	fclose(fp);	
	//plot file in png in gnuplot
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 800,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Brunt-Vaisala for %s'\n", title);
	fprintf(gnuplot, "set logscale x 10\n");
	fprintf(gnuplot, "set format x '10^{%%L}'\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xlabel 'log_{10} P\n");
	fprintf(gnuplot, "set ylabel 'log_{10} N^2 & log_{10} L_1^2 (Hz^2)\n");
	fprintf(gnuplot, "set yrange [1e-6:1e2]\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'N^2'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'L_1^2'", txtname);
	fprintf(gnuplot, "\n");
	
	//
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);	
}

#endif