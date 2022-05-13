//**************************************************************************************
//						COWLING PULSATION nlm MODE DRIVER
// Cowling ModeDriver.cpp
//		Solves the Newtonian LAWE in Dziembowski variables in the Cowling approximation
//		Does not consider perturbtations to gravitational field, making 2nd-order modes
//***************************************************************************************

#ifndef COWLINGMODEDRIVERCLASS
#define COWLINGMODEDRIVERCLASS

#include "CowlingModeDriver.h"

//constructor
CowlingModeDriver::CowlingModeDriver(Star *star, double Gamma1)
	 : ModeDriver(num_var, star), adiabatic_index(Gamma1)
{	
	len_star = star->length();
	len = (len_star%2==1? (len_star+1)/2 : len_star/2+1);
	//set the arrays that define the perturbation
	double R = star->Radius();
	r = new double[len];
	//define values of radius for the mode
	for(int X=0; X<len-1; X++){
		r[X] = star->rad(2*X)/R;
	}
	r[len-1] = 1.0;
	//initialize helper arrays
	initializeArrays();
	//initialize the boundary
	central_bc_order = BC_C;
	surface_bc_order = BC_S;
	setupBoundaries();
}

//destructor
CowlingModeDriver::~CowlingModeDriver(){
	delete[] r;
	delete[] A; delete[] U;
	delete[] C; delete[] V;
}

int CowlingModeDriver::length(){
	return len;
}

double CowlingModeDriver::Gamma1(){
	return adiabatic_index;
}

double CowlingModeDriver::rad(int X){
	return (r[X]);
}

//recalculate a11, a12, a21, a22
void CowlingModeDriver::initializeArrays(){
	A = new double[len_star]; U = new double[len_star];
	C = new double[len_star]; V = new double[len_star];
	A[0] = 0.0;
	V[0] = 0.0;
	U[0] = 3.0;
	C[0] = star->getC(0);
	int X;
	for(X=1; X<len_star-1; X++){
		A[X] = star->getAstar(X, adiabatic_index);
		U[X] = star->getU(X);
		C[X] = star->getC(X);
		V[X] = star->getVg(X, adiabatic_index);
	}
	//surface requires special treatment
	X = len_star-1;
	V[X] = 0.0;
	U[X] = 0.0;
	A[X] = 0.0;  
	C[X] = 1.0;
}

//these functions come from coupled wave equations -- see Unno et al. Ch 18
void CowlingModeDriver::getCoeff(
	double *CCI,
	int X,
	int b,
	double sig2,
	int L
)
{	
	double (*CC)[num_var] =(double (*)[num_var]) CCI;
	//avoid recalculations
	double Ax, Ux, Cx, Vx;
	Ax = A[2*X+b];
	Ux = U[2*X+b];
	Cx = C[2*X+b];
	Vx = V[2*X+b];
	//dy1/dx
	CC[0][0] = Vx - double(L+1);
	CC[0][1] = double(L*L+L)/(sig2*Cx) - Vx;
	//dy2/dx
	CC[1][0] = sig2*Cx - Ax;
	CC[1][1] = double(3-L) + Ax - Ux;
}


void CowlingModeDriver::getBoundaryMatrix(int nv, double *y0, double* ys, double **y, int* indexOrder){
	const double yy[num_var][num_var] = {
		{ 1.0, 1.0 },
		{ 1.0, 1.0 }
	};
	const int indices[num_var] = {0,0};
	for(int i=0; i<num_var; i++)
		for(int j=0; j<num_var; j++)
			y[i][j] = yy[i][j];
	for(int i=0; i<num_var; i++) indexOrder[i] = indices[i];
}


void CowlingModeDriver::setupBoundaries() {
	//set up coefficients for the central boundary
	star->getAstarCenter(Ac, central_bc_order, adiabatic_index);
	star->getVgCenter(Vc, central_bc_order, adiabatic_index);
	star->getUCenter(Uc, central_bc_order);
	star->getC1Center(cc, central_bc_order);
	cProdc[0] = 1./cc[0];
	cProdc[1] = -cc[1]*pow(cc[0],-2);
	cProdc[2] = (cc[1]*cc[1]-cc[2]*cc[0])*pow(cc[0],-3);
	
	//set up coefficients for the surface boundary
	star->getAstarSurface(As, surface_bc_order, adiabatic_index);
	star->getVgSurface(Vs, surface_bc_order, adiabatic_index);
	star->getUSurface(Us, surface_bc_order);
	star->getC1Surface(cs, surface_bc_order);
	cProds[0] =  1.0;
	cProds[1] = -cs[1];
	cProds[2] =  cs[1]*cs[1]-cs[2];
	cProds[3] = -pow(cs[1],3) + 2.*cs[1]*cs[2] - cs[3];
	cProds[4] = 0.0;
	k_surface = 0.0;
}

int CowlingModeDriver::CentralBC(double **ymode, double *y0, double omeg2, int l, int m){	
	double yy[num_var][central_bc_order/2+1];//0,2,4
	//the zero-order terms are simple
	yy[0][0] = y0[0];
	yy[1][0] = y0[0]*C[0]*omeg2/double(l);
	
	double L2 = double(l*l+l);
	double Gam1 = adiabatic_index;
	
	for(int i=0; i<num_var; i++) ymode[i][1] = ymode[i][0] = yy[i][0];
	
	double sumU, sumC, sumG, sumV, sumA;
	double trip[3]; trip[0] = yy[0][0]-yy[1][0];
	for(int j=1;j<=central_bc_order/2;j++){
		sumA = sumV = sumG = sumC = sumU = 0.0;
		for(int k=0; k<j; k++){
			sumV += Vc[j-k] *trip[k];
			sumA += Ac[j-k] *trip[k];
			sumC += cc[j-k]*yy[0][k];
			sumU += Uc[j-k]*yy[1][k];
			sumG += L2*cProdc[j-k]*yy[1][k]/omeg2;
		}
		yy[0][j] = ( double(2*j+l)*(sumG+sumV) - L2*(sumU+sumA+omeg2*sumC)/(omeg2*cc[0]) )/double(2*j*(1+2*j+2*l));
		yy[1][j] = (omeg2*cc[0]*yy[0][j] - omeg2*sumC - sumU - sumA)/double(2*j+l);

		trip[j] = yy[0][j]-yy[1][j];
	}
	//noew solve for y using expansion
	int start = 1;
	for(int X=0; X<=start; X++){
		for(int i=0; i<num_var; i++){
			ymode[i][X] = yy[i][0];
			for(int j=1;j<=(central_bc_order/2);j++) ymode[i][X] += yy[i][j]*pow(r[X],2*j);
		}
	}
	return start;
}

int CowlingModeDriver::SurfaceBC(double **ymode, double *ys, double omeg2, int l, int m){
	//specify initial conditions at surface
	double yy[num_var][surface_bc_order+1];	//coefficients y = yy[0] + yy[1]t + ... + yy[k]t^k
	double yyn[num_var]; for(int i=0;i<num_var;i++) yyn[i]=0.0;
	yy[0][0] = ys[0];
	yy[1][0] = ys[0];
	
	//constants that show up
	double L2 = double(l*l+l);
	
	int O=1;
	double trip[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	double sumV, sumA, sumU, sumC, sumG, sumUX, sumUY, sumUZ;
	for(int k=0; k<surface_bc_order; k++){
		//calculate "triples" (W[k]-X[k]+Y[k]) which appear in recursion relation
		trip[k] = yy[0][k] - yy[1][k];
		sumV=sumA=sumC=sumU=sumG=sumUX=sumUY=sumUZ = 0.0;
		for(int j=0; j<=k; j++){
			sumV += Vs[O+k-j]*trip[j];
			sumA += As[O+k-j]*trip[j];
			sumC += omeg2*cs[k-j]*yy[0][j];
			sumG += L2/omeg2*cProds[k-j]*yy[1][j];
			sumU += Us[k-j]*yy[1][j];
		}
		trip[k+1] = (sumC-sumG-sumU  - sumV-sumA
					+ double(k+l-3)*trip[k] + 4.*yy[0][k])/(double(k+1)+Vs[O-1]+As[O-1]);
		sumV += Vs[O-1]*trip[k+1];
		sumA += As[O-1]*trip[k+1];
		yy[0][k+1] = (double(k+l+1)*yy[0][k] - sumG - sumV)/double(k+1);
		yy[1][k+1] = (double(k+l-3)*yy[1][k] - sumC + sumA + sumU)/double(k+1);
	}
	//the number of terms to calculate
	int start = len-2;
	double t;
	for(int i=0; i<num_var; i++) ymode[i][len-1] = yy[i][0];
	for(int X=len-2; X>=start; X--){
		t = (1.-r[X]);
		for(int i=0; i<num_var; i++){
			ymode[i][X] = yy[i][0];
			for(int k=1; k<=surface_bc_order; k++) ymode[i][X] += yy[i][k]*pow(t,k);
		}
	}
	return start;
}

// *** this method is assuming a uniform grid ***
// *** we cannot assume a unifirm grid for all stars
double CowlingModeDriver::SSR(double omeg2, int l, ModeBase* mode){
	double checkCont=0.0;
	double checkNewt=0.0;
	double e1,e2;
	double n1,n2;
	double difxi, dDP;
	double rho,P,drhodr,dPdr,g,r,G1,freq2,L2,xi,chi,DPhi,dDPhi,Drho;
	double b3,b2,b1,a1,a2,a3, h1;
	freq2 = omeg2;//sigma2*star->mr(len_star-1)*pow(star->rad(len_star-1),-3);
	for(int X=4; X<len-4; X++){
		int XX = 2*X;
		//stellar variables, to simplify equations
		rho    = star->rho(XX);
		P      = star->P(XX);
		drhodr = star->drhodr(XX);
		dPdr   = star->dPdr(XX);
		g      = star->dPhidr(XX);
		r      = star->rad(XX);
		G1     = (Gamma1()==0.0? star->Gamma1(XX) : Gamma1());
		//mode variables	
		L2     = double(l*l+l);
		xi   = r  *mode->getY(0,X);	//r  *y1
		chi  = r*g*mode->getY(1,X);	//r*g*y2
		Drho = rho*rho/(G1*P)*(chi - g*xi) - xi*drhodr;
		
		//calculate numerical derivatives
		//double b3,b2,b1,a1,a2,a3, h1;
		h1= star->rad(XX+2)-star->rad(XX );
		//next dxi/dr, xi = r y
		b3=star->rad(XX-6)*mode->getY(0,X-3);
		b2=star->rad(XX-4)*mode->getY(0,X-2);
		b1=star->rad(XX-2)*mode->getY(0,X-1);
		a1=star->rad(XX+2)*mode->getY(0,X+1);
		a2=star->rad(XX+4)*mode->getY(0,X+2);
		a3=star->rad(XX+6)*mode->getY(0,X+3);
		difxi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
		//now dDP/dr, DP = rho*(chi-DPhi), chi=rg*y2, DPhi=rg*y3
		b3=star->rho(XX-6)*star->dPhidr(XX-6)*star->rad(XX-6)*(mode->getY(1,X-3));
		b2=star->rho(XX-4)*star->dPhidr(XX-4)*star->rad(XX-4)*(mode->getY(1,X-2));
		b1=star->rho(XX-2)*star->dPhidr(XX-2)*star->rad(XX-2)*(mode->getY(1,X-1));
		a1=star->rho(XX+2)*star->dPhidr(XX+2)*star->rad(XX+2)*(mode->getY(1,X+1));
		a2=star->rho(XX+4)*star->dPhidr(XX+4)*star->rad(XX+4)*(mode->getY(1,X+2));
		a3=star->rho(XX+6)*star->dPhidr(XX+6)*star->rad(XX+6)*(mode->getY(1,X+3));
		dDP = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
		//Now calculate residuals
		//continuity equation
		e1 = fabs(
				Drho + xi*drhodr + rho*( 2.*xi/r + difxi - L2*pow(r,-2.0)/freq2*chi )
		);
		n1 = fabs(Drho) + fabs(xi*drhodr) + fabs(2.*rho/r*xi)
				+ fabs(rho*difxi)
				+ fabs(rho*L2*chi*pow(r,-2.0)/freq2);
		//newton's equation -- the r component. The theta component defines chi = r sigma2 xiH
		e2 = fabs( rho*freq2*xi - g*Drho  - dDP );
		n2 = fabs( rho*freq2*xi ) + fabs( g*Drho ) + fabs( dDP );
		//normalize residuals
		e1 = e1/n1;
		e2 = e2/n2;
		//collect residuals
		checkCont += e1*e1;
		checkNewt += e2*e2;
	}
	return sqrt((checkCont+checkNewt)/double(2*len-7));
}


double CowlingModeDriver::tidal_overlap(ModeBase* mode){
	double omega2 = mode->getOmega2();
	int k,l,m;
	mode->modeNumbers(k,l,m);
	double L2 = double(l*(l+1));
	double Mstar=star->Mass(), Rstar=star->Radius();
	double dscale = Mstar*pow(Rstar,-3); //dimensionless, set equal to M/R^3
	double rscale = 1.0;//Rstar;		 //dimensionless, set equal to 1
	double freq2 = omega2*star->Gee()*Mstar*pow(Rstar,-3);
	double xir=0.0, xiH=0.0, rho=0.0, rlp1=0.0;
	double dx;
	int xx;
	//need to numerically integrate
	double sum1=0.0, sum2=0.0, integral = 0.0;
	double N1=0.0, N2=0.0, NN=0.0;
	for(int x=1; x<len; x++){
		xx = 2*x;
		if(x==len-1) xx = len_star-1;
		dx = (r[x]-r[x-1])*rscale;
		xir = mode->getY(0,x)*r[x]*rscale;
		xiH = mode->getY(1,x)*star->dPhidr(xx)/freq2*rscale/Rstar;
		//
		rho  = star->rho(xx)/dscale;
		rlp1 = pow(r[x]*rscale,l+1);
		//the integral using trapezoidal rule
		sum2 = sum1;
		sum1 = rho*rlp1*(double(l)*xir + L2*xiH);
		integral += 0.5*dx*(sum1+sum2);
		//integrate normalization constant
		N2 = N1;
		N1 = rho*pow(r[x]*rscale,2)*(xir*xir+L2*xiH*xiH);
		NN += 0.5*dx*(N1+N2);
	}
	return integral/sqrt(NN);
}

double CowlingModeDriver::innerproduct(ModeBase* mode1, ModeBase* mode2){
	double omega21 = mode1->getOmega2();
	double omega22 = mode2->getOmega2();
	int k1,l1,m1, k2,l2,m2;
	mode1->modeNumbers(k1,l1,m1);
	mode2->modeNumbers(k2,l2,m2);
	if(l2!=l1) return 0.0;
	double L2 = double(l1*l1+l1);
	double Mstar=star->Mass(), Rstar=star->Radius();
	double dscale = Mstar*pow(Rstar,-3);
	double rscale = 1.0;
	double freq21 = omega21*star->Gee()*Mstar*pow(Rstar,-3);
	double freq22 = omega22*star->Gee()*Mstar*pow(Rstar,-3);
	double xir1=0.0, xiH1=0.0, rho=0.0;
	double xir2=0.0, xiH2=0.0;
	double dx;
	int xx;
	//need to numerically integrate
	double sum1=0.0, sum2=0.0, integral = 0.0;
	double N1=0.0, N2=0.0, NN=0.0;
	double M1=0.0, M2=0.0, MM=0.0;
	for(int x=1; x<len; x++){
		xx = 2*x;
		if(x==len-1) xx = len_star-1;
		dx = (r[x]-r[x-1])*rscale;
		xir1 = mode1->getY(0,x)*r[x]*rscale;
		xir2 = mode2->getY(0,x)*r[x]*rscale;
		xiH1 = mode1->getY(1,x)*star->dPhidr(xx)/freq21*rscale/Rstar;
		xiH2 = mode2->getY(1,x)*star->dPhidr(xx)/freq22*rscale/Rstar;
		//
		rho  = star->rho(xx)/dscale;
		//the integral using trapezoidal rule
		sum2 = sum1;
		sum1 = rho*pow(r[x]*rscale,2)*(xir1*xir2 + L2*xiH1*xiH2);
		integral += 0.5*dx*(sum1+sum2);
		//integrate normalization constants
		N2 = N1;
		N1 = rho*pow(r[x]*rscale,2)*(xir1*xir1+L2*xiH1*xiH1);
		NN += 0.5*dx*(N1+N2);
		M2 = M1;
		M1 = rho*pow(r[x]*rscale,2)*(xir2*xir2+L2*xiH2*xiH2);
		MM += 0.5*dx*(M1+M2);
	}
	return integral/sqrt(NN*MM);
}



#endif