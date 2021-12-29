//**************************************************************************************
//						NONRADIAL PULSATION nlm MODE DRIVER
// NonradialModeDriver.cpp
//		Solves the Newtonian LAWE in Dziembowski variables
//		Includes perturbtations to gravitational field, making 4th-order modes
//***************************************************************************************

#ifndef FULLMODEDRIVERCLASS
#define FULLMODEDRIVERCLASS

#include "NonradialModeDriver.h"

//constructor
NonradialModeDriver::NonradialModeDriver(Star *star, double Gamma1)
	 : ModeDriver(num_var, star), adiabatic_index(Gamma1)
{	
	len_star = star->length();
	len = (len_star%2==1? (len_star+1)/2 : len_star/2+1);
	//set the arrays that define the perturbation
	double R = star->rad(len_star-1);
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
NonradialModeDriver::~NonradialModeDriver(){
	delete[] r;
	delete[] A; delete[] U;
	delete[] C; delete[] V;
}

int NonradialModeDriver::length(){
	return len;
}

double NonradialModeDriver::Gamma1(){
	return adiabatic_index;
}

double NonradialModeDriver::rad(int X){
	return r[X];
}

//recalculate a11, a12, a21, a22
void NonradialModeDriver::initializeArrays(){
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
void NonradialModeDriver::getCoeff(
	double *CCI,
	int X,
	int b,
	double omeg2,
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
	CC[0][1] = double(L*L+L)/(omeg2*Cx) - Vx;
	CC[0][2] = Vx;
	CC[0][3] = 0.0;
	//dy2/dx
	CC[1][0] = omeg2*Cx - Ax;
	CC[1][1] = double(3-L) + Ax - Ux;
	CC[1][2] =-Ax;
	CC[1][3] = 0.0;
	//dy3/dx
	CC[2][0] = 0.0;
	CC[2][1] = 0.0;
	CC[2][2] = double(3-L) - Ux;
	CC[2][3] = 1.0;
	//dy4/dx
	CC[3][0] = Ax*Ux;
	CC[3][1] = Vx*Ux;
	CC[3][2] = double(L*L+L) - Vx*Ux;
	CC[3][3] = double(2-L) - Ux;
}


void NonradialModeDriver::getBoundaryMatrix(int nv, double *y0, double* ys, double **y, int* indexOrder){
	double one = 1.0; //this is used only to provide visual distinction in grid
	//a basis of independent solutions to be used at center and surface
	const double yy[num_var][num_var] = {
		{ one, one, 0.0, 0.0},
		{ 0.0, 0.0, one, one},
		{ 0.0, one, one, one},
		{ one, one, 0.0, 0.0}
	};
	for(int i=0; i<nv; i++)
		for(int j=0; j<nv; j++)
			y[i][j] = yy[i][j];
	//when matching outward/inward, indexOrder matches correct coefficient to correct solution
	const int indices[num_var] = {y1,y3,y3,y1};
	for(int i=0; i<nv; i++) indexOrder[i] = indices[i];
}


void NonradialModeDriver::setupBoundaries() {
	//set up coefficients for the central boundary
	star->getAstarCenter(Ac, central_bc_order, adiabatic_index);
	star->getVgCenter(   Vc, central_bc_order, adiabatic_index);
	star->getUCenter(    Uc, central_bc_order);
	star->getC1Center(   cc, central_bc_order);
	cProdc[0] = 1./cc[0];
	cProdc[1] = -cc[1]*pow(cc[0],-2);
	cProdc[2] = (cc[1]*cc[1]-cc[2]*cc[0])*pow(cc[0],-3);
	
	//set up coefficients for the surface boundary
	star->getAstarSurface(As, surface_bc_order, adiabatic_index);
	star->getVgSurface(   Vs, surface_bc_order, adiabatic_index);
	star->getUSurface(    Us, surface_bc_order);
	star->getC1Surface(   cs, surface_bc_order);
	cProds[0] =  1.0;
	cProds[1] = -cs[1];
	cProds[2] =  cs[1]*cs[1]-cs[2];
	cProds[3] = -pow(cs[1],3) + 2.*cs[1]*cs[2] - cs[3];
	cProds[4] = 0.0;
	k_surface = 0.0;// *** This needs to be re-implemented
//	k_surface = structure_coefficients_surface[offset+indxC];
}


//************************* CENTRAL BOUNDARY CONDITION *********************************
//	Uses a recursion algorithm to solve for power series terms of y_1,...,y_4
//		based on the power series terms of A*,Vg, c1, U from stellar model
//	Limited by number of terms returned by stellar model
//	For the derivation of this algorithm, see Mathematica notebook NewtonianBCs.nb
//	The expansions must be in terms of x = r/R
int NonradialModeDriver::CentralBC(double **ymode, double *y0, double omeg2, int l, int m){
	double yy[num_var][central_bc_order/2+1];//0,2,4
	//the zero-order terms are simple
	yy[y1][0] = y0[0];
	yy[y2][0] = y0[0]*C[0]*omeg2/double(l);
	yy[y3][0] = y0[2];
	yy[y4][0] = y0[2]*double(l);
	
	double L2 = double(l*l+l);
	double Gam1 = adiabatic_index;
	
	for(int i=0; i<num_var; i++) ymode[i][1] = ymode[i][0] = yy[i][0];
	
	double sumU, sumC, sumG, sumV, sumA, sumUY, sumUZ;
	double trip[3]; trip[0] = yy[0][0]-yy[1][0]+yy[2][0];
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
		sumUY = sumUZ = sumU = sumV = 0.0;
		for(int k=0; k<j; k++){
			sumUY += Uc[j-k]*yy[2][k];
			sumUZ += Uc[j-k]*yy[3][k];
			sumV = 0.0;
			for(int m=0; m<k; m++) sumV += Ac[k-m]*yy[0][m] + Vc[k-m]*(yy[1][m]-yy[2][m]);
			sumU += Uc[j-k]*sumV;
		}
		//the sumU requires an extra step
		for(int m=0, sumV=0.0; m<j; m++) sumV += Ac[j-m]*yy[0][m] + Vc[j-m]*(yy[1][m]-yy[2][m]);
		sumU += Uc[0]*sumV;
		//now calculate the coefficients
		yy[2][j] = (double(l+2*j+1)*sumUY+sumUZ-sumU)/double(2*j*(2*l+2*j+1));
		yy[3][j] =  double(l+2*j)*yy[2][j] + sumUY;
		
		trip[j] = yy[0][j]-yy[1][j]+yy[2][j];
	}
	int start = int(1);
	for(int X=0; X<=start; X++){
		for(int i=0; i<num_var; i++){
			ymode[i][X] = yy[i][0];
			for(int j=1;j<=(central_bc_order/2);j++) ymode[i][X] += yy[i][j]*pow(r[X],2*j);
		}
	}
	return start;
}


//************************* SURFACE BOUNDARY CONDITION *********************************
//	Uses a recursion algorithm to solve for power series terms of y_1,...,y_4
//		based on the power series terms of A*,Vg, c1, U from stellar model
//	Limited by number of terms returned by stellar model
//	For the derivation of this algorithm, see Mathematica notebook NewtonianBCs.nb
//	The expansiions must be in terms of t = 1-r/R
int NonradialModeDriver::SurfaceBC(double **ymode, double *ys, double omeg2, int l, int m){
	//int surface_bc_order = 4;
	//specify initial conditions at surface
	double yy[num_var][surface_bc_order+1];	//coefficients y = yy[0] + yy[1]t + ... + yy[k]t^k
	double yyn[num_var]; for(int i=0;i<num_var;i++) yyn[i]=0.0;
	yy[0][0] = ys[0];
	yy[1][0] = ys[0] + ys[2];
	yy[2][0] = ys[2];
	yy[3][0] = -double(l+1)*ys[2];
	
	//constants that show up
	double L2 = double(l*l+l);
	
	int O=1;
	double trip[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	double sumV, sumA, sumUX, sumUY, sumUZ, sumU, sumC, sumG;
	for(int k=0; k<surface_bc_order; k++){
		//calculate "triples" (W[k]-X[k]+Y[k]) which appear in recursion relation
		trip[k] = yy[0][k] - yy[1][k] + yy[2][k];
		sumV=sumA=sumC=sumG=sumUX=sumUY=sumUZ = 0.0;
		for(int j=0; j<=k; j++){
			sumV += Vs[O+k-j]*trip[j];
			sumA += As[O+k-j]*trip[j];
			sumC += omeg2*cs[k-j]*yy[0][j];
			sumG += L2/omeg2*cProds[k-j]*yy[1][j];
			sumUX += Us[k-j]*yy[1][j];
			sumUY += Us[k-j]*yy[2][j];
			sumUZ += Us[k-j]*yy[3][j];
		}
		trip[k+1] = (sumC-sumG-sumUX+sumUY - sumV-sumA
					+ double(k+l-3)*trip[k] + 4.*yy[0][k] - yy[3][k])/(double(k+1)+Vs[O-1]+As[O-1]);
		sumV += Vs[O-1]*trip[k+1];
		sumA += As[O-1]*trip[k+1];
		yy[0][k+1] = (double(k+l+1)*yy[0][k] - sumG - sumV)/double(k+1);
		yy[1][k+1] = (double(k+l-3)*yy[1][k] - sumC + sumA + sumUX)/double(k+1);
		yy[2][k+1] = (double(k+l-3)*yy[2][k] - yy[3][k] + sumUY)/double(k+1);
		sumU=0.0;
		for(int m=-1; m<=k; m++){
			sumV = 0.0;
			for(int j=0; j<=m+1; j++){
				sumV += As[O+m-j]*yy[0][j] + Vs[O+m-j]*(yy[1][j]-yy[2][j]);
			}
			sumU += Us[k-m]*sumV;
		}
		yy[3][k+1] = (double(k+l-2)*yy[3][k] - L2*yy[2][k] - sumU + sumUZ)/double(k+1);
	}
	//if n<1, then an additional term is needed
	if(k_surface != 0.0){
		//find k using the ridiculous formulas A13,A14 from JCD-DJM
		double n=0.0;
		double Gam1 = adiabatic_index;
		double kn = pow(k_surface, n);
		double tripn = -3.*kn*ys[0]/(2.*n+1.);
		yyn[0] = -tripn/Gam1;
		yyn[1] = (3.*kn*ys[1] + (n-(n+1.)/Gam1)*tripn)/(n+1.);
		yyn[2] = (3.*kn*ys[2])/(n+1.);
		yyn[3] =  3.*kn*(ys[3]+(n+1.)/Gam1*trip[1] - n*yy[0][1])/(n+1.);
	}
	//the number of terms to calculate
	int start = len-2;
	double t;
	for(int i=0; i<num_var; i++) ymode[i][len-1] = yy[i][0];
	//calculate the y in terms of power series solution
	for(int X=len-2; X>=start; X--){
		t = (1.-r[X]);
		for(int i=0; i<num_var; i++){
			ymode[i][X] = yy[i][0];
			for(int k=1; k<=surface_bc_order; k++) ymode[i][X] += yy[i][k]*pow(t,k);
			//ymode[i][X] += yyn[i]*pow(t, n+1);
		}
	}
	return start;
}

// *** this method is assuming a uniform grid ***
// *** we cannot assume a unifirm grid for all stars
double NonradialModeDriver::SSR(double omega2, int l, ModeBase* mode){
	double checkCont=0.0;
	double checkNewt=0.0;
	double checkPois=0.0;
	double e1,e2,e3;
	double n1,n2,n3;
	double d2Phi, difxi, dDP;
	double rho,P,drhodr,dPdr,g,r,G1,freq2,L2,xi,chi,DPhi,dDPhi,Drho;
	double b3,b2,b1,a1,a2,a3, h1;
	double Gee = star->Gee();
	double Mstar  = star->Mass();
	double Rstar  = star->Radius();
	double Dscale = Mstar*pow(Rstar,-3);
	double Pscale = Gee*pow(Mstar,2)*pow(Rstar,-4);
	double Gscale = Gee*Mstar*pow(Rstar,-2);
	freq2 = omega2;
	Gee=1.0;
	for(int X=6; X<len-7; X++){
		int XX = 2*X;
		//stellar variables, to simplify equations
		rho    = star->rho(XX)   /Dscale;
		P      = star->P(XX)     /Pscale;
		drhodr = star->drhodr(XX)/Dscale*Rstar;
		dPdr   = star->dPdr(XX)  /Pscale*Rstar;
		g      = star->dPhidr(XX)/Gscale;
		r      = star->rad(XX)   /Rstar;
		G1     = (Gamma1()==0.0? star->Gamma1(XX) : Gamma1());
		//mode variables
		
		L2     = double(l*l+l);
		xi   = r  *mode->getY(y1,X);	//r  *y1
		chi  = r*g*mode->getY(y2,X);	//r*g*y2
		DPhi = r*g*mode->getY(y3,X);	//r*g*y3
		dDPhi=   g*mode->getY(y4,X);	//  g*y4
		Drho = rho*rho/(G1*P)*(chi - DPhi - g*xi) - xi*drhodr;
		
		//calculate numerical derivatives 
		h1= star->rad(XX+2)-star->rad(XX );
		// d2DPhi/dr2 = (d/dr)(dDphi/dr), dDPhi/dr = g y4
		b3=star->dPhidr(XX-6)*mode->getY(y4,X-3);
		b2=star->dPhidr(XX-4)*mode->getY(y4,X-2);
		b1=star->dPhidr(XX-2)*mode->getY(y4,X-1);
		a1=star->dPhidr(XX+2)*mode->getY(y4,X+1);
		a2=star->dPhidr(XX+4)*mode->getY(y4,X+2);
		a3=star->dPhidr(XX+6)*mode->getY(y4,X+3);
		d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1)/Gscale*Rstar;
		//next dxi/dr, xi = r y
		b3=star->rad(XX-6)*mode->getY(y1,X-3);
		b2=star->rad(XX-4)*mode->getY(y1,X-2);
		b1=star->rad(XX-2)*mode->getY(y1,X-1);
		a1=star->rad(XX+2)*mode->getY(y1,X+1);
		a2=star->rad(XX+4)*mode->getY(y1,X+2);
		a3=star->rad(XX+6)*mode->getY(y1,X+3);
		difxi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
		//now dDP/dr, DP = rho*(chi-DPhi), chi=rg*y2, DPhi=rg*y3
		b3=star->rho(XX-6)*star->dPhidr(XX-6)*star->rad(XX-6)*(mode->getY(1,X-3)-mode->getY(y3,X-3));
		b2=star->rho(XX-4)*star->dPhidr(XX-4)*star->rad(XX-4)*(mode->getY(1,X-2)-mode->getY(y3,X-2));
		b1=star->rho(XX-2)*star->dPhidr(XX-2)*star->rad(XX-2)*(mode->getY(1,X-1)-mode->getY(y3,X-1));
		a1=star->rho(XX+2)*star->dPhidr(XX+2)*star->rad(XX+2)*(mode->getY(1,X+1)-mode->getY(y3,X+1));
		a2=star->rho(XX+4)*star->dPhidr(XX+4)*star->rad(XX+4)*(mode->getY(1,X+2)-mode->getY(y3,X+2));
		a3=star->rho(XX+6)*star->dPhidr(XX+6)*star->rad(XX+6)*(mode->getY(1,X+3)-mode->getY(y3,X+3));
		dDP = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1)/Pscale*Rstar;
		
		if(isnan(d2Phi)) d2Phi = 0.0;
		if(isnan(difxi)) difxi = 0.0;
		if(isnan(dDP)) dDP = 0.0;
		
		//Now calculate residuals
		//Perturbed Poisson equation
		e1 = fabs(
				4.0*m_pi*Drho*Gee
				+ L2*pow(r,-2.0)*DPhi 
				- 2.0*dDPhi/r  - d2Phi
		);
		n1 = fabs(4.0*m_pi*Drho*Gee )
				+ fabs( L2*pow(r,-2.0)*DPhi )
				+ fabs( 2.0*dDPhi/r ) + fabs( d2Phi );
		//continuity equation
		e2 = fabs(
				Drho + xi*drhodr + rho*( 2.*xi/r + difxi - L2*pow(r,-2.0)/freq2*chi )
		);
		n2 = fabs(Drho) + fabs(xi*drhodr) + fabs(2.*rho/r*xi)
				+ fabs(rho*difxi)
				+ fabs(rho*L2*chi*pow(r,-2.0)/freq2);
		//newton's equation -- the r component. The theta component defines chi = r sigma2 xiH
		e3 = fabs( rho*freq2*xi - g*Drho - rho*dDPhi - dDP );
		n3 = fabs( rho*freq2*xi ) + fabs( g*Drho ) + fabs( rho*dDPhi ) + fabs( dDP );
		//normalize residuals
		e1 = e1/n1;
		e2 = e2/n2;
		e3 = e3/n3;
		if(isnan(e1)) e1 = 0.0;
		if(isnan(e2)) e2 = 0.0;
		if(isnan(e3)) e3 = 0.0;
		//collect residuals
		checkPois += e1*e1;
		checkCont += e2*e2;
		checkNewt += e3*e3;
	}
	return sqrt((checkCont+checkPois+checkNewt)/double(3*len-7));
}

double NonradialModeDriver::tidal_overlap(ModeBase* mode){
	double omega2 = mode->getOmega2();
	int k,l,m;
	mode->modeNumbers(k,l,m);
	double L2 = double(l*(l+1));
	double Mstar=star->Mass(), Rstar=star->Radius();
	double dscale = Mstar*pow(Rstar,-3); //dimensionless, set equal to M/R^3
	double rscale = 1.0;//Rstar;		//dimensionless, set equal to 1
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
		xir = mode->getY(y1,x)*r[x]*rscale;
		xiH = mode->getY(y2,x)*star->dPhidr(xx)/freq2*rscale/Rstar;
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

double NonradialModeDriver::innerproduct(ModeBase* mode1, ModeBase* mode2){
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
		xir1 = mode1->getY(y1,x)*r[x]*rscale;
		xir2 = mode2->getY(y1,x)*r[x]*rscale;
		xiH1 = mode1->getY(y2,x)*star->dPhidr(xx)/freq21*rscale/Rstar;
		xiH2 = mode2->getY(y2,x)*star->dPhidr(xx)/freq22*rscale/Rstar;
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