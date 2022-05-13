//**************************************************************************************
//						GR COWLING PULSATION nlm MODE DRIVER
// GRCowlingModeDriver.cpp
//		Based on the equation of Yoshida & Lee (2002), using Dziembowski-like variables,
//		See their equations (42), (43), and definitions in (37) to (41), (44)
//		The metric component has been corrected, their 2nu -> nu, their 2lambda -> lambda
//			total metric: ds^2 = -e^nu dt^2 + e^lambda dr^2 + r^2 dOmega^2
//		In their (43), the "U" should be U_1
//***************************************************************************************

#ifndef GRCOWLINGMODECLASS
#define GRCOWLINGMODECLASS

#include "GRCowlingModeDriver.h"

//constructor
GRCowlingModeDriver::GRCowlingModeDriver(GRStar *star, double Gamma1)
	 : star(star), ModeDriver(num_var, star), adiabatic_index(Gamma1)
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
	
	Gee = star->Gee();
	cee2 = star->light_speed2();
	zsurf = star->getRedshift();
	printf("REDSHIFT = %le\n", zsurf);
	
	//initialize helper arrays
	initializeArrays();
	//initialize the boundary
	central_bc_order = BC_C;
	surface_bc_order = BC_S;
	setupBoundaries();
}

//destructor
GRCowlingModeDriver::~GRCowlingModeDriver(){
	delete[] r;
	delete[] A; delete[] U1; delete[] U2;
	delete[] C; delete[] V;
}

int GRCowlingModeDriver::length(){
	return len;
}

double GRCowlingModeDriver::Gamma1(){
	return adiabatic_index;
}

double GRCowlingModeDriver::rad(int X){
	return r[X];
}

// functions of the stellar background, as defined in Yoshida & Lee 2002, eq 38-41
// note that their 2nu -> nu, their 2lambda -> lambda
// the definition of V1 in eq 38 corresponds to the usual V1 = -d(log P)/d(log r)
void GRCowlingModeDriver::initializeArrays(){
	A = new double[len_star]; U1 = new double[len_star]; U2 = new double[len_star];
	C = new double[len_star]; V  = new double[len_star];
	A[0]  = 0.0;
	V[0]  = 0.0;
	U1[0] = 3.0;
	U2[0] = 0.0;
	C[0]  = star->getC(0); //this value DOES matter, in 0-order central BC
	int X;
	for(X=1; X<len_star-1; X++){
		A[X]  = star->getAstar(X, adiabatic_index);
		V[X]  = star->getVg(X, adiabatic_index);
		C[X]  = star->getC(X);
		//
		U1[X] = 2. + star->rad(X)/star->dNudr(X)*star->d2Nudr2(X);
		U2[X] = star->rad(X)*star->dLambdadr(X)/2.;
	}
	//surface requires special treatment
	X = len_star-1;
	A[X]  = 0.0;  
	V[X]  = 0.0;
	U1[X] = 0.0;
	U2[X] = 0.0;
	C[X]  = star->getC(X);
}


//these equations come from Yoshida & Lee 2002, eq (42),(43)
// the definition of U1 has been changed to create greater agreement with classical Cowling equations
void GRCowlingModeDriver::getCoeff(
	double *CCI,
	int X,
	int b,
	double omeg2,
	int L
){	
	double (*CC)[num_var] =(double (*)[num_var]) CCI;
	//avoid recalculations
	double Ax, U1x, U2x, Cx, Vx;
	Ax = A[2*X+b];
	double expL = std::exp(star->Lambda(2*X+b));
	U1x = U1[2*X+b];
	U2x = U2[2*X+b];
	Cx = C[2*X+b];
	Vx = V[2*X+b];
	//dy1/dx
	CC[0][0] = double(2-L) + Vx - 3. - U2x;
	CC[0][1] = double(L*L+L)/(omeg2*Cx) - Vx;
	//dy2/dx
	CC[1][0] = omeg2*expL*Cx - Ax;
	CC[1][1] = double(2-L) + 1. + Ax - U1x;
}


void GRCowlingModeDriver::getBoundaryMatrix(int nv, double *y0, double* ys, double **y, int* indexOrder){
	//a basis of independent solutions to be used at center and surface
	const double yy[num_var][num_var] = {
		{ 1.0, 1.0},
		{ 1.0, 1.0}
	};
	for(int i=0; i<nv; i++)
		for(int j=0; j<nv; j++)
			y[i][j] = yy[i][j];
	//when matching outward/inward, indexOrder matches correct coefficient to correct solution
	const int indices[num_var] = {0,0};
	for(int i=0; i<nv; i++) indexOrder[i] = indices[i];
}

void GRCowlingModeDriver::setupBoundaries() {
	//set up coefficients for the central boundary
	star->getAstarCenter( Ac, central_bc_order, adiabatic_index);
	star->getVgCenter(   Vgc, central_bc_order, adiabatic_index);
	star->getC1Center(    cc, central_bc_order);
	star->getLambdaCenter(lc, central_bc_order);
	star->getNuCenter(   nuc, central_bc_order);
	switch(central_bc_order){
		case BC_C:
			U1c[2] = 4.*(3.*nuc[1]*nuc[3]-2.*nuc[2]*nuc[2])/(nuc[1]*nuc[1]);
			U2c[2] = 2.*lc[2];
		case 3:
		case 2:
			U1c[1] = 4.*nuc[4]/nuc[2];
			U2c[1] = lc[1];
		case 1:
		case 0:
			U1c[0] = 3.0;
			U2c[0] = 0.0;
			break;
	}
	
	//set up coefficients for the surface boundary
	star->getAstarSurface( As, surface_bc_order, adiabatic_index);
	star->getVgSurface(   Vgs, surface_bc_order, adiabatic_index);
	star->getC1Surface(    cs, surface_bc_order);
	star->getLambdaSurface(ls, surface_bc_order);
	star->getNuSurface(   nus, surface_bc_order);
	switch(central_bc_order){
		case BC_S:
			U1s[4] = 0.0;
			U2s[4] = 2.*ls[4] - 5.*ls[5]/2.;
		case 3:
			U1s[3] = 0.0;
			U2s[3] = 3.*ls[3]/2. - 2.*ls[4];
		case 2:
			U1s[2] = 0.0;
			U2s[2] = ls[2] - 3.*ls[3]/2.;
		case 1:
			U1s[1] = 4.*nus[2]*nus[2] + 2.*nus[1]*(nus[2]-3.*nus[3])/nus[1]/nus[1];
			U2s[1] =  0.5*ls[1]-ls[2];
		case 0:
			U1s[0] = 2. - 2.*nus[2]/nus[1];
			U2s[0] = -0.5*ls[1];
			break;
	}
	k_surface = 0.0;	
}

int GRCowlingModeDriver::CentralBC(double **ymode, double *y0, double omeg2, int l, int m){	
	double yy[num_var][central_bc_order/2+1];//0,2,4
	double L = double(l);
	double L2 = double(l*l+l);
	double Gam1 = (adiabatic_index == 0.0 ? star->Gamma1(0) : adiabatic_index);
	
		
	//the zero-order terms are simple
	y0[y2] = y0[y1]*cc[0]*omeg2/L;
	for(int i=0; i<num_var; i++) yy[i][0] = y0[i];
	
	//the matrix A0, which only appears as a prefactor
	double A0[num_var][num_var] = {
		{ -1. - L - U2c[0], L2/(omeg2*cc[0])},
		{      omeg2*cc[0], 3. - L - U1c[0]}
	};
	//the matrix A_2, the second-order component in A = A_0 + A_2*x^2
	double A2[num_var][num_var] = {
		{ Vgc[1] - U2c[1]              , -L2*cc[1]/omeg2/(cc[0]*cc[0])-Vgc[1]},
		{ -Ac[1] + omeg2*(cc[1]+cc[0]*lc[1]), Ac[1] - U1c[1]}
	};
	double A4[num_var][num_var] = {
		{ Vgc[2] - U2c[2] , L2*(cc[1]*cc[1] - cc[0]*cc[2])/omeg2/(cc[0]*cc[0]*cc[0])- Vgc[2] },
		{ -Ac[2] + omeg2*(cc[2]+cc[1]*lc[1] + cc[0]*(0.5*lc[1]*lc[1]+lc[2])), Ac[2] - U1c[2] }
	};
	
	double RHS[num_var], ynext[num_var]={0.};
	double coeff[3][num_var][num_var];
	double kA0[num_var][num_var];
	for(int i=0; i<num_var; i++){
		for(int j=0; j<num_var; j++){
			coeff[0][i][j] = A0[i][j];
			coeff[1][i][j] = A2[i][j];
			coeff[2][i][j] = A4[i][j];
		}
	}
	
	for(int k=1; k<central_bc_order/2+1; k++){
		//make the RHS array, k Y[k] - A[k-m]Y[m]
		for(int i=0; i<num_var; i++){
			RHS[i] = 0.0;
			for(int j=0; j<num_var; j++){
				//add up matrix products over lower-orders
				for(int m=0; m<k; m++){
					RHS[i] += coeff[k-m][i][j]*yy[j][m];
				}
				kA0[i][j] = -coeff[0][i][j];
			}
			kA0[i][i] += double(2*k);
		}
		invertMatrix(kA0, RHS, ynext);
		for(int i=0; i<num_var; i++) yy[i][k] = ynext[i];
	}

	int start = 2;
	for(int X=0; X<=start; X++){
		for(int i=0; i<num_var; i++){
			ymode[i][X] = yy[i][0];
			for(int j=1; j<=central_bc_order/2; j++) ymode[i][X] += yy[i][j]*pow(r[X],2*j);
		}
	}
	return start;		
}

int GRCowlingModeDriver::SurfaceBC(double **ymode, double *ys, double omeg2, int l, int m){
	//specify initial conditions at surface
	double yy[num_var][surface_bc_order+1];	//coefficients y = yy[0] + yy[1]t + ... + yy[k]t^k
	double yyn[num_var] = {0.0};
	//constants that show up
	double L  = double(l);
	double L2 = double(l*l+l);
	double Gam1 = (adiabatic_index == 0.0 ? star->Gamma1(0) : adiabatic_index);
	
	//the zero-order terms are simple
	ys[y2] = ys[y1];
	for(int i=0; i<num_var; i++) yy[i][0] = ys[i];
	
	//set up the matrix of surface coefficients -- terms with l or omeg2 are not included
	double A11, A12, A21, A22;
	int O = 1; //anchor index
	//the matrix A(-1)
	A11 =  Vgs[O-1];
	A12 = -Vgs[O-1];
	A21 = - As[O-1];
	A22 =   As[O-1];
	double Am1[num_var][num_var] = {
		{ A11,  A12 },
		{ A21,  A22 }
	};
	//the matrix A(0)
	//begin with non-product terms consistent at this order
	A11 =  Vgs[O] - (L+1.) - U2s[0];
	A12 =  L2/omeg2/cs[0]  - Vgs[O];
	A21 = -As[O] + std::exp(ls[0])*omeg2*cs[0];
	A22 =  (3.-L) + As[O]  - U1s[0];
	double A0[num_var][num_var] = {
		{  A11,  A12 },
		{  A21,  A22 }
	};
	//the matrix A(1)
	A11 =  Vgs[O+1] - U2s[1];
	A12 = -Vgs[O+1] - L2*cs[1]/omeg2/cs[0]/cs[0];
	A21 = - As[O+1] + std::exp(ls[0])*omeg2*(cs[1] + cs[0]*ls[1]);
	A22 =   As[O+1] - U1s[1];
	double A1[num_var][num_var] = {
		{  A11,  A12 },
		{  A21,  A22 }
	};
	
	
	//the matrix A(n)
	// used to capture first non-integer power in surface expansion
	A11  = Vgn - U2sn;
	A12  =-Vgn;
	A21  =-Asn;
	A22  = Asn - U1sn;
	double AN[num_var][num_var] = {
		{  A11,  A12 },
		{  A21,  A22 }
	};	
	
	
	double coeff[O+surface_bc_order][num_var][num_var];
	double RHS[num_var], ynext[num_var]={0.};
	
	for(int i=0; i<num_var; i++){
		for(int j=0; j<num_var; j++){
			coeff[O-1][i][j] = Am1[i][j];
			coeff[O+0][i][j] = A0[i][j];
			coeff[O+1][i][j] = A1[i][j];
		}
	}
	
	for(int k=0; k<surface_bc_order; k++){
		//make the RHS array, k Y[k] - A[k-m]Y[m]
		for(int i=0; i<num_var; i++){
			RHS[i] = double(k)*yy[i][k];
			for(int j=0; j<num_var; j++){
				//add up matrix products over lower-orders
				for(int m=0; m<=k; m++){
					RHS[i] -= coeff[O+k-m][i][j]*yy[j][m];
				}
				Am1[i][j] = coeff[O-1][i][j];
			}
			Am1[i][i] += double(k+1);
		}
		invertMatrix(Am1, RHS, ynext);
		for(int i=0; i<num_var; i++) yy[i][k+1] = ynext[i];
	}
	
	//the number of terms to calculate
	int start = len-2;
	double t;
	for(int X=len-1; X>=start; X--){
		t = (1.-r[X]);
		for(int i=0; i<num_var; i++){
			ymode[i][X] = yy[i][0];
			for(int k=1; k<=surface_bc_order; k++) ymode[i][X] += yy[i][k]*pow(t,k);
		}
	}
	
	return start;
}//*/

// *** this method is assuming a uniform grid ***
// *** we cannot assume a uniform grid for all stars ***
double GRCowlingModeDriver::SSR(double omega2, int l, ModeBase* mode){
	double checkCont=0.0;
	double checkNewt=0.0;
	double e1,e2;
	double n1,n2;
	double difxiR, dDP;
	double r,G1,xiR,xiH, rho,P,energ, nu,la, drhodr,dldr,dndr, Drho, DP;
	double b3,b2,b1,a1,a2,a3, h1;
	double freq2 = omega2*star->mr(len_star-1)*pow(star->rad(len_star-1),-3);
	double L2    = double(l*l+l);
	for(int X=4; X<len-4; X++){
		int XX = 2*X;
		//stellar variables, to simplify equations
		rho    = star->rho(XX);
		P      = star->P(XX);
		nu     = star->Nu(XX);
		la     = star->Lambda(XX);
		energ  = (rho + P/cee2);
		drhodr = star->drhodr(XX);
		dldr   = star->dLambdadr(XX);
		dndr   = star->dNudr(XX);
		r      = star->rad(XX);
		G1     = (Gamma1()==0.0? star->Gamma1(XX) : Gamma1());
		//mode variables			
		xiR  = r*mode->getY(0,X);	            //r*y1
		xiH  =   mode->getY(1,X)/(C[XX]*omega2);//  y2/(c1*omeg2)
		Drho = 0.0;//rho*rho/(G1*P)*(chi - g*xi) - xi*drhodr;
		DP   = 0.0;
		
		
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
		difxiR = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
		//now dDP/dr, DP = rho*(chi-DPhi), chi=rg*y2, DPhi=rg*y3
		b3=star->rho(XX-6)*star->dPhidr(XX-6)*star->rad(XX-6)*(mode->getY(1,X-3));
		b2=star->rho(XX-4)*star->dPhidr(XX-4)*star->rad(XX-4)*(mode->getY(1,X-2));
		b1=star->rho(XX-2)*star->dPhidr(XX-2)*star->rad(XX-2)*(mode->getY(1,X-1));
		a1=star->rho(XX+2)*star->dPhidr(XX+2)*star->rad(XX+2)*(mode->getY(1,X+1));
		a2=star->rho(XX+4)*star->dPhidr(XX+4)*star->rad(XX+4)*(mode->getY(1,X+2));
		a3=star->rho(XX+6)*star->dPhidr(XX+6)*star->rad(XX+6)*(mode->getY(1,X+3));
		dDP = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
		//Now calculate residuals
		// conservation time component,   D_aT^{a0} = 0
		e1 = fabs(
				Drho - L2*energ*xiH - energ*difxiR + xiR*(2./r*energ + drhodr + 0.5*energ*dldr)
		);
		n1 = fabs(Drho);
		// conservation radial component, D_aT^{ar} = 0
		e2 = fabs( dDP +  0.5*dndr*(Drho*cee2+DP) - freq2*xiR*energ*std::exp(nu-la) );
		n2 = fabs( dDP );
		//normalize residuals
		e1 = e1/n1;
		e2 = e2/n2;
		//collect residuals
		checkCont += e1*e1;
		checkNewt += e2*e2;
	}
	return sqrt((checkCont+checkNewt)/double(2*len-7));
}

double GRCowlingModeDriver::tidal_overlap(ModeBase* mode){
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

double GRCowlingModeDriver::innerproduct(ModeBase* mode1, ModeBase* mode2){
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