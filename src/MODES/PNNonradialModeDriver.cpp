//***************************************************************************************
//				POST-NEWTONIAN NONRADIAL PULSATION nlm MODE DRIVER
// PNNonradialMode.cpp
// 		Solves the 1PN LAWE in the 1PN Djiembowski variables
//		Modes are 8th order, reduced by harmonic coordinate condition
//***************************************************************************************

#ifndef PNFULLMODEDRIVERCLASS
#define PNFULLMODEDRIVERCLASS

#include "PNNonradialModeDriver.h"

//constructor
PNNonradialModeDriver::PNNonradialModeDriver(PNStar *star, double Gamma1)
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
PNNonradialModeDriver::~PNNonradialModeDriver(){
	delete[] r;
	delete[] A; delete[] U;
	delete[] C; delete[] V;
	delete[] beta;
	delete[] Phi;
	delete[] xstar;
}

int PNNonradialModeDriver::length(){
	return len;
}

double PNNonradialModeDriver::Gamma1(){
	return adiabatic_index;
}

double PNNonradialModeDriver::rad(int X){
	return r[X];
}

void PNNonradialModeDriver::initializeArrays(){
	double Phi1 = star->Phi(len_star-1);
	double Gam1 = adiabatic_index;
	double R = star->Radius();
	A = new double[len_star];  U = new double[len_star];
	C = new double[len_star];  V = new double[len_star];
	beta  = new double[len_star];
	Phi   = new double[len_star];
	xstar = new double[len_star];
	A[0]  = 0.0;
	V[0]  = 0.0;
	U[0]  = 3.0;
	C[0]  = star->getC(0);
	beta[0] = star->sound_speed2(0,Gam1)/Phi1;
	Phi[0] = star->Phi(0)/Phi1;
	xstar[0] = 0.0;
	int X;
	for(X=1; X<len_star-1; X++){
		A[X] = star->getAstar(X,Gam1);
		U[X] = star->getU(X);
		C[X] = star->getC(X);
		V[X] = star->getVg(X,Gam1);
		beta[X]  = star->sound_speed2(X,Gam1)/Phi1;
		Phi[X]   = star->Phi(X)/Phi1;
		xstar[X] = star->rad(X)/R;
	}
	//surface requires special treatment
	X = len_star-1;
	V[X]  = 0.0;//indeterminate
	U[X]  = 0.0;//either 3.0 or 0.0
	A[X]  = 0.0;//indeterminate
	C[X]  = star->getC(X);
	beta[X] = 0.0; //almost always 0
	Phi[X]  = 1.0; //1 by definition
	xstar[X] = 1.0;//1 by definition
}

void PNNonradialModeDriver::getCoeff(
	double *CCI,
	int X,
	int b,
	double omeg2,
	int L
)
{	
	double (*CC)[num_var] =(double (*)[num_var]) CCI;
	int here = 2*X+b;
	//Gamma1
	double iG1 = (adiabatic_index==0.0 ? 1.0/star->Gamma1(here) : 1.0/adiabatic_index);
	double Ax, Ux, Cx, Vx, Phix, betax;
	double x2 = xstar[here]*xstar[here];
	Ax    =    A[here];
	Ux    =    U[here];
	Cx    =    C[here];
	Vx    =    V[here];
	Phix  =  Phi[here];
	betax = beta[here];
	double L2  = double(L*L+L);
	double Lm2 = double(2-L);

	//set everything to zero, then reset the others
	for(int i=0;i<num_var;i++)
		for(int j=0;j<num_var;j++)
			CC[i][j]=0.0;
	//dy1/dx
	CC[y1][y1] = Lm2 + Vx - 3. - zsurf*(3.*Vx*betax);
	CC[y1][y2] = L2/(omeg2*Cx) - Vx - 4.*zsurf*Vx*Phix;
	CC[y1][y3] = Vx - 3.*zsurf*Vx*betax;
	CC[y1][z1] = -zsurf*Vx; //z1
	CC[y1][z5] = -zsurf*Vx; //z5
	//dy2/dx
	CC[y2][y1] = omeg2*Cx - Ax + zsurf*(4.*Ax*Phix - 4.*Ux*omeg2*x2/L2);
	CC[y2][y2] = Lm2 + 1. + Ax - Ux - 4.*zsurf*Vx*betax;
	CC[y2][y3] =-Ax + 4.*zsurf*Ax*Phix;
	CC[y2][y4] =-zsurf*4.*omeg2*x2/L2;
	CC[y2][z1] = zsurf*Ax; //z1
	CC[y2][z5] = zsurf*Ax; //z5
	//dy3/dx
	CC[y3][y3] = double(2-L) + 1.-Ux;
	CC[y3][y4] = 1.0;
	//dy4/dx
	CC[y4][y1] = Ux*Ax * (1. + 2.*zsurf*(betax*iG1-Phix) );
	CC[y4][y2] = Ux*Vx * (1. + 2.*zsurf*(betax*iG1+Phix));
	CC[y4][y3] = L2 - Ux*Vx*(1. + 2.*zsurf*(betax*iG1-Phix));
	CC[y4][y4] = Lm2 - Ux;
	CC[y4][z1] = zsurf*Ux*Vx;	//z1
	CC[y4][z5] = zsurf*Ux*Vx;	//z5
	
	//dz1/dx
	CC[z1][z1] = Lm2 + 1.-Ux;
	CC[z1][z2] = 1.0;
	//dz2/dx
	CC[z2][y1] = - 2.0*Ux*Ax*Phix;
	CC[z2][y2] = Ux*Vx*(3.*betax- 2.0*Phix);
	CC[z2][y3] = omeg2*x2  + Ux*Vx*(2.*Phix - 5.*betax);
	CC[z2][z1] = L2;		//z1
	CC[z2][z2] = Lm2 - Ux;	//z2
	
	//dz3/dx
	CC[z3][y3] =-4.*omeg2*x2;
	CC[z3][z3] = Lm2 - Ux;	//z5
	CC[z3][z5] = L2;        //z6
	//dz5/dx
	CC[z5][y1] = 4.*omeg2*Ux*x2/L2;//y1
	CC[z5][y4] = 4.*omeg2*   x2/L2;//y2
	CC[z5][z3] = 1.0;		            //z5
	CC[z5][z5] = Lm2 + 1.-Ux;	//z6			
}

void PNNonradialModeDriver::getBoundaryMatrix(int nv, double *y0, double* ys, double **y, int* indexOrder){
	double one = 1.0; //this is only used to make visual distinction in grid
	const double yy[num_var][num_var] = {	
		{ one, one, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},//a
		{ 0.0, 0.0, one, one, 0.0, 0.0, 0.0, 0.0},//b
		{ 0.0, 0.0, 0.0, 0.0, one, one, 0.0, 0.0},//c
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, one, one},//d
		//
		{ 0.0, one, one, one, 0.0, 0.0, 0.0, 0.0},//e
		{ 0.0, 0.0, 0.0, 0.0, one, one, 0.0, 0.0},//f
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, one, one},//g
		{ one, one, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} //1
	};
	const int indices[num_var] = {y1,y3,z1,z5,y3,z1,z5,y1};
	for(int i=0; i<nv; i++)
		for(int j=0; j<nv; j++)
			y[i][j] = yy[i][j];
	for(int i=0; i<nv; i++) indexOrder[i] = indices[i];
}

void PNNonradialModeDriver::setupBoundaries() {
	//set up coefficients for the central boundary
	double read[2];
	star->getAstarCenter(read, central_bc_order, adiabatic_index);
	Ac2 = read[1];	//Ac0 = 0
	star->getVgCenter(read, central_bc_order, adiabatic_index);
	Vgc2 = read[1];	//Vgc0 = 0
	star->getUCenter(read, central_bc_order);
	Uc2 = read[1];	//Uc0 = 3
	star->getC1Center(cc, central_bc_order);	//need both
	//
	star->getBetaCenter(read, central_bc_order, adiabatic_index);
	betac0 = read[0];//betac2 != 0. but does not appear in equations
	star->getPhiCenter(read, central_bc_order);
	phic0 = read[0]; //phic2 != 0. but does not appear in equations
		
	//set up coefficients for the surface boundary
	const int indxAll=BC_S+1;
	double As[indxAll], Vgs[indxAll], Us[indxAll],
		 betas[indxAll], phis[indxAll];
	double Vgn=0, Usn=0, Asn=0;
	star->getAstarSurface(As, surface_bc_order, adiabatic_index);
	star->getVgSurface(   Vgs,surface_bc_order, adiabatic_index);
	star->getUSurface(    Us, surface_bc_order);
	star->getC1Surface(   cs, surface_bc_order);
	//
	star->getBetaSurface(betas,surface_bc_order, adiabatic_index);
	star->getPhiSurface( phis, surface_bc_order);
	betas0 = betas[0];
	phis0  = phis[0];
	//fix surface values
	C[len_star-1] = cs[0];
	beta[len_star-1] = betas0;
	U[len_star-1] = Us[0];
	
	//set up the matrix of surface coefficients -- terms with l or omeg2 are not included
	double Gam1 = (adiabatic_index == 0.0 ? star->Gamma1(len_star-1) : adiabatic_index);
	double z = zsurf;
	double A11, A12, A13, A21, A22, A23, A41, A42, A43, A4z1, Az21, Az22, Az23;
	int O = 1; //anchor index
	//the matrix A(-1)
	A11 =  Vgs[O-1] - 3.*z*Vgs[O-1]*betas[0];
	A12 = -Vgs[O-1] - 4.*z*Vgs[O-1]* phis[0];
	A13 =  Vgs[O-1] - 3.*z*Vgs[O-1]*betas[0];
	A21 = - As[O-1] + 4.*z* As[O-1]* phis[0];
	A22 =   As[O-1] - 4.*z*Vgs[O-1]*betas[0];
	A23 = - As[O-1] + 4.*z* As[O-1]* phis[0];
	A41 =   As[O-1]*Us[0]*(1.+2.*z*betas[0]/Gam1-2.*z*phis[0]);
	A42 =  Vgs[O-1]*Us[0]*(1.+2.*z*betas[0]/Gam1+2.*z*phis[0]);
	A43 = -Vgs[O-1]*Us[0]*(1.+2.*z*betas[0]/Gam1-2.*z*phis[0]);
	A4z1 = Vgs[O-1]*Us[0]*z;
	Az21 = -2.* As[O-1]*Us[0]*phis[0];
	Az22 = -2.*Vgs[O-1]*Us[0]*phis[0] + 3.*Vgs[O-1]*Us[0]*betas[0];
	Az23 =  2.*Vgs[O-1]*Us[0]*phis[0] - 5.*Vgs[O-1]*Us[0]*betas[0];
	double Am1[num_var][num_var] = {
		{ A11,  A12,  A13, 0,-z*Vgs[O-1],0,0,-z*Vgs[O-1]},
		{ A21,  A22,  A23, 0, z* As[O-1],0,0, z* As[O-1]},
		{   0,    0,    0, 0,          0,0,0,          0},
		{ A41,  A42,  A43, 0,       A4z1,0,0,       A4z1},
		{   0,    0,    0, 0,          0,0,0,          0},
		{Az21, Az22, Az23, 0,          0,0,0,          0},
		{   0,    0,    0, 0,          0,0,0,          0},
		{   0,    0,    0, 0,          0,0,0,          0}
	};
	//the matrix A(0)
	//begin with non-product terms consistent at this order
	A11 =  Vgs[O] - 3.;
	A12 = -Vgs[O];
	A13 =  Vgs[O];
	A21 = -As[O];
	A22 =  1. + As[O] - Us[0];
	A23 = -As[O];
	A41 = A42 = A43 = A4z1 = 0.0;
	Az21 = Az22 = Az23 = 0.0;
	//then handle the product terms separately as sums of indices, k=0
	for(int j=0; j<=1; j++){
		A11 -= 3.*z*(Vgs[O-j]*betas[j]);
		A12 -= 4.*z*(Vgs[O-j]* phis[j]);
		A13 -= 3.*z*(Vgs[O-j]*betas[j]);
		A21 += 4.*z*( As[O-j]* phis[j]);
		A22 -= 4.*z*(Vgs[O-j]*betas[j]);
		A23 += 4.*z*( As[O-j]* phis[j]);
		A41 +=  As[O-j]*Us[j];
		A42 += Vgs[O-j]*Us[j];
		A43 -= Vgs[O-j]*Us[j];
		A4z1+= Vgs[O-j]*Us[j]*z;
		for(int m=0; m<=1-j; m++){
			A41 += 2.* As[O-j-m]*Us[j]*z*(betas[m]/Gam1-phis[m]);
			A42 += 2.*Vgs[O-j-m]*Us[j]*z*(betas[m]/Gam1+phis[m]);
			A43 -= 2.*Vgs[O-j-m]*Us[j]*z*(betas[m]/Gam1-phis[m]);
			Az21 -= 2.* As[O-j-m]*Us[j]*phis[m];
			Az22 += 3.*Vgs[O-j-m]*Us[j]*betas[m]- 2.*Vgs[O-j-m]*Us[j]*phis[m];
			Az23 += 2.*Vgs[O-j-m]*Us[j]*phis[m] - 5.*Vgs[O-j-m]*Us[j]*betas[m];
		}
	}
	double A0[num_var][num_var] = {
		{  A11,  A12,    A13,      0,-z*Vgs[O],     0,     0,-z*Vgs[O]},
		{  A21,  A22,    A23,  -4.*z,  z*As[O],     0,     0,  z*As[O]},
		{   0.,    0,1.-Us[0],    1.,        0,     0,     0,        0},
		{  A41,  A42,    A43, -Us[0],     A4z1,     0,     0,     A4z1},
		{   0.,    0,      0,      0, 1.-Us[0],    1.,     0,        0},
		{ Az21, Az22,   Az23,      0,        0,-Us[0],     0,        0},
		{   0.,    0,    -4.,      0,        0,     0,-Us[0],        0},
		{4.*Us[0], 0,      0,     4.,        0,     0,    1.,  1.-Us[0]}
	};
	//the matrix A(1)
	A11 =  Vgs[O+1];
	A12 = -Vgs[O+1];
	A13 =  Vgs[O+1];
	A21 = -As[O+1];
	A22 =  As[O+1] - Us[1];
	A23 = -As[O+1];
	A41 = A42 = A43 = A4z1 = 0.0;
	Az21 = Az22 = Az23 = 0.0;
	for(int j=0; j<=2; j++){
		A11 -= 3.*z*(Vgs[O+1-j]*betas[j]);
		A12 -= 4.*z*(Vgs[O+1-j]*phis[j]);
		A13 -= 3.*z*(Vgs[O+1-j]*betas[j]);
		A21 += 4.*z*( As[O+1-j]*phis[j]);
		A22 -= 4.*z*(Vgs[O+1-j]*phis[j]);
		A23 += 4.*z*( As[O+1-j]*phis[j]);
		A41 +=   As[O+1-j]*Us[j];
		A42 +=  Vgs[O+1-j]*Us[j];
		A43 -=  Vgs[O+1-j]*Us[j];
		A4z1 += Vgs[O+1-j]*Us[j]*z;
		for(int m=0; m<=2-j; m++){
			A41 += 2.* As[O+1-j-m]*Us[j]*z*(betas[m]/Gam1-phis[m]);
			A42 += 2.*Vgs[O+1-j-m]*Us[j]*z*(betas[m]/Gam1+phis[m]);
			A43 -= 2.*Vgs[O+1-j-m]*Us[j]*z*(betas[m]/Gam1-phis[m]);
			Az21 -= 2.* As[O+1-j-m]*Us[j]*phis[m];
			Az22 += 3.*Vgs[O+1-j-m]*Us[j]*betas[m] - 2.*Vgs[O+1-j-m]*Us[j]*phis[m];
			Az23 += 2.*Vgs[O+1-j-m]*Us[j]*phis[m]  - 5.*Vgs[O+1-j-m]*Us[j]*betas[m];
		}
	}
	double A1[num_var][num_var] = {
		{  A11,  A12,   A22,      0,-z*Vgs[O+1],     0,     0,-z*Vgs[O+1]},
		{  A21,  A22,   A23,   8.*z, z* As[O+1],     0,     0, z* As[O+1]},
		{    0,    0,-Us[1],      0,          0,     0,     0,          0},
		{  A41,  A42,   A43, -Us[1],       A4z1,     0,     0,       A4z1},
		{    0,    0,     0,      0,     -Us[1],     0,     0,          0},
		{ Az21, Az22,  Az23,      0,          0,-Us[1],     0,          0},
		{    0,    0,    8.,      0,          0,     0,-Us[1],          0},
		{4.*(Us[1]-2.*Us[0]),0,0,-8.,         0,     0,     0,     -Us[1]}
	};
	
	A11  = Vgn;
	A12  =-Vgn - 4.*z*Vgn*phis[0];
	A13  = Vgn;
	A21  =-Asn*(1.-4.*z*phis[0]);
	A22  = Asn - Usn;
	A23  =-Asn + 4.*Asn*z*phis[0];
	A41  = Usn*As[O]*(1.+betas[0]/Gam1) + Usn*As[O-1]*betas[1]/Gam1;
	A42  = Usn*( Vgs[O]*(1.+betas[0]/Gam1+4.*z*phis[0]) + Vgs[O-1]*(betas[1]/Gam1+4.*z*phis[1]) );
	A43  =-Usn*(Vgs[O]+Vgs[O]*betas[0]/Gam1+Vgs[O-1]*betas[1]/Gam1);
	A4z1 = z*Usn*Vgs[O];
	Az21 =-2.*Usn*(As[O]*phis[0]+As[O-1]*phis[1]);
	Az22 = Usn*(-2.*Vgs[O]*phis[0]-2.*Vgs[O-1]*phis[1]);
	Az23 = Usn*(+2.*Vgs[O]*phis[0]+2.*Vgs[O-1]*phis[1]);
	double AN[num_var][num_var] = {
		{  A11,  A12,  A13,    0, -z*Vgn,    0,    0,-z*Vgn},
		{  A21,  A22,  A23,    0,  z*Asn,    0,    0, z*Asn},
		{    0,    0, -Usn,    0,      0,    0,    0,     0},
		{  A41,  A42,  A43, -Usn,   A4z1,    0,    0,  A4z1},
		{    0,    0,    0,    0,   -Usn,    0,    0,     0},
		{ Az21, Az22, Az23,    0,      0, -Usn,    0,     0},
		{    0,    0,    0,    0,      0,    0, -Usn,     0},
		{4.*Usn,   0,    0,    0,      0,    0,    0,  -Usn}
	};	
	for(int i=0; i<num_var; i++){
		for(int j=0; j<num_var; j++){
			surfCoeff[O-1][i][j] = Am1[i][j];
			surfCoeff[O  ][i][j] =  A0[i][j];
			surfCoeff[O+1][i][j] =  A1[i][j];
			surfCoeffN[i][j] = AN[i][j];
		}
	}
}

void PNNonradialModeDriver::updateSurface(double (&A)[3][8][8], double omeg2, int l){
//the point of this method is to include the factors of omeg2 and L that are specific to mode
//the coefficients A12, A21, A24, A43, Az23, Az2z1, Az51, Az54 should contain L or omeg2
	double L2 = double(l*l+l);
	double L  = double(l);
	int O = 1;//anchor index -- A[O-1] = A_-1, A[O] = A_0, A[O+1] = A_1
	for(int k=-1; k<=1; k++){
		for(int i=0; i<num_var; i++){
			for(int j=0; j<num_var; j++){
				A[O+k][i][j] = surfCoeff[O+k][i][j];
			}
		}
	}
	//  A12
	A[O  ][y1][y2] += L2/(omeg2*cs[0]);
	A[O+1][y1][y2] -= L2*cs[1]/(omeg2*cs[0]*cs[0]);
	//  A21
	A[O  ][y2][y1] += omeg2*cs[0] - zsurf*A[O  ][7][0]*omeg2/L2;
	A[O+1][y2][y1] += omeg2*cs[1] - zsurf*A[O+1][7][0]*omeg2/L2;
	//  A24
	A[O  ][y2][y4] *= omeg2/L2;
	A[O+1][y2][y4] *= omeg2/L2;
	//  A43
	A[O  ][y4][y3]  += L2;
	//  AZ23
	A[O  ][z2][y3] += omeg2;
	A[O+1][z2][y3] -= 2.*omeg2;
	//  AZ2Z1
	A[O  ][z2][z1] += L2;
	//  AZ33
	A[O  ][z3][y3] *= omeg2;
	A[O+1][z3][y3] *= omeg2;
	//  Az3z5
	A[O  ][z3][z5] += L2;
	//  AZ51
	A[O  ][z5][y1] *= omeg2/L2;
	A[O+1][z5][y1] *= omeg2/L2;
	//  AZ54
	A[O  ][z5][y4] *= omeg2/L2;
	A[O+1][z5][y4] *= omeg2/L2;
	//include contribution from regularization term
	for(int i=0; i<num_var; i++) A[O][i][i] += 2.-L;
}

int PNNonradialModeDriver::CentralBC(double **ymode, double *y0, double omeg2, int l, int m){
	double L = double(l);
	double L2 = double(l*l+l);
	double Gam1 = (adiabatic_index == 0.0 ? star->Gamma1(0) : adiabatic_index);
	double z = zsurf;
	
	double yy[num_var][BC_C/2+1];
	//the zero-order terms are simple
	y0[1] = y0[0]*cc[0]*omeg2/double(l);//initial chi
	y0[3] = y0[2]*double(l);			//initial d(DPhi)/dr
	y0[5] = y0[4]*double(l);			//initial d(DPsi)/dr
	y0[6] = y0[7]*double(l);
	for(int i=0; i<num_var; i++) yy[i][0] = y0[i];
	
	//the matrix (2k - A_0), which only appears as a prefactor
	double kA0[num_var][num_var] = {
		{ 3.+L, -L2/(omeg2*cc[0]), 0,0,0,0,0,0},
		{-omeg2*cc[0], 2.+L,       0,0,0,0,0,0},
		{0,0, 2.+L, -1., 0,0,0,0},
		{0,0, -L2, 3.+L, 0,0,0,0},
		{0,0,0,0, 2.+L, -1., 0,0},
		{0,0,0,0, -L2, 3.+L, 0,0},
		{0,0,0,0,0,0,  3.+L, -L2},
		{0,0,0,0,0,0,   -1.,2.+L}
	};
	double A2[num_var][num_var] = {
		{Vgc2-3.*z*Vgc2*betac0, -L2*cc[1]/omeg2/(cc[0]*cc[0])-Vgc2-4.*z*Vgc2*phic0,
			Vgc2-3.*z*Vgc2*betac0, 0,-z*Vgc2,0,0,-z*Vgc2},
		{-Ac2+omeg2*cc[1] - 12.*z*omeg2/L2 + 4.*z*Ac2*phic0,
			Ac2 - Uc2 - 4.*z*Vgc2*betac0, -Ac2+4.*z*Ac2*phic0, -4.*z*omeg2/L2, z*Ac2,0,0,z*Ac2},
		{0,0,-Uc2,0,0,0,0,0},
		{3.*Ac2*(1.+2.*z*(betac0/Gam1-phic0)), 3.*Vgc2*(1.+2.*z*(betac0/Gam1+phic0)),
			-3.*Vgc2*(1.+2.*z*(betac0/Gam1-phic0)), -Uc2, 3.*z*Vgc2, 0,0,  3.*z*Vgc2},
		{0,0,0,0,-Uc2,0,0,0},
		{-6.*Ac2*phic0, 9.*Vgc2*betac0-6.*Vgc2*phic0, omeg2-15.*Vgc2*betac0+6.*Vgc2*phic0, 0,0,-Uc2,0,0},
		{0,0,-4.*omeg2,0,0,0,-Uc2,0},
		{12.*omeg2/L2,0,0,4.*omeg2/L2,0,0,0,-Uc2}
	};
	
	double RHS[num_var];
	for(int i=0; i<num_var; i++){
		RHS[i]=0.0;
		for(int j=0; j<num_var; j++){
			RHS[i] += A2[i][j]*yy[j][0];
		}
	}
	
	double ynext[num_var];
	invertMatrix(kA0, RHS, ynext);
	for(int i=0; i<num_var; i++) yy[i][1] = ynext[i];

	int start = 1;
	for(int X=0; X<=start; X++){
		for(int i=0; i<num_var; i++){
			ymode[i][X] = yy[i][0];
			for(int j=1; j<=central_bc_order/2; j++) ymode[i][X] += yy[i][j]*pow(r[X],2*j);
		}
	}
	return start;
}

int PNNonradialModeDriver::SurfaceBC(double **ymode, double *ys, double omeg2, int l, int m){
	//specify initial conditions at surface
	double yy[ num_var][BC_S+1];
	double yyn[num_var] = {0.0};
	
	int O=1; //anchor index
	double RHS[num_var], coeff[BC_S+1][num_var][num_var], Am1[num_var][num_var], ynext[num_var];
	//update the coefficient matrices at the surface
	updateSurface(coeff, omeg2, l);
		
	//first-order boundary conditions
	ys[y2] = (ys[y1] + ys[y3])*(1.-3.*zsurf*betas0)/(1.+4.*zsurf*phis0) 
				- zsurf*(ys[z1]+ys[z5])/(1.+4.*zsurf*phis0);//solve (A_-1)Y0 = 0 for y1
	ys[y4] = -double(l+1)*ys[y3];// y4 = -(l+1)*y3
	ys[z2] = -double(l+1)*ys[z1];// z2 = -(l+1)*z1
	ys[z3] = -double(l+1)*ys[z5] + 4.*omeg2*ys[2]/double(l);//solve for z3 in harmonic condition
	for(int i=0; i<num_var; i++) yy[i][0] = ys[i];
	
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
		if(invertMatrix(Am1, RHS, ynext)){
			printf("surface %d\n", k+1);
		}
		for(int i=0; i<num_var; i++) yy[i][k+1] = ynext[i];
	}
	
	//this handles n<2 for non-integer index polytropes
	n_surface = 3;
	if(double(int(n_surface))!=n_surface && n_surface<2.0){
		double AN[num_var][num_var];
		for(int i=0;i<num_var;i++){
			for(int j=0;j<num_var;j++){
				AN[i][j] = surfCoeffN[i][j];
				Am1[i][j] = coeff[O-1][i][j];
			}
			Am1[i][i] += n_surface+1.;
		}
		AN[z5][y1] *= omeg2/double(l*l+l);
		AN[y2][y1] -= AN[z5][y1];
		for(int i=0;i<num_var;i++){
			RHS[i] = 0.0;
			for(int j=0;j<num_var;j++){
				RHS[i] -= AN[i][j]*yy[j][0];
			}
		}
		invertMatrix(Am1, RHS, yyn);
	}
	
	//the number of terms to calculate
	int start = len-2;
	double t;
	for(int X=len-1; X>=start; X--){
		t = (1.-r[X]);
		for(int i=0; i<num_var; i++){
			ymode[i][X] = yy[i][0];
			for(int k=1; k<=surface_bc_order; k++) ymode[i][X] += yy[i][k]*pow(t,k);
			//ymode[i][X] += yyn[i]*pow(t, n_surface+1.);
		}
	}
	return start;
}

double PNNonradialModeDriver::DPfunc(int X, double freq, ModeBase* mode){
	double rho  = star->rho(2*X);
	double P    = star->P(  2*X);
	double q    = star->dPhidr(2*X) + star->dPsidr(2*X)/cee2;
	double Phi  = star->Phi(2*X);
	double Phi1 = star->Phi(len_star-1);
	double r    = star->rad(2*X);
	//mode variables
	double chi  = r*q*mode->getY(y2,X);	//r*g*y2
	double DPhi = r*q*mode->getY(y3,X);	//r*g*y3
	double DPsi = r*q*mode->getY(z1,X)*Phi1;
	double DWH  =   q*mode->getY(z5,X)*Phi1/(4.*freq);//divive by i
	double DP   = (	(P+rho*cee2)*(chi-DPhi)
						- rho*(DPsi + 2.*Phi*DPhi + 4.*r*freq*DWH +2.*Phi*chi)
					)/(cee2+2.*Phi);
	return DP;
}

double PNNonradialModeDriver::SSR(double omega2, int l, ModeBase* mode) {
	double checkCont=0.0;
	double checkNewt=0.0;
	double checkPois=0.0;
	double e1,e2,e3;
	double n1,n2,n3;
	double d2DPhi, difxi, dDP;
	double Phi1 = star->Phi(len_star-1);
	for(int X=6; X<len-7; X++){
		int XX = 2*X;
		//stellar variables, to simplify equations
		double rho    = star->rho(XX);
		double P      = star->P(XX);
		double drhodr = star->drhodr(XX);
		double dPdr   = star->dPdr(XX);
		double g      = star->dPhidr(XX);
		double gPN    = star->dPsidr(XX);
		double q      = g + gPN/cee2;
		double Phi    = star->Phi(XX);
		double r      = star->rad(XX);
		double G1     = (Gamma1()==0.0 ? star->Gamma1(XX) : Gamma1());
		double A      = star->Schwarzschild_A(XX,G1);
		//mode variables
		double freq2 = omega2*C[XX]*q/r;
		double freq  = sqrt(freq2);
		double L2    = double(l*l+l);
		double xi   = r  *mode->getY(y1,X);//r  *y1
		double chi  = r*q*mode->getY(y2,X);//r*g*y2
		double DPhi = r*q*mode->getY(y3,X);//r*g*y3
		double dDPhi=   q*mode->getY(y4,X);//  g*y4
		double DPsi = r*q*mode->getY(z1,X)*Phi1;
		double dDPsi=   q*mode->getY(z2,X)*Phi1;
		double DWR  = q*Phi1/(4.*freq)*mode->getY(z3,X);//divide by i
		double DWH  = q*Phi1/(4.*freq)*mode->getY(z5,X);//divide by i
		double Drho = 
			-(rho+P/cee2)*A*xi
			+(rho+P/cee2)*(
				chi*(   P + rho*(cee2 - 2.*Phi))
				- DPhi*(P + rho*(cee2 + 2.*Phi))
				- DPsi*rho
				-4.*r*freq*rho*DWH
			)/(P*G1*(cee2+2.*Phi));
		double DP = G1*P/(rho+P/cee2)*Drho + G1*P*xi*A;
		
		//calculate numerical derivatives
		double b3,b2,b1,a1,a2,a3, h1;
		h1= star->rad(XX+2)-star->rad(XX );
		// d2DPhi/dr2 = (d/dr)(dDphi/dr), dDPhi/dr = q y4
		b3=(star->dPhidr(XX-6)+star->dPsidr(XX-6)/cee2)*mode->getY(y4,X-3);
		b2=(star->dPhidr(XX-4)+star->dPsidr(XX-4)/cee2)*mode->getY(y4,X-2);
		b1=(star->dPhidr(XX-2)+star->dPsidr(XX-2)/cee2)*mode->getY(y4,X-1);
		a1=(star->dPhidr(XX+2)+star->dPsidr(XX+2)/cee2)*mode->getY(y4,X+1);
		a2=(star->dPhidr(XX+4)+star->dPsidr(XX+4)/cee2)*mode->getY(y4,X+2);
		a3=(star->dPhidr(XX+6)+star->dPsidr(XX+6)/cee2)*mode->getY(y4,X+3);
		d2DPhi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
		//next dxi/dr, xi = r y1
		b3=star->rad(XX-6)*mode->getY(y1,X-3);
		b2=star->rad(XX-4)*mode->getY(y1,X-2);
		b1=star->rad(XX-2)*mode->getY(y1,X-1);
		a1=star->rad(XX+2)*mode->getY(y1,X+1);
		a2=star->rad(XX+4)*mode->getY(y1,X+2);
		a3=star->rad(XX+6)*mode->getY(y1,X+3);
		difxi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
		//now dDP/dr
		b3 = DPfunc(X-3, freq, mode);
		b2 = DPfunc(X-2, freq, mode);
		b1 = DPfunc(X-1, freq, mode);
		a1 = DPfunc(X+1, freq, mode);
		a2 = DPfunc(X+2, freq, mode);
		a3 = DPfunc(X+3, freq, mode);
		dDP = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
		//Now calculate residuals
		//Perturbed Poisson equation
		e1 = fabs(
				4.0*m_pi*Gee*Drho
				+ L2*pow(r,-2)*DPhi 
				- 2.0*dDPhi/r  - d2DPhi
		);
		n1 = fabs(4.0*m_pi*Gee*Drho )
				+ fabs( L2*pow(r,-2)*DPhi )
				+ fabs( 2.0*dDPhi/r ) + fabs( d2DPhi )
		;
		//continuity equation
		e2 = fabs(
				Drho*(1.-2.*Phi/cee2) - 3.*rho*DPhi/cee2
				+ xi*(drhodr+dPdr/cee2-2.*drhodr*Phi/cee2-2.*rho*g/cee2)
				+ (rho+P/cee2-2.*rho*Phi/cee2)*( 2.*xi/r + difxi - L2*pow(r,-2)/freq2*chi )		
		);
		n2 = fabs(Drho) + fabs(2.*Drho*Phi/cee2) + fabs(3.*rho*DPhi/cee2)
				+ fabs(xi*drhodr)+fabs(xi*dPdr/cee2)+2.*fabs(xi*drhodr*Phi/cee2) + 2.*fabs(xi*rho*g/cee2)
				+ ( fabs(rho+P/cee2)+fabs(2.*rho*Phi)/cee2 )*( fabs(2.*xi/r)+fabs(difxi)+fabs(L2*pow(r,-2)/freq2*chi) )
		;
		//newton's equation -- the r component. The theta component defines chi
		e3 = fabs( -(rho+(P-4.*rho*Phi)/cee2)*freq2*xi + dDP + rho*dDPhi + g*Drho
 				+ ( rho*dDPsi + Drho*gPN + P*dDPhi + DP*g + 4.*rho*freq*DWR )/cee2 );
 		n3 = ( fabs(rho+P/cee2)+fabs(4.*rho*Phi/cee2) )*fabs(freq2*xi) + fabs( dDP ) + fabs( rho*dDPhi ) + fabs( g*Drho )
 				+ ( fabs(rho*dDPsi) + fabs(Drho*gPN) + fabs(P*dDPhi) + fabs(DP*g) + fabs(4.*rho*freq*DWR) )/cee2;
		//normalize residuals
		e1 = e1/n1;
		e2 = e2/n2;
		e3 = e3/n3;
		//collect residuals
		checkPois += e1*e1; 
		checkCont += e2*e2; 
		checkNewt += e3*e3;
	}
	return sqrt((checkCont+checkPois+checkNewt)/double(3*len));
}

double PNNonradialModeDriver::tidal_overlap(ModeBase* mode){
	double omega2 = mode->getOmega2();
	int k,l,m;
	mode->modeNumbers(k,l,m);
	double L2 = double(l*(l+1));
	double Mstar=star->Mass(), Rstar=star->Radius();
	double dscale = Mstar*pow(Rstar,-3);
	double rscale = 1.0;
	double freq2 = omega2*star->Gee()*Mstar*pow(Rstar,-3);
	double xir=0.0, xiH=0.0, rho=0.0, rlp1=0.0;
	double dx;
	int xx;
	//need to numerically integrate
	double sum1=0.0, sum2=0.0, integral = 0.0;
	double N1=0.0, N2=0.0, NN=0.0;
	for(int x=1; x<len; x++){
		xx = 2*x;
		dx = (r[x]-r[x-1])*rscale;
		xir = mode->getY(y1,x)*r[x]*Rstar;
		xiH = mode->getY(y2,x)*(star->dPhidr(xx)+star->dPsidr(xx)/cee2)/freq2;
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

double PNNonradialModeDriver::innerproduct(ModeBase* mode1, ModeBase* mode2){
	double omega21 = mode1->getOmega2();
	double omega22 = mode2->getOmega2();
	int k1,l1,m1, k2,l2,m2;
	mode1->modeNumbers(k1,l1,m1);
	mode2->modeNumbers(k2,l2,m2);
	if(l2!=l1) return 0.0;
	double L2 = double(l1*l1+l1);
	double Mstar=star->Mass(), Rstar=star->Radius();
	double dscale = 1.0;
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
		dx = (r[x]-r[x-1])*rscale;
		xir1 = mode1->getY(y1,x)*r[x]*Rstar;
		xir2 = mode2->getY(y1,x)*r[x]*Rstar;
		xiH1 = mode1->getY(y2,x)*(star->dPhidr(xx)+star->dPsidr(xx)/cee2)/freq21;
		xiH2 = mode2->getY(y2,x)*(star->dPhidr(xx)+star->dPsidr(xx)/cee2)/freq22;
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