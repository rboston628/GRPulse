//**************************************************************************************
//						POST-NEWTONIAN CHANDRASEKHAR WHITE DWARF
// PNChandrasekharWD.h
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar (1935)
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//  	Builds off of work in Tooper (1964) to implement CHWD in GR
//		This model assumes T=0, and ignores Coulombic and other effects
//		The surface is not treated in any special way
//		Updated to include a chemical profile, indicated by mu_electron
//
//		Based on my own work with perturbed 1PN equations of stellar structure
// 		We use units where G=c=1
//**************************************************************************************

#ifndef PNChandrasekharWDCLASS
#define PNChandrasekharWDCLASS

#include "PNChandrasekharWD++.h"

/*void PNChandrasekharWD::chemical_gradient(
	const double xi, 	// the radial position
	const double dx, 	// the radial step size -- used to scale total radius
	double& mu, 		// return for the chemical potential
	double& dmu			// return for the derivative of the chemical potential
){
	//if the relative core is set to 1, the chem gradient is constant
	if(acore>=1.){
		mu  = mu0;
		dmu = 0.0;
		return;
	}
	//otherwise calculate the chem gradient as a sigmoidal curve from mu0 to 1
	double xcore = dx*(acore*len);	// position of the 
	double xswap = dx*(aswap*len);	// position of sigmoid midpoint
	double EXP = exp(k*(xi-xswap));	// the expoential
	// within the core, constant distribution
	if(xi<xcore){
		mu  = mu0;
		dmu = 0.0; 
	}
	// outside the core, sigmoidal distribution
	else {
		mu = (mu0-1.) + 1.0/(1.+EXP);
		dmu= -k*EXP*pow(1.+EXP,-2);
	}
}
//*/

void PNChandrasekharWD::chemical_gradient(
	const double x, 	// the degeneracy factor
	const double dxdxi, // the derivative dx/dxi, needed to find dmue/dxi
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
	//   dmudxi = dmudx * dxdxi
	//dmu = k/x * EXP*pow(1.+EXP,-2) * dxdxi;
	// IF USING X
	//   dmudx = -k * EXP/(1+EXP)^2
	//   dmudxi = dmudx * dxdxi
	//dmu = -k *  EXP*pow(1.+EXP,-2) * dxdxi;
	// IF USING F
	//   dmudx = - k f'(x) * EXP/(1+EXP)^2 = - 8kx^4/y * EXP/(1+EXP)^2
	//   dmudxi = dmudx * dxdxi = -8kx^4/y * EXP/(1+EXP)^2 *dxdxi
	//dmu  = -8.*k*pow(x,4)/sqrt(1.+x*x) * EXP*pow(1.+EXP,-2) * dxdxi;
	// IF USING LOG F
	//   dmudx = k f'(x)/f(x) * EXP/(1+EXP)^2  = 8k x^4/y/f * EXP/(1+EXP)^2
	//   dmudxi = dmudx * dxdxi  = 8k x^4/y/f * dx/xi * EXP/(1+EXP)^2
	dmu = 8.*k*pow(x,4)/sqrt(1.+x*x)/F * EXP*pow(1.+EXP,-2) * dxdxi;
	//printf("%le %le %le %le %le\n", x, Chandrasekhar::factor_f(x), z, EXP, mu); 	
}

void PNChandrasekharWD::basic_setup(){
	//exclude unphysical values of Y0
	if(Y0 < 1.) Y0 = 1.;
	setX0Y0(Y0, setY0);

	//now set physical properties
	Rn = sqrt(2.*A0/(m_pi*GG))/B0;
	// sigma is determined by A0, B0, and the units
	sigma = A0/(B0*cee2);
	printf("sigma=%le\n", sigma);
	
	//name this star for files
	sprintf(name, "PNChandrasekharWD%1.3f", Y0);
	
	Y = new double*[len];
	for(int i=0; i<len; i++)
		Y[i] = new double[numvar];
}

void PNChandrasekharWD::setX0Y0(double xoyo, X0Y0setState set){
	switch(set){
		case setX0:
			X0 = xoyo;
			X02 = X0*X0;
			Y02 = 1.0 + X02;
			Y0 = sqrt(Y02);
		break;
		case setY0:
			Y0 = xoyo;
			Y02 = Y0*Y0;
			X02 = Y02 - 1.0;
			X0 = sqrt(X02);
		break;
		default:
			X0 = Y0 = X02 = Y02 = 0.0;
	}
}

void PNChandrasekharWD::init_arrays(){
	//calculate derivative and other useful things
	mass = new double[len];
	mue  = new double[len];
	dmue = new double[len];
	dxdxi= new double[len];
	dfdx = new double[len];
	dhdx = new double[len];
	double eggs = Y[0][x];
	for(int X=0; X<len; X++){
		eggs = Y[X][x];
		chemical_gradient(Y[X][x], dxdxi[X], mue[X], dmue[X]);
		Y[X][y]  = sqrt(1.0 + eggs*eggs);
		Y[X][f]  = Chandrasekhar::factor_f(eggs);
		Y[X][h]  = Chandrasekhar::factor_h(eggs, sigma, mue[X]);
		dxdxi[X] = -((Y[X][h] + sigma*Y[X][f])*Y[X][u] + sigma*Y[X][h]*Y[X][v])*Y[X][y]*pow(eggs,-4);
		chemical_gradient(Y[X][x], dxdxi[X], mue[X], dmue[X]);
		dfdx[X]  = 8.*pow(eggs,4)/Y[X][y];
		dhdx[X]  = 3.*eggs*eggs*(mue[X] + 8.*sigma*(Y[X][y]-1.));
		mass[X]  = 4.*m_pi*Y[X][u]*pow(Y[X][xi],2);
	}
	dxdxi[len-1] = -mue[len-1]*(Y[len-1][u] + sigma*Y[len-1][v]);
}

//initalize polytrope from index and length
// note that in this model sigma is fixed by units (mass of proton/electron)
//	and that zsurf is determined uniquely from Y0
PNChandrasekharWD::PNChandrasekharWD(
	double Y0, int Len,
	double mu0, double k, double acore, double aswap
)
	: Y0(Y0), len(Len),
	A0(Chandrasekhar::A0), B0(Chandrasekhar::B0),
	GG(G_CGS), cee2(C_CGS*C_CGS),
	//A0(6.799403e-05), B0(1.0), GG(1.0), cee2(1.0),
	mu0(mu0), k(k), acore(acore), aswap(aswap)
{
	basic_setup();
	
	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with y[len-1]=0.0
	//we will find the proper dx with a bisection search
	
	//an initial guess -- based on 0pn constant-volume model
	dxi = sqrt(6.0)/(len-1);
	phic0 = -1.0e-4/sigma;
	
	size_t const np = 2;
	double z1[np] = {phic0,     dxi};
	double dz[np] = {0.1*z1[0], dxi};
	double f1[np];
	RK4integrate(z1[0], z1[1]);
	RK4integrate(z1[0], z1[1]);
	RK4integrate(z1[0], z1[1]);

	double target[np] = {0.,0.};
	// require the two expressions for zsurf to agree at the surface, and for x to vanish at surface
	std::function<void(double[np],double[np])> match_zsurf = [this](double f[np], double z[np]){
		this->RK4integrate(z[0], z[1]);
		f[0] = this->Y[len-1][phi] + this->Y[len-1][u]*this->Y[len-1][xi];
		f[1] = this->Y[len-1][x];
	};
	// phic0 is required to be negative, while dxi cannot be negative
	std::function<bool(double[np])> sign_limit = [](double z[np])->bool{
		return (z[0]<0.0) && (z[1]>0.0);
	};
	//now run the newtonian search
	newton_search<np>(match_zsurf, z1, dz, 1.e-8, 100, sign_limit);
	phic0 = z1[0];
	dxi = z1[1];
	
	RK4integrate(len, dxi, (int)1);
	zsurf = setZsurf(Y[len-1]);
	
	printf("dxi=%le\n", dxi);
	printf("sigma = \t%le\n", sigma);
	printf("zsurf = \t%0.16le\t%0.16le\n", zsurf, 8.*sigma*Y[len-1][u]*Y[len-1][xi]);
	printf("edge  = \t%0.16le\n", Y[len-1][x]);
	
	init_arrays();
	indexFit = len/2;
	for(int X=1; X<len-1; X++){
		//as we scan through x,y,z, set matching point as where x = 0.5
		if(Y[X-1][x]>0.5 & Y[X+1][x]<=0.5) indexFit = X;
	}
	indexFit /= 2;
	printf("xfit = %d\n", indexFit);
	setupCenter();
	setupSurface();
	printf("%0.2lf\t%le\t%le\n", 1./Y02, Mass()/MSOLAR, Radius());
}


//initalize polytrope from index, length, and dx
PNChandrasekharWD::PNChandrasekharWD(
	double Y0, int Len, const double dxi,
	double mu0, double k, double acore, double aswap
)
	: Y0(Y0), len(Len), dxi(dxi),
	A0(1.0e-4), B0(1.0), GG(1.0), cee2(1.0),
	mu0(mu0), k(k), acore(acore), aswap(aswap)
{
	basic_setup();
	
	phic0 = -1.0e-4/sigma;
	double pmin=-10.0, pmax=0.0;
	std::function<double(double)> match_phi = [this](double p)->double{
		this->phic0 = p;
		this->RK4integrate(this->len, this->dxi, (int)1);
		return Y[len-1][phi]+Y[len-1][u]*Y[len-1][xi];
	};
	bisection_find_brackets_move(match_phi, phic0, pmin, pmax);
	bisection_search(match_phi, phic0, pmin, pmax);
	
	RK4integrate(len, dxi, (int)1);
	//continue adjusting the value of sigma until the surface redshift is correct
	/*int tracker = 0;
	while(fabs(Y[len-1][phi] + zsurf/8./sigma) > 1e-16  & tracker<100){
		phic0 -= (zsurf/8./sigma + Y[len-1][phi]);
		RK4integrate(len, dxi,1);
		tracker++;
	}//*/

	//calculate derivative and other useful thngs
	init_arrays();
	indexFit = 512*round(double(len)/1024.0);
	indexFit /= 2;
	printf("  indexFit  = %d\n", indexFit);
	printf("r[indexFit] = %0.32le\n", rad(indexFit));
	setupCenter();
	setupSurface();
	printf("%0.2lf\t%le\t%le\n", 1./Y02, Mass()/MSOLAR, Radius());
}


PNChandrasekharWD::~PNChandrasekharWD(){
	for(int i=0; i>len; i++)
		delete[] Y[i];
	delete[] dxdxi;
	delete[] mass;
	delete[] mue;
	delete[] dmue;
	delete[] dfdx;
	delete[] dhdx;
}

void PNChandrasekharWD::centerInit(double ycenter[numvar]){
	ycenter[xi]  = 0.0;				//radius in center
	ycenter[x]   = X0;				//central X0, determiend from Y0
	ycenter[y]   = Y0;				//central Y0, specified by model
	ycenter[f]   = Chandrasekhar::factor_f(X0);				// f(X0)
	ycenter[h]   = Chandrasekhar::factor_h(X0,sigma,mu0);	// h(X0) w/ correct sigma, mu
	ycenter[phi] = phic0;			//central value of phi, so that match at surface
	//these derivatives are zero at center
	ycenter[u]   = 0.0;
	ycenter[v]   = 0.0;
}

double PNChandrasekharWD::setZsurf(double ysurface[numvar]){
	// (1 + z)^2 = exp(la) = 1/[1 - 2Gm/c^2r]
	// zsurf = [1 - 2GM/c^2R]^(-1/2) - 1 
	// zsurf = exp(la/2) - 1
	return -8.*sigma*ysurface[phi];
}

void PNChandrasekharWD::RK4step(double dx, double yin[numvar], double yout[numvar]){
	double YC[numvar];
	double K[4][numvar];
	static const double B[4] = {0.5,0.5,1.0,0.0};
	
	double MUC = 1.0, DMUC = 0.0, DX = 0.0;
	for(int b=xi; b<=phi; b++)
		YC[b] = yin[b];
	chemical_gradient(YC[x], DX, MUC, DMUC);
	YC[y] = sqrt(1.0 + yin[x]*yin[x]);
	YC[f] = Chandrasekhar::factor_f(yin[x]);
	YC[h] = Chandrasekhar::factor_h(yin[x], sigma, MUC);
	DX = (  (YC[h] + sigma*YC[f])*YC[u] + sigma*YC[h]*YC[v] )*YC[y]*pow(YC[x],-4);
	chemical_gradient(YC[x], DX, MUC, DMUC);
	
	for(int a = 0; a<4; a++){
		// find the shift vector
		K[a][xi] = dx;
		K[a][x] =-dx*(  (YC[h] + sigma*YC[f])*YC[u] + sigma*YC[h]*YC[v] )*YC[y]*pow(YC[x],-4);
		K[a][u] = dx*(   YC[h]                     - 2.*YC[u]/YC[xi] );
		K[a][v] = dx*(3.*YC[f] - 16.*YC[h]*YC[phi] - 2.*YC[v]/YC[xi] ); 
		K[a][phi] = dx*YC[u];
		// possible divide-by-zero errors at pure middle
		if(YC[xi]==0) {
			K[0][u] = dx/3.0*YC[h];
			K[0][v] = dx/3.0*(3.*YC[f] - 16.*YC[h]*YC[phi]);
		}
		if(YC[x]==0.0) {
			K[a][x] = -dx/3.0*MUC*(YC[u]+sigma*YC[v]);
		}
		//prepare for next step
		for(int b=xi; b<=phi; b++)
			YC[b] = yin[b] + B[a]*K[a][b];
		chemical_gradient(YC[x], DX, MUC, DMUC);
		YC[y] = sqrt(1.0 + YC[x]*YC[x]);
		YC[f] = Chandrasekhar::factor_f(YC[x]);
		YC[h] = Chandrasekhar::factor_h(YC[x], sigma, MUC);
		DX = (  (YC[h] + sigma*YC[f])*YC[u] + sigma*YC[h]*YC[v] )*YC[y]*pow(YC[x],-4);
		chemical_gradient(YC[x], DX, MUC, DMUC);
	}
	yout[xi] = yin[xi] + dx;
	for(int b=x; b<=phi; b++)
		yout[b] = yin[b] + K[0][b]/6.0 + K[1][b]/3.0 + K[2][b]/3.0 + K[3][b]/6.0;
	chemical_gradient(yout[x], DX, MUC, DMUC);
	yout[y] = sqrt(1.0 + yout[x]*yout[x]);
	yout[f] = Chandrasekhar::factor_f(yout[x]);
	yout[h] = Chandrasekhar::factor_h(yout[x], sigma, MUC);
	DX = (  (yout[h] + sigma*yout[f])*yout[u] + sigma*yout[h]*yout[v] )*yout[y]*pow(yout[x],-4);
	chemical_gradient(yout[x], DX, MUC, DMUC);
}

//integrate the polytrope up to Len using RK4
double PNChandrasekharWD::RK4integrate(double p0, double &dx, bool f){
	double y1[numvar], y2[numvar];
	
	//set the initial conditions
	centerInit(y2);
	y2[phi] = p0;
	
	int X = 0;
	while(y2[x] > 0.0 && !isnan(y2[x])){
		for(int b=0; b<numvar; b++) y1[b] = y2[b];
		RK4step(dx, y1, y2);
		X++;
		if(X>2+len) break;
	}
	if(!f) printf("%d %le %le %le\n", X, dx, y1[x], setZsurf(y1));
	if(f && y1[xi] != 0){
		dx = y1[xi]/(len-1);
	}
	for(int b=0; b<numvar; b++){
		Y[len-1][b] = y1[b];
	}
	return setZsurf(y1);
}

//integrate the polytrope up to Len using RK4
double PNChandrasekharWD::RK4integrate(const int Len, double dx){

	//set our initial conditions
	centerInit(Y[0]);
	
	//integrate with RK4
	int X = 0;
	for(X = 0; X<Len-1; X++){
		RK4step(dx, Y[X], Y[X+1]);
		//sometimes values y<0 lead to nans
		//if these occur, it is safe to terminate the integration
		if(Y[X+1][x] <0.0 || isnan(Y[X+1][x])){
			return Y[X+1][x]; //dx too large should give negative
		}
	}
	return Y[Len-1][x];
}

int    PNChandrasekharWD::RK4integrate(const int Len, double dx, int grid){
	if(grid<1) grid=1;

	//set our initial conditions
	centerInit(Y[0]);
	//integrate with RK4
	int X = 0;
	for(X = 0; X<Len-2; X++){
		RK4step(dx, Y[X], Y[X+1]);
	}
	return Len;
}
	

//Here we define functions to access radius, pressure, etc.
double PNChandrasekharWD::rad(int X){
	return Rn*Y[X][xi];
}
double PNChandrasekharWD::rho(int X){
	return B0*Y[X][h];
}
double PNChandrasekharWD::drhodr(int X){
	return B0*( dhdx[X]*dxdxi[X] + pow(Y[X][x],3)*dmue[X] )/Rn;
}
double PNChandrasekharWD::P(int X){
	return A0*Y[X][f];
}
double PNChandrasekharWD::dPdr(int X){
	return A0*dfdx[X]*(dxdxi[X]/Rn);
}
double PNChandrasekharWD::Phi(int X){
	return 8.*A0/B0*Y[X][phi];
}
double PNChandrasekharWD::dPhidr(int X){
	return 8.*A0/B0*(Y[X][u]/Rn);
}
double PNChandrasekharWD::mr(int X){
	//m = 4 pi Int[r^2 rho,r] 
	//  = 4 pi Rn^3 B Int[xi^2 h(x),xi] 
	//  = 4 pi Rn^3 B xi^2 dphi/dxi
	return B0*pow(Rn,3)*mass[X];
}
double PNChandrasekharWD::Psi(int){
	return 0.0; //TO BE CONTINUED...
}
double PNChandrasekharWD::dPsidr(int X){
	return 8.*A0*A0/B0/B0*(Y[X][v]/Rn);
}
//the gravitomagnetic potential vansihes for non-rotating star
double PNChandrasekharWD::Wx(int X){return 0.0;}
double PNChandrasekharWD::dWxdr(int ){return 0.0;}
double PNChandrasekharWD::Wy(int ){return 0.0;}
double PNChandrasekharWD::dWydr(int ){return 0.0;}
double PNChandrasekharWD::Wz(int ){return 0.0;}
double PNChandrasekharWD::dWzdr(int ){return 0.0;}

double PNChandrasekharWD::Gamma1(int X){
	if(Y[X][f]==0.0) return  (Y[X-1][h]+sigma*Y[X-1][f])/Y[X-1][f] * dfdx[X-1]/dhdx[X-1];
	return (Y[X][h]+sigma*Y[X][f])/Y[X][f] * dfdx[X]/dhdx[X];
}

double PNChandrasekharWD::sound_speed2(int X, double GamPert){
	// vs2 = Gamma1* P/(rho+P/c^2) 
	//     = Gamma1* A0 * f/(B0 h + A0 f/c2)
	//     = Gamma1* (A0/B0) * f/(h + sigma f)
	if(GamPert==0.0) return Gamma1(X)*A0/B0*Y[X][f]/(Y[X][h]+sigma*Y[X][f]);
	else             return GamPert  *A0/B0*Y[X][f]/(Y[X][h]+sigma*Y[X][f]);
}

//Schwarzschild discriminant as for GR case
double PNChandrasekharWD::Schwarzschild_A(int X, double GamPert){
	double A;
	if(GamPert==0.0) A = pow(Y[X][x],3)/(Y[X][h]+sigma*Y[X][f])*dmue[X]/Rn;
	else             A = pow(Y[X][x],3)/(Y[X][h]+sigma*Y[X][f])*dmue[X]/Rn
							+ (dhdx[X]/(Y[X][h]+sigma*Y[X][f]) - dfdx[X]/GamPert/Y[X][f])*dxdxi[X]/Rn;
	return A;
}

double PNChandrasekharWD::getAstar(int X, double GamPert){
	double AS;
	if(GamPert==0.0) AS = - pow(Y[X][x],3)/(Y[X][h]+sigma*Y[X][f])*Y[X][xi]*dmue[X];
	else        	 AS = - pow(Y[X][x],3)/(Y[X][h]+sigma*Y[X][f])*Y[X][xi]*dmue[X]
						  + (dfdx[X]/(GamPert*Y[X][f]) - dhdx[X]/(Y[X][h]+sigma*Y[X][f]))*Y[X][xi]*dxdxi[X];				  
	return AS;
}


double PNChandrasekharWD::getU(int X){
	double U;
	if(X==0) U = 3.0;
	else     U = Y[X][xi]*(3.*sigma*Y[X][f] + Y[X][h]*(1. - 16.*sigma*Y[X][phi]))/(Y[X][u]+sigma*Y[X][v]);
	return U;
}

double PNChandrasekharWD::getVg(int X, double GamPert){
	double V;
	if(GamPert==0.0) V = - Y[X][xi]/Gamma1(X)/Y[X][f] *dfdx[X]*dxdxi[X];
	else			 V = - Y[X][xi]/GamPert  /Y[X][f] *dfdx[X]*dxdxi[X];
	return V;
}

double PNChandrasekharWD::getC(int X){
	if(X==0) return 0.5*Y[len-1][u]*Y[len-1][xi]/(φc[1]+sigma*ψc[1]);
	return Y[len-1][u]/(Y[X][u]+sigma*Y[X][v])  *  Y[X][xi]/Y[len-1][xi];
}

double PNChandrasekharWD::Radius(){return rad(len-1);}	//total radius
double PNChandrasekharWD::Mass(){return mr(len-1);}//total mass

double PNChandrasekharWD::SSR(){
	double checkEuler=0.0, checkPoiss=0.0, checkPsi=0.0;
		
	//sum up errors in equations
	double r = 0.0;
	double R = Radius();
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
		e1 = fabs( 4.*m_pi*GG*rho(X)*r - 2.*dPhidr(X) - d2Phi*r );
		n1 = fabs( 4.*m_pi*GG*rho(X)*r ) + fabs( 2.*dPhidr(X) ) + fabs(d2Phi*r);
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
		e3 = fabs( 12.*m_pi*GG*P(X) - 8.*m_pi*GG*rho(X)*Phi(X) - 2.*dPsidr(X)/r - d2Psi ); 
		n3 = fabs( 12.*m_pi*GG*P(X) )
				+ fabs( 8.*m_pi*GG*rho(X)*Phi(X) )
				+ fabs( 2.*dPsidr(X)/r) + fabs(d2Psi);
		e1 /= n1;
		e2 /= n2;
		e3 /= n3;
		checkPoiss += e1*e1;
		checkEuler += e2*e2;
		checkPsi   += e3*e3;	
	}
	
	return sqrt((checkPoiss + checkEuler + checkPsi)/double(3*len));
}//*/

// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
void PNChandrasekharWD::setupCenter(){
	//the central coefficients -- expanded in terms of r/R
	double z1 = dxdxi[len-1]*Y[len-1][x]/Y[len-1][y];
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


	//the central coefficients -- expanded in terms of xi
	double xi1  = Y[len-1][xi];
	double xi12 = xi1*xi1;
	//first order
	//xc[0] = X0;
	//yc[0] = Y0;
	//fc[0] = Chandrasekhar::factor_f(X0);
	hc[0] = Chandrasekhar::factor_h(X0, sigma, mu0) ;
	φc[0] = phic0;
	ψc[0] = 0.0;
	//second order
	//xc[1] = xi12*(4.*hc[0]*yc[0]*sigma*(4.*hc[0]*φc[0]-fc[0]) - hc[0]*hc[0]*yc[0])/6.*pow(X0,-4);
	//yc[1] = X0*xc[1]/Y0;
	//fc[1] = 8.*(xc[1]/yc[0] - yc[0]*xc[1] + xc[0]*xc[0]*xc[1]*yc[0]);
	hc[1] = 3.*X02*xc[1] + 24.*sigma*X02*xc[1]*(Y0-1.);
	φc[1] = xi12/6.*hc[0];
	ψc[1] = xi12/6.*(3.*fc[0] - 16.*hc[0]*φc[0]);
	//fourth order
	φc[2] = xi12/20.*hc[1];
	ψc[2] = xi12/20.*(3.*fc[1]- 16.*hc[1]*φc[0] - 16.*hc[0]*φc[1]);
	//*/
	
}

/*
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
	double x12 = Y[len-1][xi]*Y[len-1][xi];
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] =-16.*X02*X02*xc[1]/(Gam1*fc[0]*Y0);
	if(maxPow> 2) maxPow = 2; 
}
void PNChandrasekharWD::getUCenter(    double *Uc, int& maxPow){
	double x12 = Y[len-1][xi]*Y[len-1][xi];
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] = 6.*hc[1]/(5.*hc[0])
		- 2.*sigma*(8.*pow(hc[0],3)*x12 - 9.*hc[0]*fc[1] + 9.*fc[0]*hc[1])/(5.*hc[0]*hc[0]);
	if(maxPow> 2) maxPow = 2;
}
//*/
void PNChandrasekharWD::getC1Center(   double *cc, int& maxPow){
	double u1x1 = Y[len-1][u]*Y[len-1][xi];
	if(maxPow>=0) cc[0] = 0.5*u1x1/(φc[1]+sigma*ψc[1]);
	if(maxPow>=2) cc[1] =    -u1x1*(φc[2]+sigma*ψc[2])*pow(φc[1]+sigma*ψc[1],-2);
	if(maxPow> 2) maxPow = 2;
}


void PNChandrasekharWD::getBetaCenter( double *bc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double phi1 = Y[len-1][phi];
	if(maxPow>=0) bc[0] = Gam1* fc[0]/(hc[0]+sigma*fc[0])/8./phi1;
	if(maxPow>=2) bc[1] = Gam1*(fc[1]*hc[0]-fc[0]*hc[1])*pow(hc[0]+sigma*fc[0],-2)/8./phi1;
	//second order not appear in equations
	if(maxPow> 2) maxPow = 2; 
}
void PNChandrasekharWD::getPhiCenter(  double *pc, int& maxPow){
	if(maxPow>=0) pc[0] =-8.*sigma*phic0/zsurf;
	if(maxPow>=2) pc[1] =-8.*sigma*φc[2]/zsurf; //this order does not appear in equations
	if(maxPow> 2) maxPow = 2;
}



void PNChandrasekharWD::getAstarCenter(double *Ac, int& maxPow, double g){
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

void PNChandrasekharWD::getVgCenter(double *Vc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double x1 = Y[len-1][xi];
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] = -16.*X02*X0*yc[1]/Gam1/fc[0];
	if(maxPow>=4) Vc[2] =-8.*   pow(x1,2)*pow(X0,6)*(5.*fc[1]+4.*pow(x1,2)*fc[0]*X0*Y0)/(15.*Gam1*fc[0]*fc[0]);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4; 
}

void PNChandrasekharWD::getUCenter(double *Uc, int& maxPow){
	double x1 = Y[len-1][xi];
	double mu2 = mue[0]*mue[0];
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] =-0.6*mu2*x1*x1*X0*Y0;
	if(maxPow>=4) Uc[2] = 
		9.*xc[2]/X0 - pow(x1*x1*mu2,2)*X02*(0.08+0.134*X02);
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
void PNChandrasekharWD::setupSurface(){
	double xi1 =  Y[len-1][xi];
	double xs1 = -dxdxi[len-1]*xi1;
	double φs1 = - Y[len-1][u]*xi1;
	double ψs1 = - Y[len-1][v]*xi1;
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

/*
void PNChandrasekharWD::getAstarSurface(double *As, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-2) : g);
	double AV = 3. - 5./Gam1;
	int O=1;
	if(maxPow>=-1) As[O-1] = AV;
	if(maxPow>= 0) As[O  ] = AV*(xs[2]-xs[1])/xs[1];
	if(maxPow>= 1) As[O+1] = AV*( -xs[1]*xs[1]/7. - pow(xs[2]/xs[1],2) - (xs[2]-2.*xs[3])/xs[1] )
							  + xs[1]*xs[1]*(3./7.);
	if(maxPow > 1) maxPow = O+1;
}
void PNChandrasekharWD::getVgSurface(double *Vs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-2) : g);
	double VG = 5./Gam1;
	int O=1;
	if(maxPow>=-1) Vs[O-1] = VG;
	if(maxPow>= 0) Vs[O  ] = VG*(xs[2]-xs[1])/xs[1];
	if(maxPow>= 1) Vs[O+1] =-VG*( 
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
	double x1u1 = Y[len-1][xi]*Y[len-1][u];
	if(maxPow>=0) cs[0]  =  - x1u1/(φs[1]+sigma*ψs[1]);
	if(maxPow>=1) cs[1]  = 3.*x1u1/(φs[1]+sigma*ψs[1]);
	if(maxPow>=2) cs[2]  =-3.*x1u1/(φs[1]+sigma*ψs[1]);
	if(maxPow> 2) maxPow = 2;
}
//*/

void PNChandrasekharWD::getAstarSurface(double *As, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-2) : g);
	double x1 =  Y[len-1][xi];
	double a1 = -mue[len-1]*(Y[len-1][u] + sigma*Y[len-1][v])*x1;
	int O=1;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=-1) As[O-1] = 0.0;//1.5;
	if(maxPow>= 0) As[O  ] = 0.75*a1;
	if(maxPow>= 1) As[O+1] = 0.75*a1 - (8./Gam1 + 0.375)*a1*a1;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 1) maxPow = O+1;
}

void PNChandrasekharWD::getVgSurface(double *Vs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-2) : g);
	double x1 =  Y[len-1][xi];
	double a1 = -mue[len-1]*(Y[len-1][u] + sigma*Y[len-1][v])*x1;
	int O=1;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=-1) Vs[O-1] = 1.5;
	if(maxPow>= 0) Vs[O  ] = 0.75*a1;
	if(maxPow>= 1) Vs[O+1] = 0.75*a1 + 8.*a1*a1/Gam1 - 3.*a1*a1/8.;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 1) maxPow = O+1;
}

void PNChandrasekharWD::getUSurface(double *Us, int& maxPow){
	//these happen to all be zero at surface
	for(int j=0; j<=maxPow; j++){
		Us[j] = 0.0;
	}
}

void PNChandrasekharWD::getC1Surface(double *cs, int& maxPow){
	if(maxPow >= 0) cs[0] =  1.;
	if(maxPow >= 1) cs[1] = -3.;
	if(maxPow >= 2) cs[2] =  3.;
	if(maxPow >= 3) cs[3] = -1.;
	if(maxPow >= 4) cs[4] =  0.; //does not actually appear in equations
	//if more  terms than this requested, cap number of terms
	if(maxPow  > 4) maxPow = 4;
}

void PNChandrasekharWD::getBetaSurface(double *bs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-1) : g);
	double phi1 = Y[len-1][phi];
	if(maxPow>=0) bs[0]  = 0.0;
	if(maxPow>=1) bs[1]  = 0.0;
	if(maxPow>=2) bs[2]  = Gam1*xs[1]*xs[1]/5./phi1;
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
	char pathname[256];
	if(c==NULL)	sprintf(pathname, "./out/%s", name);
	else{
		sprintf(pathname, "./%s/star/", c);
	}
	
	char command[300];
	sprintf(command, "mkdir -p %s", pathname);
	
	printStar(pathname);
	printBV(pathname);
	printCoefficients(pathname);
	printDegeneracy(pathname);
	printChemicalGradient(pathname);
}
	
void PNChandrasekharWD::printDegeneracy(char *pathname){	
	char txtname[256];
	char outname[256];
	sprintf(txtname, "%s/%s.txt", pathname, "degeneracy");
	sprintf(outname, "%s/%s.png", pathname, "degeneracy");

	FILE *fp = fopen(txtname, "w");
	double R = Y[len-1][xi];
	//print results to text file
	for(int X=0; X< length(); X++){
		fprintf(fp, "%d\t%le", X, Y[X][xi]/R);
		fprintf(fp, "\t%le\t%le\t%le\t%le", Y[X][x], Y[X][y], Y[X][f], Y[X][h]);
		fprintf(fp, "\t%le", Y[X][phi]);
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
	fprintf(gnuplot, "set title 'Degeneracy functions for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'rho/rho_c, P/P_c, m/M, g/g_S'\n");
	fprintf(gnuplot, "plot '%s' u 2:3 w l t 'x'", txtname);
	fprintf(gnuplot, ", '%s' u 2:4 w l t 'y'", txtname);
	fprintf(gnuplot, ", '%s' u 2:5 w l t 'f'", txtname);
	fprintf(gnuplot, ", '%s' u 2:6 w l t 'h'", txtname);
	fprintf(gnuplot, "\n"); 
	
	//print fthe phi
	fprintf(gnuplot, "set output '%s/%s.png'\n", pathname, "phi");
	fprintf(gnuplot, "set title 'Newtonian Potential for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'Phi'\n");
	fprintf(gnuplot, "plot ");
	fprintf(gnuplot, " '%s' u 1:7 w l t 'Phi'", txtname);
	fprintf(gnuplot, "\n");
	//
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
}

void PNChandrasekharWD::printChemicalGradient(char *pathname){
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
			(1.-mr(X)/Mass()),
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