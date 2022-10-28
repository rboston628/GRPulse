//**************************************************************************************
//							SIMPLE WHITE DWARF STAR
// SimpleWD.cpp
// Uses simple C/O core with He/H envelope and atmosphere
//  User must specify total H, He, C, and O fractions within the star
//  Layer stratification is determined using the method described in 
//		- Wood 1990
//		- Tassoul, Fontaine & Winget 1990
//		- Arcoragi & Fontaine 1980
//  The interior equation of state includes pressure contributions from
//		P = P_deg + P_ions + P_coulomb + P_rad
//		the electron degeneracy is approximated as T=0 in the core
//		in atmosphere, uses partial degeneracy (see Cox and Giuli)
//	Integrates with Schwarzschild's dimensionless logarithmic variables
// **************************************************************************************/

#ifndef SIMPLEWDCLASS
#define SIMPLEWDCLASS

#include "SimpleWD.h"
#include "ChandrasekharWD++.h"

double radiative_opacity(const StellarVar& ly, const Abundance& X){
//	static const double kappa0 = 4.34e24*(X.C12+X.O16)*(1.+X.H1); //Schwarzschild 1946, Cox&Giuli, Shapiro & Teuskolsky 1983 eq 4.1.8
//	return kappa0*exp(ly[dens]-3.5*ly[temp]);// Shapiro & Teukolsky 1983 eq 4.1.3
//the below is from Hansen & Kawaler (also found in Shapiro & Teukolsky 1983 and Schwarzschild 1958)
	double meanZA = X.H1 + X.He4 + 3.*X.C12 + 4.*X.O16;
	double ke = 0.195*(1+X.H1);
	double kbf = 4.34e25*(X.C12+X.O16)*(1.+X.H1)*exp(ly[dens]-3.5*ly[temp]);
	double kff = 3.68e22*(X.H1 +X.He4 + (X.C12+X.O16)*meanZA)*(1.+X.H1)*exp(ly[dens]-3.5*ly[temp]);
//	return kbf + kff + ke;
//	if(X.H1<1e-3){ //from Iben 1975
		double t6 = ly[temp]-6.*log(10.0);
		double C  = pow(2.019 + 1e-4*exp(ly[dens]-1.7*t6),2.425);
		double Ap = 1. + C*(1.+C/24.55);
		double mu = meanZA - 1.;
		double B  =  3.86 + 0.252*sqrt(mu) + 0.018*mu;
		double A  = (1.25 + 0.488*sqrt(mu) + 0.092*mu)/0.67;
		double lr = -A + B*t6;
//		double k  = 0.67*(ly[dens]-lr) + log(Ap);
		double k = Ap*exp(0.67*(ly[dens]-lr));
		return k + ke + kbf + kff;
//	}//*/
//	else{ //if(X.H1<0.8){
//	double kff = 0.645e23*(X.mean_Z()*X.mean_Z()*X.mu_e()/X.mean_A())*exp(ly[dens]-3.5*ly[temp]);
//	double k_H  = 2.5e-31*(X.C12+X.O16)/0.02*exp(0.5*ly[dens]+9.*ly[temp]);
//	double k_K  = (4.3e25*(X.C12+X.O16) + 3.7e22*(X.H1+X.He4))*(1.+X.H1)*exp(ly[dens]-3.5*ly[temp]);
//	return ke+kff+kbf + k_K + k;
//	return k_bf+k_e+k_ff;//*/
//  the below is from .... Henyey et al 1959, or Bodenheimer 1965
/*	double lt7 = ly[temp] - 7.0*log(10.0), t7 = exp(lt7);
	double rho = exp(ly[dens]), rho1X = rho*(1.+X.H1);
	double k1 = 6.48e-3*X.H1*rho1X*exp(-4.25*lt7);
	double k2 = 186.0*rho1X*X.He4*( 1.05e-3*exp(-3.*lt7)/(1.+10*t7) + (2.-X.He4)*pow(0.15 + t7*(300.0*t7-7.),-4.));
	double k3 = 4.38e3*(X.C12+X.O16)/sqrt( (0.455*sqrt(rho1X)+1.)*(5.38*t7*t7/rho1X+1.)*(1.34e5*exp(5*lt7)/rho1X+1.)*(1.44*t7*t7 + 0.1*exp(-1.5*lt7)));
	return ( (k1+k2+ kbf + kff)*(k3+ke) + k3*ke )/(k1+k2+k3+ke + kbf + kff) + k1 + k2 + k3 + ke;
//	}//*/
//	else{
/*		double Pe = pressure_deg_partial(exp(ly[dens]),exp(ly[temp]), X);
		double V  = exp(-ly[dens]);
		double T  = exp(ly[temp] - 4.*log(10.));
		double k1 = 0.194*(1.+X.H1);//5.4e-13*V/T;
		double k2 = X.H1/(4.5*pow(T,6) + 1./T/( 4.e-3/pow(T,4)+2.e-4*pow(V,0.25) ));
		double k3 = X.He4/(1.4e3*T+pow(T,6));
		double k4 = (X.C12+X.O16)*sqrt(T)/(20.*T+5.*pow(T,4)+pow(T,5));
		return Pe*(k1+(k2+k3+k4)*exp(ly[dens]));
//	}//*/
}

double conductive_opacity(const StellarVar& ly, const Abundance& X){
//	return 4.e-8*pow(X.mu_e()*X.mean_Z(),2)/X.mean_A()*exp(2.*(ly[temp]-ly[dens]));
//	double t7 = ly[temp] - 7.0*log(10.0);
//	return 301.8*exp(-t7)/(1. + pow(0.004364*(1.+X.H1)*exp(ly[dens]-1.5*t7),5./3.) );
// below is based on Cassis 2007, which relies on Yakovlev & Urpin 1980 and Potekhin et al 1999
	double rho = exp(ly[dens]);
	double x = pow(rho/Chandrasekhar::B0/X.mu_e(), 1./3.);
	//the Coulomb factor, which also controls crystallization
	double GC = 2.275e5*(X.mean_coulomb())*exp(ly[dens]/3.-ly[temp]);
	//Coloumb logarithm from Yakovlev&Urpin 1980, eq 11, vastly simplified from Potekhin 1999
	double Coulomb_Logarithm = log(2.*m_pi*X.mean_Z()/3.)/3. + log(1.5+3./GC)/2. - 0.5*x*x/(1.+x*x); 
	// now get opacity, by combining Cassisi 2007 eq 1, 2, 8
	double opaccoeff = 512.*boltzmann_sigma*pow(electron.mass_CGS,2)*pow(electron.charge_CGS,4)
				*proton.mass_CGS/(3.*pow(planck_h_CGS,3)*pow(boltzmann_k,2));
	double opacity = opaccoeff*X.mean_A()*sqrt(1.+x*x)*exp(2.*(ly[temp]-ly[dens]))*Coulomb_Logarithm;
	return opacity;//*/
}

double SimpleWD::opacity(const StellarVar& ly, const Abundance& X){
	double k_rad = radiative_opacity(ly,X);
	double k_cond= conductive_opacity(ly,X);
//	double Pcore = core_pressure(exp(ly[dens]), exp(ly[temp]), X);
//	double Patm  =  atm_pressure(exp(ly[dens]), exp(ly[temp]), X);
//	if(exp(ly[dens]) < 1e-3) k_cond=1e24;
	return 1./(1./k_rad+1./k_cond);
//	else return k_rad;
}

//initalize SimpleWD model
SimpleWD::SimpleWD(
	double M,   //mass, in solar masses
	double Teff,//effective temperature, in kelvin
	int Len
)	: Msolar(M), Teff(Teff), Ntot(Len) 
{
	//begin by assigning values to EOS and chemical abundance based on file
	setup();
	//prepare a starting model from a ChandrasekharWD, to guess R, P0
	initFromChandrasekhar();
	
	//
	if(Ntot%2==0) Ntot = Ntot+1;
	int Ntrue = Ntot;
	Ntot  = (Ntrue+1)/2;
	
	//set indices for different regions
	//the boundary of the core
	Ncore = (int) (0.8*double(Ntot));
	//if( Ncore%2==0) Ncore = Ncore+1;
	Natm = Ntot - Ncore;
	if( Ncore%2==0) {Ncore = Ncore+1; Natm = Natm-1;}
	//the boundary of atmosphre
	//if( Natm%2==0) Natm = Natm+1;
	//Ntot = Ncore + Natm;
	
	double Mcore = 0.99;
	//initilize the log-radial grid
	logQ = new double[Ntot];	//logQ could be either r or m
	setupGrid(Mcore, Ncore);
	//initialize equilibrium structure arrays
	logY = new StellarVar[Ntot];
	dlogY= new StellarVar[Ntot];
	//initialize chemical gradient arrays
	Xelem = new Abundance[Ntot];
	dXelem= new Abundance[Ntot];
	for(int n=0; n<Ntot; n++) Xelem[n] = findAbundance(logQ[n], 0.0, dXelem[n]);
	Xtot  = massFraction();
	Xmass = Xtot;
	
	printf("star  :\t M=%le\tR=%le\tL=%le\n", Msolar, Rsolar, Lsolar);
	sprintf(name, "WD_M%0.2f_L%0.2f_X%0.2f_Y%0.2f", Msolar, Rsolar, Xtot[0], Xtot[1]);
	// find relevant scales for re-dimensionalizing the quantities
	Dscale = Mstar*pow(Rstar,-3)/(4.*m_pi);
	Pscale = G_CGS*pow(Mstar,2)*pow(Rstar,-4)/(4.*m_pi);
	Tscale = Xtot.mean_A()*proton.mass_CGS/boltzmann_k*G_CGS*Mstar/Rstar;
	printf("scales:\t M=%le\tR=%le\tL=%le\n", Mstar, Rstar, Lstar);
	printf("scales:\t D=%le\tP=%le\tT=%le\n", Dscale, Pscale, Tscale);
	printf("Xtot  :\t %0.8lf %0.8lf %0.8lf %0.8lf\n", Xtot[0],Xtot[1],Xtot[2],Xtot[3]);
	
	Yscale = StellarVar(Dscale,Rstar,Pscale,Mstar,Tscale,Lstar);
	logYscale = log(Yscale);
	Ysolar = StellarVar(  0.0 , Rsolar, 0.0 , Msolar, 0.0 , Lsolar);
	Ystar  = StellarVar(  0.0 , Rstar,  0.0 , Mstar,  0.0 , Lstar );
			
	//begin the calculation to find values for P0, T0, and Rstar = (1+qs)*R
	printf("Converging model to P0, T0, R\n");
	double P0 = Ystart0[pres]/Pscale;//1.e24/Pscale;// //dimensionless
	double T0 = 1.e3*Teff/Tscale; //dimensionless
	double qs = 1.e-3;
	//variables of the Newtonian gradient-descent method
	double x[numv] = {P0, T0, qs};
	double z[numv] = {0., 0.,-1.e0};
	double f[numv], dx[numv], F=1.0;
	double M1, M2=0.0;
	double dFdx[numv], dfdx[numv][numv];
	double x2[numv], f2[numv], varied[numv], used[numv], F2;
	//first iteration
	int count=1;
	double maxDF = 1.0;
	while( maxDF > 1.0e-10 ){
		printf("ITERATION: %d\n", count++);
		//recompute difference, and Jacobian matrix
		joinAtCenter(x, f, F);
		printf("\tX1:\t"); for(int i=0; i<numv; i++) printf("%le ", x[i]); printf("\n");
		printf("\tf1:\t"); for(int i=0; i<numv; i++) printf("%le ", f[i]); printf("\n");
		printf("\tF1:\t%le\n", F);
		
		//calculate a Jacobian matrix by varying each variable
		for(int i=0; i<numv; i++) varied[i] = 1.01*x[i];
		for(int i=0; i<numv; i++){
			for(int j=0; j<numv; j++) used[j] = x[j];
			used[i] = varied[i];
			joinAtCenter(used, f2, F2);
			for(int j=0; j<numv; j++) dfdx[j][i] = (f2[j]-f[j])/(varied[i]-x[i]);
			dFdx[i] = (F2-F)/(varied[i]-x[i]);
		}
		
		//Newton's algorithm calls for -f
		for(int i=0; i<numv; i++) f2[i] = -f[i];
		
		//save the past gradient in case matrix inversion fails
		double dxsave[numv]; for(int i=0;i<numv;i++) dxsave[i] = dx[i];
		//invert the matrix -- check for errors
		if(invertMatrix(dfdx,f2, dx)){
			//if the matrix is singular or otherwise fails
			// then just do something to try to recover
			printf("ERROR: Matrix inversion failed!\n");
			//use the last gradient
			for(int i=0; i<numv; i++) dx[i] = dxsave[i];
		}

		bool failed = false;
		for(int i=0; i<numv; i++) if(isnan(dx[i])) failed=true;
		//if(failed) for(int i=0; i<numv; i++) dx[i] = dxsave[i];
		bool allzero = true;
		for(int i=0; i<numv; i++) allzero &= (dx[i] == 0.0);
		if(allzero) for(int i=0; i<numv; i++) {
			printf("all zero:\t%le %le\n", dx[0], dx[1]);
			dx[i] = 0.1*x[i];
		}
	
		double L = 1.0;
		for(int i=0; i<numv; i++) x2[i] = x[i] + dx[i];
		bool negative = false;
		for(int i=0; i<numv; i++) if(x2[i]<z[i]) negative=true;
		while(negative){
			for(int i=0; i<numv; i++) dx[i] *= 0.1;
			for(int i=0; i<numv; i++) x2[i] = x[i]+dx[i];
			printf("\tX2:\t"); for(int i=0; i<numv; i++) printf("%le ", x2[i]); printf("\n");
			negative = false;
			for(int i=0; i<numv; i++) if(x2[i]<z[i]) negative=true;	
		}
		
		joinAtCenter(x2, f2, F2);
		int freeze = 0;
		while(F2 > F){
			if(L<1.e-3) {L = 1.e-3; break;}
			L *= 0.1;
			for(int i=0; i<numv; i++) x2[i] = x[i] + L*dx[i];
			joinAtCenter(x2, f2, F2);
		}
		printf("\tdx:\t"); for(int i=0; i<numv; i++) printf("%le ", L*dx[i]); printf("\n");
		F = F2;
		for(int i=0; i<numv; i++) f[i] = f2[i];
		for(int i=0; i<numv; i++) x[i] = x2[i];
		
		printf("\tX2:\t"); for(int i=0; i<numv; i++) printf("%le ", x[i]); printf("\n");
		printf("\tf2:\t"); for(int i=0; i<numv; i++) printf("%le ", f[i]); printf("\n");
		printf("\tF2:\t%le\n", F);
				
		// determine the max size of the differences
		maxDF = -1.0;
		for(int i=0; i<numv; i++) if(fabs(f[i]) > maxDF) maxDF = fabs(f[i]);		
		printf("MAX DIF = %le\n", maxDF);
	}
	
	printf("MODEL CONVERGED!\n");
	joinAtCenter(x,f,F);
	
	//all done!  begin post-production
	rescaleR();
	expandGrid(Ntrue);

	populateBruntVaisala();
	Ncore = int(0.8*double(Ntot));
	if(Ncore%2==1) Ncore = Ncore - 1;
	indexFit = (Ncore)/2;
	
	printf("Teff=%le\tTsurf=%le\n", Teff, Tscale*exp(logY[Ntot-1][temp]));
	printf("Xtot :\t%0.8lf %0.8lf %0.8lf %0.8lf\n", Xtot.H1 , Xtot.He4 , Xtot.C12 , Xtot.O16);
		
	setupCenter();
	setupSurface();
}

SimpleWD::~SimpleWD(){
	delete[] logQ;
	delete[] logY;
	delete[] dlogY;
	delete[] Xelem;
	//
	delete[] dXelem;
	delete[] adiabatic_1;
	delete[] nabla;
	delete[] nabla_ad;
	delete[] brunt_vaisala;
	delete[] ledoux;
	delete[] kappa;
}


void SimpleWD::setup(){
	FILE *input_file;
	if(!(input_file=fopen("swd.txt", "r"))){
		printf("ERROR: EOS file not found... using a default\n");
		//do something
		std::vector<PartialPressure> corePres{deg_zero, rad_gas, ideal, coul};
		std::vector<PartialPressure> atmPres{rad_gas, ideal};
		for(PartialPressure p : corePres) core_pressure.push_back(p);
		for(PartialPressure p : atmPres)   atm_pressure.push_back(p);
		zy = 10.0; by = 3.0; my=1.0;
		zc = 3.0;  bc = 3.0; mc=1.0;
		zo = 2.0;  bo = 2.0; mo=0.6;
		return;
	}
	size_t buffer_size = 256;
	ssize_t line_size;
	char *input_buffer = NULL, *pressure;
	std::string instring;
	EOS *pres = NULL;
	line_size = getline(&input_buffer, &buffer_size, input_file);
	printf("%s", input_buffer);fflush(stdout);
	while(line_size > 0){
		//fscanf(input_file, "%s\n", input_buffer);
		line_size = getline(&input_buffer, &buffer_size, input_file);
		if(line_size >1) printf("%s", input_buffer);
		if(     !strcmp(input_buffer, "core:\n")) pres = &core_pressure;
		else if(!strcmp(input_buffer, "atm:\n"))  pres = &atm_pressure;
		if(pres!=NULL){
			getline(&input_buffer, &buffer_size, input_file);
			pressure = strtok(input_buffer, " \t\n");
			while(pressure != NULL){
				printf("\t%s", pressure);
				if(     !strcmp(pressure, "rad"))
					pres->push_back(rad_gas);
				else if(!strcmp(pressure, "ideal"))
					pres->push_back(ideal);
				else if(!strcmp(pressure, "coul"))
					 pres->push_back(coul);
				else if(!strcmp(pressure, "deg_zero"))
					pres->push_back(deg_zero);
				else if(!strcmp(pressure, "deg_partial"))
					pres->push_back(deg_partial);
				else if(!strcmp(pressure, "deg_finite"))
					pres->push_back(deg_finite);
				else if(!strcmp(pressure, "deg_trap"))
					pres->push_back(deg_trap);
				//else if(!strcmp(pressure, "deg_large"))
				//	pres->push_back(pressure_deg_large);
				else printf("ERROR: partial pressure term unrecognized!\n");
				pressure = strtok(NULL, " \t\n");
			}
			printf("\n");
		}
		if(!strcmp(input_buffer, "# chemical parameters\n")) {
			fscanf(input_file, "%*[^0123456789] %lf %lf %lf\n", &zy, &by, &my);
			fscanf(input_file, "%*[^0123456789] %lf %lf %lf\n", &zc, &bc, &mc);
			fscanf(input_file, "%*[^0123456789] %lf %lf %lf\n", &zo, &bo, &mo);
		}
		pres = NULL;
	}
	printf("\the\t%lf\t%lf\t%lf\n", zy, by, my);
	printf("\tc \t%lf\t%lf\t%lf\n", zc, bc, mc);
	printf("\to \t%lf\t%lf\t%lf\n", zo, bo, mo);
}

void SimpleWD::initFromChandrasekhar(){
	printf("Preparing starting values from Chandrasekhar model\n");
	Mstar = Msolar*MSOLAR;
	int Ntest = 500;
	double y0 = 1.58,  ymin = 1.0, ymax = 2.0;
	double Mtry = 0.0, Mmin = 0.0, Mmax = 2.02;
	ChandrasekharWD *testStar = new ChandrasekharWD(y0, Ntest, 2.,1.,1.,1.);
	Mtry = testStar->Mass()/MSOLAR-Msolar;
	Mmin = Mmin - Msolar;
	Mmax = Mmax - Msolar;
	
	while(fabs(ymin-ymax) > 1.0e-6){
		if(Mtry*Mmax > 0.0 ){
			ymax = y0;
			Mmax = Mtry;
		}
		else if(Mtry*Mmin > 0.0){
			ymin = y0;
			Mmin = Mtry;
		}
		y0 = 0.5*(ymin+ymax);
		delete testStar;
		testStar = new ChandrasekharWD(y0, Ntest, 2.,1.,1.,1.);
		Mtry = testStar->Mass()/MSOLAR-Msolar;
	}
	
	Rstar = testStar->Radius();
	Rsolar = Rstar/REARTH;
	
	Lstar = 4.*m_pi*boltzmann_sigma*pow(Rstar,2)*pow(Teff,4);
	Lsolar = Lstar/LSOLAR;
	
	//make initial guesses for core pressure based on this model
	// we must OVER-estimate, to avoid radiation pressure dominating
	Ystart0[dens] = 1.5*testStar->rho(0);
	Ystart0[pres] = 1.5*testStar->P(0);
	//guesses of surface -- not very good, don't use
	YstartS[dens] = testStar->rho(Ntest-2);
	YstartS[pres] = testStar->P(Ntest-2);
	
	printf("Radius in CM: %le\n", Rstar);
	printf("Radius in RE: %le\n", Rsolar);
	
	printf("Core guess: %le \t %le\n", Ystart0[pres], Ystart0[dens]);
	printf("Surf guess: %le \t %le\n", YstartS[pres], YstartS[dens]);
	
	//use the Chandrasekhar model to handle the central mode BCs
//	int CBC = 4;
//	testStar->getAstarCenter(A0, CBC, 0);
//	testStar->getVgCenter(   V0, CBC, 0); 
//	testStar->getUCenter(    U0, CBC); 
//	testStar->getC1Center(   c0, CBC); 
	return;
}

//**************************************************************************************/
//  The Simple WD Grid
//    we will use a simple core+envelope grid, similar to that used by Hansen & Kawaler
//        in their ZAMS.f program
//    the core is considered to comprise the inner 99% of the star's mass
//        within the core, we assume a constant step size
//    both the atmospher (He4) and envelope (H1) use a decreasing step size
//        the step size is modelled after a geometric series, and slowly approaches 0.9999
//**************************************************************************************/
void SimpleWD::setupGrid(double Qcore, int Ncenter){
	//in the core, just use constant step size
/*	double dq = Qcore/(Ncenter-1);
	double logQcore = log(Qcore);
		
	//now we define the core with constant step sizes
	double q = 0.0; //the center
	int n=0;
	for( ; n<Ncenter; n++){
		logQ[n] = log(q);
		q = q + dq;
	}
	
	//in the envelope, use a decreasing step size based on geometric series
	double dq1 = 1. - dq/(1.-Qcore);
	// fixed-point iteration to ensure self-consistency of dx1
	for(int a=0;a<20;a++) dq1 = 1.0 - dq/(1.-Qcore)*(1.-pow(dq1, (Ntot-Ncenter-1)/2));
	//now form the grid, starting at end of core
	for( ; n<=Ntot-1; n++){
		logQ[n] = log(q);
		if(n%2==0) dq = dq*dq1;
		q = q + dq;
	}
	
	// FOR SOME REASON, THIS HAS SLIGHTLY WRONG LIMIT
	for(n=0; n<Ntot; n++){
		logQ[n] = logQ[n]-logQ[Ntot-1];
	}
	*/
	
	int n=0;
	double q = 0.0; //the center
	
	double Ndenom = double(Ncenter*(Ncenter-1.)*(2.*Ncenter-1.));
	double dx = pow(6.*Qcore/Ndenom,1./3.);
	double dx3 = pow(dx,3);
	double dq = dx3;
	double logq = log(q);
	logQ[0] = logq;
	q = q + 0.5*dx3;
	logQ[1] = log(q);
	for(n=2 ; n<Ncenter; n++){
		dq = n*n*dx3;
		q = q + dq;
		logQ[n] = log(q);
	}
			
	logq = log(1.-q); //q = 1 - m, start at surface
	dq = (logq-log(1e-14))/(Ntot-n);
	logq = log(1e-14);
	for(n=Ntot-1;n>=Ncenter;n--){
		logQ[n] = log(1.-exp(logq));
		logq = logq + dq;
	}
	//logQ[Ntot-1] = 0.0;//*/
				
	FILE *fp = fopen("grid.txt", "w");
	for(int n=0; n<Ntot; n++){
		fprintf(fp, "%d\t%0.16le\t%0.16le\n", n, exp(logQ[n]), (1.-exp(logQ[n])));
	}
	fclose(fp);


	indexFit = Ncenter;
}

void SimpleWD::expandGrid(int Ntrue){
	//The RK4 method for finding eigenmodes requires stellar variables specified on a half-grid
	// the below will expand the WD grid to include the half-grid points, and populate them

	//first, replace each variable array with a longer one
	//	replace logQ
	double* tlogq = new double[Ntot];
	for(int x=0; x<Ntot; x++) tlogq[x] = logQ[x];
	delete[] logQ;
	logQ = new double[Ntrue];
	for(int x=0; x<Ntot; x++) logQ[2*x] = tlogq[x];
	delete[] tlogq;
	// replace logY
	StellarVar* tlogy = new StellarVar[Ntot];
	for(int x=0; x<Ntot; x++) tlogy[x] = logY[x];
	delete[] logY;
	logY = new StellarVar[Ntrue];
	for(int x=0; x<Ntot; x++) logY[2*x] = tlogy[x];
	delete[] tlogy;
	// replace dlogy
	StellarVar* tdlogy = new StellarVar[Ntot];
	for(int x=0; x<Ntot; x++) tdlogy[x] = dlogY[x];
	delete[] dlogY;
	dlogY = new StellarVar[Ntrue];
	for(int x=0; x<Ntot; x++) dlogY[2*x] = tdlogy[x];
	delete[] tdlogy;
	// replace Xelem
	Abundance* tx = new Abundance[Ntot];
	for(int x=0; x<Ntot; x++) tx[x] = Xelem[x];
	delete[] Xelem;
	Xelem = new Abundance[Ntrue];
	for(int x=0; x<Ntot; x++) Xelem[2*x] = tx[x];
	delete[] tx;
	// replace dXelem
	Abundance* tdx = new Abundance[Ntot];
	for(int x=0; x<Ntot; x++) tdx[x] = dXelem[x];
	delete[] dXelem;
	dXelem = new Abundance[Ntrue];
	for(int x=0; x<Ntot; x++) dXelem[2*x] = tdx[x];
	delete[] tdx;
	
	//now fill in the logQ, logY, dlogY arrays at the half-points
	//  for logY, we make two guesses of midpoint, integrating foward and backward, and average
	double dlogr=0., dr=0., r=0.;
	double nabla = 0.;
	StellarVar logy1, logy2;
	
	double rholast = Ystart0[dens];
	//variables to use in RK4
	StellarVar YC, logYC, Y[3];
	Y[0] = exp(logY[0]);
	Abundance chemC(0.,0.,0.,1.);
	StellarVar K[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	//begin at the center
	//setting the initial step size
	dr = 0.5*exp(logY[2][radi]);//
	//now begin integrations
	for(int x = 0; x<1; x++){
		YC = Y[x];
		chemC = findAbundance(logY[x+2][mass], logY[x+2][radi], dXelem[x+1]);
		for(int a = 0; a<4; a++){
			nabla = energyTransport(log(YC), chemC);
			K[a] = dYdR(YC, nabla)*dr;
			//at very center, account for r=0
			if(YC[radi]==0.0) K[a][pres] = K[a][temp] = 0.0;
			//calculate "corrected" positions using previous shift vectors
			YC = Y[x] + K[a]*B[a];
			chemC = findAbundance(log(YC)[mass], log(YC)[radi], dXelem[x+1]);
			YC[dens] = exp(equationOfState(log(YC), chemC, rholast));
			K[a][dens] = YC[dens];
		}
		//update the structure variables
		Y[x+1] = Y[x] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		Y[x+1][dens] = (K[0][dens]+2.*K[1][dens]+2.*K[2][dens]+K[3][dens])/6.;
		//calculate the log versions
		logY[x+1] = log(Y[x+1]);
		//logY[x+1][dens] = equationOfState(logY[x+1], Xelem[x+1], rholast);
		//Y[x+1][dens] = exp(logY[x+1][dens]);
		//fill in other properties at half-grid
		logQ[x+1] = logY[x+1][mass];
		Xelem[x+1] = findAbundance(logY[x+1][mass], logY[x+1][radi], dXelem[x+1]);
		nabla = energyTransport(logY[x+1], Xelem[x+1]);
		dlogY[x+1] = dlogYdlogR(logY[x+1], nabla);
		dlogY[x+1][dens] = (logY[x+2][dens]-logY[x][dens])/(logY[x+2][radi]-logY[x][radi]);
	}
	for(int x=3; x<Ntrue-3; x+=2){
		r = exp(logY[x-1][radi]);
		dr = (exp(logY[x+1][radi])-exp(logY[x-1][radi]));		//linear midpoint
		//integrate forward from x-1
		dlogr = log(1. + 0.5*dr/r);		//logairthmic midpoint from x-1
		nabla = energyTransport(logY[x-1], Xelem[x-1]);
		logy1 = logY[x-1] + dlogYdlogR(logY[x-1],nabla)*dlogr; 
		logy1[dens] = equationOfState(logy1, Xelem[x-1], rholast);
		//integrate backward from x+1
		r = exp(logY[x+1][radi]);
		dlogr =-log(1. - 0.5*dr/r);		//logarithmic midpoint from x+1
		nabla = energyTransport(logY[x+1], Xelem[x+1]);
		logy2 = logY[x+1] - dlogYdlogR(logY[x+1],nabla)*dlogr;
		logy2[dens] = equationOfState(logy2, Xelem[x+1], rholast);
		//average guesses
		logY[x] = (logy1+logy2)*0.5;//(logY[x-1]+logY[x+1])*0.5; //
		//fill in other properties at half-grid
		logQ[x]  = logY[x][mass];
		Xelem[x] = findAbundance(logY[x][mass], logY[x][radi], dXelem[x]);
		nabla = energyTransport(logY[x], Xelem[x]);
		dlogY[x] = dlogYdlogR(logY[x],nabla);
		dlogY[x][dens] = (logY[x+1][dens]-logY[x-1][dens])/(logY[x+1][radi]-logY[x-1][radi]);
	}
	int x=Ntrue-2;
	logY[x] = (logY[x+1]+logY[x-1])*0.5;
	logQ[x] = logY[x][mass];
	Xelem[x] = findAbundance(logY[x][mass], logY[x][radi], dXelem[x]);
	dlogY[x] = dlogYdlogR(logY[x],nabla);
	dlogY[x][dens] = (logY[x+1][dens]-logY[x-1][dens])/(logY[x+1][radi]-logY[x-1][radi]);
	Ntot = Ntrue;
}

void SimpleWD::rescaleR(){
	double logr1 = logY[Ntot-1][radi];
	for(int x=0; x<Ntot; x++){
		logY[x][radi] = logY[x][radi] - logr1;
		logY[x][dens] = logY[x][dens] + 3.*logr1;
		logY[x][pres] = logY[x][pres] + 4.*logr1;
		logY[x][temp] = logY[x][temp] + logr1;
	}
	Rstar = Rstar*exp(logr1);
	Dscale = Mstar*pow(Rstar,-3)/(4.*m_pi);
	Pscale = G_CGS*pow(Mstar,2)*pow(Rstar,-4)/(4.*m_pi);
	Tscale = Xtot.mean_A()*proton.mass_CGS/boltzmann_k*G_CGS*Mstar/Rstar;
	
	Yscale = StellarVar(Dscale,Rstar,Pscale,Mstar,Tscale,Lstar);
	logYscale = log(Yscale);
	Rsolar = Rstar/REARTH;
	Lstar = 4.*m_pi*boltzmann_sigma*pow(Rstar,2)*pow(Teff,4);
	Lsolar = Lstar/LSOLAR;
	Ysolar = StellarVar(  0.0 , Rsolar, 0.0 , Msolar, 0.0 , Lsolar);
	Ystar  = StellarVar(  0.0 , Rstar,  0.0 , Mstar,  0.0 , Lstar );
}

//**************************************************************************************/
//  The SimpleWD Core
//    the EOS will assume completely degenerate electron gas at zero temperature
//    temperature in core is finite, allowed to vary through conduction to surface
//    luminosity in core is constant
//    the opacity is taken to be some constant small value, entirely due to conduction
//    the composition is a mixture of Carbon and Oxygen
//**************************************************************************************/
double SimpleWD::calculateCore(const double x[numv], int Nmax){
	//prepare variables for RK4
	double dlogQ;
	StellarVar logYC;
	Abundance chemC(0.,0.,0.,0.);
	double nabla;
	double rholast = 10.0;
	
	//the array K contains the corrections to logY
	StellarVar K[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	
	//the first step must be handled specially
	int X = firstCoreStep(x, rholast, Nmax);
	//now begin integrations
	for(; X<Nmax; X++){
		dlogQ = logQ[X+1]-logQ[X];
		chemC = Xelem[X];
		logYC = logY[X];
		for(int a = 0; a<4; a++){
			nabla = energyTransport(logYC, chemC);
			//now from these, calculate next shift
			//use radius as independent variable
			K[a] = dlogYdlogM(logYC, nabla)*dlogQ;
			//calculate "corrected" positions using previous shift vectors
			logYC = logY[X] + K[a]*B[a];
			chemC = findAbundance(logYC[mass], logYC[radi], dXelem[X+1]);
			logYC[dens] = equationOfState(logYC, chemC, rholast);
		}
		//update the structure variables
		logY[X+1] = logY[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		//update the chemical abundances and calculate the density from EOS
		Xelem[X+1] = findAbundance(logY[X+1][mass], logY[X+1][radi], dXelem[X+1]);
		logY[X+1][dens] = equationOfState(logY[X+1], Xelem[X+1],rholast);
		//save the derivatives
		nabla = energyTransport(logY[X+1], Xelem[X+1]);
		dlogY[X+1] = dlogYdlogR(logY[X+1], nabla);
		dlogY[X][dens] = (logY[X+1][dens]-logY[X][dens])/(logY[X+1][radi]-logY[X][radi]);
	}
	return rholast;
}

//**************************************************************************************/
//  First Step of Core
//    this method performs the first step in the integration from the core
//    this must be done using non-log variables, beause log variables are singular
//    this step must be integrated in radius, because all equations in mass are singular
//**************************************************************************************/
int SimpleWD::firstCoreStep(const double x[numv], double& rholast, int Nmax){
	//prepare variables
	double P0 = (x[0]), T0 = (x[1]);

	//variables to use in RK4
	double dQ;
	StellarVar YC, logYC;
	Abundance chemC(0.,0.,0.,1.);
	double nabla;
	//the arrays K and L contain the corrections to y and z to be applied
	// when generating the next step of the correction
	StellarVar K[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};

	//set initial conditions
	//printf("Core dens: %le\t%le\t%le\n", Ystart0[dens], sqrt(P0*Pscale), 3.9308846494269581E+06);
	//printf("Core pres: %le\t%le\t%le\n", x[0]*Pscale,  Ystart0[pres], 2.0782717264828863E+23);
	rholast = Ystart0[dens]; //sqrt(P0*Pscale); //an initial guess
	StellarVar Ycenter(rholast/Dscale, 0.0, P0, 0.0, T0, 0.0);
	StellarVar logYcenter = log(Ycenter);
	Xelem[0] = findAbundance(logQ[1], logYcenter[1], dXelem[0]);
	logYcenter[dens] = equationOfState(logYcenter, Xelem[0], rholast);
	Ycenter[dens] = exp(logYcenter[dens]);
	Ystart0[dens] = rholast;
	
	//begin RK4 at center
	const int start = 1;	//number of points to compute in this way
	StellarVar Y[start+1];		//the array of points computed this way
	//begin at the center
	Y[0]     = Ycenter;			//the central values
	logY[0]  = logYcenter;		//the central log-values
	dlogY[0] = StellarVar();	//the central log-derivative is undefined
	//setting the initial step size
	dQ = exp((log(3.) + logQ[1] - logY[0][dens])/3.);//
	//now begin integrations
	for(int X = 0; X<start; X++){
		YC = Y[X];
		chemC = findAbundance(logQ[1], logY[1][radi], dXelem[X+1]);
		for(int a = 0; a<4; a++){
			nabla = energyTransport(log(YC), chemC);
			K[a] = dYdR(YC, nabla)*dQ;
			//at very center, account for r=0
			if(YC[radi]==0.0) K[a][pres] = K[a][temp] = 0.0;
			//calculate "corrected" positions using previous shift vectors
			YC = Y[X] + K[a]*B[a];
			chemC = findAbundance(log(YC)[mass], log(YC)[radi], dXelem[X+1]);
			YC[dens] = exp(equationOfState(log(YC), chemC, rholast));
			K[a][dens] = YC[dens];
		}
		//update the structure variables
		Y[X+1] = Y[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		Y[X+1][dens] = (K[0][dens]+2.*K[1][dens]+2.*K[2][dens]+K[3][dens])/6.;
		//calculate the log versions
		logY[X+1] = log(Y[X+1]);
		logY[X+1][mass] = logQ[X+1];
		Xelem[X+1] = findAbundance(logY[X+1][mass], logY[X+1][radi], dXelem[X+1]);
		//logY[X+1][dens] = equationOfState(logY[X+1], Xelem[X+1], rholast);
		//Y[X+1][dens] = exp(logY[X+1][dens]);
		//calculte derivatives at this location
		nabla = energyTransport(logY[X+1], Xelem[X+1]);
		dlogY[X+1] = dlogYdlogR(logY[X+1], nabla);
		dlogY[X+1][dens] = (logY[X+1][dens]-logY[X][dens])/(logY[X+1][radi]-logY[X][radi]);
	}
	return start;
}

//**************************************************************************************/
//  The SimpleWD Atmosphere
//    the EOS incudes electrons, ions, Coulomb corrections, partial ionization, and radiation
//    either radiative or convective heat transport allowed to take place
//    the composition is a mixture of Helium, Carbon, and Hydrogen
//	Must integrate from surface (photosphere) to the core
//**************************************************************************************/
void SimpleWD::calculateAtmosphere(const double x[numv]){	
	//prepare variables
	double dlogQ;
	StellarVar logYC;
	Abundance chemC(0.,0.,0.,0.);
	double nabla;
	double rholast = 1.0;
	
	//the array K contains the corrections to logY to be applied
	//    when generating the next step of the correction
	StellarVar K[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	
	//set our initial conditions
	int X = firstAtmosphereStep(x, rholast);
	for(; X > Ncore ; X--){
		dlogQ = logQ[X-1] - logQ[X];
		chemC = Xelem[X];
		logYC = logY[X];		
		for(int a = 0; a<4; a++){
			//find energy production, heat transport
			nabla = energyTransport(logYC, chemC);
			//now from these, calculate next shift
			//use radius as independent variable
			K[a] = dlogYdlogM(logYC, nabla)*dlogQ;
			//calculate "corrected" positions using previous shift vectors
			logYC = logY[X] + K[a]*B[a];
			chemC = findAbundance(logYC[mass], logYC[radi], dXelem[X-1]);
			logYC[dens] = equationOfState(logYC, chemC, rholast);
		}
		//update the chemical composition and the density
		logY[X-1] = logY[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		//update the chemical abundances and calculate the density from EOS
		Xelem[X-1] = findAbundance(logY[X-1][mass], logY[X-1][radi], dXelem[X-1]);
		logY[X-1][dens] = equationOfState(logY[X-1], Xelem[X-1], rholast);
		//calculate the derivatives at this location
		nabla = energyTransport(logY[X-1], Xelem[X-1]);
		dlogY[X-1] = dlogYdlogR(logY[X-1], nabla);
		dlogY[X-1][dens] = (logY[X-1][dens]-logY[X][dens])/(logY[X-1][radi]-logY[X][radi]);
	}
	return;
}

//**************************************************************************************/
//  First Step of Atmosphere
//    this method performs the first step in the integration from the surface
//    this must be done using non-log variables, beause log variables are singular
//**************************************************************************************/
int SimpleWD::firstAtmosphereStep(const double x[numv], double& rholast){
	//prepare variables
	double rs = (1.0 + x[2]);
	double ms = (1.0);
	double LS = (1.0);//pow(1.0+x[2],2);
	
	//values of true surface, photosphere, and first inward grid point
	StellarVar Ysurf, Yphoto, Yfirst;
	StellarVar lys, lyp, lyf;
	
	//chemical composition in the outer atmosphere
	Xelem[Ntot-1] = findAbundance(0.0, 0.0, dXelem[Ntot-1]);
	Abundance Xsurf = Xelem[Ntot-1];	
	
	double taup = 2./3.; // optical depth at photosphere, in Eddington approximation
	double taus = 1.e-5; // optical depth at surface - a small number
	double tauf, mf;     // optical dpeth at first grid point
	
	//begin by finding the photosphere, using Eddington approximation
	// this sets the value of opacity used in gray atmosphere
	//initialize a guess for photosphere, then iterate to perfect
	Yphoto[dens] = 1e-6;     //a first guess
	Yphoto[temp] = Teff;	 //photosphere defined by effective temperature
	Yphoto[pres] = atm_pressure(Yphoto[dens], Yphoto[temp], Xsurf);
	Yphoto[radi] = Rstar*rs; //surface value
	Yphoto[mass] = Mstar*ms; //surface value
	Yphoto[lumi] = Lstar*LS; //surface value
	lyp = log(Yphoto); //dimensional log version	
	double gs = G_CGS*exp(lyp[mass]-2.*lyp[radi]); //the surface gravity
	double a = gs*taup;
	double b = radiation_a*exp(4.*lyp[temp])/6.;
	double kp = radiative_opacity(lyp, Xsurf);
	lyp[pres] = log(a/kp+b);	// refine guess for PS
	//prepare a Newton method to find photosphere
	double f1 = atm_pressure(Yphoto[dens], Yphoto[temp], Xsurf)-a/kp-b, f2=1.0e2;
	double x1=1e2, x2 = Yphoto[dens], dx;
	while(fabs(x1-x2)/x1 > 1e-10){
		x1 = x2;
		x2 = 1.01*x1;
	//  perturb
		lyp[dens] = log(x2);
		kp = radiative_opacity(lyp,Xsurf);
		f2 = atm_pressure(x2,Yphoto[temp],Xsurf) - a/kp-b;
	//  correct
		dx = -f1*(x2-x1)/(f2-f1);
		while(x1+dx < 0.0){
			x2 *= 0.1;
			lyp[dens] = log(x2);
			kp = radiative_opacity(lyp,Xsurf);
			f2 = atm_pressure(x2,Yphoto[temp],Xsurf) - a/kp-b;
			dx = -f1*(x2-x1)/(f2-f1);
		}
		x2 = x1 + dx;
	//  find new
		lyp[dens] = log(x2);
		kp = radiative_opacity(lyp,Xsurf);
		f1 = atm_pressure(x2,Yphoto[temp],Xsurf) - a/kp-b;
	}	
	Yphoto[dens] = x2;
	Yphoto[pres] = a/kp+b;
	lyp = log(Yphoto);
	kp = radiative_opacity(lyp, Xsurf);
//	printf("first guess photosphere : %le %le\n", Yphoto[dens], Yphoto[pres]);
		
	//now find values of P, rho, T at the true surface, tau=0.0
	Ysurf[temp] = pow(3./4.*(taus + 2./3.),0.25)*Teff; // surface temperature from Eddington
	Ysurf[pres] = gs*taus/kp + b;
	Ysurf[dens] = atm_pressure.invert(Yphoto[dens], Ysurf[pres], Ysurf[temp], Xsurf);
	Ysurf[radi] = Rstar*rs;
	Ysurf[mass] = Mstar*ms;
	Ysurf[lumi] = Lstar*LS;
	lys = log(Ysurf);
	double ks = radiative_opacity(lys, Xsurf);
//	printf("first guess surface     : %le %le\n", Ysurf[dens], Ysurf[pres]);
		
	mf = exp(logQ[Ntot-1]);
	tauf = kp/(4.*m_pi)*exp(lys[mass]-2.*lyp[radi])*(1.-mf);
//	printf("tauf=%le\t", tauf);
//	printf("kp =%le\tks = %le\n", kp, ks);
	StellarVar Y1, Y2, dYdt, ly1;
	double tau1, tau2, dtau;
	double g, kappa, nabla;
	tau1 = taup;
	dtau = (tauf - tau1)/100.0;
	Y1 = Yphoto;
//	printf("integrating first grid:\t%le %le %le %le %le %le\n", tau1, Y1[dens], Y1[pres], Y1[temp], Y1[radi], 1.-mf);
	//while(tau1 < tauf){
	while(Y1[mass] > mf*Mstar){
		ly1 = log(Y1);
		g = G_CGS*exp(ly1[mass]-2.*ly1[radi]);
		kappa = radiative_opacity(ly1, Xsurf);
		nabla = atm_pressure.nabla_ad(Y1[dens],Y1[temp],Xsurf);
	//	printf("g=%le\tk=%le\tD=%le\n", g, kappa, nabla);
		//
		tau2 = tau1 + dtau;
		dYdt[pres] = g/kappa;
		dYdt[temp] = g/kappa*exp(ly1[temp]-ly1[pres])*nabla;
		dYdt[mass] = -4.*m_pi*exp(2.*ly1[radi])/kappa;
		dYdt[radi] = -1./(Y1[dens]*kappa);
		//
		Y2 = Y1 + dYdt*dtau;
		Y2[temp] = pow(3./4.*(tau2 + 2./3.),0.25)*Teff;
		Y2[dens] = atm_pressure.invert(Y1[dens], Y2[pres],Y2[temp],Xsurf);
		//
		Y1 = Y2;
		tau1 = tau2;		
	//	printf("integrating first grid:\t%le %le %le %le %le %le\n", tau1, Y1[dens], Y1[pres], Y1[temp], Y1[radi], 1.-Y1[mass]/Mstar);
	}
//	printf("integrating first grid:\t%le %le %le %le %le %le\n", tau1, Y1[dens], Y1[pres], Y1[temp], Y1[radi], 1.-Y1[mass]/Mstar);
	Yfirst = Y1;
	Yfirst[mass] = Mstar*ms*mf;
	Yfirst[lumi] = Lstar*LS*mf;
	//printf("integrating first grid:\t%le %le %le %le %le %le\n", tau1, Yfirst[dens], Yfirst[pres], Yfirst[temp], Yfirst[radi], 1.-Yfirst[mass]/Mstar);
	
	Ysurf = StellarVar(Yfirst[dens],Rstar*rs,Yfirst[pres],Mstar*ms,Ysurf[temp],Lstar*LS);
	
	//int start = Ntot-1, first = Ntot-1;
	//StellarVar Y[(Ntot-start)];
	
	
	const int npoints = 2, start = Ntot-npoints+1, first = Ntot-1;
	StellarVar Y[npoints];
	Y[0] = Ysurf;
	Y[1] = Yfirst;
	rholast = Y[1][dens];
	logY[Ntot-1] = log(Y[0]) - logYscale;
	logY[Ntot-2] = log(Y[1]) - logYscale;
	Y[0] = exp(logY[Ntot-1]);
	Y[1] = exp(logY[Ntot-2]);
	//calculate surface derivatives
	nabla = atm_pressure.nabla_ad( Y[1][dens], Y[1][temp], Xsurf);
	dlogY[Ntot-1] = dlogYdlogR(logY[Ntot-1], nabla); 
	dlogY[Ntot-2] = dlogYdlogR(logY[Ntot-2], nabla);
	Xelem[Ntot-1] = findAbundance(logY[Ntot-1][mass], logY[Ntot-1][radi], dXelem[Ntot-1]);
	Xelem[Ntot-2] = findAbundance(logY[Ntot-2][mass], logY[Ntot-2][radi], dXelem[Ntot-2]);

	//now integrate whatever remains
	double dQ;
	StellarVar YC, logYC;
	Abundance chemC;
	//the array K contains the corrections to Y to be applied
	//    when generating the next step of the correction
	StellarVar K[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	for(int X = first; X > start; X--){
		YC = Y[Ntot-1-X];
		chemC = Xelem[X];
		dQ = exp(logQ[X-1]) - exp(logQ[X]);
		for(int a = 0; a<4; a++){
			nabla = energyTransport(log(YC), chemC);
			//now from these, calculate next shift
			K[a] = dYdM(YC, nabla)*dQ;
			//calculate "corrected" positions using previous shift vectors
			YC = Y[Ntot-1-X] + K[a]*B[a];
			chemC = findAbundance(log(YC)[mass], log(YC)[radi], dXelem[X-1]);
			YC[dens] = exp(equationOfState(log(YC), chemC, rholast));
		}
		//update the chemical composition and the density
		Y[Ntot-X] = Y[Ntot-X-1] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		//calculate the log versions
		logY[X-1] = log(Y[Ntot-X]);
		Xelem[X-1] = findAbundance(logY[X-1][mass], logY[X-1][radi], dXelem[X-1]);
		logY[X-1][dens] = equationOfState(logY[X-1], Xelem[X-1], rholast);
		Y[Ntot-X][dens] = exp(logY[X-1][dens]);
		//calculate the derivatives at this location
		nabla = energyTransport(logY[X-1], Xelem[X-1]);
		dlogY[X-1] = dlogYdlogR(logY[X-1], nabla);
		dlogY[X-1][dens] = (logY[X-1][dens]-logY[X][dens])/(logY[X-1][radi]-logY[X][radi]);
	}
	return start;//*/

/*//this will create "correct" convection zone for Ntot=1000... save in case needed

	//prepare variables
	double rs = (1.0 + x[2]);
	double ms = (1.0);
	double LS = (1.0);
	
	//set our initial conditions
	//set the surface values of logY
	Xelem[Ntot-1] = findAbundance(0.0, 0.0, dXelem[Ntot-1]);
	Abundance Xsurf = Xelem[Ntot-1];
		
	//begin by finding values of P, rho, T at the true surface, tau=0.0
	rholast = 1e-10;    // a first guess
	double To = pow(0.5,0.25)*Teff; //surface temperature from Eddington equation
	double PS = atm_pressure(rholast, To, Xsurf);	// a first guess
	StellarVar Ysurf(rholast/Dscale, rs,  PS, ms, To/Tscale, LS);
	StellarVar ly = log(Ysurf)+logYscale; //dimensional log version
	double tau = 1.e-5;//
	double gs = G_CGS*Mstar*pow(Rstar*rs,-2); //the surface gravity
	double A = Xsurf.mean_A()/(N_Avogadro*boltzmann_k*To);
	double a = gs*tau;
	double b = radiation_a*pow(Teff,4)/6.;
	double kp = radiative_opacity(ly, Xsurf);
	ly[pres] = log(a/kp+b);	// refine guess for PS
	//prepare a Newton method to find surface value of rho
	double f1 = atm_pressure(rholast, To, Xsurf)-a/kp-b, f2=1.0e2;
	double x1=1e2, x2 = rholast;//A*(a/kp - b);
	printf("rho guess: %le %le %le\n", x1, x2, f1);
	while(fabs(x1-x2)/x1 > 1e-10){
		x1 = x2;
		x2 = 1.01*x1;
	//  perturb
		ly[dens] = log(x2);
		kp = radiative_opacity(ly,Xsurf);
		f2 = atm_pressure(x2,To,Xsurf) - a/kp-b;
	//  correct
		x2 = x1 - f1*(x2-x1)/(f2-f1);
	//  find new
		ly[dens] = log(x2);
		kp = radiative_opacity(ly,Xsurf);
		f1 = atm_pressure(x2,To,Xsurf) - a/kp-b;
	//
		printf("rho guess: %le %le %le\n", x1, x2, f1);
	}	
	rholast = x2;
	PS = a/kp+b;
	printf("P   guess: %le\n", PS);
	
	if(isnan(PS)) return 0;
	
	double mp = exp(logQ[Ntot-2]);
	double taup = kp*Mstar/4./m_pi*pow(Rstar*rs,-2)*(1.-mp);
	printf("taup = %le\n", taup);
	
	//now integrate from the true surface to the photosphere
	double T1, P1, rho1, T2, P2, rho2;
	double dtau = 1e-2, tau2, tau1;
	//double taup = 20.;//2./3.; //the photospheric value, where integration stops
	tau1 = 0.0;
	P1 = gs*tau1/kp + b;
	T1 = pow(3./4.*(tau1 + 2./3.),0.25)*Teff;
	rho1 = rholast;
	printf("integrating photosphere:\t%le %le %le %le\n", tau1, rho1, P1, T1);
	while(tau2 < taup){
		ly[dens] = log(rho1);
		ly[temp] = log(T1);
		ly[pres] = log(P1);
	//	kp = radiative_opacity(ly, Xsurf);
		tau2 = tau1 + dtau;
		P2 = P1 + gs/kp*dtau;
		T2 = pow(3./4.*(tau1 + 2./3.),0.25)*Teff;
		rho2 = atm_pressure.invert(rho1, P2,T2,Xsurf);
		tau1 = tau2;
		P1 = P2;
		rho1 = rho2;
		T1 = T2;
	}
	printf("integrating photosphere:\t%le %le %le %le\n", tau2, rho2, P2, T2);
	PS = P1;
	rholast = rho1;
	
	//oldie but goodie
	//PS = (1.e9);			//an initial guess
	//rholast = sqrt(PS);   //an initial guess
				
	//set surface variables
	Ysurf = StellarVar( rholast/Dscale, rs, PS/Pscale, ms, Teff/Tscale, LS );
	
	logY[Ntot-1] = log(Ysurf);
	logY[Ntot-1][dens] = equationOfState(logY[Ntot-1], Xsurf, rholast);
	Ysurf[dens] = exp(logY[Ntot-1][dens]);
	//calculate surface derivative
	double nabla = energyTransport(log(Ysurf), Xsurf);
	dlogY[Ntot-1] = dlogYdlogR(logY[Ntot-1], nabla);

	//set up RK4 integration
	double dQ;
	StellarVar YC, logYC;
	Abundance chemC;
	//the array K contains the corrections to logY to be applied
	//    when generating the next step of the correction
	StellarVar K[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};

	int start = Ntot-2;
	StellarVar Y[(Ntot-start)];
	Y[0] = Ysurf;
	for(int X = Ntot-1; X > start; X--){
		YC = Y[Ntot-1-X];
		chemC = Xelem[X];
		dQ = exp(logQ[X-1]) - exp(logQ[X]);
		for(int a = 0; a<4; a++){
			nabla = energyTransport(log(YC), chemC);
			//now from these, calculate next shift
			K[a] = dYdM(YC, nabla)*dQ;
			//calculate "corrected" positions using previous shift vectors
			YC = Y[Ntot-1-X] + K[a]*B[a];
			chemC = findAbundance(log(YC)[mass], log(YC)[radi], dXelem[X-1]);
			YC[dens] = exp(equationOfState(log(YC), chemC, rholast));
			if(isnan(YC[pres])) printf("%d %d NAN pres\n", X, a);
			if(isnan(YC[dens])) printf("%d %d NAN dens\n", X, a);
			if(isnan(YC[mass])) printf("%d %d NAN mass\n", X, a);
			if(isnan(YC[temp])) printf("%d %d NAN temp\n", X, a);
			if(isnan(YC[lumi])) printf("%d %d NAN lumi\n", X, a);
		}
		//update the chemical composition and the density
		Y[Ntot-X] = Y[Ntot-X-1] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		//calculate the log versions
		logY[X-1] = log(Y[Ntot-X]);
		Xelem[X-1] = findAbundance(logY[X-1][mass], logY[X-1][radi], dXelem[X-1]);
		logY[X-1][dens] = equationOfState(logY[X-1], Xelem[X-1], rholast);
		Y[Ntot-X][dens] = exp(logY[X-1][dens]);
		//calculate the derivatives at this location
		nabla = energyTransport(logY[X-1], Xelem[X-1]);
		dlogY[X-1] = dlogYdlogR(logY[X-1], nabla);
		dlogY[X-1][dens] = (logY[X-1][dens]-logY[X][dens])/(logY[X-1][radi]-logY[X][radi]);
	}
	return start;//*/
}


//**************************************************************************************/
//  JOIN AT CENTER -- Double-Shooting Algorithm 
//    Here we begin integration from core and from surface, joining at central location
//    Will return the difference in two integrations at the core
//    Will also return Jacobian matrix, how differences change with dependent parameters
//**************************************************************************************/
void SimpleWD::joinAtCenter(double x[numv], double f[numv], double& F){
	StellarVar Y1, Y2;
	//calculate first using the current best guesses
	calculateCore(x,Ncore);
	Y1 = exp(logY[Ncore]);
	calculateAtmosphere(x);
	Y2 = exp(logY[Ncore]);
	Abundance Xm1 = massFraction();
	
	f[0] = (Y1[pres]-Y2[pres]);
	f[1] = (Y1[temp]-Y2[temp]);
	f[2] = (Y1[radi]-Y2[radi]);
	
	F = 0.0;
	for(int i=0; i<numv; i++) F += 0.5*( f[i]*f[i] );

	return;
}

// differential equations
StellarVar SimpleWD::dYdR(const StellarVar& y, const double& nabla, const double epsilon){
	StellarVar deriv;
	deriv[dens] = 0.0;
	deriv[pres] =-y[mass]*y[dens]*pow(y[radi],-2);
	deriv[temp] =-y[mass]*y[dens]*pow(y[radi],-2)*y[temp]/y[pres]*nabla;
	deriv[mass] = y[dens]*pow(y[radi],2);
	deriv[lumi] = y[dens]*pow(y[radi],2)*epsilon;
	deriv[radi] = 1.0;
	return deriv;
}

// logarithmic derivatives of variables
StellarVar SimpleWD::dlogYdlogR(  const StellarVar& logy, const double& nabla, const double epsilon){
	StellarVar deriv;
	deriv[dens] =  0.0;
	deriv[pres] = -exp(logy[mass] + logy[dens] - logy[pres] - logy[radi]);
	deriv[temp] = -exp(logy[mass] + logy[dens] - logy[pres] - logy[radi])*nabla;
	deriv[mass] =  exp(logy[dens] - logy[mass]          + 3.0*logy[radi]);
	deriv[lumi] =  exp(logy[dens] - logy[lumi]          + 3.0*logy[radi])*epsilon;
	deriv[radi] =  exp(0.0);
	return deriv;
}

// differential equations using mass
StellarVar SimpleWD::dYdM(const StellarVar& y, const double& nabla, const double epsilon){
	StellarVar deriv;
	deriv[dens] = 0.0;
	deriv[pres] =-y[mass]*pow(y[radi],-4);
	deriv[temp] =-y[mass]*pow(y[radi],-4)*y[temp]/y[pres]*nabla;
	deriv[radi] = pow(y[radi],-2)/y[dens];
	deriv[lumi] = epsilon;
	deriv[mass] = 1.0;
	return deriv;
}

// logarithmic derivatives of variables
StellarVar SimpleWD::dlogYdlogM(  const StellarVar& logy, const double& nabla, const double epsilon){
	StellarVar deriv;
	deriv[dens] =  0.0;
	deriv[pres] = -exp(2.*logy[mass] - logy[pres] - 4.*logy[radi]);
	deriv[temp] = -exp(2.*logy[mass] - logy[pres] - 4.*logy[radi])*nabla;
	deriv[radi] =  exp(   logy[mass] - logy[dens] - 3.*logy[radi]);
	deriv[lumi] =  exp(   logy[mass] - logy[lumi])*epsilon;
	deriv[mass] =  exp(0.0);
	return deriv;
}


//**************************************************************************************/
//  Constitutive Physics
//**************************************************************************************/

// EPSILON  the source of luminosity
double SimpleWD::energyProduction(const StellarVar& logy, const Abundance& X){
	//in our setting, dm/dr = rho*r^2, with BCs 0 and 1
	//but also,       dL/dr = rho*r^2*eps, with BCs 0 and 1 and eps constant
	//the only way for both of these to be satisfied is if eps = 1.0
	//if a non-constant eps is needed, obviously we will need to modify this
	return 1.0;
}

// NABLA    the heat transportation, calculting from radiative/conductive transport
double SimpleWD::energyTransport( const StellarVar& logy, const Abundance& X)
{
	//rescale variables to CGS
	StellarVar ly = logy + logYscale;
	double kappa = opacity(ly, X);
	double D_coeff = 3./(16.*m_pi*radiation_a*C_CGS*G_CGS);
	double delta = D_coeff*kappa*exp(ly[lumi]+ly[pres]-ly[mass]-4.*ly[temp]);
	if(isnan(delta)) delta = D_coeff*kappa*exp(ly[pres]-4.*ly[temp])*Lstar/Mstar;
	
	double delta_ad = getEOS(exp(ly), X)->nabla_ad(exp(ly[dens]), exp(ly[temp]), X);
		
	if(delta_ad < delta) {
		return delta_ad;
	}
	else if(delta < delta_ad) {
		return delta;
	}
	else return delta;
}

// EQUATION OF STATE  relates pressure and density
EOS* SimpleWD::getEOS(const StellarVar &y, const Abundance& X){
	double Patm  =  atm_pressure(y[dens], y[temp], X);
	double Pcore = core_pressure(y[dens], y[temp], X);
	if(X.H1 < 0.5) return &core_pressure;
	else if (Pcore <  Patm) return &atm_pressure;
	else if (Pcore >= Patm) return &core_pressure;
	else return &core_pressure; //*/
}

double SimpleWD::equationOfState(const StellarVar& logy, const Abundance& chem, double& rho_last){
	StellarVar y = exp(logy+logYscale);
	//find the radiation pressure
	double Prad, rho=0.0;
	Prad = radiation_a/3.*pow(y[temp],4);
	if(Prad > y[pres]) {
		printf("\nERROR: RADIATION PRESSURE TOO HIGH!\n");
	//	y[temp] = 0.0;
	}
	
	rho = getEOS(y, chem)->invert(rho_last, y[pres], y[temp], chem);
	
	if(rho<0.0) rho = -rho;
	if(isnan(rho)) printf("RHO IS NAN!\n");
	if(isnan(log(rho))){
		printf("%le %le %le %le %le %le\n", exp(logy[radi]), y[pres], y[temp], rho, rho/Dscale, log(rho/Dscale));
	}
	
	//save the density	
	rho_last = rho;
	//return the log density
	return log(rho/Dscale);
}

//**************************************************************************************/
//  Chemical Composition
//**************************************************************************************/
Abundance SimpleWD::findAbundance(
	const double logm,
	const double logr,
	Abundance& dX
)
{
	double tol = 0.001; //the limit for element to be conidered to have "vanished"

	Abundance XX(0.0,0.0,0.0,0.0);
	Abundance Xc(0.0,0.0,0.0,0.0);
	dX = Abundance(0.,0.,0.,0.);
		
	double zee = -log(1.-exp(logm));	
	Xc.He4 = my/(1.+exp(by*(zee-zy)));
	Xc.C12 = mc/(1.+exp(bc*(zee-zc)));
	Xc.O16 = mo/(1.+exp(bo*(zee-zo)));//*/
	
	//any values are are negligibly small, set to zero
	for(int i=0; i<4; i++) if(Xc[i]<=tol) Xc[i]=0.0;
	
	//adjust values to ensure all abundances always sum to zero
	XX.O16 = Xc.O16;
	XX.C12 = Xc.C12 - Xc.O16;
	XX.He4 = Xc.He4 - Xc.C12;
	XX.H1  = 1.0 - XX.He4 - XX.C12 - XX.O16;
	
	//calculate the derivatives -- with respect to zee
	dX.O16 = -bo*Xc.O16*exp(bo*(zee-zo))/(1.+exp(bo*(zee-zo)))*exp(zee);
	dX.C12 = -bc*Xc.C12*exp(bc*(zee-zc))/(1.+exp(bc*(zee-zc)))*exp(zee);
	dX.He4 = -by*Xc.He4*exp(by*(zee-zy))/(1.+exp(by*(zee-zy)))*exp(zee);//*/
	// change derivatives -- with respect to P, so that dX holds dlogX/dlogP
	dX = dX*(-exp(4.*logr-logm)/Pscale);
	
	//adjust derivatives ensuring abundances always sum to zero
	dX.H1  =-dX.He4 - dX.C12 - dX.O16;
	dX.He4 = dX.He4 - dX.C12;
	dX.C12 = dX.C12 - dX.O16;
	
	XX.enforce();
	XX.e = chemical::h;
//	XX.e = (XX.He4 > 0.5 | XX.H1 >= 0.5 ? chemical::h : chemical::c);
		
	return XX;
}

Abundance SimpleWD::massFraction(){
	Abundance massTot(0.,0.,0.,0.);
	Abundance Xint1=Xelem[0], Xint2=Xelem[1];
	
	//first integrate the core
	massTot = massTot + (Xint1+Xint2)*(exp(logQ[1]))*0.5;
	Xint2 = Xelem[1]*exp(logQ[1]);
	for(int x=1; x<Ntot-1; x++){
		if(exp(logQ[x+1]) >1.0) break;
		Xint1 = Xint2;
		Xint2 = Xelem[x+1]*exp(logQ[x+1]);
		massTot = massTot + (Xint1+Xint2)*(logQ[x+1]-logQ[x])*0.5;
	}
	massTot.enforce();
	return massTot;
}

//**************************************************************************************/
//  ACCESSORS
//    Here we define functions to access radius, pressure, etc.
//**************************************************************************************/
double SimpleWD::rad(int X){
	return Rstar*exp(logY[X][radi]);
}
double SimpleWD::rho(int X){
	return Dscale*exp(logY[X][dens]);
}
double SimpleWD::drhodr(int X){
	 // dlogD/dlogr = dD/dr*r/D --> dD/dr = D/r * dlogD/dlogr
	return Dscale/Rstar*exp(logY[X][dens]-logY[X][radi])*dlogY[X][dens];
}
double SimpleWD::P(int X){
	return Pscale*exp(logY[X][pres]);
}
double SimpleWD::dPdr(int X){
	return Pscale/Rstar*exp(logY[X][pres]-logY[X][radi])*dlogY[X][pres];
}
double SimpleWD::Phi(int X){
	//zeroed to join exterior solution at surface, where Phi->0 at infty
	return 0.0;
}
double SimpleWD::dPhidr(int X){
	return G_CGS*Mstar*pow(Rstar,-2)*exp(logY[X][mass] - 2.*logY[X][radi]);
}
double SimpleWD::mr(int X){
	return Mstar*exp(logY[X][mass]);
}

double SimpleWD::Schwarzschild_A(int X, double GamPert){
	if(GamPert==0.0) return -brunt_vaisala[X]*exp(2.*logY[X][radi]-logY[X][mass])/Rstar;
	else        	 return (dlogY[X][dens] - dlogY[X][pres]/GamPert)*exp(-logY[X][radi])/Rstar;
}

double SimpleWD::getAstar(int X, double GamPert){
	if(GamPert==0.0) return brunt_vaisala[X]*exp(3.*logY[X][radi]-logY[X][mass]);
	else        	 return (dlogY[X][pres]/GamPert - dlogY[X][dens]);
}

//the Ledoux part of the Schwarzschild discriminant
double SimpleWD::Ledoux(int X, double GamPert){
	if(GamPert==0.0) return ledoux[X];
	else             return 0.0;
}

double SimpleWD::getU(int X){
	if(X==0) return 3.0;
	else     return dlogY[X][mass];	//this is already dlog(m)/dlog(r)
}

double SimpleWD::getVg(int X, double GamPert){
	if(GamPert==0.0) return -dlogY[X][pres]/Gamma1(X);//we already have dlog(P)/dlog(r)
	else			 return -dlogY[X][pres]/GamPert;
}

double SimpleWD::getC(int X){
	if(X==0) return 3.*exp(-logY[0][dens]);
	else     return exp(3.*logY[X][radi] - logY[X][mass]);
}

double SimpleWD::Gamma1(int X){
	return adiabatic_1[X];
}

double SimpleWD::sound_speed2(int X, double GamPert){
	if(GamPert == 0.0) return Gamma1(X) * exp(logY[X][pres]-logY[X][dens])*Pscale/Dscale;
	else               return GamPert   * exp(logY[X][pres]-logY[X][dens])*Pscale/Dscale;
}

void SimpleWD::populateBruntVaisala(){
	adiabatic_1 = new double[Ntot];
	nabla = new double[Ntot];
	nabla_ad = new double[Ntot];
	brunt_vaisala = new double[Ntot];
	ledoux = new double[Ntot];
	kappa = new double[Ntot];
	
	StellarVar YY = exp(logY[0]+logYscale);
	StellarVar ly;
	double g, chiT, chir, nbad;
	double P, rho, Pcore, Patm;
	double kappa_rad, kappa_cond, kappa0, D_coeff;
	EOS* myEOS;
	D_coeff = 3./(16.*m_pi*radiation_a*C_CGS*G_CGS);
	for(int X=0; X<Ntot;X++){
		ly = logY[X] + logYscale;
		YY = exp(ly);
			
		kappa[X] = opacity(ly, Xelem[X]);
		
		nabla[X] = D_coeff*kappa[X]*exp(ly[pres]+ly[lumi]-ly[mass]-4.*ly[temp]);
		if(isnan(nabla[X])) nabla[X] = D_coeff*kappa[0]*exp(ly[pres] - 4.*ly[temp])*Lstar/Mstar;
		myEOS = getEOS(YY, Xelem[X]);

		adiabatic_1[X] = myEOS->Gamma1(  YY[dens],YY[temp],Xelem[X]);
		nabla_ad[X]    = myEOS->nabla_ad(YY[dens],YY[temp],Xelem[X]);
		ledoux[X]      = myEOS->Ledoux(  YY[dens],YY[temp],Xelem[X],dXelem[X]);
		chiT           = myEOS->chiT(    YY[dens],YY[temp],Xelem[X]);
		chir           = myEOS->chiRho(  YY[dens],YY[temp],Xelem[X]);
		if(nabla_ad[X] < nabla[X]) nabla[X] = nabla_ad[X];
		if(isnan(ledoux[X])) ledoux[X]=0.0;
		brunt_vaisala[X] = exp(2.*logY[X][mass]-4.*logY[X][radi]+logY[X][dens]-logY[X][pres])
				*chiT/chir*(nabla_ad[X]-nabla[X]+ledoux[X]);
	}
}

double SimpleWD::Radius(){return Rstar;}	//total radius
double SimpleWD::Mass(){  return Mstar;}	//total mass
double SimpleWD::Gee(){   return G_CGS;}
//in Newtonian, light speed is infinity...
double SimpleWD::light_speed2(){return C_CGS*C_CGS;}

// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R
void SimpleWD::setupCenter(){
	nc = 1./(Gamma1(1)-1.);
	double rho0 = Dscale*exp(logY[0][dens]);
	double P0   = Pscale*exp(logY[0][pres]);
	double T0   = Tscale*exp(logY[0][temp]);
	double dPdd = core_pressure.partialRho(rho0,T0,Xelem[0])*Dscale/Pscale;
	double dPdt = core_pressure.partialT(  rho0,T0,Xelem[0])*Tscale/Pscale;
	//second derivative terms required for x^4 part... just ignore for now
	// in future, EOS can be modified to also return these derivatives, if needed
	double dPddd = Gamma1(1)*(Gamma1(1)-1.)*exp(logY[0][pres]-2.*logY[0][dens]); //second partial rho
	double dPdtt = 0.0; //second partial T
	double dPddt = (Gamma1(1)-1.)*dPdt*exp(-logY[0][dens]); //the mixed partial derivative
		
	double term = m_pi*rho0*rho0/P0/(nc+1.);
	ac[0] = 1.0;
	ac[2] = -2./3.*term;;
	ac[4] = 2.*nc/15.*pow(term,2);
	ac[6] = -4.*nc*(8.*nc-5.)/(945.)*pow(term,3);
	ac[1] = ac[3] = ac[5] = ac[7] = 0.0;
	
	//terms x^0
	dc[0] = exp(logY[0][dens]);
	pc[0] = exp(logY[0][pres]);
	tc[0] = exp(logY[0][temp]);
	double del = nabla[0];
	printf("Establishing center:\n");
	printf("x^0:\t%le %le %le %le\n", dc[0], pc[0], tc[0], del);
	
	//terms x^2
	tc[1] = -del*tc[0]*dc[0]*dc[0]/6./pc[0];
	pc[1] = -dc[0]*dc[0]/6.;
	dc[1] = (pc[1] - tc[1]*dPdt)/dPdd;
	printf("x^2:\t%le %le %le\n", dc[1], pc[1], tc[1]);
	
	//terms x^4
	tc[2] = (5.*pc[1]*tc[0]*dc[0]*dc[0]-5.*pc[0]*tc[1]*dc[0]*dc[0]-8.*pc[0]*tc[0]*dc[0]*dc[1]);
	tc[2] = del*tc[2]/(60.*pc[0]*pc[0]);
	pc[2] = -dc[0]*dc[1]*2./15.;
	dc[2] = (2.*pc[2]-2.*tc[2]*dPdt
			-tc[1]*tc[1]*dPdtt-2.*tc[1]*dc[1]*dPddt-dc[1]*dc[1]*dPddd)/(2.*dPdd);
	//dc[2] = dc[0]*(nc*ac[4] + 0.5*nc*(nc-1)*ac[2]*ac[2]);
	printf("x^4:\t%le %le %le\n", dc[2], pc[2], tc[2]);
	
	FILE *fp = fopen("./development_tests/testSimpleWD/testCenter/center.txt", "w");
	fprintf(fp, "dens:\t%0.16le\t%0.16le\t%0.16le\n", dc[0],dc[1],dc[2]);
	fprintf(fp, "pres:\t%0.16le\t%0.16le\t%0.16le\n", pc[0],pc[1],pc[2]);
	fprintf(fp, "temp:\t%0.16le\t%0.16le\t%0.16le\n", tc[0],tc[1],tc[2]);
	
	//now make the actual terms
	// c
	c0[0] = 3./dc[0];
	c0[1] = -9.*dc[1]/(5.*dc[0]*dc[0]);
	c0[2] = 27.*dc[1]*dc[1]/(25.*dc[0]*dc[0]*dc[0]) - 9.*dc[2]/(7.*dc[0]*dc[0]);
	// U
	U0[0] = 3.;
	U0[1] = 6.*dc[1]/(5.*dc[0]);
	U0[2] = 12.*dc[2]/(7.*dc[0]) - 18.*dc[1]*dc[1]/(25.*dc[0]*dc[0]);
	// A
	A0[0] = 0.0;
	A0[1] = -2.*dc[1]/dc[0];
	A0[2] =  2.*dc[1]*dc[1]/dc[0]/dc[0] - 4.*dc[2]/dc[0];
	// Vg
	V0[0] = 0.0;
	V0[1] = -2.*pc[1]/pc[0];
	V0[2] =  2.*pc[1]*pc[1]/pc[0]/pc[0] - 4.*pc[2]/pc[0];	
	//
	fprintf(fp, "A*:\t%0.16le\t%0.16le\t%0.16le\n", A0[0],A0[1],A0[2]);
	fprintf(fp, "U :\t%0.16le\t%0.16le\t%0.16le\n", U0[0],U0[1],U0[2]);
	fprintf(fp, "Vg:\t%0.16le\t%0.16le\t%0.16le\n", V0[0],V0[1],V0[2]);
	fprintf(fp, "c1:\t%0.16le\t%0.16le\t%0.16le\n", c0[0],c0[1],c0[2]);
	fclose(fp);
}

void SimpleWD::getAstarCenter(double *Ac, int& maxPow, double g){
	double gam1 = (g==0.0 ? Gamma1(0) : g);
	if(maxPow>=0) Ac[0] = A0[0] - V0[0]/gam1;
	if(maxPow>=2) Ac[1] = A0[1] - V0[1]/gam1;
	if(maxPow>=4) Ac[2] = A0[2] - V0[2]/gam1;
	if(maxPow> 4) maxPow = 4;
}

void SimpleWD::getVgCenter(double *Vc, int& maxPow, double g){
	double gam1 = (g==0.0 ? Gamma1(0) : g);
	if(maxPow>=0) Vc[0] = V0[0]/gam1;
	if(maxPow>=2) Vc[1] = V0[1]/gam1;
	if(maxPow>=4) Vc[2] = V0[2]/gam1;
	if(maxPow> 4) maxPow = 4;
}

void SimpleWD::getUCenter(double *Uc, int& maxPow){
	if(maxPow>=0) Uc[0] = U0[0];
	if(maxPow>=2) Uc[1] = U0[1];
	if(maxPow>=4) Uc[2] = U0[2];
	if(maxPow> 4) maxPow = 4;
}

void SimpleWD::getC1Center(double *cc, int& maxPow){
	if(maxPow>=0) cc[0] = c0[0];
	if(maxPow>=2) cc[1] = c0[1];
	if(maxPow>=4) cc[2] = c0[2];
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
void SimpleWD::setupSurface(){
	
	//for all derivatives, we assume the EOS is ideal gas + radiation
	double Cideal = N_Avogadro*boltzmann_k/Xelem[Ntot-2].mean_A()*Dscale*Tscale/Pscale;
	double Crad   = radiation_a/3.*pow(Tscale,4)/Pscale;
	
	// the x^0 part
	ds[0] = exp(logY[Ntot-1][dens]);
	ts[0] = exp(logY[Ntot-1][temp]);
	ps[0] = exp(logY[Ntot-1][pres]);
	double del = nabla[Ntot-2];

	// the x^1 part
	ts[1] = del*ts[0]*ds[0]/ps[0];
	ps[1] = ds[0];
	ds[1] = (ps[1]-4.*Crad*pow(ts[0],3)*ts[1]-Cideal*ts[1]*ds[0])/(Cideal*ts[0]);
	
	
	// the x^2 part
	double tempprod=1.0;
	int mp = 2;
	ts[2] =(2.*ps[0]*ts[0]*ds[0]-ps[1]*ts[0]*ds[0]+ps[0]*ts[1]*ds[0]-ps[0]*ts[0]*ds[0]*ds[0]+ps[0]*ts[0]*ds[1]);
	ts[2] = del*ts[2]/(2.*ps[0]*ps[0]);
	ps[2] = 0.5*(2.*ds[0] - ds[0]*ds[0] + ds[1]);
	ds[2] = (ps[2]
		 - Crad*(6.*ts[0]*ts[0]*ts[1]*ts[1] + 4.*ts[0]*ts[0]*ts[0]*ts[2])
		 - Cideal*(ts[2]*ds[0] + ts[1]*ds[1]))/(Cideal*ts[0]);
	
	//the x^3 part
	mp = 3;
	ts[3] = -(del/(6.*pow(ps[0],3)))*(
		- 6.0*ps[0]*ps[0]*ts[0]*ds[0] + 4.*ps[0]*ps[1]*ts[0]*ds[0] - 2.*ps[1]*ps[1]*ts[0]*ds[0]
		+ 2.0*ps[0]*ps[2]*ts[0]*ds[0] - 4.*ps[0]*ps[0]*ts[1]*ds[0] + 2.*ps[0]*ps[1]*ts[1]*ds[0]
		- 2.0*ps[0]*ps[0]*ts[2]*ds[0] + 2.*ps[0]*ps[0]*ts[0]*ds[0]*ds[0]
		- 4.0*ps[0]*ps[0]*ts[0]*ds[1] + 2.*ps[0]*ps[1]*ts[0]*ds[1]
		- 2.0*ps[0]*ps[0]*ts[1]*ds[0]*ds[0] - 4.0*ps[0]*ps[0]*ts[0]*ds[1]
		+ 2.0*ps[0]*ps[1]*ts[0]*ds[1] - 2.*ps[0]*ps[0]*ts[1]*ds[1]
		+ 3.0*ps[0]*ps[0]*ts[0]*ds[0]*ds[1] - 2.0*ps[0]*ps[0]*ts[0]*ds[2]);
		//2.*del*ps[1]*ts[1]*ds[1]-del*ps[2]*ts[1]*ds[1]+del*ps[1]*ts[2]*ds[1] + del*ps[1]*ts[1]*ds[2])/(3.*ps[1]*ps[1]);
	ps[3] = ds[0] - ds[0]*ds[0]/3. + 2.*ds[1]/3. - 0.5*ds[0]*ds[1] + ds[2]/3.;
	ds[3] = (ps[3]
		- Crad*(4.*ts[0]*ts[1]*ts[1]*ts[1]+12.*ts[0]*ts[0]*ts[1]*ts[2]+4.*ts[0]*ts[0]*ts[0]*ts[3])
		- Cideal*(ts[3]*ds[0]+ts[2]*ds[1]+ts[1]*ds[2]))/(Cideal*ts[0]);
	
	// the x^4 part
	mp = 4;
	ts[4] = (del/(24.*pow(ps[0],4) ))*
 	(-6.*pow(ps[1],3)*ts[0]*ds[0] + 6.*ps[0]*ps[1]*(2.*ps[2]*ts[0]*ds[0]
 		 + ps[1]*(ts[1]*ds[0] + ts[0]*(-(-2. + ds[0])*ds[0] + ds[1]))) + 
  		3.*ps[0]*ps[0]*(-2*ps[3]*ts[0]*ds[0] + 2.*ps[2]*((-ts[1] + ts[0]*(-2 + ds[0]))*ds[0] - ts[0]*ds[1])
  		+ ps[1]*(2*(-ts[2] + ts[1]*(-2 + ds[0]))*ds[0] - 2.*ts[1]*ds[1]
  		+ ts[0]*(ds[0]*(-6. + 2.*ds[0] + 3.*ds[1]) - 2.*(2.*ds[1] + ds[2]))))
  		+ pow(ps[0],3)*(6.*ts[3]*ds[0] + 6.*ts[2]*(-(-2. + ds[0])*ds[0] + ds[1])
  		+ 3.*ts[1]*(ds[0]*(6. - 2.*ds[0] - 3.*ds[1]) + 2.*(2.*ds[1] + ds[2])) 
  		+ ts[0]*(-8.*ds[0]*ds[0] - 3.*(-6. + ds[1])*ds[1] - 8.*ds[0]*(-3. + ds[1] + ds[2]) 
  		+ 6.*(2.*ds[2] + ds[3])))); //magical formula from mathematica
	ps[4] = (24.*ds[0] - 8.*ds[0]*ds[0] + 18.*ds[1] - 8.*ds[0]*ds[1] - 3.*ds[1]*ds[1] 
			+ 12.*ds[2] - 8.*ds[0]*ds[2] + 6.*ds[3])/24.;
	ds[4]= ps[4];
	tempprod=0.0;
	for(int a=0; a<=mp; a++){
		for(int b=0; b<=mp-a; b++){
			for(int c=0; c<=mp-a-b; c++){
				tempprod += ts[a]*ts[b]*ts[c]*ts[mp-a-b-c];
			}
		}
	}
	ds[mp] -= Crad*tempprod;
	for(int a=0; a<mp; a++){
		ds[mp] -= Cideal*ts[mp-a]*ds[a];
	}
	ds[mp] /= (Cideal*ts[0]);
		
	//now set the structure variables
	// c
	c1[0] =  1.;
	c1[1] = -3. + ds[0];
	c1[2] =  3. - 4.*ds[0] + ds[0]*ds[0] + 0.5*ds[1];
	c1[3] = -1. + 19.*ds[0]/3. - 5.*ds[0]*ds[0] + ds[0]*ds[0]*ds[0] 
			- 13.*ds[1]/6. + ds[0]*ds[1] + ds[2]/3.;
	c1[4] = 0.; // does not atually appear in equations
	// U
	U1[0] = ds[0];
	U1[1] = ds[1] + ds[0]*ds[0] - 3.*ds[0];
	U1[2] = 3.*ds[0] - 4.*ds[0]*ds[0] + ds[0]*ds[0]*ds[0] - 3.*ds[1] + 1.5*ds[0]*ds[1] + ds[2];
	U1[3] = ds[3] + 4./3.*ds[0]*ds[2] - 3.*ds[2] + ds[1]*ds[1]/2. + 2.*ds[0]*ds[0]*ds[1] - 37./6.*ds[0]*ds[1]
			+ 3.*ds[1] + pow(ds[0],4) - 5.*pow(ds[0],3) + 19.*ds[0]*ds[0]/3. - ds[0];
	U1[4] = ds[4] + 1.25*ds[0]*ds[3] - 3.*ds[3] + 5./6.*ds[1]*ds[2] + 5./3.*ds[0]*ds[0]*ds[2]
			- 5.5*ds[0]*ds[2] + 3.*ds[2] + 1.25*ds[0]*ds[1]*ds[1] - 13./6.*ds[1]*ds[1] + 2.5*pow(ds[0],3)*ds[1]
			- 31./3.*ds[0]*ds[0]*ds[1] + 121./12.*ds[0]*ds[1] - ds[1] + pow(ds[0],5) - 6.*pow(ds[0],4)
			+ 32./3.*pow(ds[0],3) - 5.*ds[0]*ds[0];
	int O = 1;
	// Vg
	ps[1] /= ps[0];
	ps[2] /= ps[0];
	ps[3] /= ps[0];
	ps[4] /= ps[0];
	V1[O-1] = 0.0;
	V1[O+0] = ps[1];
	V1[O+1] = 2.*ps[2] - ps[1]*ps[1] - ps[1];
	V1[O+2] =  ps[1]*ps[1]*(1.+ps[1]) - 2.*ps[2] - 3.*ps[1]*ps[2] + 3.*ps[3];
	V1[O+3] =  3.*ps[1]*ps[2] - pow(ps[1],3)*(1.+ps[1]) + 4.*ps[1]*ps[1]*ps[2]
				-2.*ps[2]*ps[2] - 3.*ps[3] - 4.*ps[1]*ps[3] + 4.*ps[4];
	// A*
	ds[1] /= ds[0];
	ds[2] /= ds[0];
	ds[3] /= ds[0];
	ds[4] /= ds[0];
	A1[O-1] = 0.0;
	A1[O+0] = ds[1];
	A1[O+1] = 2.*ds[2] - ds[1]*ds[1] - ds[1];
	A1[O+2] = ds[1]*ds[1]*(1.+ds[1]) - 2.*ds[2] - 3.*ds[1]*ds[2] + 3.*ds[3];
	A1[O+3] = 3.*ds[1]*ds[2] - pow(ds[1],3)*(1.+ds[1]) + 4.*ds[1]*ds[1]*ds[2]
				 - 2.*ds[2]*ds[2] - 3.*ds[3] - 4.*ds[1]*ds[3] + 4.*ds[4];
	//using the Brassard relation instead
	A1[O-1] = 0.0;
	A1[O+0] = brunt_vaisala[Ntot-2];
	A1[O+1] = A1[O+0]*(ds[0]-3.);
	A1[O+2] = A1[O+0]*(3. - 4.*ds[0] + ds[0]*ds[0] + 0.5*ds[1]);
	A1[O+3] = A1[O+0]*(-1. + 19./3.*ds[0] - 5.*ds[0]*ds[0] + ds[0]*ds[0]*ds[0] - 13./6.*ds[1] + ds[0]*ds[1] + ds[2]/3.);	
	
	for(int i=1; i<=4; i++) { ds[i] *= ds[0]; ps[i]*=ps[0];}
}

void SimpleWD::getAstarSurface(double *As, int& maxPow, double g){
	double gam1 = (g==0.0 ? Gamma1(Ntot-1) : g);
	int O=1;
	if(maxPow>= -1) As[O-1] = A1[O-1];// - V1[O-1]/gam1;
	if(maxPow>=  0) As[O  ] = A1[O  ];// - V1[O  ]/gam1;
	if(maxPow>=  1) As[O+1] = A1[O+1];// - V1[O+1]/gam1;
	if(maxPow>=  2) As[O+2] = A1[O+2];// - V1[O+2]/gam1;
	if(maxPow>=  3) As[O+3] = A1[O+3];// - V1[O+3]/gam1;
	//if more  terms than this requested, cap number of terms
	if(maxPow > 3) maxPow = O+3;
}

void SimpleWD::getVgSurface(double *Vs, int& maxPow, double g){
	double gam1 = (g==0.0 ? Gamma1(Ntot-1) : g);
	int O=1;
	if(maxPow>= -1) Vs[O-1] = V1[O-1]/gam1;
	if(maxPow>=  0) Vs[O  ] = V1[O  ]/gam1;
	if(maxPow>=  1) Vs[O+1] = V1[O+1]/gam1;
	if(maxPow>=  2) Vs[O+2] = V1[O+2]/gam1;
	if(maxPow>=  3) Vs[O+3] = V1[O+3]/gam1;
	//if more  terms than this requested, cap number of terms
	if(maxPow > 3) maxPow = O+3;
}

void SimpleWD::getUSurface(double *Us, int& maxPow){
// coefficients of U must extend up to order maxPow	
	if(maxPow>= 0) Us[0] = U1[0];
	if(maxPow>= 1) Us[1] = U1[1];
	if(maxPow>= 2) Us[2] = U1[2];
	if(maxPow>= 3) Us[3] = U1[3];
	if(maxPow>= 4) Us[4] = U1[4];
	//if more  terms than this requested, cap number of terms
	if(maxPow > 4) maxPow = 4;
}

void SimpleWD::getC1Surface(double *cs, int& maxPow){
// coefficients of c1 are only needed up to order maxPow-1
	if(maxPow>=0) cs[0] = c1[0];// 1.;
	if(maxPow>=1) cs[1] = c1[1];//-3.;
	if(maxPow>=2) cs[2] = c1[2];// 3.;
	if(maxPow>=3) cs[3] = c1[3];//-1.;
	if(maxPow>=4) cs[4] = c1[4];// 0.; //does not actually appear in equations
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;	
}

// **************************  WRITE STAR  ***************************************
// will print the following graphs of the star, using gnuplot
// 1) graph of T,m,rho vs. radius
// 2) graph of chemical composition vs. inverse mass fraction
// 3) graph of chemical composition vs. radius
// 4) graph of Brunt-Vaisala and Lamb frequency vs. logP (to compare to Goldreich & Wu)
// optional argument c[] is calculation directory where files should be written
// if no argument, files written to ./out/[name], where name is shown in constructor
void SimpleWD::printChem(char *c){
	char filename[256];
	char rootname[256];
	char txtname[256];
	char outname[256];
	sprintf(filename, "%s", c);
	char title[256]; graph_title(title);
	
	//print the chemical gradient
	sprintf(txtname, "%s/chemical.txt", filename);
	sprintf(outname, "%s/chemical.png", filename);
	FILE* fp  = fopen(txtname, "w");
	double H,He;
	double MM = logY[Ntot-1][mass];
	fprintf(fp, "num\t1-r\t1-m\tm\tH\tHe\tC\tO\n");
	for(int X=Ntot-1; X >= 0; X--){
	//	if(isnan(Xelem[X].H1)) {MM=logY[X-1][mass]; continue;};
		fprintf(fp, "%d", X);									//col 1
		fprintf(fp, "\t%0.24le", (1.-exp(logY[X][radi])));	//col 2
		fprintf(fp, "\t%0.24le", (1.-exp(logY[X][mass])));	//col 3
		fprintf(fp, "\t%0.16le", exp(logY[X][mass]));			//col 4
		fprintf(fp, "\t%0.16le\t%0.16le\t%0.16le\t%0.16le",
			Xelem[X].H1, Xelem[X].He4, Xelem[X].C12, Xelem[X].O16);	//col 5-8
		fprintf(fp, "\t%d\n", Xelem[X].e);
	}
	fclose(fp);	
	//plot file in png in gnuplot
	FILE* gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Chemical composition for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'log_{10}(1-m/M)'\n");
	fprintf(gnuplot, "set ylabel 'mass fraction\n");
	//fprintf(gnuplot, "set xlabel 'log_{10} r/R'\n");
	fprintf(gnuplot, "set logscale x 10\n");
	fprintf(gnuplot, "set format x '%%L'\n");
	fprintf(gnuplot, "set xrange [1e-8:1e0]\n");
	//fprintf(gnuplot, "set xtics 1.0\n set ytics 0.1\n");
	//fprintf(gnuplot, "set hidden3d noundefined\n");
	fprintf(gnuplot, "set yrange [-0.01: 1.01]\n");
	fprintf(gnuplot, "plot '%s' u 3:5 w l t 'X (H1)'", txtname);
	fprintf(gnuplot, ",    '%s' u 3:6 w l t 'Y (He4)'", txtname);
	fprintf(gnuplot, ",    '%s' u 3:7 w l t 'Z (C12)'", txtname);
	fprintf(gnuplot, ",    '%s' u 3:8 w l t 'Z (O16)'", txtname);
	//fprintf(gnuplot, ",    '%s' u 3:($9/4) w l t 'excl'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "set output '%s/chem_log.png'\n", filename);
	fprintf(gnuplot, "set xlabel 'log_{10}(1-m/M)'\n");
	fprintf(gnuplot, "set ylabel 'log_{10} X_i'\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set logscale y 10\n set logscale x 10\n");
	fprintf(gnuplot, "set xrange [1e-8:1e0]\n");
	fprintf(gnuplot, "set yrange [1e-8:1.01]\n");
	fprintf(gnuplot, "set ytics 10\n");
	//fprintf(gnuplot, "set xtics 0.1\n");
	fprintf(gnuplot, "plot '%s' u 3:5 w l t 'X (H1)'", txtname);
	fprintf(gnuplot, ",    '%s' u 3:6 w l t 'Y (He4)'", txtname);
	fprintf(gnuplot, ",    '%s' u 3:7 w l t 'Z (C12)'", txtname);
	fprintf(gnuplot, ",    '%s' u 3:8 w l t 'Z (O16)'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);	
}
void SimpleWD::printBV(char *c){
	char filename[256];
	char rootname[256];
	char txtname[256];
	char outname[256];
	sprintf(filename, "%s", c);
	char title[256]; graph_title(title);
	
	//print the Brunt-Vaisala frequency
	sprintf(txtname, "%s/BruntVaisala.txt", filename);
	sprintf(outname, "%s/BruntVaisala.png", filename);
	FILE* fp  = fopen(txtname, "w");
	double N2 = -1.0, f0 = G_CGS*Mstar*pow(Rstar,-3);
	fprintf(fp, "1-m\tN2\tL1\n");
	for(int X=0; X< Ntot; X++){
		N2 = brunt_vaisala[X]*f0;
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\n",
			(1.-exp(logY[X][mass])),
			N2,
			2.*sound_speed2(X,0.)*exp(-2.*logY[X][radi])*pow(Rstar,-2));
	}
	fclose(fp);	
	//plot file in png in gnuplot
	FILE* gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	//fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Brunt-Vaisala and Lamb for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'log_{10}(1-m/M)'\n");
	fprintf(gnuplot, "set ylabel 'log_{10} N^2 & log_{10} L_1^2 (Hz^2)\n");
	fprintf(gnuplot, "set logscale x 10\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "set format x '%%L'\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
//	fprintf(gnuplot, "set xrange [1e-14:1e0]\n");
//	fprintf(gnuplot, "set yrange [1e-8:1e2]\n");
	fprintf(gnuplot, "set ytics 10\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'N^2'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'L_1^2'", txtname);
	//fprintf(gnuplot, ",    '%s' u 1:(-$2) w l t '-N^2'", txtname);
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
}
void SimpleWD::printOpacity(char *c){
	char filename[256];
	char rootname[256];
	char txtname[256];
	char outname[256];
	sprintf(filename, "%s", c);
	char title[256]; graph_title(title);
	
	//make a graph of T, nabla, opacity
	sprintf(txtname, "%s/opacity.txt", filename);
	sprintf(outname, "%s/opacity.png", filename);
	FILE* fp  = fopen(txtname, "w");
	Abundance X;
	StellarVar ly;
	double logT=logY[0][temp];
	double GC = 0.0;
	fprintf(fp, "1-m\tT/T1\tdel\tdel_ad\tkappa\tH\tHe\tGamma1\n");
	for(int x=0; x < Ntot; x++){
		X  = Xelem[x];
		ly = logY[x];
		fprintf(fp, "%0.24le\t%0.16le", (1.-exp(ly[mass])), exp(ly[temp]-logT));	//col1, col2
		//get the layer, for using in energyTransport()
		fprintf(fp, "\t%0.16le\t%0.16le", nabla[x], nabla_ad[x]); //col3, col4
		//form the approximation for opacity
		ly = logY[x] + logYscale;
		GC = 2.275e5*(X.mean_coulomb())*exp(ly[dens]/3.-ly[temp]);
		fprintf(fp, "\t%0.16le", kappa[x]); //col5
		fprintf(fp, "\t%0.16le\t%0.16le", X.H1, X.He4);  //col6, col7
		fprintf(fp, "\t%0.16le", adiabatic_1[x]); //col8
		fprintf(fp, "\t%0.16le", GC);
		fprintf(fp, "\n");
	}
	fclose(fp);	
	FILE* gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'opacity, temperature for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'log(1-m/M)'\n");
	fprintf(gnuplot, "set y2label 'temperature'\n");
	fprintf(gnuplot, "set logscale x 10\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "set xrange [1e-14:1e0]\n");
//	fprintf(gnuplot, "set yrange [1e-5:1e2]\n");
	fprintf(gnuplot, "set y2range [0:1]\n");
	fprintf(gnuplot, "set ytics 10 nomirror\n");
	fprintf(gnuplot, "set y2tics 0.1 nomirror\n");
	fprintf(gnuplot, "set format x '%%L'\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'T' axes x1y2", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'nabla'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:4 w l t 'nabla_{ad}'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:5 w l t 'opacity'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:6 w l t 'H' lc 7 axes x1y2", txtname);
	fprintf(gnuplot, ",    '%s' u 1:8 w l t 'Gamma1' lc 10", txtname);
	fprintf(gnuplot, ",    '%s' u 1:9 w l t 'Gamma_C' lc 11", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
}
void SimpleWD::printCoefficients(char *c){
	char filename[256];
	char rootname[256];
	char txtname[256];
	char outname[256];

	char command[300];
	sprintf(filename, "%s/wave_coefficient", c);
	sprintf(command, "mkdir %s", filename);
	system(command);
	
	char title[256]; graph_title(title);
	
	//print the coefficients of the center and surface, for series analysis
	sprintf(txtname, "%s/center.txt", filename);
	FILE *fp = fopen(txtname, "w");
	double gam1 = adiabatic_1[0];
	fprintf(fp, "dens:\t%0.16le\t%0.16le\t%0.16le\n", dc[0],dc[1],dc[2]);
	fprintf(fp, "pres:\t%0.16le\t%0.16le\t%0.16le\n", pc[0],pc[1],pc[2]);
	fprintf(fp, "temp:\t%0.16le\t%0.16le\t%0.16le\n", tc[0],tc[1],tc[2]);
	fprintf(fp, "A*:\t%0.16le\t%0.16le\t%0.16le\n", (A0[0]-V0[0]/gam1),(A0[1]-V0[1]/gam1),(A0[2]-V0[2]/gam1));
	fprintf(fp, "U :\t%0.16le\t%0.16le\t%0.16le\n", U0[0],U0[1],U0[2]);
	fprintf(fp, "Vg:\t%0.16le\t%0.16le\t%0.16le\n", V0[0]/gam1,V0[1]/gam1,V0[2]/gam1);
	fprintf(fp, "c1:\t%0.16le\t%0.16le\t%0.16le\n", c0[0],c0[1],c0[2]);
	fclose(fp);
	sprintf(txtname, "%s/surface.txt", filename);
	fp = fopen(txtname, "w");
	gam1 = adiabatic_1[Ntot-1];
	fprintf(fp, "dens:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", ds[0],ds[1],ds[2],ds[3],ds[4]);
	fprintf(fp, "pres:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", ps[0],ps[1],ps[2],ps[3],ps[4]);
	fprintf(fp, "temp:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", ts[0],ts[1],ts[2],ts[3],ts[4]);
	fprintf(fp, "A*:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", A1[0],A1[1],A1[2],A1[3],A1[4]);
	fprintf(fp, "U :\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", U1[0],U1[1],U1[2],U1[3],U1[4]);
	fprintf(fp, "Vg:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", V1[0]/gam1,V1[1]/gam1,V1[2]/gam1,V1[3]/gam1,V1[4]/gam1);
	fprintf(fp, "c1:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", c1[0],c1[1],c1[2],c1[3],c1[4]);
	fclose(fp);
	
	//print fits to those coefficients at center and surface
	int NC=15, NS=15;
	sprintf(txtname, "%s/centerfit.txt", filename);
	fp = fopen(txtname, "w");
	double x2 = exp(2.*logY[0][radi]);
	gam1 = adiabatic_1[0];
	fprintf(fp, "x\trho\trho_fit\tP\tP_fit\tT\tT_fit\tA*\tA*_fit\tU\tU_fit\tVg\tVg_fit\tc1\tc1_fit\n");
	for(int X=0; X<NC; X++){
		x2 = exp(2.*logY[X][radi]);
		fprintf(fp, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
			exp(logY[X][radi]),
			exp(logY[X][dens]),
			dc[0]+dc[1]*x2+dc[2]*x2*x2,
			exp(logY[X][pres]),
			pc[0]+pc[1]*x2+pc[2]*x2*x2,
			exp(logY[X][temp]),
			tc[0]+tc[1]*x2+tc[2]*x2*x2,
			brunt_vaisala[X]*exp(3.*logY[X][radi]-logY[X][mass]),
			(A0[0]-V0[0]/gam1) + (A0[1]-V0[1]/gam1)*x2 + (A0[2]-V0[2]/gam1)*x2*x2,
			getU(X),
			U0[0] + U0[1]*x2 + U0[2]*x2*x2,
			getVg(X,gam1),
			V0[0]/gam1 + V0[1]/gam1*x2 + V0[2]/gam1*x2*x2,
			getC(X),
			c0[0] + c0[1]*x2 + c0[2]*x2*x2
		);
	}
	fclose(fp);
	sprintf(txtname, "%s/surfacefit.txt", filename);
	fp = fopen(txtname, "w");
	double t = 1.-exp(logY[0][radi]);
	gam1 = adiabatic_1[Ntot-1];
	fprintf(fp, "x\trho\trho_fit\tP\tP_fit\tT\tT_fit\tA*\tA*_fit\tU\tU_fit\tVg\tVg_fit\tc1\tc1_fit\n");
	for(int X=Ntot-1; X>=Ntot-NS-1; X--){
		t = 1. - exp(logY[X][radi]);
		fprintf(fp, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
			exp(logY[X][radi]),
			exp(logY[X][dens]),
			ds[0]+ds[1]*t+ds[2]*t*t+ds[3]*t*t*t+ds[4]*t*t*t*t,
			exp(logY[X][pres]),
			ps[0]+ps[1]*t+ps[2]*t*t+ps[3]*t*t*t+ps[4]*t*t*t*t,
			exp(logY[X][temp]),
			ts[0]+ts[1]*t+ts[2]*t*t+ts[3]*t*t*t+ts[4]*t*t*t*t,
			brunt_vaisala[X]*exp(3.*logY[X][radi]-logY[X][mass]),
			(A1[0])/t + (A1[1]) + (A1[2])*t + (A1[3])*t*t + (A1[4])*t*t*t,
			dlogY[X][mass],
			U1[0] + U1[1]*t + U1[2]*t*t + U1[3]*t*t*t + U1[4]*t*t*t*t,
			-dlogY[X][pres]/gam1,
			(V1[0]/gam1)/t + (V1[1]/gam1) + (V1[2]/gam1)*t + (V1[3]/gam1)*t*t + (V1[4]/gam1)*t*t*t,
			exp(3.*logY[X][radi]-logY[X][mass]),
			c1[0] + c1[1]*t + c1[2]*t*t + c1[3]*t*t*t + c1[4]*t*t*t*t
		);
	}
	fclose(fp);
	
	//print the pulsation coeffcients
	sprintf(txtname, "%s/coefficients.txt", filename);
	sprintf(outname, "%s/coefficients.png", filename);
	fp  = fopen(txtname, "w");
	fprintf(fp, "m\tA*\tU\tVg\tc1\n");
	for(int X=0; X<Ntot; X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			exp(logQ[X]),
			-exp(logY[X][radi])*Schwarzschild_A(X, 0.)*Rstar,
			getU(X),
			getVg(X, 0.),
			getC(X)
		);
	}
	fclose(fp);	
	//plot file in png in gnuplot
	//gnuplot = popen("gnuplot -persist", "w");
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Pulsation Coefficients for %s'\n", title);
	//fprintf(gnuplot, "set xlabel 'log_{10} r/R'\n");
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
	sprintf(txtname, "%s/centerfit.txt", filename);
	sprintf(outname, "%s/centerfit.png", filename);
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Central Fitting by Power Series for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'difference'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 100\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xrange [0:%le]\n", 1.01*exp(logY[NC][radi]));
	fprintf(gnuplot, "plot '%s' u 1:(abs($2-$3)/abs($2)) w lp t 'rho'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($4-$5)/abs($4)) w lp t 'P'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($6-$7)/abs($6)) w lp t 'T'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($8-$9)/abs($8)) w lp t 'A*'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($10-$11)/abs($10)) w lp t 'U'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($12-$13)/abs($12)) w lp t 'Vg'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($14-$15)/abs($14)) w lp t 'c1'", txtname);
	fprintf(gnuplot, "\n");
	sprintf(txtname, "%s/surfacefit.txt", filename);
	sprintf(outname, "%s/surfacefit.png", filename);
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Surface Fitting by Power Series for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'difference'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 100\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xrange [%le:1]\n", exp(logY[Ntot-NS-1][radi]));
	fprintf(gnuplot, "plot '%s' u 1:(abs($2-$3)/abs($2)) w lp t 'rho'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($4-$5)/abs($4)) w lp t 'P'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($6-$7)/abs($6)) w lp t 'T'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($8-$9)/abs($8)) w lp t 'A*'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($10-$11)/abs($10)) w lp t 'U'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($12-$13)/abs($12)) w lp t 'Vg'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:(abs($14-$15)/abs($14)) w lp t 'c1'", txtname);
	fprintf(gnuplot, "\n");
	
	//now leave gnuplot
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);	
}
void SimpleWD::printBigASCII(char *c){
	char filename[256];
	char rootname[256];
	char txtname[256];
	char outname[256];
	sprintf(filename, "%s", c);
	char title[256]; graph_title(title);

	sprintf(txtname, "%s/M%0.3le_R%0.2le_Teff%4.0lf.txt", filename, Msolar, Rsolar, Teff);
	FILE* fp = fopen(txtname, "w");
	fprintf(fp, "## SimpleWD output file\n");
	fprintf(fp, "# Mass \t%le (g)\t%0.3le (Msun) \n", Mstar, Msolar);
	fprintf(fp, "# Radius \t%le (cm)\t%1.2le (Rearth) \n", Rstar, Rsolar);
	fprintf(fp, "# Teff \t%le\n", Teff);
	fprintf(fp, "# Layer masses:\t H: %le\tHe: %le\tC: %le\tO: %le\n", Xmass.H1, Xmass.He4, Xmass.C12, Xmass.O16);
	fprintf(fp, "#######################################################\n");
	//below we print a title for each column, centered on the column
	fprintf(fp, "# N   \t");
	fprintf(fp, "log(1-m/M)   \t");
	fprintf(fp, "log(1-r/R)   \t");
	fprintf(fp, "m(r)         \t");
	fprintf(fp, "radius       \t");
	fprintf(fp, "log g        \t");
	fprintf(fp, "log rho      \t");
	fprintf(fp, "log T        \t");
	fprintf(fp, "log P        \t");
	fprintf(fp, "log P_ideal  \t");
	fprintf(fp, "log P_rad    \t");
	fprintf(fp, "log P_deg    \t");
	fprintf(fp, "log(-P_coul) \t");
	fprintf(fp, "log U        \t");
	fprintf(fp, "log U_ideal  \t");
	fprintf(fp, "log U_rad    \t");
	fprintf(fp, "log U_deg    \t");
	fprintf(fp, "log U_coul   \t");
	fprintf(fp, "x (p_F/mec)  \t");
	fprintf(fp, "eta=mu/kT    \t");
	fprintf(fp, "log(epsilon) \t");
	fprintf(fp, "log(X_H)     \t");
	fprintf(fp, "log(X_He)    \t");
	fprintf(fp, "log(X_C)     \t");
	fprintf(fp, "log(X_O)     \t");
	fprintf(fp, "kappa        \t");
	fprintf(fp, "kappa_rad    \t");
	fprintf(fp, "kappa_cond   \t");
	fprintf(fp, "Gamma1       \t");
	fprintf(fp, "Gamma3       \t");
	fprintf(fp, "del          \t");
	fprintf(fp, "del_ad       \t");
	fprintf(fp, "Ledoux B     \t");
	fprintf(fp, "log(N^2)     \t");
	fprintf(fp, "log(L1^2)    \t");
	fprintf(fp, "chi_T        \t");
	fprintf(fp, "chi_rho      \t");
	fprintf(fp, "chi_H        \t");
	fprintf(fp, "chi_He       \t");
	fprintf(fp, "chi_C        \t");
	fprintf(fp, "chi_O        \n");
	double log10e = log10(exp(1.0));
	for(int X=0; X<Ntot; X++){
		fprintf(fp, "%6d\t", X+1);
		//the independent variables
		fprintf(fp, "%+13le\t%+13le\t%+13le\t%+13le\t", 
			log10(1.-exp(logY[X][mass])), 
			log10(1.-exp(logY[X][radi])),
			Mstar*exp(logY[X][mass]), 
			Rstar*exp(logY[X][radi])
		);
		//gravitational field
		if(X==0) fprintf(fp, "%+13le\t", 0.0);
		else 
			fprintf(fp, "%+13le\t", (logY[X][mass]-2.*logY[X][radi]+log(G_CGS)+logYscale[mass]-2.*logYscale[radi])*log10e);
		//log density
		fprintf(fp, "%+13le\t", (logY[X][dens]+logYscale[dens])*log10e);
		//log temperature
		fprintf(fp, "%+13le\t", (logY[X][temp]+logYscale[temp])*log10e);	
		//logs of all partial pressures
		EOS* eos = getEOS(exp(logY[X]+logYscale),Xelem[X]);
		bool incore = &core_pressure==eos;
		double rho=Dscale*exp(logY[X][dens]), T=Tscale*exp(logY[X][temp]);
		fprintf(fp, "%+13le\t%+13le\t%+13le\t%+13le\t%+13le\t", (logY[X][pres]+logYscale[pres])*log10e, 
			log10(pressure_ideal(rho,T, Xelem[X])),
			log10(pressure_rad(rho,T, Xelem[X])),
			incore ? log10(pressure_deg_zero(rho,T, Xelem[X])) : 0.0,
			incore ? log10(-pressure_coul(rho,T, Xelem[X])) : 0.0
		);
		// the log internal energy density U goes here
		fprintf(fp, "%+13le\t%+13le\t%+13le\t%+13le\t%+13le\t",
			log10(eos->U(rho,T,Xelem[X])),
			log10(energy_ideal(rho,T,Xelem[X])),
			log10(energy_rad(rho,T,Xelem[X])),
			incore ? log10(energy_deg_zero(rho,T,Xelem[X])) : 0.0,
			incore ? log10(-energy_coul(rho,T,Xelem[X])) : 0.0
		);
		// the fermi-momentum goes here
		fprintf(fp, "%+le\t", pow(rho/Chandrasekhar::B0/Xelem[X].mu_e(), 1./3.));
		double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
		double eta = findEta(rho, beta, Xelem[X].mu_e());
		// the degeneracy parameter goes here
		fprintf(fp, "%+13le\t", eta);
		//log effective energy generation
		fprintf(fp, "%+13le\t", log10(Lstar/Mstar));
		//log chemical abundance
		Abundance logX( log10(Xelem[X][0]), log10(Xelem[X][1]), log10(Xelem[X][2]),log10(Xelem[X][3]) );
		fprintf(fp, "%+13le\t%+13le\t%+13le\t%+13le\t", logX[0],logX[1],logX[2],logX[3]);
		//opacity, kappa_rad, kappa_cond
		fprintf(fp, "%+13le\t%+13le\t%+13le\t", kappa[X], 
			radiative_opacity(logY[X]+logYscale,Xelem[X]),
			conductive_opacity(logY[X]+logYscale,Xelem[X])
		);
		//Gamma1, Gamma_3, nabla, nabla_ad
		fprintf(fp, "%+13le\t%+13le\t%+13le\t%+13le\t", adiabatic_1[X], eos->Gamma3(rho,T,Xelem[X]), nabla[X],nabla_ad[X]);
		fprintf(fp, "%+13le\t", ledoux[X]);
		//log BV and Lamb
		double logWscale = log(G_CGS)+logYscale[mass]-3.*logYscale[radi];
		fprintf(fp, "%+13le\t%+13le\t", 
			(log(brunt_vaisala[X])+logWscale)*log10e, 
			(log(2.*adiabatic_1[X])+logY[X][pres]-logY[X][dens]-2.*logY[X][radi]+logWscale)*log10e
		);
		//chi_T, chi_rho
		fprintf(fp, "%+13le\t%+13le\t", eos->chiT(rho,T,Xelem[X]), eos->chiRho(rho,T,Xelem[X]));
		fprintf(fp, "%+13le\t%+13le\t%+13le\t%+13le\n",
			eos->chiY(chemical::h,rho,T,Xelem[X]),
			eos->chiY(chemical::he,rho,T,Xelem[X]),
			eos->chiY(chemical::c,rho,T,Xelem[X]),
			eos->chiY(chemical::o,rho,T,Xelem[X])
		);
	}
	fclose(fp);
}
void SimpleWD::writeStar(char *c){
	printf("PRINTING STAR\n");
	//create names for files to be opened
	char filename[256];
	char rootname[256];
	char txtname[256];
	char outname[256];
	if(c==NULL)	sprintf(filename, "./out/%s", name);
	else{
		sprintf(filename, "./%s/star", c);
	}
	sprintf(txtname, "%s/star.txt", filename);
	sprintf(outname, "%s/star.png", filename);

	char title[256]; graph_title(title);

	//prepare the output directory, making sure it exists
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
	double T0 = exp(logY[0][temp]);
	double RS = exp(logY[Ntot-1][radi])*Rstar;
	double logT = logY[0][temp];
	double logL = logY[Ntot-1][lumi];
	double logM = logY[Ntot-1][mass];
	double logR = logY[Ntot-1][radi];
	fprintf(fp, "num\tq\trho\tradius\tP\tm\tT\tL\tm\n");
	for(int X=0; X< Ntot; X++){
		fprintf(fp, "%d\t", X);
		fprintf(fp, "%0.16le\t", (1.-exp(logQ[X])));	//col 2
		fprintf(fp, "%0.16le\t", exp(logY[X][dens]));	//col 3
		fprintf(fp, "%0.16le\t", exp(logY[X][radi]));	//col 4
		fprintf(fp, "%0.16le\t", exp(logY[X][pres]));	//col 5
		fprintf(fp, "%0.16le\t", exp(logY[X][mass]));	//col 6
		fprintf(fp, "%0.16le\t", exp(logY[X][temp]));	//col 7
		fprintf(fp, "%0.16le\t", exp(logY[X][lumi]));	//col 8
		fprintf(fp, "%0.16le\n", exp(logQ[X]));			//col 9
		//fflush(fp);
	}
	fclose(fp);	
	//plot file in png in gnuplot, and open png
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Profile for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'log_{10}(1-m/M)'\n");
	fprintf(gnuplot, "set ylabel 'rho, P'\n");
	fprintf(gnuplot, "set logscale y 10\n set ytics 10 nomirror\nset format y '10^{%%L}'\n");
	fprintf(gnuplot, "set logscale x 10\n set format x '%%L'\n");
	//fprintf(gnuplot, "set yrange [1e0:1e24]\n");
	//fprintf(gnuplot, "set xrange [1e-8:1e0]\n");
	fprintf(gnuplot, "set y2range [-0.01:1.01]\nset y2tics 0.2 nomirror\n");
	fprintf(gnuplot, "set y2label 'normalized temperature, radius'\n");
	fprintf(gnuplot, "plot ");
	fprintf(gnuplot, "     '%s' u 2:($7/%le) w l t 'T' axes x1y2", txtname, T0);
	fprintf(gnuplot, ",    '%s' u 2:3 w l t 'rho'", txtname);
	fprintf(gnuplot, ",    '%s' u 2:5 w l t 'P'", txtname);
	fprintf(gnuplot, ",    '%s' u 2:($4/%le) w l t 'r' axes x1y2", txtname, 1.0);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	
	//call files to print other properties
	printChem(filename);
	printBV(filename);
	printOpacity(filename);
	printCoefficients(filename);	
	//printBigASCII(filename);
}

//calculte a scaled sum square residual from backsubstitution into physical equations
double SimpleWD::SSR(){	
	double checkEuler;
	double checkPoiss;
	int len = Ntot;
			
	//sum up errors in equations
	checkEuler = 0.0;
	checkPoiss = 0.0;
	double d2Phi = 0.0;
	double e1, e2, n1, n2;
	FILE *fp = fopen("SSR.txt", "w");
	for(int X=4; X<len-4; X++){
		//Euler equation
		e1 = fabs(dPdr(X) + rho(X)*dPhidr(X) );
		n1 = fabs(dPdr(X)) + fabs(rho(X)*dPhidr(X));
		//Poisson equation
		//calculate numerical derivatives
		d2Phi = Gee()*dlogY[X][mass]*pow(rad(X),-3)*mr(X) - 2.*Gee()*mr(X)*pow(rad(X),-3);
		e2 = fabs( 4.0*Gee()*m_pi*rho(X)*rad(X)  -      2.0*dPhidr(X)  -      d2Phi*rad(X) );
		n2 = fabs( 4.0*Gee()*m_pi*rho(X)*rad(X)) + fabs(2.0*dPhidr(X)) + fabs(d2Phi*rad(X) );
		//add absolute error
		e1 = e1/n1;
		e2 = e2/n2;
		fprintf(fp, "%d\t%le\t%le\n", X, e1, e2);
		checkEuler += e1*e1;
		checkPoiss += e2*e2;
	}
	fclose(fp);
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", "SSR.png");
	char title[256]; graph_title(title);
	fprintf(gnuplot, "set title 'error for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'gird point'\n");
	fprintf(gnuplot, "set ylabel 'error'\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set arrow 1 from %d, graph 0 to %d, graph 1 lc rgb 'red' nohead\n", Ncore, Ncore);
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'Euler error'", "SSR.txt");
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'Poisson error'", "SSR.txt");
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	return sqrt((checkPoiss+checkEuler)/double(2*len));
}

#endif