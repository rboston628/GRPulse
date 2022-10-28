#ifndef STELLARHELPSCPP
#define STELLARHELPSCPP

#include "stellar.h"

double chemical::partial_mean_A(elem i, Abundance const &X){
	return -pow(X.mean_A(),2)*
			(1./chemical::A[i] - 1./chemical::A[X.e]);
}

double chemical::partial_mean_Z(elem i, Abundance const &X){
	return chemical::Z[i] - chemical::Z[X.e];
}

double chemical::partial_mu_e(elem i, Abundance const &X){
	return -pow(X.mu_e(),2)*
			(chemical::Z[i]/chemical::A[i] - chemical::Z[X.e]/chemical::A[X.e]);
}

double chemical::partial_mean_coul(elem i , Abundance const &X){
	return pow(chemical::Z[i],2)*pow(chemical::A[i],-1./3.) 
			- pow(chemical::Z[X.e],2)*pow(chemical::A[X.e],-1./3.);
}

//exponent and logarithm
StellarVar exp(const StellarVar &x){
	StellarVar val;
	for(int j=0; j<num_var; j++) val[j] = std::exp(x[j]);
	return val;
}
StellarVar log(const StellarVar &x){
	StellarVar val;
	for(int j=0; j<num_var; j++) val[j] = std::log(x[j]);
	return val;
}

/*//some opacities
double radiative_opacity(StellarVar ly, Abundance const &X){
	//Schwarzschild 1946, Cox&Giuli, Shapiro & Teuskolsky 1983 eq 4.1.8
	const double kappa0 = 4.34e24*(X.C12+X.O16)*(1.+X.H1); 
	return kappa0*exp(ly[dens]-3.5*ly[temp]);// Shapiro & Teukolsky 1983 eq 4.1.3
}
double conductive_opacity(StellarVar ly, Abundance const &X){
	return 4.e-8*pow(X.mu_e()*X.mean_Z(),2)/X.mean_A()*exp(2.*(ly[temp]-ly[dens]));
}*/


int EOS::rando(39);
double EOS::U(double rho, double T, Abundance const &chem){
	double U = 0.0;
	for(PartialPressure p : pressure) U += p.U(rho,T,chem);
	return U;
}
double EOS::Gamma1(double rho, double T, Abundance const &chem){
	double P=0.0, U=0.0, dPdT = 0.0, dUdT = 0.0, dPdr = 0.0, dUdr = 0.0;
	for(PartialPressure p : pressure){
		P += p(rho,T,chem);
		U += p.U(rho,T,chem);
		dPdT += p.partialT(rho,T,chem);
		dUdT += p.UpartialT(rho,T,chem);
		dPdr += p.partialRho(rho,T,chem);
		dUdr += p.UpartialRho(rho,T,chem);
	}
	return rho/P*dPdr + U/P*dPdT/dUdT*(1. + P/U - rho/U*dUdr);
}
double EOS::Gamma3(double rho, double T, Abundance const &chem){
	double dPdT = 0.0, dUdT = 0.0;
	for(PartialPressure p : pressure) dPdT += p.partialT(rho,T,chem);
	for(PartialPressure p : pressure) dUdT += p.UpartialT(rho,T,chem);
	return dPdT/dUdT + 1.;
}
double EOS::nabla_ad(double rho, double T, Abundance const &chem){
	return (Gamma3(rho,T,chem)-1.)/Gamma1(rho,T,chem);
}
double EOS::chiRho(double rho, double T, Abundance const &chem){
	double P=0.0, dPdr = 0.0;
	for(PartialPressure p : pressure) {
		P += p(rho,T,chem);
		dPdr += p.partialRho(rho,T,chem);
	}
	return rho/P*dPdr;
}
double EOS::chiT(double rho, double T, Abundance const &chem){
	double P=0.0, dPdT = 0.0;
	for(PartialPressure p : pressure) {
		P += p(rho,T,chem);
		dPdT += p.partialT(rho,T,chem);
	}
	return T/P*dPdT;
}
double EOS::chiY(chemical::elem i, double rho, double T, Abundance const &chem){
	double P=0.0, dPdY = 0.0;
	for(PartialPressure p : pressure) {
		P += p(rho,T,chem);
		dPdY += p.partialX(i, rho,T,chem);
	}
	return chem[i]/P*dPdY;
}
//calculate the Ledoux term -- see Brassard etc 1991, eqn (13, 14)
// expanding out the logarithmic derivatives
// B = -1/chiT*Sum[ (delP/delX)_rho,T * dX/dP ]
double EOS::Ledoux(double rho, double T, Abundance const &chem, Abundance const &dchem){
	double P=0.0, dPdr = 0.0, dPdT = 0.0, dPdY = 0.0, U=0.0, dUdr = 0.0, dUdT = 0.0;
	double chiT, chiRho, chiY;
	for(PartialPressure p : pressure){
		P += p(rho,T,chem);
		dPdT += p.partialT(rho,T,chem);
		dPdY += p.partialX(chemical::h, rho,T,chem)*dchem[chemical::h];
		dPdY += p.partialX(chemical::he,rho,T,chem)*dchem[chemical::he];
		dPdY += p.partialX(chemical::c, rho,T,chem)*dchem[chemical::c];
		dPdY += p.partialX(chemical::o, rho,T,chem)*dchem[chemical::o];
	}
	chiT = T/P*dPdT;
	return -dPdY/chiT;
}

double EOS::partialRho(double rho, double T, Abundance const &chem){
	double dPdr = 0.0;
	for(PartialPressure p : pressure) dPdr += p.partialRho(rho,T,chem);
	return dPdr;
}
	
double EOS::partialT(double rho, double T, Abundance const &chem){
	double dPdT = 0.0;
	for(PartialPressure p : pressure) dPdT += p.partialT(rho,T,chem);
	return dPdT;
}


//retrieve a particular partial-pressure term
PartialPressure EOS::operator[](int n){
	return pressure[n];
}

// Invert the usual EOS P = P(rho,T,X), to instead solve for rho given P,T,X
// uses a combination of Newton method and bisection search
double EOS::invert(double rho_last, double P, double T, Abundance const &chem){
	double x1=rho_last, x2 = rho_last, x, x3, dx;
	double y1 = P - (*this)(x1,T,chem), y2 = y1, y, y3;
	y2=y1;
	//if the two are on same side adjust x1 until y1 is on other side
	if(y1*y2 > 0.0){
		x3 = 1.01*x2;
		y3 = P - (*this)(x3,T,chem);
		dx = -y2*(x3-x2)/(y3-y2);
		if(!isnan(dx)) dx = 0.0;
		x1 = x2 + dx;
		if(x1 < 0.0) x1 = 0.0;
		y1 = P - (*this)(x1,T,chem);
	}
	//now keep adjusting x1
	int term = 0;
	while(y1*y2 > 0.0){	
		//compute numerical derivative
		x3 = 1.01*x1;					//use x3 = x1 + dx
		y3 = P - (*this)(x3,T,chem);	//find y3 = f(x3) = f(x1+dx)
		dx = -y1*(x3-x1)/(y3-y1);		//this is the needed shift to x1 to approach zero
		if(y3 == y1) dx = 0.0;
		x3 = x1 + dx;					//now update
		if(x3  > 10.0*x1) x3 = x1 + 0.01*dx;
		if(x3  < 0.01*x1) x3 = x1 + 0.01*dx;
		//Newton's method can fail at this if it only approaches the zero from one direction
		//In such a case, slightly broaden the bracket to get to other side of zero
		if( x3 == x1 ) { 				//if the new (x3) is same as old (x1), bracket stagnant
			if( x3 > x2 ) x3 *= 1.01;	//if x3 should be bigger than x2, make even bigger
			if( x3 < x2 ) x3 *= 0.99;	//if x3 should be smaller than x2, make even smaller
		}
		//update and repeat
		x1 = x3;					    //
		y1 = P - (*this)(x1,T,chem);	//
		if(++term>1000) {
			return x1;
		}
	}
	
	//if brackets are backward, swap
	if(x2 < x1) {					//if y2 (the presumed max bracket) is smaller than y1
		double temp = y2;			//save value of y2
		y2 = y1;					//replace y2 with y1
		y1 = temp;					//replace y1 with saved y2
		temp = x2;					//same value of x2
		x2 = x1;					//replace x2 with x1
		x1 = temp;					//replace x1 with samed x2
	}
	//now close brackets
	term = 0;
	while(fabs(x2-x1)/fabs(rho_last)>1.0e-10){
		x = 0.5*(x1+x2);			//pick halfway point
		y = P - (*this)(x,T,chem);	//evalute function at that point
		if(y==y1 | y==y2){//if the problem becomes stuck in a loop
			rando = (a*rando+b)%r; //generates a psuedo-random integer in (0,r)
			x = x1 + (double(rando)/double(r))*fabs(x2-x1);
			y = P - (*this)(x,T,chem);
		}
		if(y*y1 > 0.0) {			//if on same side as min bracket
			if(x > x1){				//and if it is closer than min bracket
				y1 = y;				//replace values
				x1 = x;
			}
		}
		else if(y*y2 > 0.0) {		//if on same side as max bracket
			if(x < x2){				//and if it is closer than max bracket
				y2 = y;				//replace values
				x2 = x;
			}
		}
		//if(term>1000) printf("STUX %le\n", P);
		if(++term>1000) break;
	}
	x = 0.5*(x1+x2);				//pick the halfway point again
	rho_last =x;
	return x;
}

// Invert the usual EOS P = P(rho,T,X), to instead solve for rho given P,T,X
// uses nothing but Newton's method
double EOS::invertNewton(double rho_last, double P, double T, Abundance const &chem){
	double x1=rho_last, x2 = rho_last, x, x3, dx;
	double y1 = P - (*this)(x1,T,chem), y2 = y1, y, y3;
	int term = 0;
	while(fabs(y1) > 1.e-10){	
		//compute numerical derivative
		x3 = 1.01*x1;					//use x3 = x2 + dx
		y3 = P - (*this)(x3,T,chem);	//find y3 = f(x3) = f(x2+dx)
		dx = -y1*(x3-x1)/(y3-y1);		//this is the needed shift to x1 to approach zero
		if(y3 == y1) break;
		x3 = x1 + dx;					//now update
		//limit amount of change allowed in single step
		if(x3 > 1.01*x1) x3 = x1 + 0.01*dx;	//if we increased, don't increase too much
		if(x3 < 0.99*x1) x3 = x1 + 0.01*dx;	//if we decreased, don't decrease too much
		//update and repeat
		x1 = x3;					    //
		y1 = P - (*this)(x1,T,chem);	//
		if(++term>1000) break;
	}
	rho_last = x1;
	return x1;
}

//  Adding new terms to the EOS will add them to the vector of function pointers
void EOS::push_back(PartialPressure p){
	pressure.push_back(p);
}


//************************************************************************************
//	PARTIAL PRESSURE FUNCTIONS for EQUATION OF STATE
//		Each must take arguments double rho, double T, Abundance const &X
//		Output must be a double, the pressure 
//		For references, see:
//		*	Chandrasekhar 1939
//		*	Cox & Giulu 1980
//		*	Pichon 1989
//		*	Tassoul, Fontaine, Winget, 1990
//		*	Van Horn 1968
//		*	Koester 1976
//************************************************************************************
PartialPressure ideal = {
	pressure_ideal, partialRho_ideal, partialT_ideal, partialX_ideal,
	energy_ideal, UpartialRho_ideal, UpartialT_ideal, UpartialX_ideal};
//pressure of an ideal gas
double pressure_ideal(double rho, double T, Abundance const &X){
	return N_Avogadro*boltzmann_k/X.mean_A()*rho*T;
}
double partialRho_ideal(double rho, double T, Abundance const &X){
	return N_Avogadro*boltzmann_k/X.mean_A()*T;
}
double partialT_ideal(double rho, double T, Abundance const &X){
	return N_Avogadro*boltzmann_k/X.mean_A()*rho;
}
double partialX_ideal(chemical::elem i, double rho, double T, Abundance const &X){
	return -N_Avogadro*boltzmann_k/pow(X.mean_A(),2)*rho*T*chemical::partial_mean_A(i,X);
}
double energy_ideal(double rho, double T, Abundance const &X){
	return 1.5*pressure_ideal(rho,T,X);
}
double UpartialRho_ideal(double rho, double T, Abundance const &X){
	return 1.5*partialRho_ideal(rho,T,X);
}
double UpartialT_ideal(double rho, double T, Abundance const &X){
	return 1.5*partialT_ideal(rho,T,X);
}
double UpartialX_ideal(chemical::elem i, double rho, double T, Abundance const &X){
	return 1.5*partialX_ideal(i, rho,T,X);
}

//RADIATION GAS
//standard radiation pressure
PartialPressure rad_gas = {
	pressure_rad, partialRho_rad, partialT_rad, partialX_rad,
	energy_rad, UpartialRho_rad, UpartialT_rad, UpartialX_rad};
double pressure_rad(double rho, double T, Abundance const &X){
	return radiation_a/3.*pow(T,4);
}
double partialRho_rad(double rho, double T, Abundance const &X){
	return 0.0;
}
double partialT_rad(double rho, double T, Abundance const &X){
	return 4./3.*radiation_a*pow(T,3);
}
double partialX_rad(chemical::elem i, double rho, double T, Abundance const &X){
	return 0.0;
}
double energy_rad(double rho, double T, Abundance const &X){
	return radiation_a*pow(T,4);
}
double UpartialRho_rad(double rho, double T, Abundance const &X){
	return 0.0;
}
double UpartialT_rad(double rho, double T, Abundance const &X){
	return 4.*radiation_a*pow(T,3);
}
double UpartialX_rad(chemical::elem i, double rho, double T, Abundance const &X){
	return 0.0;
}


//COULOMB PRESSURE
//pressure contribution of ideal gas due to coloumb  correction
//  following EOS of Tassoul, Fontaine, Winget 1990
//  see also Van Horn 1968, Koester 1976
//  According to Huebner 2014, eq 9.31, GC should have a power of 4 and not a 5...
//  According to Cassisi et al 2007, eq 7, power should be 5
PartialPressure coul = {
	pressure_coul, partialRho_coul, partialT_coul, partialX_coul,
	energy_coul, UpartialRho_coul, UpartialT_coul, UpartialX_coul};
double pressure_coul(double rho, double T, Abundance const &X){
	double GC = 2.275e5*(X.mean_coulomb())*pow(rho,1./3.)/T;
	if(T<=0.0) GC = 0.0;
	double coul = -0.3*pow(GC,1.5)/(pow(GC,0.5)+1.03921);
	return coul*pressure_ideal(rho,T,X);
}
double partialRho_coul(double rho, double T, Abundance const &X){
	double GC = 2.275e5*(X.mean_coulomb())*pow(rho,1./3.)/T;
	if(T<=0.0) GC = 0.0;
	double denom = (pow(GC,0.5)+1.03921);
	double coul = -0.3*pow(GC,1.5)/denom;
	double dGCdr = 2.275e5*(X.mean_coulomb())*(1./3.)*pow(rho,-2./3.)/T;
	return coul*partialRho_ideal(rho,T,X)
			 + pressure_ideal(rho,T,X)*( -0.45*pow(GC,0.5)/denom + 0.15*GC/denom/denom )*dGCdr;
}
double partialT_coul(double rho, double T, Abundance const &X){
	double GC = 2.275e5*(X.mean_coulomb())*pow(rho,1./3.)/T;
	if(T<=0.0) GC = 0.0;
	double denom = (pow(GC,0.5)+1.03921);
	double coul = -0.3*pow(GC,1.5)/denom;
	double dGCdt = -2.275e5*(X.mean_coulomb())*pow(rho,1./3.)/T/T;
	return coul*partialT_ideal(rho,T,X)
			 + pressure_ideal(rho,T,X)*( -0.45*pow(GC,0.5)/denom + 0.15*GC/denom/denom  )*dGCdt;
}
double partialX_coul(chemical::elem i, double rho, double T, Abundance const &X){
	double GC = 2.275e5*(X.mean_coulomb())*pow(rho,1./3.)/T;
	if(T<=0.0) GC = 0.0;
	double denom = (pow(GC,1./2.)+1.03921);
	double coul = -0.3*pow(GC,3./2.)/denom;
	double dGCdX = 2.275e5*pow(rho,1./3.)/T*chemical::partial_mean_coul(i, X);
	return coul*partialX_ideal(i, rho,T,X)
			 + pressure_ideal(rho,T,X)*( -0.45*pow(GC,0.5)/denom + 0.15*GC/denom/denom  )*dGCdX;
}
double energy_coul(double rho, double T, Abundance const &X){
	return 3.*pressure_coul(rho,T,X);
}
double UpartialRho_coul(double rho, double T, Abundance const &X){
	return 3.*partialRho_coul(rho,T,X);
}
double UpartialT_coul(double rho, double T, Abundance const &X){
	return 3.*partialT_coul(rho,T,X);
}
double UpartialX_coul(chemical::elem i, double rho, double T, Abundance const &X){
	return 3.*partialX_coul(i, rho,T,X);
}


//ELECTRON DEGENERACY PRESSURE -- ZERO TEMP
//pressure of completely degenerate electron gas (T=0)
//	following Chandrasekhar 1939
PartialPressure deg_zero = {
	pressure_deg_zero, partialRho_deg_zero, partialT_deg_zero, partialX_deg_zero,
	energy_deg_zero, UpartialRho_deg_zero, UpartialT_deg_zero, UpartialX_deg_zero};
double pressure_deg_zero(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	double X0 = pow(rho/B0/X.mu_e(), 1./3.);
	return A0*factor_f(X0);
}
double partialRho_deg_zero(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	double X0 = pow(rho/B0/X.mu_e(), 1./3.);
	return A0*8.*pow(X0,4)/sqrt(1.+X0*X0) * (1./(3.*B0*X.mu_e()*X0*X0));
}
double partialT_deg_zero(double rho, double T, Abundance const &X){
	return 0.0;
}
double partialX_deg_zero(chemical::elem i, double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	double X0 = pow(rho/B0/X.mu_e(), 1./3.);
	return A0*8.*pow(X0,4)/sqrt(1.+X0*X0) 
		* pow(rho/B0,1./3.)*(-1./3.*pow(X.mu_e(),-4./3.))*chemical::partial_mu_e(i,X);
}
double energy_deg_zero(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	double X0 = pow(rho/B0/X.mu_e(), 1./3.);
	return A0*factor_g(X0);
}
double UpartialRho_deg_zero(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	double X0 = pow(rho/B0/X.mu_e(), 1./3.);
	return A0*24.*X0*X0*(sqrt(1.+X0*X0)-1.) * (1./(3.*B0*X.mu_e()*X0*X0));
}
double UpartialT_deg_zero(double rho, double T, Abundance const &X){
	return 0.0;
}
double UpartialX_deg_zero(chemical::elem i, double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	double X0 = pow(rho/B0/X.mu_e(), 1./3.);
	return A0*24.*X0*X0*(sqrt(1.+X0*X0)-1.)
		* pow(rho/B0,1./3.)*(-1./3.*pow(X.mu_e(),-4./3.))*chemical::partial_mu_e(i,X);
}


//ELECTRON DEGENERACY PRESSURE -- FINITE TEMP
PartialPressure deg_finite = {
	pressure_deg_finite, partialRho_deg_finite, partialT_deg_finite, partialX_deg_finite,
	energy_deg_finite, UpartialRho_deg_finite, UpartialT_deg_finite, UpartialX_deg_finite};
double findEta(double rho, double beta, double mue){
	using namespace Chandrasekhar;
	using namespace FermiDirac;
	//approximate eta as in Cox and Giuli eq 24.29', pg 793 (but using 24.327b)
	static const double eta_coeff = pow(8.*B0*B0*mue*mue, -1./3.);
	static const double rho_coeff = 3.*sqrt(2.)*B0*mue*pow(beta, 1.5);
	double eta1 = eta_coeff*pow(rho, 2./3.)/beta;
	//  find eta by inverting rho = rho0*( F1/2 + beta*F3/2)
	double rho1 = rho_coeff*(FermiDirac1Half(eta1,beta) + beta*FermiDirac3Half(eta1,beta));
	double eta2 = 1.01*eta1;
	double rho2 = rho_coeff*(FermiDirac1Half(eta2,beta) + beta*FermiDirac3Half(eta2,beta));
	int halt=1000, count=0;
	while((rho2-rho)/rho>1.0e-12){
		eta2 = eta2 + (rho-rho2)*(eta2-eta1)/(rho2-rho1);
		rho2 = rho_coeff*(FermiDirac1Half(eta2,beta) + beta*FermiDirac3Half(eta2,beta));
		if(++count>halt) break;
	}
	return eta2;
}
//pressure of electron gas under partial degeneracy (T finite)
//	following Chandrasekhar 1939, Cox & Giuli 1980 (especially pg 851)
//  the Fermi-Dirac functions calculated with 24-point quadrature due to Sagar 1991a&b
double pressure_deg_finite(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;	
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	double mue = X.mu_e();
	//approximate eta as in Cox and Giuli eq 24.29', pg 793 (but using 24.327b)
	static const double eta_coeff = pow(8.*B0*B0*mue*mue, -1./3.);
	static const double rho_coeff = 3.*sqrt(2.)*B0*mue*pow(beta, 1.5);
	double eta1 = eta_coeff*pow(rho, 2./3.)/beta;
	//  find eta by inverting rho = rho0*( F1/2 + beta*F3/2)
	double rho1 = rho_coeff*(FermiDirac1Half(eta1,beta) + beta*FermiDirac3Half(eta1,beta));
	double eta2 = 1.01*eta1;
	double rho2 = rho_coeff*(FermiDirac1Half(eta2,beta) + beta*FermiDirac3Half(eta2,beta));
	int halt=10000, count=0;
	while((rho2-rho)/rho>1.0e-12){
		eta2 = eta2 + (rho-rho2)*(eta2-eta1)/(rho2-rho1);
		rho2 = rho_coeff*(FermiDirac1Half(eta2,beta) + beta*FermiDirac3Half(eta2,beta));
		if(++count>halt) break;
	}
	//return the pressure, as in Cox and Giuli eq 24.326b, pg 850
	double P =  16.*sqrt(2.)*A0*pow(beta,2.5)*(
		FermiDirac3Half(eta2,beta) + 0.5*beta*FermiDirac5Half(eta2,beta));
	return P;
}
double partialRho_deg_finite(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;	
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	double mue = X.mu_e();
	double eta = findEta(rho, beta, mue);
	//return the partial pressure derivative
	double dPdrho =  16.*A0/(3.*B0)*beta/X.mu_e()*(
		FermiDirac3HalfdelEta(eta,beta) + 0.5*beta*FermiDirac5HalfdelEta(eta,beta))/(
		FermiDirac1HalfdelEta(eta,beta) +     beta*FermiDirac3HalfdelEta(eta,beta));
	return dPdrho;
}
double partialT_deg_finite(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;	
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	double mue = X.mu_e();
	double eta = findEta(rho, beta, mue);
	//return the partial pressure derivative
	double dPdT =  16.*sqrt(2.)*A0*(beta/T)*(
		2.5*pow(beta,1.5)*(	FermiDirac3Half(eta,beta) + 0.5*beta*FermiDirac5Half(eta,beta) )
		+ 0.5*pow(beta,2.5)*FermiDirac5Half(eta,beta)
		+ pow(beta,2.5)*(FermiDirac3HalfdelBeta(eta,beta)+0.5*beta*FermiDirac5HalfdelBeta(eta,beta))
	);
	return dPdT;
}
double partialX_deg_finite(chemical::elem i, double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;	
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	double mue = X.mu_e();
	double eta = findEta(rho, beta, mue);
	//return the partial pressure derivative
	double dPdX =  -16.*sqrt(2.)*A0*pow(beta,2.5)
		*(FermiDirac1Half(eta,beta)+beta*FermiDirac3Half(eta,beta))
		*(FermiDirac3HalfdelEta(eta,beta)+0.5*beta*FermiDirac5HalfdelEta(eta,beta))
		/(FermiDirac1HalfdelEta(eta,beta)+    beta*FermiDirac3HalfdelEta(eta,beta))
		*chemical::partial_mu_e(i,X)/mue;
	return dPdX;
}
double energy_deg_finite(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;	
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	double mue = X.mu_e();
	double eta = findEta(rho, beta, mue);
	//return the pressure, as in Cox and Giuli eq 24.326b, pg 850
	double U =  24.*sqrt(2.)*A0*pow(beta,2.5)*(
		FermiDirac3Half(eta,beta) + beta*FermiDirac5Half(eta,beta));
	return U;
}
double UpartialRho_deg_finite(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;	
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	double mue = X.mu_e();
	double eta = findEta(rho, beta, mue);
	//return the partial pressure derivative
	double dUdrho =  24.*A0/(3.*B0)*beta/X.mu_e()*(
		FermiDirac3HalfdelEta(eta,beta) + beta*FermiDirac5HalfdelEta(eta,beta))/(
		FermiDirac1HalfdelEta(eta,beta) + beta*FermiDirac3HalfdelEta(eta,beta));
	return dUdrho;
}
double UpartialT_deg_finite(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;	
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	double mue = X.mu_e();
	double eta = findEta(rho, beta, mue);
	//return the partial pressure derivative
	double dUdT =  24.*sqrt(2.)*A0*(beta/T)*(
		2.5*pow(beta,1.5)*(	FermiDirac3Half(eta,beta) + beta*FermiDirac5Half(eta,beta) )
		+ pow(beta,2.5)*FermiDirac5Half(eta,beta)
		+ pow(beta,2.5)*(FermiDirac3HalfdelBeta(eta,beta)+beta*FermiDirac5HalfdelBeta(eta,beta))
	);
	return dUdT;
}
double UpartialX_deg_finite(chemical::elem i, double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;	
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	double mue = X.mu_e();
	double eta = findEta(rho, beta, mue);
	//return the partial pressure derivative
	double dUdX =  -24.*sqrt(2.)*A0*pow(beta,2.5)
		*(FermiDirac1Half(eta,beta) + beta*FermiDirac3Half(eta,beta))
		*(FermiDirac3HalfdelEta(eta,beta) + beta*FermiDirac5HalfdelEta(eta,beta))
		/(FermiDirac1HalfdelEta(eta,beta) + beta*FermiDirac3HalfdelEta(eta,beta))
		*chemical::partial_mu_e(i,X)/mue;
	return dUdX;
}


//ELECTRON DEGENERACY PRESSURE -- FINITE TEMP
//  find eta by inverting rho = rho0*( F1/2 + beta*F3/2)
PartialPressure deg_trap = {
	pressure_deg_trap,  partialRho_deg_finite,  partialT_deg_finite,  partialX_deg_finite,
	  energy_deg_trap, UpartialRho_deg_finite, UpartialT_deg_finite, UpartialX_deg_finite};
double density_trap(double eta, double beta, double mue){
	using namespace Chandrasekhar;
	using namespace FermiDirac;
	static const double rho_coeff = 3.*sqrt(2.)*B0*mue*pow(beta, 1.5);
	return rho_coeff*(FermiDirac::FermiDirac(0.5,eta,beta) + beta*FermiDirac::FermiDirac(1.5,eta,beta));
}
double findEta_trap(double rho, double beta, double mue){
	using namespace Chandrasekhar;
	using namespace FermiDirac;
	//approximate eta as in Cox and Giuli eq 24.29', pg 793 (but using 24.327b)
	static const double eta_coeff = pow(8.*B0*B0*mue*mue, -1./3.);
	static const double rho_coeff = 3.*sqrt(2.)*B0*mue*pow(beta, 1.5);
	double x1 = eta_coeff*pow(rho, 2./3.)/beta;
	double y1 = rho-density_trap(x1, beta, mue);
	double x3 = 1.01*x1;
	double y3 = rho-density_trap(x3,beta,mue);
	double x2 = x1 - y1*(x3-x1)/(y3-y1);
	double y2 = rho-density_trap(x2, beta, mue);
	
	double x, dx;
	double y, dydx;
	//now keep adjusting x1
	int term = 0;
	while(y1*y2 > 0.0){	
		//compute numerical derivative
		x3 = 1.01*x1;					//use x3 = x1 + dx
		y3 = rho - density_trap(x3,beta,mue);	//find y3 = f(x3) = f(x1+dx)
		dx = -y1*(x3-x1)/(y3-y1);		//this is the needed shift to x1 to approach zero
		if(y3 == y1) dx = 0.0;
		x3 = x1 + dx;
		//if(x3  > 10.0*x1) x3 = x1 + 0.001*dx;
		//if(x3  < 0.01*x1) x3 = x1 + 0.001*dx;
		//Newton's method can fail at this if it only approaches the zero from one direction
		//In such a case, slightly broaden the bracket to get to other side of zero
		if( x3 == x1 ) { 				//if the new (x3) is same as old (x1), bracket stagnant
			if( x3 > x2 ) x3 *= 1.01;	//if x3 should be bigger than x2, make even bigger
			if( x3 < x2 ) x3 *= 0.99;	//if x3 should be smaller than x2, make even smaller
		}
		//update and repeat
		x1 = x3;					    //
		y1 = rho - density_trap(x1, beta, mue);//
		if(++term>1000) {
			return x1;
		}
	}
	//if brackets are backward, swap
	if(x2 < x1) {					//if y2 (the presumed max bracket) is smaller than y1
		double temp = y2;			//save value of y2
		y2 = y1;					//replace y2 with y1
		y1 = temp;					//replace y1 with saved y2
		temp = x2;					//same value of x2
		x2 = x1;					//replace x2 with x1
		x1 = temp;					//replace x1 with samed x2
	}
	//now close brackets
	x = 0.5*(x1+x2);
	term = 0;
	while(fabs(x2-x1)/fabs(x)>1.0e-14){
		x = 0.5*(x1+x2);			//pick halfway point
		y = rho - density_trap(x, beta, mue);	//evalute function at that point
		if(y*y1 > 0.0) {			//if on same side as min bracket
			if(x > x1){				//and if it is closer than min bracket
				y1 = y;				//replace values
				x1 = x;
			}
		}
		else if(y*y2 > 0.0) {		//if on same side as max bracket
			if(x < x2){				//and if it is closer than max bracket
				y2 = y;				//replace values
				x2 = x;
			}
		}
		if(++term>1000) break;
	}
	x = 0.5*(x1+x2);				//pick the halfway point again
	return x;
}
//pressure of electron gas under partial degeneracy (T finite)
//	following Chandrasekhar 1939, Cox & Giuli 1980 (especially pg 851)
//  the Fermi-Dirac functions calculated with trapezoidal integration
double pressure_deg_trap(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;	
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	double mue = X.mu_e();
	//approximate eta as in Cox and Giuli eq 24.29', pg 793 (but using 24.327b)
	double eta = findEta_trap(rho, beta, mue);
	//return the pressure, as in Cox and Giuli eq 24.326b, pg 850
	double P =  16.*sqrt(2.)*A0*pow(beta,2.5)*(
		FermiDirac::FermiDirac(1.5,eta,beta) + 0.5*beta*FermiDirac::FermiDirac(2.5,eta,beta));
	return P;
}
double energy_deg_trap(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;	
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	double mue = X.mu_e();
	double eta = findEta_trap(rho, beta, mue);
	//return the pressure, as in Cox and Giuli eq 24.326b, pg 850
	double U =  24.*sqrt(2.)*A0*pow(beta,2.5)*(
		FermiDirac::FermiDirac(1.5,eta,beta) + beta*FermiDirac::FermiDirac(2.5,eta,beta));
	return U;
}



//ELECTRON DEGENERACY PRESSURE -- PARTIAL DEGENERACY
//pressure of electron gas under partial degeneracy (T finite)
//  uses a combination of methods depending on degeneracy
//	following Cox and Giuli 1980 -- see equation 24.51, the definition of w
PartialPressure deg_partial = {
	pressure_deg_partial, partialRho_deg_partial, partialT_deg_partial, partialX_deg_partial,
	energy_deg_partial, UpartialRho_deg_partial, UpartialT_deg_partial, UpartialX_deg_partial};
double pressure_deg_partial(double rho, double T, Abundance const &X){
	static const double 
		w_coeff = 0.5*pow(planck_h_CGS,3)/pow(2.*m_pi*electron.mass_CGS*boltzmann_k,1.5);
	double ne = rho/proton.mass_CGS/X.mu_e();
	double w = w_coeff*ne*pow(T, -1.5);
	if      (w >10.0) return pressure_deg_zero(rho,T,X);
	//else if (w > 20.0) return pressure_deg_large(rho,T,X);
	else if (w >= 0.1) return pressure_deg_finite(rho,T,X);
	else if (w <  0.1) return ne*boltzmann_k*T;
	else return 0.0;
}
double partialRho_deg_partial(double rho, double T, Abundance const &X){
	static const double 
		w_coeff = 0.5*pow(planck_h_CGS,3)/pow(2.*m_pi*electron.mass_CGS*boltzmann_k,1.5);
	double ne = rho/proton.mass_CGS/X.mu_e();
	double w = w_coeff*ne*pow(T, -1.5);
	if      (w >10.0) return partialRho_deg_zero(rho,T,X);
	//else if (w > 20.0) return pressure_deg_large(rho,T,X);
	else if (w >= 0.1) return partialRho_deg_finite(rho,T,X);
	else if (w <  0.1) return boltzmann_k*T/proton.mass_CGS/X.mu_e();
	else return 0.0;
}
double partialT_deg_partial(double rho, double T, Abundance const &X){
	static const double 
		w_coeff = 0.5*pow(planck_h_CGS,3)/pow(2.*m_pi*electron.mass_CGS*boltzmann_k,1.5);
	double ne = rho/proton.mass_CGS/X.mu_e();
	double w = w_coeff*ne*pow(T, -1.5);
	if      (w >10.0) return partialT_deg_zero(rho,T,X);
	//else if (w > 20.0) return pressure_deg_large(rho,T,X);
	else if (w >= 0.1) return partialT_deg_finite(rho,T,X);
	else if (w <  0.1) return ne*boltzmann_k;
	else return 0.0;
}
double partialX_deg_partial(chemical::elem i, double rho, double T, Abundance const &X){
	static const double 
		w_coeff = 0.5*pow(planck_h_CGS,3)/pow(2.*m_pi*electron.mass_CGS*boltzmann_k,1.5);
	double ne = rho/proton.mass_CGS/X.mu_e();
	double w = w_coeff*ne*pow(T, -1.5);
	if      (w >10.0) return partialX_deg_zero(i,rho,T,X);
	//else if (w > 20.0) return pressure_deg_large(rho,T,X);
	else if (w >= 0.1) return partialX_deg_finite(i,rho,T,X);
	else if (w <  0.1) return rho*boltzmann_k*T*pow(X.mu_e(),-2)*chemical::partial_mu_e(i,X);
	else return 0.0;
}
double energy_deg_partial(double rho, double T, Abundance const &X){
	static const double 
		w_coeff = 0.5*pow(planck_h_CGS,3)/pow(2.*m_pi*electron.mass_CGS*boltzmann_k,1.5);
	double ne = rho/proton.mass_CGS/X.mu_e();
	double w = w_coeff*ne*pow(T, -1.5);
	if      (w >10.0) return energy_deg_zero(rho,T,X);
	//else if (w > 20.0) return pressure_deg_large(rho,T,X);
	else if (w >= 0.1) return energy_deg_finite(rho,T,X);
	else if (w <  0.1) return 1.5*ne*boltzmann_k*T;
	else return 0.0;
}
double UpartialRho_deg_partial(double rho, double T, Abundance const &X){
	static const double 
		w_coeff = 0.5*pow(planck_h_CGS,3)/pow(2.*m_pi*electron.mass_CGS*boltzmann_k,1.5);
	double ne = rho/proton.mass_CGS/X.mu_e();
	double w = w_coeff*ne*pow(T, -1.5);
	if      (w >10.0) return UpartialRho_deg_zero(rho,T,X);
	//else if (w > 20.0) return pressure_deg_large(rho,T,X);
	else if (w >= 0.1) return UpartialRho_deg_finite(rho,T,X);
	else if (w <  0.1) return 1.5*boltzmann_k*T/proton.mass_CGS/X.mu_e();
	else return 0.0;
}
double UpartialT_deg_partial(double rho, double T, Abundance const &X){
	static const double 
		w_coeff = 0.5*pow(planck_h_CGS,3)/pow(2.*m_pi*electron.mass_CGS*boltzmann_k,1.5);
	double ne = rho/proton.mass_CGS/X.mu_e();
	double w = w_coeff*ne*pow(T, -1.5);
	if      (w >10.0) return UpartialT_deg_zero(rho,T,X);
	//else if (w > 20.0) return pressure_deg_large(rho,T,X);
	else if (w >= 0.1) return UpartialT_deg_finite(rho,T,X);
	else if (w <  0.1) return 1.5*ne*boltzmann_k;
	else return 0.0;
}
double UpartialX_deg_partial(chemical::elem i, double rho, double T, Abundance const &X){
	static const double 
		w_coeff = 0.5*pow(planck_h_CGS,3)/pow(2.*m_pi*electron.mass_CGS*boltzmann_k,1.5);
	double ne = rho/proton.mass_CGS/X.mu_e();
	double w = w_coeff*ne*pow(T, -1.5);
	if      (w >10.0) return UpartialX_deg_zero(i,rho,T,X);
	//else if (w > 20.0) return pressure_deg_large(rho,T,X);
	else if (w >= 0.1) return UpartialX_deg_finite(i,rho,T,X);
	else if (w <  0.1) return 1.5*rho*boltzmann_k*T*pow(X.mu_e(),-2)*chemical::partial_mu_e(i,X);
	else return 0.0;
}


//pressure of electron gas, where degeneracy can be considered large
// following Pichon 1989, eq 27-29, based on Sommerfeld approximation
// can also be seen in Cox & Giuli, sec 24.2
// *  this partial pressure is less stable, due to high inaccuracies for very high eta
// *  does not smoothly limit to complete degeneracy case
double pressure_deg_large(double rho, double T, Abundance const &X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;
	double beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	//approximate eta as in Cox and Giuli eq 24.29', pg 793 (but using 24.327b)
	static const double eta_coeff = pow(8.*B0*B0, -1./3.);
	static const double rho_coeff = 3.*sqrt(2.)*B0*pow(beta,1.5);
	double eta1 = eta_coeff*pow(rho, 2./3.)/beta;
	//  finds eta by inverting rho = rho0*( F1/2(eta) + beta*F3/2(eta))
	double rho1 = rho_coeff*(FermiDiracEtaLarge1Half(eta1,beta)
					 + beta*FermiDiracEtaLarge3Half(eta1,beta));
	double eta2 = 1.01*eta1;
	double rho2 = rho_coeff*(FermiDiracEtaLarge1Half(eta2,beta)
					 + beta*FermiDiracEtaLarge3Half(eta2,beta));
	int halt=10000, count=0;
	while((rho2-rho)/rho>1.0e-12){
		eta2 = eta2 + (rho-rho2)*(eta2-eta1)/(rho2-rho1);
		rho2 = rho_coeff*(FermiDiracEtaLarge1Half(eta2,beta)
					 + beta*FermiDiracEtaLarge3Half(eta2,beta));
		if(++count>halt) break;
	}	
	//return the pressure, as in Cox and Giuli eq 24.326b, pg 850
	double P =  16.*sqrt(2.)*A0*pow(beta,2.5)*(
		FermiDiracEtaLarge3Half(eta2,beta) + 0.5*beta*FermiDiracEtaLarge5Half(eta2,beta));
	return P;
}

//pressure of electon gas with small degeneracy...
//  approximated as Pdeg if larger than Pideal, ignored otherwise
//  not a great method, really, but combined with deg_zero does lead to decent models
double pressure_deg_switch(double rho, double T, Abundance const &X){
	double eta, X0;
	double Pdeg = pressure_deg_zero(rho,T,X);
	double Pideal = pressure_ideal(rho,T,X);
	if(Pdeg > Pideal) return Pdeg;
	else return 0.0; 
}




// ****************************************************************************
//  NEUTRON STAR PRESSURES
//  These pressures are those used by Richardson et al (1982)
//    or by Gudmundsson et al (1983) to construct the neutron star models
//    used by McDermott et al (1988) and subsequently by Yoshia & Lee (2001)
// ****************************************************************************

// the core -- composed of protons, electrons, neutrons
//   1. zero-temperature component from Baym, Bethe, Pethic (1971)
//   2. electrons are degenerate and relativistic
//   3. neutrons  are degenerate and nonrelativistic
//   4. protons   are degenerate and nonrelativistic
// both protons and neutrons can be superfluid

// the inner crust -- free electrons, free neutrons, and heavy ions
//   1. zero-temperature component from Negele (1974), Negel & Vautherin (1973)
//   2. electrons are degenerate and relativistic
//   3. neutrons  are degenerate and nonrelativistic
//   4. ions      are nonrelativistic in BCC lattice
// over part of density, neutrons can be superfluid

//the outer crust -- fully ionized plasma
// is this the same as the "ocean"?

// the ocean -- liquid heavy ions of catalyzed matter
// based on Gudmundsson et al 1983, used ub McDermott et al 1988
//   1. catalyzed natter (iron or heavier nuclei)
//   2. electrons are degenerate and relativistic
//   3. nuclei follow ideal gas law plus Coulomb repulsion
//   4. why not radiaiton pressure?


/*//NEUTRON or PROTON DEGENERACY PRESSURE
PartialPressure press_fermi = {
	pressure_fermi, partialRho_fermi, partialT_fermi, partialX_fermi,
	energy_fermi, UpartialRho_fermi, UpartialT_fermi, UpartialX_fermi};
double pressure_fermi(double rho, double T, Abundance const & X);
double partialRho_fermi(double rho, double T, Abundance const & X);
double partialT_fermi(double rho, double T, Abundance const & X);
double partialX_fermi(chemical::elem i, double rho, double T, Abundance const & X);
double energy_fermi(double rho, double T, Abundance const & X);
double UpartialRho_fermi(double rho, double T, Abundance const & X);
double UpartialT_fermi(double rho, double T, Abundance const & X);
double UpartialX_fermi(chemical::elem i, double rho, double T, Abundance const & X);

//COLD CATALYZED MATTER PRESSURE
//  Taken from Baym, Bethe, Pethick (1971) Nuclear Physics vol A175:225-271
double pressure_BBP(double rho, double T, Abundance const & X);
double partialRho_BBP(double rho, double T, Abundance const & X);
double partialT_BBP(double rho, double T, Abundance const & X);
double partialX_BBP(chemical::elem i, double rho, double T, Abundance const & X);
double energy_BBP(double rho, double T, Abundance const & X);
double UpartialRho_BBP(double rho, double T, Abundance const & X);
double UpartialT_BBP(double rho, double T, Abundance const & X);
double UpartialX_BBP(chemical::elem i, double rho, double T, Abundance const & X);
*/





#endif

