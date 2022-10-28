#ifndef STELLARHELPSH
#define STELLARHELPSH

#include "../src/constants.h"
#include "chandra.h"

struct Abundance;

namespace chemical {
	enum elem {h=0, he, c, o};
	const double Z[] = {1., 2.,  6.,  8.};
	const double A[] = {1., 4., 12., 16.};
	
	double partial_mean_A(elem i, Abundance const & X);
	double partial_mean_Z(elem i, Abundance const & X);
	double partial_mu_e(  elem i, Abundance const & X);
	double partial_mean_coul(elem i, Abundance const & X);
}

struct Abundance {
	double H1, He4, C12, O16;
	int e;
	double mean_A() const {
		using namespace chemical;
	//	printf("\tmuA:\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", H1, He4, C12, O16);
		return 1./(H1 + He4/chemical::A[he] + C12/chemical::A[c] + O16/chemical::A[o]);
	}
	double mean_Z() const {
		using namespace chemical;
		return H1 + He4*chemical::Z[he] + C12*chemical::Z[c] + O16*chemical::Z[o];
	}
	double mu_e(double y=1.0) const {
		using namespace chemical;
		return 1./(H1*y + 0.5*(He4+C12+O16)*y);
	}
	double mean_weight(double y=1.0) const {
		using namespace chemical;
		return 1./(1./mu_e(y) + 1./mean_A());
	}
	double mean_coulomb() const{
		using namespace chemical;
		return H1 
			+ He4*pow(chemical::Z[he],2)*pow(chemical::A[he],-1./3.) 
			+ C12*pow(chemical::Z[ c],2)*pow(chemical::A[ c],-1./3.) 
			+ O16*pow(chemical::Z[ o],2)*pow(chemical::A[ o],-1./3.) ;
	}
	void print() const {
		printf("%0.2f\t%0.2f\t%0.2f\t%0.2f\n", H1, He4, C12, O16);
	}
	Abundance() :H1(0), He4(0), C12(0), O16(0), e(0) {};
	Abundance(double x1, double x2, double x3, double x4) 
		: H1(x1),He4(x2),C12(x3),O16(x4), e(0){};
	
	Abundance operator*(const double& x) {
		return Abundance(H1*x, He4*x, C12*x, O16*x);
	}
	Abundance operator+(const Abundance &x){
		return Abundance(H1+x.H1, He4+x.He4, C12+x.C12, O16+x.O16);
	}
	Abundance operator-(const Abundance &x){
		return Abundance(H1-x.H1, He4-x.He4, C12-x.C12, O16-x.O16);
	}
	void operator=(const Abundance &x){
		H1 = x.H1; He4 = x.He4; C12 = x.C12; O16 = x.O16;
		e = x.e;
	}
	double  operator[](int n) const {
		switch(n){
			case 0: return H1;  break;
			case 1: return He4; break;
			case 2: return C12; break;
			case 3: return O16; break;
			default: return 0.0;
		}
	}
	double& operator[](int n){
		switch(n){
			case 0: return H1;  break;
			case 1: return He4; break;
			case 2: return C12; break;
			case 3: return O16; break;
			default: return nothing; 
		}
	}
	void enforce(){
		if(H1 <0.0) H1  = 0.0;
		if(He4<0.0) He4 = 0.0;
		if(C12<0.0) C12 = 0.0;
		if(O16<0.0) O16 = 0.0;
		double tot = H1+He4+C12+O16;
		H1 /= tot;
		He4/= tot;
		C12/= tot;
		O16/= tot;
	}
	private:
		double nothing;
};

enum VarName {dens=0, radi, pres, mass, temp, lumi, num_var};

class StellarVar {
private:
	double var[num_var];
public:
	//access the elements
	double& operator[](const int n){ return (var[n]);};
	double  operator[](const int n) const { return (var[n]);};
	double& operator[](const VarName elem){ return (var[elem]);};
	double  operator[](const VarName elem) const { return var[elem];};
	//constructors
	StellarVar(const double x[num_var]){
		for(int j=0; j<num_var; j++) var[j] = x[j];
	};
	StellarVar(const double& x1, const double& x2, const double& x3, const double& x4, const double& x5, const double& x6){
		var[dens]=x1; var[radi]=x2; var[pres]=x3; var[mass]=x4; var[temp]=x5; var[lumi]=x6;
	};
	StellarVar(const double& x2, const double& x3, const double& x4, const double& x5, const double& x6){
		var[dens]=0.; var[radi]=x2; var[pres]=x3; var[mass]=x4; var[temp]=x5; var[lumi]=x6;
	}
	StellarVar() {
		for(int j=0; j<num_var; j++) var[j] = 0.0;\
	};
	
	//mathematical operations
	StellarVar operator*(const double &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j]*x;
		return val;
	}
	StellarVar operator/(const double &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j]/x;
		return val;
	}
	StellarVar operator+(const StellarVar &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j] + x[j];
		return val;
	}
	StellarVar operator-(const StellarVar &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j] - x[j];
		return val;
	}
	StellarVar operator*(const StellarVar &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j]*x[j];
		return val;
	}
	StellarVar operator/(const StellarVar &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j]/x[j];
		return val;
	}
		
	//assignment
	void operator=(const double x[num_var]){
		for(int j=0; j<num_var; j++) var[j] = x[j];
	}
	void operator=(const StellarVar& x){
		for(int j=0; j<num_var; j++) var[j] = x[j];
	}
};

//exponent and logarithm
StellarVar exp(const StellarVar &x);
StellarVar log(const StellarVar &x);

//some opacities
//double radiative_opacity(StellarVar ly, Abundance X);
//double conductive_opacity(StellarVar ly, Abundance X);


struct PartialPressure {
	typedef double (*funcptr)(double,double,Abundance const &);
	funcptr P, partialRho, partialT;
	double (*partialX)(chemical::elem, double,double,Abundance const &);
	funcptr U, UpartialRho, UpartialT;
	double (*UpartialX)(chemical::elem, double,double,Abundance const &);	
	double operator()(double rho, double T, Abundance const & chem) {return P(rho,T,chem);};
};


class EOS {
public:
	EOS() {rando++;};
	EOS(const std::vector<PartialPressure> p) : pressure(p) {rando++;};
	double invert(double rho_last, double P, double T, Abundance const & chem);
	double invertNewton(double rho_last, double P, double T, Abundance const & chem);
	double operator()(double rho, double T, Abundance const & chem){
		double P = 0.0;
		for(PartialPressure p : pressure) P += p(rho,T,chem);
		return P;
	}
	~EOS(){pressure.clear();std::vector<PartialPressure>().swap(pressure);};	//destructor
	//other stuff
	double U(double rho, double T, Abundance const & chem);
	double Gamma1(double rho, double T, Abundance const & chem);
	double Gamma3(double rho, double T, Abundance const & chem);
	double nabla_ad(double rho, double T, Abundance const & chem);
	double chiRho(double rho, double T, Abundance const & chem);
	double chiT(double rho, double T, Abundance const & chem);
	double chiY(chemical::elem i, double rho, double T, Abundance const & chem);
	//calculate the Ledoux term -- see Brassard etc 1991, eqn (13, 14)
	double Ledoux(double rho, double T, Abundance const & chem, Abundance const & dchem);
	
	//partial derivatives
	double partialRho(double rho, double T, Abundance const & chem);
	double partialT(double rho, double T, Abundance const & chem);
	
	PartialPressure operator[](int n);
	void push_back(PartialPressure);
private:
	//the partial pressures are stored as function pointers in a vector object
	std::vector<PartialPressure> pressure;
	//for pseudorandom positioning -- used in bisection search of "invert"
	//initialize a simple "cat map" PRNG (aka Lehmer PRNG)
	static const int a=53, b=122, r=17737;
	static int rando;
};

//************************************************************************************
//	PARTIAL PRESSURE FUNCTIONS for EQUATION OF STATE
//		Each must take arguments double rho, double T, Abundance X
//		Output must be a double, the pressure 
//		For references, see:
//		*	Chandrasekhar 1939
//		*	Cox & Giulu 1980
//		*	Pichon 1989
//		*	Tassoul, Fontaine, Winget, 1990
//		*	Van Horn 1968
//		*	Koester 1976
//************************************************************************************

//IDEAL GAS
double pressure_ideal(double rho, double T, Abundance const & X);
double partialRho_ideal(double rho, double T, Abundance const & X);
double partialT_ideal(double rho, double T, Abundance const & X);
double partialX_ideal(chemical::elem, double rho, double T, Abundance const & X);
double energy_ideal(double rho, double T, Abundance const & X);
double UpartialRho_ideal(double rho, double T, Abundance const & X);
double UpartialT_ideal(double rho, double T, Abundance const & X);
double UpartialX_ideal(chemical::elem, double rho, double T, Abundance const & X);


//RADIATION GAS
double pressure_rad(double rho, double T, Abundance const & X);
double partialRho_rad(double rho, double T, Abundance const & X);
double partialT_rad(double rho, double T, Abundance const & X);
double partialX_rad(chemical::elem, double rho, double T, Abundance const & X);
double energy_rad(double rho, double T, Abundance const & X);
double UpartialRho_rad(double rho, double T, Abundance const & X);
double UpartialT_rad(double rho, double T, Abundance const & X);
double UpartialX_rad(chemical::elem, double rho, double T, Abundance const & X);


//COULOMB PRESSURE
double pressure_coul(double rho, double T, Abundance const & X);
double partialRho_coul(double rho, double T, Abundance const & X);
double partialT_coul(double rho, double T, Abundance const & X);
double partialX_coul(chemical::elem, double rho, double T, Abundance const & X);
double energy_coul(double rho, double T, Abundance const & X);
double UpartialRho_coul(double rho, double T, Abundance const & X);
double UpartialT_coul(double rho, double T, Abundance const & X);
double UpartialX_coul(chemical::elem, double rho, double T, Abundance const & X);


//ELECTRON DEGENERACY PRESSURE -- ZERO TEMP
double pressure_deg_zero(double rho, double T, Abundance const & X);
double partialRho_deg_zero(double rho, double T, Abundance const & X);
double partialT_deg_zero(double rho, double T, Abundance const & X);
double partialX_deg_zero(chemical::elem, double rho, double T, Abundance const & X);
double energy_deg_zero(double rho, double T, Abundance const & X);
double UpartialRho_deg_zero(double rho, double T, Abundance const & X);
double UpartialT_deg_zero(double rho, double T, Abundance const & X);
double UpartialX_deg_zero(chemical::elem, double rho, double T, Abundance const & X);



//ELECTRON DEGENERACY PRESSURE -- FINITE TEMP
double findEta(double rho, double beta, double mue);
double pressure_deg_finite(double rho, double T, Abundance const & X);
double partialRho_deg_finite(double rho, double T, Abundance const & X);
double partialT_deg_finite(double rho, double T, Abundance const & X);
double partialX_deg_finite(chemical::elem i, double rho, double T, Abundance const & X);
double energy_deg_finite(double rho, double T, Abundance const & X);
double UpartialRho_deg_finite(double rho, double T, Abundance const & X);
double UpartialT_deg_finite(double rho, double T, Abundance const & X);
double UpartialX_deg_finite(chemical::elem i, double rho, double T, Abundance const & X);



//ELECTRON DEGENERACY PRESSURE -- FINITE TEMP
//  find eta by inverting rho = rho0*( F1/2 + beta*F3/2)
double density_trap(double eta, double beta, double mue);
double findEta_trap(double rho, double beta, double mue);
double pressure_deg_trap(double rho, double T, Abundance const & X);
double energy_deg_trap(double rho, double T, Abundance const & X);



//ELECTRON DEGENERACY PRESSURE -- PARTIAL DEGENERACY
//pressure of electron gas under partial degeneracy (T finite)
//  uses a combination of methods depending on degeneracy
//	following Cox and Giuli 1980 -- see equation 24.51, the definition of w
double pressure_deg_partial(double rho, double T, Abundance const & X);
double partialRho_deg_partial(double rho, double T, Abundance const & X);
double partialT_deg_partial(double rho, double T, Abundance const & X);
double partialX_deg_partial(chemical::elem i, double rho, double T, Abundance const & X);
double energy_deg_partial(double rho, double T, Abundance const & X);
double UpartialRho_deg_partial(double rho, double T, Abundance const & X);
double UpartialT_deg_partial(double rho, double T, Abundance const & X);
double UpartialX_deg_partial(chemical::elem i, double rho, double T, Abundance const & X);

/*//ELECTRON DEGENERACY PRESSURE -- the Eggleton, Faulkner, Flanner (1973) approximation
//pressure of electron gas under partial degeneracy (T finite)
//  uses a combination of methods depending on degeneracy
double pressure_deg_EFF(double rho, double T, Abundance const & X);
double partialRho_deg_EFF(double rho, double T, Abundance const & X);
double partialT_deg_EFF(double rho, double T, Abundance const & X);
double partialX_deg_EFF(chemical::elem i, double rho, double T, Abundance const & X);
double energy_deg_EFF(double rho, double T, Abundance const & X);
double UpartialRho_deg_EFF(double rho, double T, Abundance const & X);
double UpartialT_deg_EFF(double rho, double T, Abundance const & X);
double UpartialX_deg_EFF(chemical::elem i, double rho, double T, Abundance const & X);//*/



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


//NEUTRON or PROTON DEGENERACY PRESSURE
/*double pressure_fermi(double rho, double T, Abundance const & X);
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

