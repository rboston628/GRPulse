//

#ifndef CHANDRACONSTCPP
#define CHANDRACONSTCPP

#include "chandra.h"

double Chandrasekhar::factor_f(double x){
	return x*(2.*x*x-3.)*sqrt(1.+x*x) + 3.*asinh(x);
}
	
double Chandrasekhar::factor_g(double x){
		return 8.*x*x*x*(sqrt(1.+x*x) - 1.0) - factor_f(x);
}
	
//methods usign the weights and abscissae of Sagar 1991
//  these methods have tested accuracy for -1 < eta < 10
//  for regions eta > 10, Sagar suggests using asymptotic approximation
double FermiDirac::FermiDirac1Half(double eta, double beta){
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = xFermiDirac1Half[i];
		w = wFermiDirac1Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))*sqrt(1.+0.5*beta*x);
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDirac3Half(double eta, double beta){
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = xFermiDirac3Half[i];
		w = wFermiDirac3Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))*sqrt(1.+0.5*beta*x);
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDirac5Half(double eta, double beta){
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = xFermiDirac5Half[i];
		w = wFermiDirac5Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))*sqrt(1.+0.5*beta*x);
		sum += w*f;
	}
	return sum;
}


//use trapezoidal method to integrate the Fermi-Dirac functon F_k(eta)
double FermiDirac::FermiDirac(double k, double eta){
	return FermiDirac(k, eta, 0.0);
}

//use trapezoidal method to integrate the Fermi-Dirac functon F_k(eta, beta)
//  here beta  = k*T/(m c^2)
//  see Tooper 1969, ApJ 156
double FermiDirac::FermiDirac(double k, double eta, double beta) {
	double integral=0.0, int1=0.0, int2=0.0;
	double x = 0.0, dx = eta/500.0;
	while(x < 2.*eta){
		x += dx;
		int1 = int2;
		int2 = pow(x,k)*sqrt(1.+0.5*beta*x)/(1. + exp(x-eta));
		integral += 0.5*(int1+int2)*dx;
	}
	while(int2 > 1e-14){
		dx = (x/100.0);
		x += dx;
		int1 = int2;
		int2 = pow(x,k)*sqrt(1.+0.5*beta*x)/(1. + exp(x-eta));
		integral += 0.5*(int1+int2)*dx;
	}
	x += dx;
	int1 = int2;
	int2 = pow(x,k)*sqrt(1.+0.5*beta*x)/(1. + exp(x-eta));
	integral += 0.5*(int1+int2)*dx;
	return integral;
}
//partial derivatives
double FermiDirac::FermiDiracdelEta(double k, double eta, double beta) {
	double integral=0.0, int1=0.0, int2=0.0;
	double x = 0.0, dx = eta/500.0;
	while(x < 2.*eta){
		x += dx;
		int1 = int2;
		int2 = pow(x,k)*sqrt(1.+0.5*beta*x)/pow(1. + exp(x-eta),2.)*(-exp(x-eta));
		integral += 0.5*(int1+int2)*dx;
	}
	while(int2 > 1e-14){
		dx = (x/100.0);
		x += dx;
		int1 = int2;
		int2 = pow(x,k)*sqrt(1.+0.5*beta*x)/pow(1. + exp(x-eta),2.)*(-exp(x-eta));
		integral += 0.5*(int1+int2)*dx;
	}
	x += dx;
	int1 = int2;
	int2 = pow(x,k)*sqrt(1.+0.5*beta*x)/pow(1. + exp(x-eta),2.)*(-exp(x-eta));
	integral += 0.5*(int1+int2)*dx;
	return integral;
}
double FermiDirac::FermiDiracdelBeta(double k, double eta, double beta) {
	double integral=0.0, int1=0.0, int2=0.0;
	double x = 0.0, dx = eta/500.0;
	while(x < 2.*eta){
		x += dx;
		int1 = int2;
		int2 = pow(x,k)/(1. + exp(x-eta))*0.25*x/sqrt(1.+0.5*beta*x);
		integral += 0.5*(int1+int2)*dx;
	}
	while(int2 > 1e-14){
		dx = (x/100.0);
		x += dx;
		int1 = int2;
		int2 = pow(x,k)/(1. + exp(x-eta))*0.25*x/sqrt(1.+0.5*beta*x);
		integral += 0.5*(int1+int2)*dx;
	}
	x += dx;
	int1 = int2;
	int2 = pow(x,k)/(1. + exp(x-eta))*0.25*x/sqrt(1.+0.5*beta*x);
	integral += 0.5*(int1+int2)*dx;
	return integral;
}

//PARTIAL DERIVATIVES
//methods usign the weights and abscissae of Sagar 1991
//  used to calculate the partial derivatives
double FermiDirac::FermiDirac1HalfdelEta(double eta, double beta){
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = xFermiDirac1Half[i];
		w = wFermiDirac1Half[i];
		f = (1.+exp(x))/pow(1.+exp(x-eta),2)*sqrt(1.+0.5*beta*x)*exp(x-eta);
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDirac3HalfdelEta(double eta, double beta){
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = xFermiDirac3Half[i];
		w = wFermiDirac3Half[i];
		f = (1.+exp(x))/pow(1.+exp(x-eta),2)*sqrt(1.+0.5*beta*x)*exp(x-eta);
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDirac5HalfdelEta(double eta, double beta){
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = xFermiDirac5Half[i];
		w = wFermiDirac5Half[i];
		f = (1.+exp(x))/pow(1.+exp(x-eta),2)*sqrt(1.+0.5*beta*x)*exp(x-eta);
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDirac1HalfdelBeta(double eta, double beta){
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = xFermiDirac1Half[i];
		w = wFermiDirac1Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))/sqrt(1.+0.5*beta*x)*x;
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDirac3HalfdelBeta(double eta, double beta){
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = xFermiDirac3Half[i];
		w = wFermiDirac3Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))/sqrt(1.+0.5*beta*x)*x;
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDirac5HalfdelBeta(double eta, double beta){
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = xFermiDirac5Half[i];
		w = wFermiDirac5Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))/sqrt(1.+0.5*beta*x)*x;
		sum += w*f;
	}
	return sum;
}


//following Cox and Giuli, use Sommerfeld's asymptotic expansion
double FermiDirac::FermiDiracEtaLarge1Half(double eta){
	return (2./3.)*pow(eta, 1.5);
}	
double FermiDirac::FermiDiracEtaLarge3Half(double eta){
	return (2./5.)*pow(eta, 2.5);
}	
double FermiDirac::FermiDiracEtaLarge5Half(double eta){
	return (2./7.)*pow(eta, 3.5);
}	

//see Pichon 1989, eq 27-28
double FermiDirac::FermiDiracEtaLarge1Half(double eta, double beta){
	if(beta == 0.0) beta = 1.0;
	double y = sqrt( pow(1. + eta*beta, 2.) - 1.);
	double C1 = pow(M_PI*beta,2)/6./y;
	double C2 = pow(M_PI*beta/y/y,2)*C1*7./20.;
	double C3 = pow(M_PI*beta/y/y,4)*C1*31./168.;
	return (FermiDiracLittleF1Half(y) 
			+ (1.+eta*beta)*(C1 + C2 + C3*(4.*y*y+7.)))/sqrt(2.*beta*beta*beta);
}
double FermiDirac::FermiDiracEtaLarge3Half(double eta, double beta){
	if(beta == 0.0) beta = 1.0;
	double y = sqrt( pow(1. + eta*beta, 2.) - 1.);
	double C1 = pow(M_PI*beta,2)/6./y;
	double C2 = pow(M_PI*beta/y/y,2)*C1*7./20.;
	double C3 = pow(M_PI*beta/y/y,4)*C1*31./168.;
	return (FermiDiracLittleF3Half(y) 
			+ C1*(3.+2*eta*beta) - C2 - C3*( (4.*eta*beta+6.)*eta*beta + 3.))/sqrt(2.*pow(beta,5));
}
double FermiDirac::FermiDiracEtaLarge5Half(double eta, double beta){
	if(beta == 0.0) beta = 1.0;
	double y = sqrt( pow(1. + eta*beta, 2.) - 1.);
	double C1 = pow(M_PI*beta,2)/6./y;
	double C2 = pow(M_PI*beta/y/y,2)*C1*7./20.;
	double C3 = pow(M_PI*beta/y/y,4)*C1*31./168.;
	return (FermiDiracLittleF5Half(y) 
			+ C1*(5.+eta*beta) - C2*((((2.*eta*beta+10.)*eta*beta + 15.)*eta*beta+15.)*eta*beta+5.)
			+ C3*(3.+5.*eta*beta))/sqrt(2.*pow(beta,7));
}
// these functions appear in asymptotic expansion (Pichon 1989, Cox&Giuli 1980)
double FermiDirac::FermiDiracLittleF1Half(double y){
	return 0.5*(y*sqrt(1.+y*y) - asinh(y));
}
double FermiDirac::FermiDiracLittleF3Half(double y){
	return 1./3.*y*y*y - FermiDiracLittleF1Half(y);
}
double FermiDirac::FermiDiracLittleF5Half(double y){
	return 5./8.*y*(1.+0.4*y*y)*sqrt(1.+y*y) - 2./3*y*y*y;
}

#endif