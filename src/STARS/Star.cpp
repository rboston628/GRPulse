//**************************************************************************************
//							STAR ABSTRACT CLASS
// Star.cpp
//		Implements the output and sum-square residual functions for equilibrium stars
//		These methods are general enough to work with any sub-class of star
//**************************************************************************************
//EQUILIBRIUM STAR

#ifndef STARCLASS
#define STARCLASS

#include "Star.h"

//method to print pertinent values of star to .txt, and plot them in gnuplot
void Star::writeStar(char *c){
	//create names for files to be opened
	char pathname[256];
	char rootname[256];
	char txtname[256];
	char outname[256];
	if(c==NULL)	sprintf(pathname, "./out/%s", name);
	else{
		sprintf(pathname, "./%s/star/", c);
	}
	sprintf(txtname, "%s/%s.txt", pathname, name);
	sprintf(outname, "%s/%s.png", pathname, name);

	FILE *fp;
	if(!(fp = fopen(txtname, "w")) ){
		char command[256];
		sprintf(command, "mkdir -p %s", pathname);
		system(command);
		fp = fopen(txtname, "w");
	}
	//print results to text file
	// radius rho pressure gravity
	double irc=1./this->rho(0), ipc=1./P(0), R=Radius(), ig=1./dPhidr(length()-1);
	for(int X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			rad(X)/R, this->rho(X)*irc, -drhodr(X)*irc*R,
			P(X)*ipc, -dPdr(X)*ipc*R,
			mr(X)/Mass(), dPhidr(X)*ig);
		fflush(fp);
	}
	fclose(fp);	
	//plot file in png in gnuplot
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	char title[256]; graph_title(title);
	fprintf(gnuplot, "set title 'Profile for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'rho/rho_c, P/P_c, m/M, g/g_S'\n");
	fprintf(gnuplot, "set yrange [0:1]\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'rho'", txtname);
	fprintf(gnuplot, ", '%s' u 1:4 w l t 'P'", txtname);
	fprintf(gnuplot, ", '%s' u 1:6 w l t 'm'", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
}

//write the stellar profile, and display it on screen for easy viewing
void Star::printStar(char *c){
	writeStar(c);
	char outname[256];
	if(c==NULL) sprintf(outname, "./out/%s/%s.png", name, name);
	else {
		sprintf(outname, "./%s/star/%s.png", c, name);
	}
	char openmyplot[260];
	sprintf(openmyplot, "open %s", outname);
	system(openmyplot);
}

//calculte a root-mean-square scaled residual from backsubstitution into physical equations
//  note the name SSR (sum-square residual) is a misnomer
double Star::SSR(){	
	double checkEuler;
	double checkPoiss;
	int len = length();
			
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
		double a3,a2,a1,b1,b2,b3, h;
		// d2Phi/dr2 = (d/dr)(dPhi/dr)
		b3=dPhidr(X-3);
		b2=dPhidr(X-2);
		b1=dPhidr(X-1);
		a1=dPhidr(X+1);
		a2=dPhidr(X+2);
		a3=dPhidr(X+3);
		h=rad(X+1)-rad(X);
		d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h);
		e2 = fabs( 4.0*Gee()*m_pi*rho(X)*rad(X) - 2.0*dPhidr(X) - d2Phi*rad(X) );
		n2 = fabs( 4.0*Gee()*m_pi*rho(X)*rad(X)) + fabs(2.0*dPhidr(X)) + fabs(d2Phi*rad(X) );
		//add absolute error
		e1 = e1/n1;
		e2 = e2/n2;
		fprintf(fp, "%d\t%le\t%le\n", X, e1, e2);
		checkEuler += e1*e1;
		checkPoiss += e2*e2;
	}
	fclose(fp);
	//the default behavior will get print and graph this in the main directory
	//  better behavior would be to print it in the specific star's output directory
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", "SSR.png");
	char title[256]; graph_title(title);
	fprintf(gnuplot, "set title 'error for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'error'\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'Euler error'", "SSR.txt");
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'Poisson error'", "SSR.txt");
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	return sqrt((checkPoiss+checkEuler)/double(2*len-8));
}


#endif
