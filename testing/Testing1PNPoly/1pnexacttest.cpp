// *************************************************************************************
//					POST-NEWTONIAN POLYTROPE EXACT TEST
// 1pnexacttest.cpp
//		This tests the 1PN polytrope against known analytic solutions in other regimes.
//		In the Newtonian case, there are three exact solutions to the Lane-Emden equations:
//			n = 0 (the uniform, or isopycnic, model)
//			n = 1 (often a good approximation to a neutron star)
//			n = 5 (a "star" with an infinite envelope)
//		When these are compared to 1PN solution, the difference should scale with redshift z
//		In the GR case, there is one exact solution (see Tooper 1964, eq 3.3)
//			n = 0 (the uniform relativistic model)
//		When compared to the 1PN solution, the difference should scale like z^2
//		To make better, should use xproper instead of x... to do...
// *************************************************************************************

#include "../../src/constants.h"
#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/PNPolytrope.cpp"
#include "../../lib/Splinor.cpp" //it radii between 1pn (harmonic) and GR (schwarzschild)

void printCompare(int, int, double*, double*, PNPolytrope*);
void printCompareGR(int, int, double*, double*, PNPolytrope*);

//these are the three Newtonian polytropes that have exact solutions
double n[3] = {0.0, 1.0, 5.0};

int main(){
	//data for stars
	double sigma;
	
	//initialize a series of polytropes
	double *theta, *xi;
	PNPolytrope *star;
	
	int L = 10000;
	double dx= 0.0;
	double xi2analytic[10], xi2numeric[10];
	double *thetaGR;
	
	FILE *fp;
	system("mkdir ./1pnexacttest");
	
	//for each polytrope index
	for(int i = 0; i<3; i++){
		//for several magnitudes of sigma
		if(n[i]==0.0) fp = fopen("./1pnexacttest/testofradii.txt", "w");
		for(int s=2; s<6; s++){
			//calculate GR polytrope
			sigma = pow(10,-s);
			if(n[i]==5.0) star = new PNPolytrope(n[i],sigma,L,dx);
			else {
				star = new PNPolytrope(n[i], sigma, L);
				dx = star->getX(2)-star->getX(1);	
			}
			theta = new double[L];
			xi = new double[L];
		
			//compute analytic solution
			for(int x=0; x<L; x++){
				xi[x] = star->getX(x);
				//these are the three exact solutions
				if(n[i]==0.0)       theta[x] = 1.0 - xi[x]*xi[x]/6.0;
				else if (n[i]==1.0) theta[x] = sin(xi[x])/xi[x];
				else if (n[i]==5.0) theta[x] = pow( (3.0 + xi[x]*xi[x])/3.0, -0.5);
			}
			
			//print the two for comparison
			printCompare(i, s, xi, theta, star);
		
			if(n[i]==0.0){
				delete[] xi;
				xi = new double[L];
				thetaGR = new double[L];
				for(int x=0; x<L; x++){
					xi[x] = star->getX(x);
					thetaGR[x] = ( (1.+3.*sigma)*sqrt(1.-2./3.*sigma*xi[x]*xi[x]) - (1.+sigma) );
					thetaGR[x] /= ( 3.*(1.+sigma) - (1.+3.*sigma)*sqrt(1.-2./3.*sigma*xi[x]*xi[x]) );
					thetaGR[x] /= sigma;
				}
				printCompareGR(i, s, xi, thetaGR, star);				
				
				//compare to Eq 3.4 in Tooper for terminal radius
				xi2analytic[i] = 6.*(1.+2.*sigma)/((1.+3.*sigma)*(1.+3.*sigma));
				xi2numeric[i]  = xi[L-1]*xi[L-1];
				fprintf(fp, "sigma=%0.2le\t exct=%0.8le\t", sigma, xi2analytic[i]);
				fprintf(fp, "num =%0.8le\t diff = %0.2le\n", xi2numeric[i], fabs(xi2analytic[i]-xi2numeric[i]));	
				delete[] thetaGR;
			}
			delete star;
			delete[] theta;
			delete[] xi;
		}
		
		if(n[i]==0.0){
			fclose(fp);
			printf("radial comparisons\n");
			for(int s=2; s<6; s++){
				//compare to Eq 3.4 in Tooper for terminal radius
				sigma = pow(10,-s);
				printf("\tsigma=%0.1le\t exct=%0.8le\t", sigma, xi2analytic[s]);
				printf("num =%0.8le\t diff = %0.8le\n", xi2numeric[s], fabs(xi2analytic[s]-xi2numeric[s]));
			}
		}

		//plot everything in single graph, for simplicity
		char filename[91];
		char outname[91];
		sprintf(outname, "./1pnexacttest/logdif_n%1.1f.png", n[i]);
		FILE *gnuplot = popen("gnuplot -persist", "w");
		fprintf(gnuplot, "reset\n");
		fprintf(gnuplot, "set term png size 1000,1000\n");
		fprintf(gnuplot, "set output '%s'\n", outname);
		fprintf(gnuplot, "set encoding utf8\n");
		fprintf(gnuplot, "set title 'Comparison of 1pn model to 0pn exact solution, n=%1.0f'\n",n[i]);
		fprintf(gnuplot, "set xlabel 'x/X'\n");
		fprintf(gnuplot, "set logscale y 10\n");
		fprintf(gnuplot, "set format y '10^{%%L}'\n");
		fprintf(gnuplot, "set ytics 10\n");
		fprintf(gnuplot, "set ylabel 'log10|θ_{1pn}-θ_{ex}|' enhanced\n");
		sprintf(filename, "./1pnexacttest/n%1.1f_sigma%0.1e.txt", n[i], pow(10,-2));
		fprintf(gnuplot, "plot '%s' u 2:5 w l t 'sigma=%0.1le'", filename, pow(10,-2));
		for(int s=3; s<6; s++){
			sprintf(filename, "./1pnexacttest/n%1.1f_sigma%0.1e.txt", n[i], pow(10,-s));
			fprintf(gnuplot, ", '%s' u 2:5 w l t 'sigma=%0.1le'", filename, pow(10,-s));
		}
		fprintf(gnuplot, "\n");
		pclose(gnuplot);
	}
	
	//at the end, pot the GR comparison
	char filename[91];
	char outname[91];
	sprintf(outname, "./1pnexacttest/logdif_n%1.1f_GR.png", 0.0);
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,1000\n");
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set encoding utf8\n");
	fprintf(gnuplot, "set title 'Comparison of 1pn model to GR exact solution, n=%1.0f'\n",0.0);
	fprintf(gnuplot, "set xlabel 'x/X'\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set ytics 10\n");
	fprintf(gnuplot, "set ylabel 'log10|θ_{1pn}-θ_{ex}|' enhanced\n");
	sprintf(filename, "./1pnexacttest/n%1.1f_sigma%0.1e_GR.txt", 0.0, pow(10,-2));
	fprintf(gnuplot, "plot '%s' u 2:5 w l t 'sigma=%0.1le'", filename, pow(10,-2));
	for(int s=3; s<6; s++){
		sprintf(filename, "./1pnexacttest/n%1.1f_sigma%0.1e_GR.txt", 0.0, pow(10,-s));
		fprintf(gnuplot, ", '%s' u 2:5 w l t 'sigma=%0.1le'", filename, pow(10,-s));
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	
	return 0;
}

void printCompare(int i, int s, double *xi, double *theta, PNPolytrope *star){
	int L = star->length();
	
	//create names for files to be opened
	char filename[91];
	char outname[91];
	char openmyplot[45];
	//save data to folder to avoid clutter - make sure folder exists
	sprintf(filename, "./1pnexacttest/n%1.1f_sigma%0.1e.txt", n[i], pow(10.0,-s));
	FILE *fp;
	if(!(fp = fopen(filename, "w")) ){
		system("mkdir ./1pnexacttest");
		fp = fopen(filename, "w");
	}
	//print results to text file
	double XI = xi[L-1];
	for(int x=0; x<L; x++){
		fprintf(fp, "%d\t%0.16f\t%0.16f\t%0.16f\t%0.16f\n", x, xi[x]/XI, theta[x], star->getY(x), fabs(theta[x]-star->getY(x)) );
	}
	fclose(fp);	
}

void printCompareGR(int i, int s, double *xi, double *theta, PNPolytrope *star){
	int L = star->length();
	
	//create names for files to be opened
	char filename[91];
	char outname[91];
	char openmyplot[45];
	//save data to folder to avoid clutter - make sure folder exists
	sprintf(filename, "./1pnexacttest/n%1.1f_sigma%0.1e_GR.txt", n[i], pow(10.0,-s));
	FILE *fp;
	if(!(fp = fopen(filename, "w")) ){
		system("mkdir ./1pnexacttest");
		fp = fopen(filename, "w");
	}
	//print results to text file
	double XI = xi[L-1];
	for(int x=0; x<L; x++){
		fprintf(fp, "%d\t%0.16f\t%0.16f\t%0.16f\t%0.16f\n", x, xi[x]/XI, theta[x], star->getY(x), fabs(theta[x]-star->getY(x)) );
	}
	fclose(fp);	
}






