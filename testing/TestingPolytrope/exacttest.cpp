// *************************************************************************************
//					NEWTONIAN POLYTROPE EXACT TEST
// exacttest.cpp
//		In the Newtonian case, there are three exact solutions to the Lane-Emden equations:
//			n = 0 (the uniform, or isopycnic, model)
//			n = 1 (often a good approximation to a neutron star)
//			n = 5 (a "star" with an infinite envelope)
//		We numerically compute simple polytropes of all three cases, and 
//			calculate their error against the known solutions.
// *************************************************************************************

#include "../../src/constants.h"
#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"

void printCompare(int, double*, double*, Polytrope*);

//these are the three polytropes that have exact solutions
double n[3] = {0.0, 1.0, 5.0};

int main(){
	//set the length of polytropes to use
	int length = 10000;
	
	//initialize a series of polytropes
	double *theta, *xi;
	Polytrope *star;
	
	for(int i=0; i<3; i++){
		//fp = fopen("./exacttest/testof.txt", "w");
		//calculate GR polytrope
		star = new Polytrope(n[i], length);
		//star->printStar();
		int L = star->length();
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
		printCompare(i, xi, theta, star);
		
		delete star;
		delete[] theta;
		delete[] xi;
	}
	//plot everything in single graph, for simplicity
	char filename[91];
	char  outname[91];
	sprintf(outname, "./exacttest/all_dif.png");
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 2000,1000\n");
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Errors in Analytic Polytropes, N_{star} = %d'\n", length);
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
	fprintf(gnuplot, "set ylabel 'log_{10} |θ_{exact}-θ_{num}|\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	sprintf(filename, "./exacttest/test_%1.1f.txt", n[0]);
	fprintf(gnuplot, "plot '%s' u 2:5 w l t 'n=%1.0f'", filename, n[0]);
	for(int i=1; i<3; i++){
		sprintf(filename, "./exacttest/test_%1.1f.txt", n[i]);
		fprintf(gnuplot, ", '%s' u 2:5 w l t 'n=%1.0f'", filename, n[i]);
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	
	return 0;
}

void printCompare(int i, double *xi, double *theta, Polytrope *star){
	int L = star->length();
	
	//create names for files to be opened
	char filename[91];
	char outname[91];
	char openmyplot[45];
	//save data to folder to avoid clutter - make sure folder exists
	sprintf(filename, "./exacttest/test_%1.1f.txt", n[i]);
	FILE *fp;
	if(!(fp = fopen(filename, "w")) ){
		system("mkdir ./exacttest");
		fp = fopen(filename, "w");
	}
	if(i==2) L = L-1;
	double XX = star->getX(L-1);
	//print results to text file
	double sum = 0.0;
	for(int x=1; x<L-1; x++){
		if( !isnan(theta[x]) ) sum += fabs(theta[x]-star->getY(x));
		fprintf(fp, "%d\t%0.30f\t%0.30f\t%0.30f\t%0.30f\n", x, xi[x]/XX, theta[x], star->getY(x), fabs(theta[x]-star->getY(x)) );
	}
	fclose(fp);
	printf("Cumm. Error in N=%1.1f:\t10^(%lf)\n", n[i], log10(sum) );
	//plot file in png in gnuplot, and open png
	sprintf(outname, "./exacttest/test_%1.1f.png",n[i]);
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1200,1000\n");
	fprintf(gnuplot, "set samples %d\n", L);
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'analytic test, n=%f'\n", n[i]);
	fprintf(gnuplot, "set xlabel 'scaled radius r/R'\n");
	fprintf(gnuplot, "set ylabel 'log_{10} |θ_{exact}-θ_{num}|\n");
	fprintf(gnuplot, "set logscale y 10\nset format y '10^{%%L}'\n");
	fprintf(gnuplot, "plot '%s' u 2:5 w l t 'n=%0.0f'", filename, n[i]);
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
}








