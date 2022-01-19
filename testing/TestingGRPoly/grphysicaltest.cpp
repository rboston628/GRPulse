// THIS will perform a physical test, recast x, theta into variables and put in equation

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/GRPolytrope.cpp"


//data for stars
double n[10] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};

int main(){
	//ask user for tolerance of calculation
	int length = 2048;
	double zsurf = 1.0e-4;
	
	//initialize a series of polytropes
	GRPolytrope *star;
	
	FILE *fp;
	
	double checkTOV1[ 10];
	double checkTOV2[ 10];
	double checkEuler[10];
	
	
	for(int i=0; i<10; i++){
		//fp = fopen("./exacttest/testof.txt", "w");
		//calculate GR polytrope
		star = new GRPolytrope(n[i], zsurf, length);
		//star->printStar();
		int L = star->length();
		double R = star->Radius();
		
		char filename[100];
		sprintf(filename, "./physicaltest/test_%1.1f.txt", n[i]);
		if(!(fp = fopen(filename,"w"))){
			system("mkdir ./physicaltest");
		}
		
		//sum up errors in equations
		checkTOV1[i] = 0.0;
		checkTOV2[i] = 0.0;
		checkEuler[i] = 0.0;
		double e1, e2, e3;
		double n1, n2, n3;
		double r = 0.0;
		double eL = 0.0;
		for(int X=1; X<L-1; X++){
			r = star->rad(X);
			eL = exp(-star->Lambda(X));
			//TOV1 -- equation 2.3 of Tooper1
			e1 = fabs( eL*(r*star->dNudr(X)     +1.0) - 1.0 - 8.0*m_pi*star->P(X)*pow(r,2) );
			n1 = fabs( eL*r*star->dNudr(X) ) + fabs(eL) + 1.0 + fabs(8.0*m_pi*star->P(X)*pow(r,2));
			//TOV2 -- equation 2.4 of Tooper1
			e2 = fabs( eL*(r*star->dLambdadr(X) - 1.0) + 1.0 - 8.0*m_pi*star->rho(X)*pow(r,2) );
			n2 = fabs( eL*r*star->dLambdadr(X)) +fabs(eL) + 1.0 + fabs(8.0*m_pi*star->rho(X)*pow(r,2));
			//Euler equation -- equation 2.6 of Tooper1
			e3 = fabs( star->dPdr(X) + 0.5*(star->rho(X)+star->P(X))*star->dNudr(X) );
			n3 = fabs(star->dPdr(X)) + fabs( 0.5*(star->rho(X)+star->P(X))*star->dNudr(X) );
		
			checkTOV1[i] += e1/n1;
			checkTOV1[i] += e2/n2;
			checkEuler[i]+= e3/n3;
			
			fprintf(fp, "%le\t%0.30le\t%0.30le\t%0.30le\n", r/R, e1, e2, e3);
		}
		
		printf("-----------------------------------------------------\n");
		printf("Sum of Errors in TOV Equation 1:\t%le\n", checkTOV1[i]);
		printf("Sum of Errors in TOV Equation 2:\t%le\n", checkTOV2[i]);
		printf("Sum of Errors in Euler Equation:\t%le\n", checkEuler[i]);
		
		delete star;
		fclose(fp);
		
		char outname[100];
		sprintf(outname, "./physicaltest/test_%1.1f.png", n[i]);
		FILE *gnuplot = popen("gnuplot -persist", "w");
		fprintf(gnuplot, "reset\n");
		fprintf(gnuplot, "set term png size 1000,1000\n");
		fprintf(gnuplot, "set samples %d\n", L);
		fprintf(gnuplot, "set output '%s'\n", outname);
		fprintf(gnuplot, "set xlabel 'radius'\n");
		fprintf(gnuplot, "set ylabel 'log10 error'\n");
		fprintf(gnuplot, "set title 'Errors in Physical Equations, n=%f'\n", n[i]);
		fprintf(gnuplot, "plot '%s' u 1:2 w l t 'TOV1'", filename);
		fprintf(gnuplot, ",    '%s' u 1:3 w l t 'TOV2'", filename);
		fprintf(gnuplot, ",    '%s' u 1:4 w l t 'Euler'", filename);
		fprintf(gnuplot, "\n");
		pclose(gnuplot);//*/
		
	}
	//plot everything in single graph, for simplicity
	char filename[91];
	char outname[91];
	int thou = length/1000;
	double nn[4] = {0.0, 1.0, 1.5, 3.0};
	sprintf(outname, "./physicaltest/test_all_%dK.png", thou);
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,1000\n");
	fprintf(gnuplot, "set samples %d\n", length);
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Residuals in Background Equations for Select Relativistic Polytropes, N_{star} = %d'\n", length);
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
	fprintf(gnuplot, "set yrange [1.e-24: 1.e-12]\n");
	fprintf(gnuplot, "set ylabel 'Normalized Residuals'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set key opaque\n");
	sprintf(filename, "./physicaltest/test_%1.1f.txt", nn[0]);
	fprintf(gnuplot, "plot '%s' u 1:2 w p t 'TOV 1, n=%1.1f'", filename, nn[0]);
	fprintf(gnuplot, ",    '%s' u 1:3 w p t 'TOV 2, n=%1.1f'", filename, nn[0]);
	fprintf(gnuplot, ",    '%s' u 1:4 w p t 'Euler, n=%1.1f'", filename, nn[0]);
	for(int i=1; i<4; i++){
		sprintf(filename, "./physicaltest/test_%1.1f.txt", nn[i]);
		fprintf(gnuplot, ", '%s' u 1:2 w p t 'TOV 1, n=%1.1f'", filename, nn[i]);
		fprintf(gnuplot, ", '%s' u 1:3 w p t 'TOV 2, n=%1.1f'", filename, nn[i]);
		fprintf(gnuplot, ", '%s' u 1:4 w p t 'Euler, n=%1.1f'", filename, nn[i]);
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	return 0;
}