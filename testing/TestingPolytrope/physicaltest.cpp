// *************************************************************************************
//					NEWTONIAN POLYTROPE RMSR TEST
// physicaltest.cpp
//		The Lane-Emden equation solves for the variable θ, related to ρ,P
//		This takes the solution of θ, calculates ρ,P, and puts into hydrostatic equations
//			and computes a residual of those equations with these variables
//		This same process is performed by the SSR() method of each Star,
//			However, this program calculates and graphs the residual at each grid point
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"


//data for stars
double n[10] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};

int main(){
	//ask user for tolerance of calculation
	int length = 10000;
	
	//initialize a series of polytropes
	Polytrope *star;
	
	double checkEuler[10];
	double checkPoiss[10];
	
	FILE *fp;
	
	for(int i=0; i<10; i++){
		//calculate GR polytrope
		star = new Polytrope(n[i], length);
		//star->printStar();
		int L = star->length();
		double R = star->Radius();
		
		char filename[100];
		sprintf(filename, "./physicaltest/test_%1.1f.txt", n[i]);
		FILE *fp;
		if(!(fp = fopen(filename, "w")) ){
			system("mkdir ./physicaltest");
			fp = fopen(filename, "w");
		}
		
		//sum up errors in equations
		checkEuler[i] = 0.0;
		checkPoiss[i] = 0.0;
		double d2Phi = 0.0;
		double e1, e2, n1, n2;
		for(int X=4; X<L-4; X++){
			//Euler equation
			e1 = fabs( star->dPdr(X) + star->rho(X)*star->dPhidr(X) );
			n1 = fabs( star->dPdr(X)) + fabs(star->rho(X)*star->dPhidr(X));
			//Poisson equation -- requires numerical derivative for d^2Phi/dr^2
			//calculate numerical derivatives
			double a3,a2,a1,b1,b2,b3, h;
			// d2DPhi/dr2 = (d/dr)(dDphi/dr), dDPhi/dr = g y4
			b3=star->dPhidr(X-3);
			b2=star->dPhidr(X-2);
			b1=star->dPhidr(X-1);
			a1=star->dPhidr(X+1);
			a2=star->dPhidr(X+2);
			a3=star->dPhidr(X+3);
			h=star->rad(X+1)-star->rad(X);
			//d2Phi = 4.0*m_pi*star->rho(X) - 2.0*star->mr(X)*pow(star->rad(X),-3);
			d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h);
			e2 = fabs( 4.0*m_pi*star->rho(X)*star->rad(X) - 2.0*star->dPhidr(X) - d2Phi*star->rad(X) );
			n2 = fabs(4.0*m_pi*star->rho(X)*star->rad(X)) + fabs(2.0*star->dPhidr(X)) + fabs(d2Phi*star->rad(X));
			//add absolute error
			e1 = e1/n1;
			e2 = e2/n2;
			checkEuler[i] += e1;
			checkPoiss[i] += e2;
			fprintf(fp, "%le\t%0.30le\t%0.30le\n", star->rad(X)/R, e1, e2);	
		}
		
		printf("-----------------------------------------------------\n");
		printf("Sum of Errors in Euler Equation:\t%le\n", checkEuler[i]);
		printf("Sum of Errors in Poisson Equation:\t%le\n", checkPoiss[i]);
		fclose(fp);
		delete star;
		
		char outname[100];
		sprintf(outname, "./physicaltest/test_%1.1f.png", n[i]);
		FILE *gnuplot = popen("gnuplot -persist", "w");
		fprintf(gnuplot, "reset\n");
		fprintf(gnuplot, "set term png size 1000,1000\n");
		fprintf(gnuplot, "set samples %d\n", L);
		fprintf(gnuplot, "set output '%s'\n", outname);
		fprintf(gnuplot, "set xlabel 'r/R'\n");
		fprintf(gnuplot, "set ylabel 'residual'\n");
		fprintf(gnuplot, "set logscale y\n");
		fprintf(gnuplot, "set format y '10^{%%L}'\n");
		fprintf(gnuplot, "set title 'Residuals in Background Equations, Polytrope n=%1.1f, N_{star}= %d'\n", n[i], length);
		fprintf(gnuplot, "plot '%s' u 1:2 w p t 'Euler'", filename);
		fprintf(gnuplot, ",    '%s' u 1:3 w p t 'Poisson'", filename);
		fprintf(gnuplot, "\n");
		pclose(gnuplot);//*/	
	}
	//plot everything in single graph, for simplicity
	//only going to use these three in the published graphs
	char filename[91];
	char outname[91];
	int thou = length/1000;
	double nn[3] = {1.0, 1.5, 3.0};
	sprintf(outname, "./physicaltest/test_all_%dK.png", thou);
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 2000,1000\n");
	fprintf(gnuplot, "set samples %d\n", length);
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Residuals in Background Equations for Select Polytropes, N_{star} = %d'\n", length);
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
	fprintf(gnuplot, "set ylabel 'Normalized Residuals'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set key opaque\n");
	sprintf(filename, "./physicaltest/test_%1.1f.txt", nn[0]);
	fprintf(gnuplot, "plot '%s' u 1:2 w p t 'Euler, n=%1.1f'", filename, nn[0]);
	fprintf(gnuplot, ",    '%s' u 1:3 w p t 'Poisson, n=%1.1f'", filename, nn[0]);
	for(int i=1; i<3; i++){
		sprintf(filename, "./physicaltest/test_%1.1f.txt", nn[i]);
		fprintf(gnuplot, ",  '%s' u 1:2 w p t 'Euler, n=%1.1f'", filename, nn[i]);
		fprintf(gnuplot, ",  '%s' u 1:3 w p t 'Poisson, n=%1.1f'", filename, nn[i]);
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	
	
	return 0;
}





