// *************************************************************************************
//					RELATIVISTIC POLYTROPE RMSR SCALE TEST
// grscale_physical.cpp
//		This program checks that the RMSR errors scale correctly with grid size
//		The Lane-Emden equation solves for the variable θ, related to ρ,P
//		This takes the solution of θ, calculates ρ,P, and puts into hydrostatic equations
//			and computes a residual of those equations with these variables
//		We expect this RMSR to scale like N^{-4}
// IMPORTANT NOTE:
//		However, because RMSR is already at machine precision, there is no scaling
//		The ratio will cluster around 1 (2^0)
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/GRPolytrope.cpp"

//data for stars
double n[10] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};


int main(){
	//ask user for tolerance of calculation
	int fac = 2;
	int length1 = 1024;
	int length2 = fac*length1-(fac-1);
	
	double redshift = 1e-4;
	
	//initialize a series of polytropes
	//Polytrope *star;
	GRPolytrope *star1, *star2;
		
	char mkdir[100];
	sprintf(mkdir, "rm -r ./scale_physical");
	system(mkdir);
	sprintf(mkdir, "mkdir ./scale_physical");
	system(mkdir);
		
	FILE *fp;
	for(int i=0; i<10; i++){
		char filename[100];
		sprintf(filename, "./scale_physical/poly%1.1f.txt", n[i]);
		fp = fopen(filename,"w");
		//initialize the stars
		star1 = new GRPolytrope(n[i],redshift, length1);
		star2 = new GRPolytrope(n[i],redshift, length2);
		double R = star1->Radius();

		double e1, e2, e3;
		double n1, n2, n3;
		double r = 0.0;
		double eL = 0.0;
		for(int X=4; X<length1-4; X++){
			r = star1->rad(X);
			eL = exp(-star1->Lambda(X));
			//TOV1 -- equation 2.3 of Tooper1
			e1 = fabs( eL*(r*star1->dNudr(X)     +1.0) - 1.0 - 8.0*m_pi*star1->P(X)*pow(r,2) );
			//TOV2 -- equation 2.4 of Tooper1
			e2 = fabs( eL*(r*star1->dLambdadr(X) - 1.0) + 1.0 - 8.0*m_pi*star1->rho(X)*pow(r,2) );
			//Euler equation -- equation 2.6 of Tooper1
			e3 = fabs( star1->dPdr(X) + 0.5*(star1->rho(X)+star1->P(X))*star1->dNudr(X) );
			
			//repeat everything for star2
			int fX = fac*X;
			r = star2->rad(fX);				
			eL = exp(-star2->Lambda(fX));
			//TOV1 -- equation 2.3 of Tooper1
			n1 = fabs( eL*(r*star2->dNudr(fX)     +1.0) - 1.0 - 8.0*m_pi*star2->P(fX)*pow(r,2) );
			//TOV2 -- equation 2.4 of Tooper1
			n2 = fabs( eL*(r*star2->dLambdadr(fX) - 1.0) + 1.0 - 8.0*m_pi*star2->rho(fX)*pow(r,2) );
			//Euler equation -- equation 2.6 of Tooper1
			n3 = fabs( star2->dPdr(fX) + 0.5*(star2->rho(fX)+star2->P(fX))*star2->dNudr(fX) );
			
			//this is the ratio of the RMSR between star1 and star2, which should scale like fac^4
			e1 = e1/n1;
			e2 = e2/n2;
			e3 = e3/n3;
			fprintf(fp, "%le\t%0.30le\t%0.30le\t%0.30le\n", star1->rad(X)/R, e1, e2, e3);	
		}
		fclose(fp);
		delete star1;
		delete star2;
	}
	
	//plot everything in single graph, for simplicity
	//only going to use these three in the published graphs
	char filename[91];
	char outname[91];
	int thou = int(length1/1024);
	double nn[3] = {1.0, 1.5, 3.0};
	sprintf(outname, "./scale_physical/physscale_all_%dK.png", thou);
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 2000,1000\n");
	fprintf(gnuplot, "set samples %d\n", length1);
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Scaling of Residuals in Background Equations, scale factor = %d'\n", fac);
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
	//fprintf(gnuplot, "set yrange [5.e-17: 1.e-10]\n");
	fprintf(gnuplot, "set ylabel 'Ratio |r_{1K}/r_{2K}|'\n");
	fprintf(gnuplot, "set logscale y %d\n", fac);
	fprintf(gnuplot, "set format y '%d^{%%L}'\n", fac);
	fprintf(gnuplot, "set ytics %d\n", fac);
	sprintf(filename, "./scale_physical/poly%1.1f.txt", nn[0]);
	fprintf(gnuplot, "set arrow 1 from 0.0,%le to 1.01, %le lc rgb 'red' nohead\n", pow(fac,4), pow(fac,4));
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'Euler, n=%1.1f'", filename, nn[0]);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'Poisson, n=%1.1f'", filename, nn[0]);
	for(int i=1; i<3; i++){
		sprintf(filename, "./scale_physical/poly%1.1f.txt", nn[i]);
		fprintf(gnuplot, ",  '%s' u 1:2 w l t 'TOV_1, n=%1.1f'", filename, nn[i]);
		fprintf(gnuplot, ",  '%s' u 1:3 w l t 'TOV_2, n=%1.1f'", filename, nn[i]);
		fprintf(gnuplot, ",  '%s' u 1:4 w l t 'Euler, n=%1.1f'", filename, nn[i]);
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);


	return 0;
}