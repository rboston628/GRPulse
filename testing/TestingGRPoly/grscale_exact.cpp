// *************************************************************************************
//					RELATIVISTIC POLYTROPE EXACT TEST
// grscale_exact.cpp
//		In the GR case, there is one exact solution (see Tooper 1964, eq 3.3)
//			n = 0 (the uniform relativistic model)
//		We test how errors from exact case scale with grid size
//		For the RK4 method, the errors should scale like N^{-4}.
//		For doubling grid, error scaling should be at 2^4
// IMPORTANT NOTE:
//		If this does not scale, check the individual errors, 
// 		If these are very small, then proper scaling not expected
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/GRPolytrope.cpp"

//data for stars
double zsurf[4] = {1e-2, 1e-3, 1e-4, 1e-5};
double sigma[4];

int main(){
	//the scale factor
	int fac0 = 1;
	int fac1 = 2;	//try changing to 4
	//set two lengths, differing by a scale factor
	int length0 = 1024;
	int length1 = fac1*length0-(fac1-1);
	
	//initialize a series of polytropes
	double n = 0.0;
	GRPolytrope *star0, *star1;
	
	FILE *fp;
	
	for(int i=0; i<4; i++){
		printf("--------------------------------------\n");
		printf("Beginning polytrope n=0.0, zsurf=%f\n", zsurf[i]);
		star1 = new GRPolytrope(n, zsurf[i], length1);
		sigma[i] = star1->getSigma();
		double dx = star1->getX(2)-star1->getX(1);
		
		delete star1;
		printf("make anew\n");
		star0 = new GRPolytrope(n, sigma[i], length0, fac1*dx);
		star1 = new GRPolytrope(n, sigma[i], length1, fac0*dx);

		double R = star0->getX(length0-1);
				
		char filename[100];
		system("mkdir -p scale_exact");
		sprintf(filename, "./scale_exact/scale_%1.1le.txt", zsurf[i]);
		fp = fopen(filename,"w");
		double exact, x, scale, rootx2, sp1;
		for(int X=0; X<length0; X++){
			if(star0->getX(X) != star1->getX(fac1*X)) continue;
			x = star0->getX(X);
			rootx2 = (1.+3.*sigma[i])*sqrt(1.-2./3.*sigma[i]*x*x);
			sp1 = 1.+sigma[i];
			exact  = ( rootx2 - sp1 )/( 3.*sp1 - rootx2 )/sigma[i];	
			scale = fabs(  ( star0->getY(fac0*X)-exact )/( star1->getY(fac1*X)-exact )  );
			fprintf(fp, "%0.32le\t%0.32le\t%0.32le\t%0.32le\n", x/R, scale, fabs(star0->getY(fac0*X)-exact), fabs(star1->getY(fac1*X)-exact));
		}
		fclose(fp);		
		
		char outname[100];
		sprintf(outname, "./scale_exact/scale_%1.1le.png", zsurf[i]);
		FILE *gnuplot = popen("gnuplot -persist", "w");
		fprintf(gnuplot, "reset\n");
		fprintf(gnuplot, "set term png size 2000,1000\n");
		fprintf(gnuplot, "set samples %d\n", length1);
		fprintf(gnuplot, "set output '%s'\n", outname);
		fprintf(gnuplot, "set xlabel 'r/R'\n");
		int thou = int(length0/1024);
		fprintf(gnuplot, "set ylabel 'diff'\n");
		fprintf(gnuplot, "set logscale y %d\n", 10);
		fprintf(gnuplot, "set ytics %d\n", 10);
		fprintf(gnuplot, "set format y '%d^{%%L}'\n", 10);
		fprintf(gnuplot, "set title 'Scaling Grid by Factor %d, polytrope n=0, z=%1.1e'\n", fac1, zsurf[i]);
		fprintf(gnuplot, "plot '%s' u 1:3 w l t '%dK'", filename, thou);
		fprintf(gnuplot, ",    '%s' u 1:4 w l t '%dK'", filename, fac1*thou);
		fprintf(gnuplot, "\n");
		pclose(gnuplot);
		
		delete star0;
		delete star1;
	}
	
	//plot everything in single graph, for simplicity
	char filename[91];
	char outname[91];
	sprintf(outname, "./scale_exact/scale_exact_all.png");
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 2000,1000\n");
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Scaling of Errors for Analytic GR Polytropes n=0.0, scale factor = %d'\n", fac1);
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
	//fprintf(gnuplot, "set yrange [1e-16: 1e-12]\n");
	int thou = int(length0/1024);
	fprintf(gnuplot, "set ylabel 'Ratio |θ_{%dK}/θ_{%dK}|'\n", thou, fac1*thou);
	fprintf(gnuplot, "set logscale y %d\n", fac1);
	fprintf(gnuplot, "set format y '%d^{%%L}'\n", fac1);
	fprintf(gnuplot, "set ytics %d\n", fac1);
	fprintf(gnuplot, "set yrange [2**-1:2**5]\n");
	sprintf(filename, "./scale_exact/scale_%1.1le.txt", zsurf[0]);
	fprintf(gnuplot, "set arrow 1 from 0.0,%le to 1.01, %le lc rgb 'red' nohead\n", pow(fac1,4), pow(fac1,4));
	fprintf(gnuplot, "plot '%s' u 1:2 w lp t 'z_{surf}=%1.1le'", filename, zsurf[0]);
	for(int i=1; i<4; i++){
		sprintf(filename, "./scale_exact/scale_%1.1le.txt", zsurf[i]);
		fprintf(gnuplot, ", '%s' u 1:2 w lp t 'z_{surf}=%1.1le'", filename, zsurf[i]);
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	
	system("rm ./scale_exact/scale*.txt");
	
	return 0;
}