// *************************************************************************************
//					NEWTONIAN POLYTROPE EXACT ERROR SCALE TEST
// scale_exact.cpp
//		In the Newtonian case, there are three exact solutions to the Lane-Emden equations:
//			n = 0 (the uniform, or isopycnic, model)
//			n = 1 (often a good approximation to a neutron star)
//			n = 5 (a "star" with an infinite envelope)
//		We test how errors from exact case scale with grid size
//		For the RK4 method, the errors should scale like N^{-4}.
//		For doubling grid, error scaling should be at 2^4
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"


//data for stars
double n[3] = {0.0, 1.0, 5.0};

int main(){
	//the scale factor
	int fac0 = 1;
	int fac1 = 2;	//try changing to 4
	//set two lengths, differing by a scale factor
	int length0 = 1024;
	int length1 = fac1*length0-(fac1-1);
	
	//for each n, we will create two stars with these different lengths
	Polytrope *star0, *star1;
		
	FILE *fp;
	for(int i=0; i<3; i++){
		printf("--------------------------------------\n");
		printf("Beginning polytrope n=%1.1f\n", n[i]);
		
		//create a starto discover the correct step size to use
		star1 = new Polytrope(n[i], length1);
		double dx = star1->getX(2)-star1->getX(1);
		delete star1;
		
		//now create two polytropes with differing grid sizes scaled by fac1
		star0 = new Polytrope(n[i], length0, fac1*dx);
		star1 = new Polytrope(n[i], length1, fac0*dx);

		//the surface
		double R = star0->getX(length0-1);
		
		char filename[100];
		sprintf(filename, "./scale_exact/scale_%1.1f.txt", n[i]);
		if(!(fp = fopen(filename, "w")) ){
			system("mkdir ./scale_exact");
			fp = fopen(filename, "w");
		}
		double exact, x, scale;
		for(int X=0; X<length0; X++){
			//only include points where the two models have the same x-value
			if(fabs(star0->getX(X) - star1->getX(fac1*X)) > 1e-4) continue;
			x = star0->getX(X);
			switch(i){
				case 0: exact = 1.0 - x*x/6.0; break;
				case 1: exact = sin(x)/x; break;
				case 2: exact = pow( (3.0 + x*x)/3.0, -0.5); break;
			}
			//caculate the ration of the errors, showing how they scale
			scale = fabs(  ( star0->getY(fac0*X)-exact )/( star1->getY(fac1*X)-exact )  );
			fprintf(fp, "%0.32le\t%0.32le\t%0.32le\t%0.32le\n", x/R, scale, fabs(star0->getY(fac0*X)-exact), fabs(star1->getY(fac1*X)-exact));
		}
		fclose(fp);
					
		char outname[100];
		sprintf(outname, "./scale_exact/scale_%1.1f.png", n[i]);
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
		fprintf(gnuplot, "set title 'Scaling Grid by Factor %d, polytrope n=%1.1f'\n", fac1, n[i]);
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
	fprintf(gnuplot, "set title 'Scaling of Errors for Analytic Polytropes, scale factor = %d'\n", fac1);
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
	//fprintf(gnuplot, "set yrange [1e-16: 1e-12]\n");
	int thou = int(length0/1024);
	fprintf(gnuplot, "set ylabel 'Ratio |θ_{%dK}/θ_{%dK}|'\n", thou, fac1*thou);
	fprintf(gnuplot, "set logscale y %d\n", fac1);
	fprintf(gnuplot, "set format y '%d^{%%L}'\n", fac1);
	fprintf(gnuplot, "set ytics %d\n", fac1);
	sprintf(filename, "./scale_exact/scale_%1.1f.txt", n[0]);
	fprintf(gnuplot, "set arrow 1 from 0.0,%le to 1.01, %le lc rgb 'red' nohead\n", pow(fac1,4), pow(fac1,4));
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'n=%1.1f'", filename, n[0]);
	for(int i=1; i<3; i++){
		sprintf(filename, "./scale_exact/scale_%1.1f.txt", n[i]);
		fprintf(gnuplot, ", '%s' u 1:2 w l t 'n=%1.1f'", filename, n[i]);
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	
	
	return 0;
}