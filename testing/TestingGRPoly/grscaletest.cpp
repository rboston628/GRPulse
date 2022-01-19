// *************************************************************************************
//					RELATIVISTIC POLYTROPE ERROR SCALE TEST
// grscaletest.cpp
//		This program checks that the errors scale correctly with grid size
//		Suppose for given grid size N, the numerical solution from RK4 is 
//			y(N) = y(infty) + E(N^4)
//		Then y(2N)-y(N) = E(2^4N^4) - E(N^4) ~ 2^4, up to proportionality constant
//		If we construct fraction [y(2N)-y(N)]/[y(4N)-y(2N)] = 2^4
//			where the proportionality constant cancels.
//		However, we consistently get scaling at 2^2, showing our method is like N^-2
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/GRPolytrope.cpp"


//data for stars
double n[10] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};

int main(){
	//the scale factor
	int fac0 = 1;
	int fac1 = 2;
	int fac2 = fac1*fac1;
	//set three lengths
	int length0 = 1024;
	int length1 = fac1*length0-(fac1-1);
	int length2 = fac2*length0-(fac2-1);
	
	//initialize a series of polytropes
	GRPolytrope *star0, *star1, *star2;
		
	FILE *fp;
	system("rm -r ./scaletest");
	system("mkdir ./scaletest");
	
	double zsurf = 1.e-2, sigma;
	
	for(int i=0; i<9; i++){
		printf("--------------------------------------\n");
		printf("Beginning polytrope n=%1.1f\n", n[i]);
		
		star2 = new GRPolytrope(n[i], zsurf, length2);
		sigma = star2->getSigma();
		double dx = star2->getX(2)-star2->getX(1);
		delete star2;
		printf("make anew\n");
		star0 = new GRPolytrope(n[i], sigma, length0, fac2*dx);
		star1 = new GRPolytrope(n[i], sigma, length1, fac1*dx);
		star2 = new GRPolytrope(n[i], sigma, length2, fac0*dx);

		
		//at an arbitrary point, print out the radius to check they all match
		int test = 873;
		printf("   dx0 =\t %0.30le\n",            (star0->getX(     test+1)-star0->getX(     test)));
		printf("%2d*dx1 =\t %0.30le\n", fac1,fac1*(star1->getX(fac1*test+1)-star1->getX(fac1*test)));
		printf("%2d*dx2 =\t %0.30le\n", fac2,fac2*(star2->getX(fac2*test+1)-star2->getX(fac2*test)));
		
		printf("r0[%d] =\t %0.30le\n", test     , star0->rad(test));
		printf("r1[%d] =\t %0.30le\n", test*fac1, star1->rad(test*fac1));
		printf("r2[%d] =\t %0.30le\n", test*fac2, star2->rad(test*fac2));
		
		//also check total stellar radius to make sure they match
		printf("R0 =\t %0.30le\n", star0->rad(length0-1));
		printf("R1 =\t %0.30le\n", star1->rad(length1-1));
		printf("R2 =\t %0.30le\n", star2->rad(length2-1));
		//also check the lengths have scaled
		printf("L0 =\t %d\n", star0->length());
		printf("L1 =\t %d\n", star1->length());
		printf("L2 =\t %d\n", star2->length());
		
		
		//printf("sigma0 =\t %0.30le\n", star0->getSigma());
		//printf("sigma1 =\t %0.30le\n", star1->getSigma());
		//printf("sigma2 =\t %0.30le\n", star2->getSigma());
		printf("z0 =\t %0.30le\n", star0->getRedshift());
		printf("z1 =\t %0.30le\n", star1->getRedshift());
		printf("z2 =\t %0.30le\n", star2->getRedshift());
		printf("ZZ =\t %0.30le\n", (star1->getRedshift()-star0->getRedshift())/(star2->getRedshift()-star1->getRedshift()));		
				
		double R = star0->Radius();
		
		//print their comparison to file for plotting
		char filename[100];
		sprintf(filename, "./scaletest/scale_%1.1f.txt", n[i]);
		fp = fopen(filename,"w");
		double scale1, scale2;
		for(int X=0; X<length0; X++){
			scale1 = fabs(star1->getY(fac1*X)-star0->getY(fac0*X));
			scale2 = fabs(star2->getY(fac2*X)-star1->getY(fac1*X));
			fprintf(fp, "%0.32le\t%0.32le\t%0.32le\t%0.32le\n", star0->rad(X)/R, scale1/scale2, scale1, scale2);
		}
		fclose(fp);
		
		delete star0;
		delete star1;
		delete star2;
	}
	
	//plot everything in single graph, for simplicity
	char filename[91];
	char outname[91];
	int thou = int(length0/1024);
	sprintf(outname, "./scaletest/scale_all_%dK.png", thou);
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 2000,1000\n");
	fprintf(gnuplot, "set samples %d\n", length1);
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Scaling of Lane-Emden solution, for scale factor = %d'\n", fac1);
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
	fprintf(gnuplot, "set yrange [%lf: %lf]\n", pow(2,-1), pow(2,6) );
	fprintf(gnuplot, "set ylabel 'Ratio |[θ_{%dK}-θ_{%dK}]/[θ_{%dK}-θ_{%dK}]|'\n", fac1*thou,thou,fac2*thou,fac1*thou);
	fprintf(gnuplot, "set logscale y %d\n", fac1);
	fprintf(gnuplot, "set format y '%d^{%%L}'\n", fac1);
	fprintf(gnuplot, "set ytics %d\n", fac1);
	sprintf(filename, "./scaletest/scale_%1.1f.txt", n[0]);
	fprintf(gnuplot, "set arrow 1 from 0.0,%le to 1.01, %le lc rgb 'red' nohead\n", pow(fac1,4), pow(fac1,4));
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'n=%1.1f'", filename, n[0]);
	for(int i=1; i<9; i++){
		sprintf(filename, "./scaletest/scale_%1.1f.txt", n[i]);
		fprintf(gnuplot, ",  '%s' u 1:2 w l t 'n=%1.1f'", filename, n[i]);
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	
	
	
	return 0;
}