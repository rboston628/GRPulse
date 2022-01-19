// *************************************************************************************
//					NEWTONIAN POLYTROPE RMSR SCALE TEST
// scale_physical.cpp
//		This program checks that the RMSR errors scale correctly with grid size
//		The Lane-Emden equation solves for the variable θ, related to ρ,P
//		This takes the solution of θ, calculates ρ,P, and puts into hydrostatic equations
//			and computes a residual of those equations with these variables
//		We expect this RMSR to scale like N^{-4}
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"


//data for stars
double n[10] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};


int main(){
	//ask user for tolerance of calculation
	int fac = 2;
	int length1 = 1024;
	int length2 = fac*length1-(fac-1);
	
	//initialize a series of polytropes
	//Polytrope *star;
	Polytrope *star1, *star2;
		
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
		star1 = new Polytrope(n[i],length1);
		star2 = new Polytrope(n[i],length2);

		double d2Phi = 0.0;
		double e1, e2, n1, n2;
		double R = star1->Radius();
		for(int X=4; X<length1-4; X++){
			//Euler equation
			e1 = fabs( star1->dPdr(X) + star1->rho(X)*star1->dPhidr(X) );
			//Poisson equation
			//calculate numerical derivatives
			double a3,a2,a1,b1,b2,b3, h;
			// d2DPhi/dr2 = (d/dr)(dDphi/dr), dDPhi/dr = g y4
			b3=star1->dPhidr(X-3);
			b2=star1->dPhidr(X-2);
			b1=star1->dPhidr(X-1);
			a1=star1->dPhidr(X+1);
			a2=star1->dPhidr(X+2);
			a3=star1->dPhidr(X+3);
			h=star1->rad(X+1)-star1->rad(X);
			//d2Phi = 4.0*m_pi*star->rho(X) - 2.0*star->mr(X)*pow(star->rad(X),-3);
			d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h);
			e2 = fabs( 4.0*m_pi*star1->rho(X)*star1->rad(X) - 2.0*star1->dPhidr(X) - d2Phi*star1->rad(X) );
			
			//repeat everything for star2
			int fX = fac*X;
			n1 = fabs( star2->dPdr(fX) + star2->rho(fX)*star2->dPhidr(fX) );
			b3=star2->dPhidr(fX-fac*3);
			b2=star2->dPhidr(fX-fac*2);
			b1=star2->dPhidr(fX-fac*1);
			a1=star2->dPhidr(fX+fac*1);
			a2=star2->dPhidr(fX+fac*2);
			a3=star2->dPhidr(fX+fac*3);
			h=star2->rad(fX+fac*1)-star2->rad(fX);
			//d2Phi = 4.0*m_pi*star->rho(X) - 2.0*star->mr(X)*pow(star->rad(X),-3);
			d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h);
			n2 = fabs( 4.0*m_pi*star2->rho(fX)*star2->rad(fX) - 2.0*star2->dPhidr(fX) - d2Phi*star2->rad(fX) );
			
			//this is the ratio of the RMSR between star1 and star2, which should scale like fac^4
			e1 = e1/n1;
			e2 = e2/n2;
			fprintf(fp, "%le\t%0.30le\t%0.30le\n", star1->rad(X)/R, e1, e2);	
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
		fprintf(gnuplot, ",  '%s' u 1:2 w l t 'Euler, n=%1.1f'", filename, nn[i]);
		fprintf(gnuplot, ",  '%s' u 1:3 w l t 'Poisson, n=%1.1f'", filename, nn[i]);
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);


	return 0;
}