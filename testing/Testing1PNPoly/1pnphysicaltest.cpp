// *************************************************************************************
//					POST-NEWTONIAN POLYTROPE RMSR TEST
// 1pnphysicaltest.cpp
//		The 1PN Lane-Emden equation solves for the variable θ, related to ρ,P
//		This takes the solution of θ, calculates ρ,P, and puts into hydrostatic equations
//			and computes a residual of those equations with these variables
//		This same process is performed by the SSR() method of each Star,
//			However, this program calculates and graphs the residual at each grid point
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/PNPolytrope.cpp"


//data for stars
double n[10] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};

int main(){
	//ask user for tolerance of calculation
	int length = 10000;
	
	//initialize a series of polytropes
	PNPolytrope *star;
	
	double checkEuler[10];
	double checkPoiss[10];
	double checkPsi[10];
	
	FILE *fp;
	
	for(int i=0; i<10; i++){
		//calculate GR polytrope
		star = new PNPolytrope(n[i], 1.e-4, length);
		int L = star->length();
		double R = star->Radius();
		
		//open file for printing
		char filename[100];
		sprintf(filename, "./1pnphysicaltest/test_%1.1f.txt", n[i]);
		if(!(fp = fopen(filename, "w")) ){
			system("mkdir ./1pnphysicaltest");
			fp = fopen(filename, "w");
		}
		
		//sum up errors in equations
		checkPoiss[i] = 0.0;
		checkEuler[i] = 0.0;
		checkPsi[i] = 0.0;
		double r = 0.0;
		double d2Phi;
		double d2Psi;
		double e1, e2, e3, n1, n2, n3;
		for(int X=4; X<L-5; X++){
			r = star->rad(X);
			//Poisson equation for Phi
			double  h=star->rad(X+1)-star->rad(X);
			double 
					b3=star->dPhidr(X-3),
					b2=star->dPhidr(X-2),
					b1=star->dPhidr(X-1),
					a1=star->dPhidr(X+1),
					a2=star->dPhidr(X+2),
					a3=star->dPhidr(X+3);
			d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h);
			e1 = fabs( 4.*m_pi*star->rho(X)*r - 2.*star->dPhidr(X) - d2Phi*r );
			n1 = fabs( 4.*m_pi*star->rho(X)*r ) + fabs( 2.*star->dPhidr(X) ) + fabs(d2Phi*r);
			//Euler equation
			e2  = fabs( star->dPdr(X)
						+ star->dPhidr(X)*(star->rho(X)+star->P(X))
						+ star->rho(X)*star->dPsidr(X) );
			n2 = fabs(star->dPdr(X))
					+ fabs(star->dPhidr(X)*(star->rho(X)+star->P(X)))
					+ fabs(star->rho(X)*star->dPsidr(X));
			// Equation for Psi
			b3=star->dPsidr(X-3);
			b2=star->dPsidr(X-2);
			b1=star->dPsidr(X-1);
			a1=star->dPsidr(X+1);
			a2=star->dPsidr(X+2);
			a3=star->dPsidr(X+3);
			d2Psi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h);
			e3 = fabs( 12.*m_pi*star->P(X) - 8.*m_pi*star->rho(X)*star->Phi(X) - 2.*star->dPsidr(X)/r - d2Psi ); 
			n3 = fabs( 12.*m_pi*star->P(X))
					+ fabs(8.*m_pi*star->rho(X)*star->Phi(X))
					+ fabs(2.*star->dPsidr(X)/r) + fabs(d2Psi);
			e1 /= n1;
			e2 /= n2;
			e3 /= n3;
			checkPoiss[i] += e1;
			checkEuler[i] += e2;
			checkPsi[i]   += e3;
			
			fprintf(fp, "%le\t%le\t%le\t%le\n", star->rad(X)/R, log10(e1), log10(e2), log10(e3));		
		}
		printf("-----------------------------------------------------\n");
		printf("Sum of Errors in Poisson Equation:\t%le\n", checkPoiss[i]);
		printf("Sum of Errors in Psi Source Equation:\t%le\n", checkPsi[i]);
		printf("Sum of Errors in Euler Equation:\t%le\n", checkEuler[i]);
		
		fclose(fp);
		delete star;
		
		char outname[100];
		sprintf(outname, "./1pnphysicaltest/test_%1.1f.png", n[i]);
		FILE *gnuplot = popen("gnuplot -persist", "w");
		fprintf(gnuplot, "reset\n");
		fprintf(gnuplot, "set term png size 1000,1000\n");
		fprintf(gnuplot, "set samples %d\n", L);
		fprintf(gnuplot, "set output '%s'\n", outname);
		fprintf(gnuplot, "set xlabel 'radius'\n");
		fprintf(gnuplot, "set ylabel 'log10 error'\n");
		fprintf(gnuplot, "set title 'Errors in Physical Equations, n=%f'\n", n[i]);
		fprintf(gnuplot, "plot '%s' u 1:2 w l t 'Poisson'", filename);
		fprintf(gnuplot, ",    '%s' u 1:3 w l t 'Euler'", filename);
		fprintf(gnuplot, ",    '%s' u 1:4 w l t 'Psi'", filename);
		fprintf(gnuplot, "\n");
		pclose(gnuplot);		
	}
	return 0;
}