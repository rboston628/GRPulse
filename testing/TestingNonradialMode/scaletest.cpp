// *************************************************************************************
//					NEWTONIAN MODE ERROR SCALE TEST
// scaletest.cpp
//		This program checks that the errors scale correctly with grid size
//		Suppose for given grid size N, the numerical solution from RK4 is 
//			y(N) = y(infty) + E(N^4)
//		Then y(2N)-y(N) = E(2^4N^4) - E(N^4) ~ 2^4, up to proportionality constant
//		If we construct fraction [y(2N)-y(N)]/[y(4N)-y(2N)] = 2^4
//			where the proportionality constant cancels.
//		The graphs produced look crazy, but do tend around the line 2^4
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"
#include "../../src/MODES/NonradialModeDriver.cpp"
#include "../../src/MODES/Mode.cpp"


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
	
	//adiabatic exponent
	double Gamma1 = 5./3.;
	const int Lmin = 2, Lmax = 3;
	const int Kmin = 1, Kmax = 9;
	
	//initialize a series of polytropes
	Polytrope *star0, *star1, *star2;
	NonradialModeDriver  *driv0, *driv1, *driv2;
	const unsigned int num_var = NonradialModeDriver::num_var;
	Mode<num_var> *mode0, *mode1, *mode2;
		
	FILE *fp;
	char filename[256];	
	system("mkdir ./scaletest/");
	
	for(int i=0; i<7; i++){
		printf("--------------------------------------\n");
		printf("Beginning polytrope n=%1.1f\n", n[i]);
		
		//begin by finding three polytropes whose grids are suitably scaled
		star2 = new Polytrope(n[i], length2);
		double dx = star2->getX(2)-star2->getX(1);
		delete star2;
		star0 = new Polytrope(n[i], length0, fac2*dx);
		star1 = new Polytrope(n[i], length1, fac1*dx);
		star2 = new Polytrope(n[i], length2, dx);
		driv0 = new NonradialModeDriver(star0, Gamma1);
		driv1 = new NonradialModeDriver(star1, Gamma1);
		driv2 = new NonradialModeDriver(star2, Gamma1);
				int Len = driv0->length();
		char mkdir[100];
		sprintf(mkdir, "rm -r ./scaletest/polytrope_%1.1f/", n[i]);
		system(mkdir);
		sprintf(mkdir, "mkdir -p ./scaletest/polytrope_%1.1f/", n[i]);
		system(mkdir);
						
		//no creator 
		for(int l=Lmin; l<=Lmax; l++){
			printf("L=%d\n", l);
			for(int k=Kmin; k<=Kmax; k++){
				mode0 = new Mode<num_var>(k,l,l, driv0);
				mode1 = new Mode<num_var>(k,l,l, driv1);
				mode2 = new Mode<num_var>(k,l,l, driv2);
				int K0 = mode0->modeOrder(), K1=mode1->modeOrder(), K2=mode2->modeOrder();
				printf("\tK0 %d\tK1 %d\tK2 %d\n", K0,K1,K2);
				printf(" w0 %le\tw1 %le\tw2 %le\t %le\n", mode0->getOmega2(), mode1->getOmega2(), mode2->getOmega2(),
						(mode1->getOmega2()-mode0->getOmega2())/(mode2->getOmega2()-mode1->getOmega2())	);
				if(K1==K2 & K2==K0 & K0>0){
					sprintf(filename, "./scaletest/polytrope_%1.1f/scale_%1.1f_L%d_K%d.txt", n[i],n[i],l,K0);
					fp = fopen(filename,"w");

					double scale1,scale2,scale3,scale4,scaleup,scaledown;
					for(int X=0; X<Len-2; X++){
						scale1 = fabs(mode1->getYtilde(0,fac1*X)-mode0->getYtilde(0,fac0*X))/fabs(mode2->getYtilde(0,fac2*X)-mode1->getYtilde(0,fac1*X));
						scale2 = fabs(mode1->getYtilde(1,fac1*X)-mode0->getYtilde(1,fac0*X))/fabs(mode2->getYtilde(1,fac2*X)-mode1->getYtilde(1,fac1*X));
						scale3 = fabs(mode1->getYtilde(2,fac1*X)-mode0->getYtilde(2,fac0*X))/fabs(mode2->getYtilde(2,fac2*X)-mode1->getYtilde(2,fac1*X));
						scale4 = fabs(mode1->getYtilde(3,fac1*X)-mode0->getYtilde(3,fac0*X))/fabs(mode2->getYtilde(3,fac2*X)-mode1->getYtilde(3,fac1*X));
						scaleup = fabs( mode1->getY(0,fac1*X)-mode0->getY(0,fac0*X) );
						scaledown=fabs( mode2->getY(0,fac2*X)-mode1->getY(0,fac1*X) );
						fprintf(fp, "%0.32le\t%0.32le\t%0.32le\t%0.32le\t%0.32le\t%0.32le\t%0.32le\n", mode0->getRad(X),scale1, scale2, scale3, scale4, scaleup, scaledown);
					}
					fclose(fp);
				}
				delete mode0;
				delete mode1;
				delete mode2;	
			}
		}
		delete star0;
		delete star1;
		delete star2;
		
		char outname[100];
		sprintf(outname, "./scaletest/polytrope_%1.1f/scale_%1.1f.png", n[i],n[i]);
		FILE *gnuplot = popen("gnuplot -persist", "w");
		fprintf(gnuplot, "reset\n");
		fprintf(gnuplot, "set term png size 2000,1000\n");
		fprintf(gnuplot, "set samples %d\n", Len);
		fprintf(gnuplot, "set output '%s'\n", outname);
		fprintf(gnuplot, "set title 'Scaling Grid by Factor %d in polytrope n=%1.1f'\n", fac1, n[i]);
		int thou = int(length0/1024);
		fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
		fprintf(gnuplot, "set xrange [0:1.01]\n");
		fprintf(gnuplot, "set ylabel 'Ratio |[y_{%dK}-y_{%dK}]/[y_{%dK}-y_{%dK}]|'\n", fac1*thou, thou, fac2*thou, fac1*thou);
		fprintf(gnuplot, "set logscale y %d\n", fac1);
		fprintf(gnuplot, "set ytics %d\n", fac1);
		fprintf(gnuplot, "set format y '%d^{%%L}'\n", fac1);
		fprintf(gnuplot, "set yrange [%le:%le]\n", pow(fac1,-10), pow(fac1,6));
		fprintf(gnuplot, "set arrow 1 from 0.0,%le to 1.01, %le lw 4 lc rgb 'red' nohead\n", pow(fac1,4), pow(fac1,4));
		fprintf(gnuplot, "plot ");
		for(int l=Lmin; l<=Lmax; l++){
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<num_var; j++){
					if(l!=Lmin | k!=Kmin | j!=0) fprintf(gnuplot, ", ");
					sprintf(filename, "./scaletest/polytrope_%1.1f/scale_%1.1f_L%d_K%d.txt", n[i], n[i], l,k);
					fprintf(gnuplot, "'%s' u 1:%d w l t 'y_{%d}^{%d,%d}'", filename,  j+2, j+1, l, k);
				}
			}	
		}
		fprintf(gnuplot, "\n");
		pclose(gnuplot);
		

	}
	return 0;
}