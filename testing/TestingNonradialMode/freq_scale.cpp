// *************************************************************************************
//					MODE FREQUENCY ERROR SCALE TEST
// freq_scale.cpp
//		Here we observe how the errors in a eigenfrequency scale with N_{star}
//		We first compute w2 using a very large grid, considered the "true" value
//		We then calculate on gradually smaller grids and compare
//		We perform this test for several different polytropes
//		THIS WILL TAKE A LONG TIME TO RUN
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"
#include "../../src/MODES/NonradialModeDriver.cpp"
#include "../../src/MODES/Mode.cpp"

// Formula for frequencies of constant-density star, due to Pekeris (1938)
double pekeris_formula(int k, int l, double Gamma1){
	long double Dkl = Gamma1*double(k*k+k*l) + 0.5*Gamma1*double(k) - 2.0;
	return Dkl + sqrt( Dkl*Dkl + double(l*l+l) );
}

int main(){
	int LEN[10] = {1000000, 500000, 100000, 50000, 10000, 5000, 1000, 500, 100, 50};
	double nn[4] = {1.0, 1.5, 2.0, 3.0};
	
	//initialize a polytrope to serve as equilibrium model
	Polytrope *star;
//	for(int i=0; i<10; i++) star[i] = new Isopycnic(LEN[i]);
	
	//now initialize mode drivers based to create modes for each background
	double Gamma1 = 5./3.;
	NonradialModeDriver *drv ;
	const unsigned int num_var = NonradialModeDriver::num_var;
	Mode<num_var>* mode;
	
	
	FILE *fp;
	char filename[100];
	char outname[100];
	system("rm -r ./freq_scale");
	system("mkdir ./freq_scale");


	int K=0;
	double w2eqn;

	//we will calculate frequencies and compare against the Pekeris formula
	// these arrays will hold the errors
	const int Lmin=1, Lmax=3;
	const int Kmin=1, Kmax=9;
	
	//for frequencies for using fewer grid points
	double omeg[4][(Lmax-Lmin)+1][(Kmax-Kmin)+1][10] = {0.0};
	double errors[4][(Lmax-Lmin)+1][(Kmax-Kmin)+1][10] = {0.0};
	
	//first get the frequency calculated with many grid points, taken as "true" value
	double omeg0[4][(Lmax-Lmin)+1][(Kmax-Kmin)+1];
	printf("------------------------------------------\n");
	printf("  THIS WILL TAKE A VERY LONG TIME \n");
	printf("  GO GET A COFFEE OR SOMETHING\n");
	for(int m=0; m<4; m++){
		printf("------------------------------------------\n");
		printf("Initial frequencies for polytrope N=%1.1f\n", nn[m]);
		star = new Polytrope(nn[m], LEN[0]);
		drv  = new NonradialModeDriver(star, Gamma1);
		for(int l=0; l<=(Lmax-Lmin); l++){
			for(int k=0; k<=(Kmax-Kmin); k++){
				printf("\tL=%d,K=%d\t", l+Lmin, k+Kmin);fflush(stdout);
				mode = new Mode<num_var>(k+Kmin, l+Lmin, 0, drv);
				omeg0[m][l][k] = mode->getOmega2();
				omeg[m][k][k][0] = omeg0[m][l][k];
				printf("w2=%le\n", omeg0[m][l][k]);
				delete mode;
			}
		}
		delete star;
		delete drv;
	}
	
	

	
	for(int m=0; m<4; m++){		
		for(int i=1; i<10; i++){
			printf("------------------------------------------\n");
			printf("grid = %d\t", LEN[i]); fflush(stdout);
			star = new Polytrope(nn[m], LEN[i]);
			drv  = new NonradialModeDriver(star, Gamma1);
			for(int l=0; l<=(Lmax-Lmin); l++){
				for(int k=0; k<=(Kmax-Kmin); k++){	
					printf("\tL=%d, K=%d\t", l+Lmin,k+Kmin);
					mode = new Mode<num_var>(omeg0[m][l][k],l+Lmin,0,drv);
					K = mode->modeOrder();
					printf("found K=%d", K);
					if(K==k+Kmin){
						omeg[m][l][k][i] = mode->getOmega2();
						errors[m][l][k][i] = fabs(omeg[m][l][k][i]-omeg0[m][l][k])/omeg0[m][l][k];
						printf("\tw=%le %le %le\n",omeg[m][l][k][i],omeg0[m][l][k],errors[m][l][k][i]);
					}
					else printf("\tNOPE!\n");
					delete mode;
				}
			}
			delete star;
			delete drv;
		}
	}
	for(int m=0; m<4; m++){
		char command[256];
		sprintf(command, "mkdir ./freq_scale/polytrope_%1.1f", nn[m]);
		system(command);
		for(int l=0; l<=(Lmax-Lmin); l++){
			sprintf(filename, "./freq_scale/polytrope_%1.1f/scale_L%d.txt", nn[m], l+Lmin);
			sprintf(outname,  "./freq_scale/polytrope_%1.1f/scale_L%d.png", nn[m], l+Lmin);
			fp = fopen(filename, "w");
			for(int i=1; i<10; i++){
				fprintf(fp, "%d", LEN[i]);
				for(int k=0; k<=(Kmax-Kmin); k++){
					//ignore outliers -- caused by mode misidentification
					if(errors[m][l][k][i] >1.0) errors[m][l][k][i] = nan("");
					fprintf(fp, "\t%0.32le", errors[m][l][k][i]);
					fflush(fp);
				}
				fprintf(fp,"\n");
			}
			fclose(fp);
			//now plot the scaling of different K for a given L
			FILE *gnuplot = popen("gnuplot -persist", "w");
			fprintf(gnuplot, "reset\n");
			fprintf(gnuplot, "set term png size 1200,1000\n");
			fprintf(gnuplot, "set output '%s'\n", outname);
			fprintf(gnuplot, "set title 'Scaling of Errors in Polytrope n=%1.1f Frequencies, L=%d'\n", nn[m], l+Lmin);
			//fix the x-axis
			fprintf(gnuplot, "set xlabel 'N_{star}'\n");
			fprintf(gnuplot, "set logscale x\n");
			fprintf(gnuplot, "set format x '10^{%%L}'\n");
			//fix the y-axis
			fprintf(gnuplot, "set ylabel 'relative error in ω'\n");
			fprintf(gnuplot, "set logscale y\n");
			fprintf(gnuplot, "set format y '10^{%%L}'\n");
			//draw a line for the resired slope -4
			fprintf(gnuplot, "set arrow from 1e2,1e-3 to 1e5,1e-15 nohead lw 3\n");
			//make the filename to read
			fprintf(gnuplot, "plot ");
			for(int k=0; k<=(Kmax-Kmin); k++){
				fprintf(gnuplot, "%c '%s' u 1:%d w lp t 'k=%d'", (k==0?' ':','), filename, k+2, k+Kmin);
			}
			fprintf(gnuplot, "\n");
			pclose(gnuplot);		
		}
	}
		
	return 0;
}
