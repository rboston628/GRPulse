// *************************************************************************************
//					PEKERIS MODE FREQUENCY EXACT SCALE TEST
// pekeris_scale.cpp
//		In the Newtonian case of a uniform star (n=0 polytrope)
//		We can find:
//			1. an exact background solution, theta = 1 - x^2
//			2. an exact formula for the dimensionless frequency, due to Pekeris
//			3. an asymptotic series solution for eigenmodes, due to Pekeris
//		See Pekeris (1938), or Shapiro & Teukoslky (1983) sec. 6.5, Cox (1980) sec. 17.7
//		We use the special stellar model Isopycnic.cpp as a background and test
//		We are testing how these errors scale with N_{star}, size of stellar grid
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"
#include "../../src/STARS/Isopycnic.cpp"
#include "../../src/MODES/NonradialModeDriver.cpp"
#include "../../src/MODES/Mode.cpp"

// Formula for frequencies of constant-density star, due to Pekeris (1938)
double pekeris_formula(int k, int l, double Gamma1){
	long double Dkl = Gamma1*double(k*k+k*l) + 0.5*Gamma1*double(k) - 2.0;
	return Dkl + sqrt( Dkl*Dkl + double(l*l+l) );
}

int main(){
	int LEN[10] = {50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000};
	
	//initialize a polytrope to serve as equilibrium model
	Isopycnic *unif_star;
//	for(int i=0; i<10; i++) star[i] = new Isopycnic(LEN[i]);
	
	//now initialize mode drivers based to create modes for each background
	double Gamma1 = 5./3.;
	NonradialModeDriver *unif_drv ;
	const unsigned int num_var = NonradialModeDriver::num_var;
	Mode<num_var>* unif_mode;
	
	
	FILE *fp;
	char filename[100];
	char outname[100];
	system("rm -r ./pekeris_scale");
	system("mkdir ./pekeris_scale");


	int K=0;
	double w2eqn;

	
	//we will calculate frequencies and compare against the Pekeris formula
	// these arrays will hold the errors
	const int Lmin=1, Lmax=3;
	const int Kmin=1, Kmax=9;
	
	double unif_omeg[(Lmax-Lmin)+1][(Kmax-Kmin)+1][10] = {0.0};
	double unif_errors[(Lmax-Lmin)+1][(Kmax-Kmin)+1][10] = {0.0};
	
	for(int i=0; i<10; i++){
		printf("------------------------------------------\n");
		printf("grid = %d\t", LEN[i]);fflush(stdout);
		unif_star = new Isopycnic(LEN[i]);
		unif_drv  = new NonradialModeDriver(unif_star, Gamma1);
		for(int l=0; l<=(Lmax-Lmin); l++){
			for(int k=0; k<=(Kmax-Kmin); k++){	
				printf("\tL=%d, K=%d\t", l+Lmin,k+Kmin);
				w2eqn = pekeris_formula(k+Kmin, l+Lmin, Gamma1);
				unif_mode = new Mode<num_var>(w2eqn,l+Lmin,0,unif_drv);
				K = unif_mode->modeOrder();
				printf("found K=%d", K);
				if(Kmin<=K && K<=Kmax){
					unif_omeg[l][K-Kmin][i] = unif_mode->getOmega2();
					w2eqn = pekeris_formula(K, l+Lmin, Gamma1);
					unif_errors[l][K-Kmin][i] = fabs(unif_omeg[l][K-Kmin][i]-w2eqn)/w2eqn;
					printf("\tw=%le %le %le\n",unif_omeg[l][K-Kmin][i],w2eqn,unif_errors[l][K-Kmin][i]);
				}
				else printf("\tNOPE!\n");
				delete unif_mode;
			}
		}
		delete unif_star;
		delete unif_drv;
	}
	for(int l=0; l<=(Lmax-Lmin); l++){
		sprintf(filename, "./pekeris_scale/scale_L%d.txt", l+Lmin);
		sprintf(outname,  "./pekeris_scale/scale_L%d.png", l+Lmin);
		fp = fopen(filename, "w");
		for(int i=0; i<10; i++){
			fprintf(fp, "%d", LEN[i]);
			for(int k=0; k<=(Kmax-Kmin); k++){
				//ignore outliers -- caused by mode misidentification
				if(unif_errors[l][k][i] >1.0) unif_errors[l][k][i] = nan("");
				fprintf(fp, "\t%0.32le", unif_errors[l][k][i]);
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
		fprintf(gnuplot, "set title 'Scaling of Errors in Isopycnic Frequencies, L=%d'\n", l+Lmin);
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
		
	return 0;
}
