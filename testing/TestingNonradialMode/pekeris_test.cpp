// *************************************************************************************
//					PEKERIS MODE FREQUENCY EXACT TEST
// pekeris_test.cpp
//		In the Newtonian case of a uniform star (n=0 polytrope)
//		We can find:
//			1. an exact background solution, theta = 1 - x^2
//			2. an exact formula for the dimensionless frequency, due to Pekeris
//			3. an asymptotic series solution for eigenmodes, due to Pekeris
//		We use the special stellar model Isopycnic.cpp as a background and test
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"
#include "../../src/STARS/Isopycnic.cpp"
#include "../../src/MODES/NonradialModeDriver.cpp"
#include "../../src/MODES/Mode.cpp"


// Formula for frequencies of constant-density star, due to Pekeris (1938)
//  See also Shapiro & Teukolksi for a more modern treatment
double pekeris_formula(int k, int l, double Gamma1){
	long double Dkl = Gamma1*double(k*k+k*l) + 0.5*Gamma1*double(k) - 2.0;
	return Dkl + sqrt( Dkl*Dkl + double(l*l+l) );
}

int main(){
	//initialize the equilibrium background models
	// we will use a polytrope with n=0
	// we will also use a special uniform-density model using the exact n=0 solution
	int length = 100000;
	Isopycnic *unif_star = new Isopycnic(length);
	Polytrope *poly_star = new Polytrope(0.0, length);
	
	//now initialize mode drivers based to create modes for each background
	double Gamma1 = 5./3.;
	NonradialModeDriver *unif_drv = new NonradialModeDriver(unif_star, Gamma1);
	NonradialModeDriver *poly_drv = new NonradialModeDriver(poly_star, Gamma1);
	const unsigned int num_var = NonradialModeDriver::num_var;
	Mode<num_var>* unif_mode;
	Mode<num_var>* poly_mode;
	
	int K=0;
	double w2eqn;
	
	//for both models, we will calculate frequencies and compare against the Pekeris formula
	// these arrays will hold the errors
	const int Lmin=1, Lmax=3;
	const int Kmin=1, Kmax=15;
	double unif_errors[(Lmax-Lmin)+1][(Kmax-Kmin)+1] = {0.0};
	double poly_errors[(Lmax-Lmin)+1][(Kmax-Kmin)+1] = {0.0};
	double unif_omeg[(Lmax-Lmin)+1][(Kmax-Kmin)+1] = {0.0};
	double poly_omeg[(Lmax-Lmin)+1][(Kmax-Kmin)+1] = {0.0};
	
	for(int l=0; l<=(Lmax-Lmin); l++) for(int k=0;k<=(Kmax-Kmin);k++) unif_errors[l][k]=poly_errors[l][k] = nan("");
	printf("beginning\n");
	for(int l=0; l<=(Lmax-Lmin); l++){
		for(int k=0; k<=(Kmax-Kmin); k++){
			printf("Looking for mode l=%d, k=%d\n", l+Lmin, k+Kmin);
			printf("\tuniform star mode:\n");
			w2eqn = pekeris_formula(k+Kmin, l+Lmin, Gamma1);
			unif_mode = new Mode<num_var>(w2eqn, l+Lmin,0, unif_drv);
			K = unif_mode->modeOrder();
			printf("\tfound mode l=%d,k=%d, w = %le", l+Lmin, K, unif_mode->getOmega2());
			if(Kmin <= K && K <= Kmax){
				unif_omeg[l][K-Kmin] = unif_mode->getOmega2();
				w2eqn = pekeris_formula(K, l+Lmin, Gamma1);
				printf(" %le", w2eqn);
				unif_errors[l][K-Kmin] = fabs(unif_omeg[l][K-Kmin]-w2eqn)/w2eqn;
			}
			printf("\n");
			
			printf("\tpolytropic star mode:\n");
			poly_mode = new Mode<num_var>(w2eqn, l+Lmin,0, poly_drv);
			K = poly_mode->modeOrder();
			printf("\tfound mode l=%d,k=%d, w = %le", l+Lmin, K, poly_mode->getOmega2());
			if(Kmin <= K && K <= Kmax){
				poly_omeg[l][K-Kmin] = poly_mode->getOmega2();
				w2eqn = pekeris_formula(K, l+Lmin, Gamma1);
				printf(" %le", w2eqn);
				poly_errors[l][K-Kmin] = fabs(poly_omeg[l][K-Kmin]-w2eqn)/w2eqn;
			}
			printf("\n");
			
			//remove outliers -- sometimes mode identification is wrong
			if(unif_errors[l][K-Kmin] > 1.0) unif_errors[l][K-Kmin] = nan("");
			if(poly_errors[l][K-Kmin] > 1.0) poly_errors[l][K-Kmin] = nan("");
			
			delete unif_mode;
			delete poly_mode;
		}
	}
	delete unif_star;
	delete poly_star;
	delete unif_drv;
	delete poly_drv;
	
	FILE *fp;
	char filename[91];
	sprintf(filename, "./pekeris_test/pekeris_test.txt");
	if(!(fp = fopen(filename, "w")) ){
		system("mkdir ./pekeris_test");
		fp = fopen(filename, "w");
	}
	for(int k=0; k<=(Kmax-Kmin); k++){
		fprintf(fp, "%d", k+Kmin);
		for(int l=0; l<=(Lmax-Lmin); l++)
			fprintf(fp, "\t%12le",unif_errors[l][k]);
		for(int l=0; l<=(Lmax-Lmin); l++)
			fprintf(fp, "\t%12le",poly_errors[l][k]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 2000,1000\n");
	fprintf(gnuplot, "set output './pekeris_test/pekeris_test.png'\n");
	fprintf(gnuplot, "set xlabel 'K'\n");
	int thou = int(length/1000);
	fprintf(gnuplot, "set ylabel 'Relative Error for Mode L,K'\n");
	fprintf(gnuplot, "set logscale y %d\n", 10);
	fprintf(gnuplot, "set format y '%d^{%%L}'\n", 10);
	fprintf(gnuplot, "set ytics %d\n", 10);

	fprintf(gnuplot, "set title 'Relative Errors in Frequencies for Various L,K, using N_{star}=%dK'\n", thou);
	fprintf(gnuplot, "plot ");
	for(int l=0; l<=(Lmax-Lmin); l++)
		fprintf(gnuplot, "%c '%s' u 1:%d w lp lw 2 t 'exact, L=%d'", (l==0?' ':','),filename, l+2, l+Lmin);
	for(int l=0; l<=(Lmax-Lmin); l++)
		fprintf(gnuplot, ", '%s' u 1:%d w lp t 'polytrope, L=%d'", filename, l+3+(Lmax-Lmin), l+Lmin);
	fprintf(gnuplot, "\n");
	pclose(gnuplot);//*/			

	return 0;
}
