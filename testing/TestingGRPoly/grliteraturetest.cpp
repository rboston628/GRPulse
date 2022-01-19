// *************************************************************************************
//					RELATIVISTIC POLYTROPE LITERATURE TEST
// grliteraturetest.cpp
//		In Tooper (1964) there are published tables of values for a GR polytrope 
//			for several n and sigma, showing radii, central densities, etc.,
//			along with several graphs of results
//		These tables were partly extended by Bludman (1973)
//		This program will generate several relativistic polytropes, and compare to 
//			these published literature results.
//		Output will be:
//			1. a table mimicking Table 1 of Tooper
//			2. as above, but the entires are differences
//			3. a graph of the differences
//			3. a set of graphs corresponding to Figs 1-5 of Tooper
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/GRPolytrope.cpp"

void printGraph(GRPolytrope*, double n, double s);

int main(){
	//Tooper computed tabulated values for polytropes n=1,...,3, sigma = 0,...,sigmax
	const int nn = 5;
	//the number of sigmas to test
	const int sn = 9;
	const double n[nn] = {1.0, 1.5, 2.0, 2.5, 3.0};
	//stars with these particular values of sigma are given in plots
	const double ngraph[nn] = {1.0, 1.5, 2.0, 2.5, 3.0};
	const double sgraph[nn] = {0.2, 0.2, 0.2, 0.3, 0.3};
	//tolerance of calculation
	const int D = 3;
	const double tol = pow(10.0, -D);
	
	//copy of table 1 from Tooper (1964)
	//	each array is a different variable
	//  each row a different polytrope
	//  each column correspponds to a different value of sigma
	const double xi1Tooper[nn][sn] = {
		{3.1416, 2.599, 2.277, 2.064, 1.913, 1.801, 0.000, 0.000, 0.000},//n=1.0
		{3.6538, 3.088, 2.699, 2.493, 2.361, 2.275, 2.219, 0.000 ,0.000},//n=1.5
		{4.3529, 3.699, 3.398, 3.271, 3.248, 3.296, 3.399, 3.493, 0.000},//n=2.0
		{5.3553, 4.782, 4.721, 4.986, 5.545, 6.434, 7.728, 9.523, 9.985},//n=2.5
		{6.8968, 6.826, 7.951, 10.83, 17.83, 37.21, 91.08, 162.6, 180.5} //n=3.0
	};
	const double xibar1Tooper[nn][sn] = {
		{3.1416, 2.853, 2.657, 2.517, 2.412, 2.330, 0.000, 0.000, 0.000},//n=1.0
		{3.6538, 3.381, 3.219, 3.120, 3.061, 3.029, 3.018, 0.000, 0.000},//n=1.5
		{4.3529, 4.150, 4.093, 4.128, 4.228, 4.381, 4.577, 4.730, 0.000},//n=2.0
		{5.3553, 5.372, 5.657, 6.183, 6.976, 8.095, 9.632, 11.69, 12.03},//n=2.5
		{6.8968, 7.606, 9.261, 12.65, 20.26, 40.55, 95.83, 169.0, 187.6} //n=3.0
	};
	const double v1Tooper[nn][sn] = {
		{3.1416, 1.751, 1.143, .8192, .6249, .4981, 0.000, 0.000, 0.000},//n=1.0
		{2.7141, 1.482, .9604, .6883, .5270, .4226, .3505, 0.000, 0.000},//n=1.5
		{2.4110, 1.299, .8402, .6055, .4680, .3800, .3201, .2901, 0.000},//n=2.0
		{2.1872, 1.169, .7606, .5556, .4386, .3664, .3202, .2904, .2872},//n=2.5
		{2.0182, 1.078, .7130, .5386, .4516, .4214, .4493, .5266, .5654} //n=3.0
	};
	const double rhoc2avgTooper[nn][sn] = {
		{3.2899, 3.342, 3.443, 3.579, 3.736, 3.909, 0.000, 0.000, 0.000},//n=1.0
		{5.9907, 6.310, 6.827, 7.504, 8.326, 9.286, 10.394, 0.000, 0.000},//n=1.5
		{11.403, 12.99, 15.57, 19.27, 24.41, 31.42, 40.90, 48.87, 0.000},//n=2.0
		{23.406, 31.18, 46.10, 74.35, 129.6, 242.3, 480.4, 991.2, 1101.},//n=2.5
		{54.18, 98.35, 235.0, 786.8, 4177., 40750., 560600., 2722000., 3465000.} //n=3.0
	};
	const double OmegaTooper[nn][sn] = {
		{0.0000, .1240, .2108, .2744, .3228, .3605, 0.000, 0.000, 0.000},//n=1.0
		{0.0000, .1307, .2207, .2857, .3344, .3719, .4014, 0.000, 0.000},//n=1.5
		{0.0000, .1355, .2274, .2926, .3405, .3764, .4036, .4181, 0.000},//n=2.0
		{0.0000, .1388, .2312, .2950, .3396, .3702, .3900, .4008, .4018},//n=2.5
		{0.0000, .1408, .2318, .2900, .3222, .3273, .3000, .2582, .2433} //n=3.0
	};
		
	//initialize a series of polytropes
	double xi1[nn][sn], xibar1[nn][sn], v1[nn][sn], rhoc2avg[nn][sn];
	int LEN = int(512*sqrt(6)); //approximately the grid size used by Tooper (see Sec IVa)
	double dx = 1./double(512); //the grid spacing used by Tooper
	GRPolytrope *star;
	
	FILE *fp1, *fp2, *fp3;
	if(!(fp1 = fopen("./littest/duplicate_table1.txt", "w"))){
		system("mkdir -p littest");
	}	
	fp2 = fopen("./littest/difs_table1.txt", "w");
	char filename[256];
	for(int k=0; k<nn; k++){
		double ntemp = n[k];				//current index n
		double gitemp = ntemp/(ntemp+1.0);	//inverse of Gamma of current index n
		//printf("n=%f\t\txi1\txibar1\tv1\trhoc/rhoavg\n", ntemp);
		fprintf(fp1, "n=%0.2f\t        \txi1     \txibar1  \tv1      \trhoc/rhoavg\n", ntemp);
		fprintf(fp2, "n=%0.2f\t        \txi1     \txibar1  \tv1      \trhoc/rhoavg\n", ntemp);
		sprintf(filename, "./littest/poly_%1.1f.txt", ntemp);
		fp3 = fopen(filename, "w");
		//for stars with this index, test each sigma in range [0, 1/Gamma]
		int s = 0;
		for(double sigma=0.0; sigma < gitemp; sigma += 0.1){
			//make a new star
			LEN = int(512*xi1Tooper[k][s])+1;
			star = new GRPolytrope(ntemp, sigma, LEN, dx);
			xi1[k][s] = star->getX(LEN-1);
			xibar1[k][s] = star->getXproper(LEN-1);
			v1[k][s] = star->getV(LEN-1);
			rhoc2avg[k][s]  = pow(xi1[k][s],3)/( 3.0*v1[k][s] );
			
			//for values of sigma with graph in Tooper, print a graph for comparison
			if(int(sigma*10) == int(sgraph[k]*10)){
				//printf("%lf\t%lf\n", sigma, sgraph[k]);
				printGraph(star, n[k], sigma);	
			}
			//print results  to file
			//duplicate table
			fprintf(fp1, "\tsig=%.2f \t%4.4le\t%4.4le", sigma, xi1[k][s], xibar1[k][s]);
			fprintf(fp1, "\t%4.4le\t%4.4le\n", v1[k][s], rhoc2avg[k][s]);
			//differences		
			fprintf(fp2, "\tsig=%.2f ", sigma);
			fprintf(fp2, "\t%le", fabs(xi1[k][s]-xi1Tooper[k][s])/xi1Tooper[k][s]);
			fprintf(fp2, "\t%le", fabs(xibar1[k][s]-xibar1Tooper[k][s])/xibar1Tooper[k][s]);
			fprintf(fp2, "\t%le", fabs(v1[k][s]-v1Tooper[k][s])/v1Tooper[k][s]);
			fprintf(fp2, "\t%le", fabs(rhoc2avg[k][s]-rhoc2avgTooper[k][s])/rhoc2avgTooper[k][s]);
			fprintf(fp2, "\n");
			//for graphing
			fprintf(fp3, "%.2f ", sigma);
			fprintf(fp3, "\t%le", fabs(xi1[k][s]-xi1Tooper[k][s])/xi1Tooper[k][s]);
			fprintf(fp3, "\t%le", fabs(xibar1[k][s]-xibar1Tooper[k][s])/xibar1Tooper[k][s]);
			fprintf(fp3, "\t%le", fabs(v1[k][s]-v1Tooper[k][s])/v1Tooper[k][s]);
			fprintf(fp3, "\t%le", fabs(rhoc2avg[k][s]-rhoc2avgTooper[k][s])/rhoc2avgTooper[k][s]);
			fprintf(fp3, "\n");
			delete star;
			s++;
		}
		//in each table Tooper lists a row for the maximum allowed sigma
		LEN = int(512*xi1Tooper[k][s])+1;
		star = new GRPolytrope(ntemp, gitemp, LEN, dx);
		xi1[k][s] = star->getX(LEN-1);
		xibar1[k][s] = star->getXproper(LEN-1);
		v1[k][s] = star->getV(LEN-1);
		rhoc2avg[k][s]  = pow(xi1[k][s],3)/(3.0*v1[k][s]);
		//print results to file
		//duplicate table
		fprintf(fp1, "\tsig=%.2f \t%4.4le\t%4.4le", gitemp, xi1[k][s], xibar1[k][s]);
		fprintf(fp1, "\t%4.4le\t%4.4le\n", v1[k][s], rhoc2avg[k][s]);
		//differences		
		fprintf(fp2, "\tsig=%.2f ", gitemp);
		fprintf(fp2, "\t%le", fabs(xi1[k][s]-xi1Tooper[k][s])/xi1Tooper[k][s]);
		fprintf(fp2, "\t%le", fabs(xibar1[k][s]-xibar1Tooper[k][s])/xibar1Tooper[k][s]);
		fprintf(fp2, "\t%le", fabs(v1[k][s]-v1Tooper[k][s])/v1Tooper[k][s]);
		fprintf(fp2, "\t%le", fabs(rhoc2avg[k][s]-rhoc2avgTooper[k][s])/rhoc2avgTooper[k][s]);
		fprintf(fp2, "\n");
		//for graphing
		fprintf(fp3, "%.2f ", gitemp);
		fprintf(fp3, "\t%le", fabs(xi1[k][s]-xi1Tooper[k][s])/xi1Tooper[k][s]);
		fprintf(fp3, "\t%le", fabs(xibar1[k][s]-xibar1Tooper[k][s])/xibar1Tooper[k][s]);
		fprintf(fp3, "\t%le", fabs(v1[k][s]-v1Tooper[k][s])/v1Tooper[k][s]);
		fprintf(fp3, "\t%le", fabs(rhoc2avg[k][s]-rhoc2avgTooper[k][s])/rhoc2avgTooper[k][s]);
		fprintf(fp3, "\n");
		fclose(fp3);
		delete star;
	}
	fclose(fp1);
	fclose(fp2);
	
	
	//print all errors
	FILE* gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1200,1000\n");
	fprintf(gnuplot, "set output '%s'\n", "./littest/differences.png");
	fprintf(gnuplot, "set title 'Relative Differences from Tooper for Relativistic Polytropes'\n");
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'relativity parameter sigma'\n");
	fprintf(gnuplot, "set yrange [1e-5: 1e0]\n");
	fprintf(gnuplot, "set ylabel 'Relative difference from Tooper (1964)'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set label 'Note: Tooper (1964) lists only four digits' at graph 0.02,0.98\n");
	fprintf(gnuplot, "plot ");
	for(int i=0; i<nn; i++){
		sprintf(filename, "./littest/poly_%1.1f.txt", n[i]);
		fprintf(gnuplot, "%c '%s' u 1:2 w l t 'xi1, n=%1.1f'", (i==0?' ':','), filename, n[i]);
		fprintf(gnuplot, ",  '%s' u 1:3 w l t 'xibar1, n=%1.1f'",     filename, n[i]);
		fprintf(gnuplot, ",  '%s' u 1:4 w l t 'v1, n=%1.1f'",         filename, n[i]);
		fprintf(gnuplot, ",  '%s' u 1:5 w l t 'rhoc/<rho>, n=%1.1f'", filename, n[i]);
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	system("rm littest/poly*");
	
	return 0;
}

void printGraph(GRPolytrope *star, double n, double sigma){
	int L = star->length();
	
	//create names for files to be opened
	char filename[91];
	char outname[91];
	char openmyplot[45];
	//save data to folder to avoid clutter - make sure folder exists
	sprintf(filename, "./littest/star%.1f_%.1f.txt", n, sigma);
	FILE *fp;
	if(!(fp = fopen(filename, "w")) ){
		system("mkdir ./littest");
		fp = fopen(filename, "w");
	}
	//variables to be calculated
	double r, rho, P, m, enu, elam;
	//constants used in calculation
	double rhoc = star->rho(0), Pc = star->P(0), R = star->Radius(), M = star->Mass();
	double np1 = n+1.0, enuc = 1.0 - 2*sigma*np1*star->getV(L-1)/star->getX(L-1);
	//print all variabels to text file
	R = star->getXproper(L-1);
	for(int x=0; x<L; x++){
		r = star->getXproper(x)/R;
		m = star->mr(x)/M;
		rho = star->rho(x)/rhoc;
		P = star->P(x)/Pc;
		elam = std::exp(star->Lambda(x));//1.0/( 1 - 2.0*sigma*np1*star->getV(x)/star->getX(x) );
		enu  = enuc*pow(1. + sigma*star->getY(x), -2*np1);
		fprintf(fp, "%0.16lf\t%0.16lf\t%0.16lf\t%0.16lf\t%0.16lf\t%0.16lf\n", r, rho, P, m, enu, elam);
	}
	//for the purposes of the comparison, we extend the result beyond the stellar radius
	double v1 = star->getV(L-1);
	double dr = (star->getXproper(L-1)-star->getXproper(L-2))/R;
	double dx = star->getX(L-1)-star->getX(L-2);
	double Rn = np1*Pc/(4.*m_pi*rhoc*rhoc);
	double x=star->getX(L-1);
	while(r<1.5){
		r = r + dr;
		x = x + dx;
		elam = 1./(1. - 2.*sigma*np1*v1/x);
		enu  = (1. - 2.*sigma*np1*v1/x);
		fprintf(fp, "%0.16lf\t%0.16lf\t%0.16lf\t%0.16lf\t%0.16lf\t%0.16lf\n", r, 0.0, 0.0, m, enu, elam);
	}
	fclose(fp);	
	//plot file in png in gnuplot
	sprintf(outname, "./littest/star%0.1f_%0.1f.png", n, sigma);
	FILE *gnuplot = popen("gnuplot -persist", "w");
	//-------
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,1000\n");
	fprintf(gnuplot, "set samples %d\n", L);
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set encoding utf8\n");
	fprintf(gnuplot, "set title 'n=%0.1f, sigma=%0.1f'\n", n, sigma);
	fprintf(gnuplot, "set xrange [0:1.5]\n");
	fprintf(gnuplot, "set xlabel 'Relative radius rbar/Rbar'\n");
	fprintf(gnuplot, "set arrow 1 from 1.0, graph 0 to 1.0, graph 1 lc rgb 'black' lw 2 nohead\n");
	//--------
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'ρ/ρ_c' enhanced", filename);
	fprintf(gnuplot, ", '%s' u 1:3 w l t 'P/P_c'", filename);
	fprintf(gnuplot, ", '%s' u 1:4 w l t 'm(r)/M'", filename);
	fprintf(gnuplot, ", '%s' u 1:5 w l t 'e^{ν}' enhanced", filename);
	fprintf(gnuplot, ", '%s' u 1:6 w l t 'e^{λ}' enhanced", filename);
	//--------
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	system("rm littest/star*.txt");
}//*/








