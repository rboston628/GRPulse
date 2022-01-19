// *************************************************************************************
//					RELATIVISTIC POLYTROPE EXACT TEST
// grexacttest.cpp
//		In the GR case, there is one exact solution (see Tooper 1964, eq 3.3)
//			n = 0 (the uniform relativistic model)
//		We numerically compute simple GR polytropes for n=0 and 
//			calculate the error against the known solution, in logdif_4K.png
//		We also compare radii to exact solution in Tooper eq. 3.4, in testofradii.txt
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/GRPolytrope.cpp"

void printCompare(int, double*, double*, GRPolytrope*);

int main(){
	//data for stars
	double n=0.0;
	double redshift[10] = {1.e-1,8e-2,6e-2,4e-2,2e-2,1.e-2,1.e-3,1.e-4,1.e-5,1.e-6};
	double sigma[10];
	
	//initialize a series of polytropes
	double *theta, *xi;	
	GRPolytrope *star;
	
	int length = 4096;
	int Len[10];
	double xi2analytic[10], xi2numeric[10];
	
	FILE *fp;
	if(!(fp = fopen("./exacttest/testofradii.txt", "w"))){
		system("mkdir exactest");
		fp = fopen("./exacttest/testofradii.txt", "w");
	}
	for(int i=0; i<10; i++){
		//calculate GR polytrope
		star = new GRPolytrope(n, redshift[i], length);
		sigma[i] = star->getSigma();
		printf("sigma=%le\tredshift=%le\n", sigma[i], redshift[i]);
		int L = Len[i] = star->length();
		theta = new double[L];
		xi = new double[L];
		
		double rootx2=0.0, sp1 = 0.0;
		//compute analytic solution
		for(int x=0; x<L; x++){
			xi[x] = star->getX(x);
			rootx2 = (1.+3.*sigma[i])*sqrt(1.-2./3.*sigma[i]*xi[x]*xi[x]);
			sp1 = 1.+sigma[i];
			theta[x]  = ( rootx2 - sp1 )/sigma[i]/( 3.*sp1 - rootx2 );
		}
		
		//print the two for comparison
		printCompare(i, xi, theta, star);
		
		//compare to Eq 3.4 in Tooper for terminal radius
		xi2analytic[i] = 6.*(1.+2.*sigma[i])/((1.+3.*sigma[i])*(1.+3.*sigma[i]));
		xi2numeric[i]  = xi[L-1]*xi[L-1];
		fprintf(fp, "sigma=%0.1lf\t exct=%0.16le\t", sigma[i], xi2analytic[i]);
		fprintf(fp, "num =%0.16le\t diff = %0.16le\n", xi2numeric[i], fabs(xi2analytic[i]-xi2numeric[i]));
		
		delete star;
		delete[] theta;
		delete[] xi;
	}
	fclose(fp);
	
	for(int i=0; i<10; i++){
		//compare to Eq 3.4 in Tooper for terminal radius
		printf("\tsigma=%0.1lf\t exct=%0.16le\t", sigma[i], xi2analytic[i]);
		printf("num =%0.16le\t diff = %0.16le\n", xi2numeric[i], fabs(xi2analytic[i]-xi2numeric[i]));
		fprintf(fp, "sigma=%0.1le\t exct=%0.16le\t", sigma[i], xi2analytic[i]);
		fprintf(fp, "num =%0.16le\t diff = %0.16le\n", xi2numeric[i], fabs(xi2analytic[i]-xi2numeric[i]));
	}
	
	
	//plot everything in single graph, for simplicity
	char filename[91];
	char outname[91];
	int thou = length/1000;
	sprintf(outname, "./exacttest/logdif_%dK.png", thou);
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,1000\n");
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Errors in Relativistic Polytrope n=0, N_{star} = %d'\n", length);
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
	//fprintf(gnuplot, "set yrange [1e-16: 1e-8]\n");
	fprintf(gnuplot, "set ylabel 'Error in θ'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	sprintf(filename, "./exacttest/sigma%0.1e.txt", sigma[0]);
	fprintf(gnuplot, "plot '%s' u 2:5 w l t 'sigma=%0.1e'", filename, sigma[0]);
	for(int i=1; i<10; i++){
		sprintf(filename, "./exacttest/sigma%0.1e.txt", sigma[i]);
		fprintf(gnuplot, ", '%s' u 2:5 w l t 'sigma=%0.1e'", filename, sigma[i]);
	}
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	
	system("rm ./exacttest/sigma*.txt");
	
	return 0;
}

void printCompare(int i, double *xi, double *theta, GRPolytrope *star){
	int L = star->length();
	
	//create names for files to be opened
	char filename[91];
	char outname[91];
	char openmyplot[45];
	//save data to folder to avoid clutter - make sure folder exists
	sprintf(filename, "./exacttest/sigma%0.1e.txt", star->getSigma());
	FILE *fp;
	if(!(fp = fopen(filename, "w")) ){
		system("mkdir ./exacttest");
		fp = fopen(filename, "w");
	}
	//print results to text file
	double XI = xi[L-1];
	for(int x=0; x<L; x++){
		fprintf(fp, "%d\t%le\t%0.30le\t%0.30le\t%0.30le\n", x, xi[x]/XI, theta[x], star->getY(x), fabs(theta[x]-star->getY(x)) );
	}
	fclose(fp);	
}








