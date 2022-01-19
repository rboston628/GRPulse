//meant to test different interpolation techniques

#include "../../src/constants.h"
#include "../../lib/Splinor.cpp"


int main(){
	int N;
	double *x, *y;
	
	//open a file and read the length of the data from first line
	char name[6][20] = {"test_step", "test_cosine", "test_cubic", "test_log", "test_runge", "test_window"};
	char data_table[255];
	
	FILE *fp;
	char mkdir[100];
	sprintf(mkdir, "rm -r ./splinetest");
	system(mkdir);
	sprintf(mkdir, "mkdir ./splinetest");
	system(mkdir);
	
	for(int i=0; i<6; i++){
		//open the data table to be fit
		sprintf(data_table, "%s.dat", name[i]);
		FILE *data_file = fopen(data_table, "r");
		fscanf(data_file, "%d %*[^\n]", &N);
		printf("%d\n", N);
		
		//create new arrays to hold data
		x = new double[N+2];
		y = new double[N+2];
	
		//now read in the data from the file
		for(int i=0; i<N; i++){
			fscanf(data_file, "%lf\t%lf\n", &x[i], &y[i]);
		}
		fclose(data_file);
		
		y[N] = y[N+1] = 0.0;
		x[N] = x[N+1] = 0.0;
		int N_sub = 10;
		printf("N=%d\tN_sub=%d\tN_sub*(N-1)=%d\n", N, N_sub, N_sub*(N-1));
		int NN = N_sub*(N-1)+1;
		printf("NN=%d\n", NN);
	
		//now make the splinor
		Splinor spline(x, y, N);
	
		double xtest[NN];
		xtest[0]    = x[0];
		xtest[NN-1] = x[N-1];
		double dx = (xtest[NN-1]-xtest[0])/(NN-2);
		for(int k=1; k<NN; k++){
			xtest[k] = xtest[k-1]+dx;
		}
	
	
		char inname[] = "compare.txt";
		FILE* infile = fopen(inname, "w");
		for(int k=0; k<NN-1; k++){
			fprintf(infile, "%lf\t%lf\n", xtest[k], spline(xtest[k]));
		}
		fclose(infile);
	
		//first print all to one file
		FILE *gnuplot = popen("gnuplot -persist", "w");
		fprintf(gnuplot, "reset\n");
		fprintf(gnuplot, "set term png size 1600,800\n");
		fprintf(gnuplot, "set output './splinetest/%s_spline.png'\n", name[i]);
		fprintf(gnuplot, "set title 'interpolation of data with cubic Splinor'\n");
		fprintf(gnuplot, "set xlabel 'x'\n");
		//plot interpolations with lines
		fprintf(gnuplot, "plot '%s' u 1:2 w l t 'Splinor'", inname);
		//plot the data with points
		fprintf(gnuplot, "   , '%s.dat' u 1:2 w p ls 4 t 'data'", name[i]);
		fprintf(gnuplot, "\n");	
		pclose(gnuplot);
	
		delete[] x;
		delete[] y;
	}
	
	system("rm compare.txt");
}



