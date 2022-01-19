// *************************************************************************************
//					NEWTONIAN POLYTROPE LITERATURE TEST
// littest.cpp
//		A comparison to literature values, tables from 
//			Mohan & Al-Bayaty, Astrophys Space Sci vol 73, 227–239 (1980)
//			https://doi.org/10.1007/BF00642378
//		Values of theta for n = 0.5, 1.0, ..., 4.9
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"

int main(){
	//character buffer for filenames
	char filenameRK4[40];
	//char filenameRKF[40];
	char filenameOut[40];
	
	//file pointers for opening files for reading: four total
	FILE *fileLit;	//table of literature values
	FILE *fp;
		
	char mkdir[100];
	sprintf(mkdir, "rm -r ./littest");
	system(mkdir);
	sprintf(mkdir, "mkdir ./littest");
	system(mkdir);
	
 	char bufferRK4[255];
	
	//read in literature values from table
	int L  = 21;
	double literature[L][10];
	double litXi[L];
	fileLit = fopen("LiteratureTable.txt", "r");
	if(!fileLit) printf("Something went wrong...");
	for(int i=0; i<L; i++){
		fscanf(fileLit, "%lf,", &litXi[i]);
		printf("%lf\t", litXi[i]);
		for(int j = 0; j<9; j++){
			fscanf(fileLit, "%lf,", &literature[i][j]);
			printf("%lf,", literature[i][j]);
		}
		fscanf(fileLit, "%lf\n", &literature[i][9]);
		printf("%lf\n", literature[i][9]);
	}
	printf("done reading... ");
	fclose(fileLit);
	printf("file closed\n");
		

	Polytrope *star;
	int length = 10000;
	double xi1, xi2, xiF;
	double numerical[L][10];
	
	int ll = 0;
	int lj = 0;
	int l  = 0;
	
	for(float n = 0.5; n <= 5.0; n += 0.5){
		if(n==5.0) n = 4.9;
		length = 10000;
		printf("POLYTROPE INDEX N = %f\n", n);
		//create polytrope of given index
		star = new Polytrope(n, length);
		length = star->length();
		xiF = star->getX(length-1);
		
		ll = 0;
		lj = int(n*2 - 1);
		if(n==4.9) lj = 9;
		
		//run through solution to find values where xi corresponds to tabulated values
		
		numerical[ll][lj] = star->getY(0);
		ll = 1;
		for(int x=1; x<length-1 & ll<L; x++){
			xi1 = star->getX(x-1)/xiF;
			xi2 = star->getX(x+1)/xiF;
			if(xi1 < litXi[ll] & litXi[ll] < xi2) {
				numerical[ll][lj] = star->getY(x);
				ll++;
			}
		}
		numerical[L-1][lj] = star->getY(length-1);
	}

	fp = fopen("./littest/littest.txt", "w");
	for(int i=0; i<L; i++){
		fprintf(fp, "%0.16lf", litXi[i]);
		printf("%lf\t", litXi[i]);
		for(int j=0; j<10; j++){
			fprintf(fp, "\t%0.16lf", fabs(numerical[i][j]-literature[i][j]));
			printf("%lf, ", numerical[i][j]);
		}
		fprintf(fp, "\n");
		printf("\n");
	}
	fclose(fp);
	
	//plot everything in single graph, for simplicity
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 2000,1000\n");
	fprintf(gnuplot, "set output './littest/littest.png'\n");
	fprintf(gnuplot, "set title 'Comparison to Mohan & Al-Bayaty (1980)\n");
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
	fprintf(gnuplot, "set ylabel 'log_{10}|num-lit|'\n");
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set term png size 2400,1600\n");
	fprintf(gnuplot, "set xrange [0.0 to 1.0]\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "plot ");
	for(int i=0; i<10; i++){
		fprintf(gnuplot, "%c './littest/littest.txt' u 1:%d w l title 'n=%1.1f'", (i==0?' ':','),(i+2), (i!=9?0.5*(i+1):4.9));
	}	
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	
	

	return 0;
}