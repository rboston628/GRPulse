//meant to test the Splinor being used on MESA data

#include <math.h>
#include <stdio.h>
#include "../../src/constants.h"
#include "../../lib/Splinor.cpp"

int main(){
	printf("Beginning read-in of MESA data.\n");
	
	char filename[] = "../../mesa_co_wd_hot.dat";
	int len = 3000;
	
	FILE *infile;
	if( !(infile = fopen(filename, "r")) ) {
		printf("No such file\n");
		exit(1);
	} 
	int Ntot;
	double Mtot, Rtot, Ltot;
	fscanf(infile, "  %d     %le     %le     %le     %*d\n", &Ntot, &Mtot, &Rtot, &Ltot);
	printf("\tMass = %0.2lg\tRadis = %0.2lg\tLuminosity = %0.2lg\n", Mtot, Rtot, Ltot);
	printf("Number of grid points in MESA calculation %d\n", Ntot);
	
	//set scale quantities
	double Dscale = Mtot*pow(Rtot,-3);
	double Pscale = G_CGS*pow(Mtot,2)*pow(Rtot,-4);
	double Gscale = G_CGS*Mtot*pow(Rtot,-2);
			
	//create arrays to hold values read from file
	double *rt = new double[Ntot]; //radius
	double *dt = new double[Ntot]; //density
	double *pt = new double[Ntot]; //pressure
	double *mt = new double[Ntot]; //mass
	double *gt = new double[Ntot]; //gravitational field
	double *Nt = new double[Ntot]; //Brunt-Vaisala frequency
	double *Gt = new double[Ntot]; //Gamma1
		
	//read from the file
	for(int k=0; k<Ntot; k++){
		fscanf(infile, "     %*d     %lg     %lg     %*lg", &rt[k], &mt[k]);
		fscanf(infile, "     %lg     %*lg    %lg     %*lg", &pt[k], &dt[k]);
		fscanf(infile, "     %lg     %lg     %*lg    %*lg", &Nt[k], &Gt[k]);
		fscanf(infile, "%*[^\n]");
		//dis-dimensionalize the variables
		rt[k] /= Rtot;
		mt[k] /= Mtot;
		pt[k] /= Pscale;
		dt[k] /= Dscale;
		if(rt[k]>0)	gt[k] = mt[k]*pow(rt[k],-2);
		else gt[k] = 0.0;
		Nt[k] /= (G_CGS*Dscale);
	}
	printf("Done reading in data!\n");
	fclose(infile);
	
	//divide intervals up into binary fractions
	//must use binary fractions to mesh with RK4 method requiring half-points
	double n = log2(len-1) - log2(Ntot-1);
	if(n < 1) n = 1; //must at least divide each interval in half
	len = pow(2,int(n))*(Ntot-1) + 1;
	int subgrid = pow(2,int(n));
	printf("number of grid points to be used %d\n", len);
			
	double* radi = new double[len];
	int kk=0;
	for(int k=0; k<Ntot-1; k++){
		double dr = (rt[k+1]-rt[k])/subgrid;
		for(int ki=0; ki<subgrid; ki++,kk++) radi[kk] = rt[k]+double(ki)*dr;
	}
	radi[len-1] = rt[Ntot-1];
	
	//now splne fit all quantitites -- remove arrays to conserve space
	Splinor dens(rt, dt, Ntot);
	delete[] dt;
	Splinor pres(rt, pt, Ntot);
	delete[] pt;
	Splinor mass(rt, mt, Ntot);
	delete[] mt;
	Splinor grav(rt, gt, Ntot);
	delete[] gt;
	Splinor BVfq(rt, Nt, Ntot);	
	delete[] Nt;
	
	char mkdir[100];
	sprintf(mkdir, "rm -r ./MESA_spline");
	system(mkdir);
	sprintf(mkdir, "mkdir ./MESA_spline");
	system(mkdir);
	char txtname[] = "./MESA_spline/splinetest.txt";
	infile = fopen(txtname, "w");
	for(int k=0; k<len; k++){
		fprintf(infile, "%0.16lf\t", radi[k]);			//col1
		fprintf(infile, "%0.16lf\t", dens(radi[k]));	//col2
		fprintf(infile, "%0.16lf\t", pres(radi[k]));	//col3
		fprintf(infile, "%0.16lf\t", mass(radi[k]));	//col4
		fprintf(infile, "%0.16lf\t", grav(radi[k]));	//col5
		fprintf(infile, "%0.16lf\t", BVfq(radi[k]));	//col6
		fprintf(infile, "\n");
	}
	fclose(infile);
	
	// PLOT IN GNUPLOT
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set xlabel 'radius'\n");
	fprintf(gnuplot, "set xrange [0:1]\n");
	//density
	fprintf(gnuplot, "set output './MESA_spline/%s_spline.png'\n", "dens");
	fprintf(gnuplot, "set title 'interpolation of MESA density with Splinor'\n");
	//			plot interpolations with lines
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'Splinor'", txtname);
	//			plot the data with points
	fprintf(gnuplot, "   , '%s' u ($2/%le):($7/%le) w p ls 4 t 'data'", filename, Rtot, Dscale);
	fprintf(gnuplot, "\n");	
	//pressure
	fprintf(gnuplot, "set output './MESA_spline/%s_spline.png'\n", "pres");
	fprintf(gnuplot, "set title 'interpolation of MESA pressure with Splinor'\n");
	//			plot interpolations with lines
	fprintf(gnuplot, "plot '%s' u 1:3 w l t 'Splinor'", txtname);
	//			plot the data with points
	fprintf(gnuplot, "   , '%s' u ($2/%le):($5/%le) w p ls 4 t 'data'", filename, Rtot, Pscale);
	fprintf(gnuplot, "\n");	
	//mass
	fprintf(gnuplot, "set output './MESA_spline/%s_spline.png'\n", "mass");
	fprintf(gnuplot, "set title 'interpolation of MESA mass with Splinor'\n");
	//			plot interpolations with lines
	fprintf(gnuplot, "plot '%s' u 1:4 w l t 'Splinor'", txtname);
	//			plot the data with points
	fprintf(gnuplot, "   , '%s' u ($2/%le):($3/%le) w p ls 4 t 'data'", filename, Rtot, Mtot);
	fprintf(gnuplot, "\n");	
	//gravitational field
	fprintf(gnuplot, "set output './MESA_spline/%s_spline.png'\n", "grav");
	fprintf(gnuplot, "set title 'interpolation of MESA gravitational field with Splinor'\n");
	//			plot interpolations with lines
	fprintf(gnuplot, "plot '%s' u 1:5 w l t 'Splinor'", txtname);
	//			plot the data with points
	fprintf(gnuplot, "   , '%s' u ($2/%le):(6.6725985e-8*$3/$2**2/%le) w p ls 4 t 'data'", filename, Rtot, Gscale);
	fprintf(gnuplot, "\n");	
	//Brunt-Vaisala frequency
	fprintf(gnuplot, "set output './MESA_spline/%s_spline.png'\n", "BVfq");
	fprintf(gnuplot, "set title 'interpolation of MESA BV frequency with Splinor'\n");
	fprintf(gnuplot, "set xrange [0:0.99]\n");
	//			plot interpolations with lines
	fprintf(gnuplot, "plot '%s' u 1:6 w l t 'Splinor'", txtname);
	//			plot the data with points
	fprintf(gnuplot, "   , '%s' u ($2/%le):($9/%le) w p ls 4 t 'data'", filename, Rtot, G_CGS*Mtot/Rtot/Rtot/Rtot);
	fprintf(gnuplot, "\n");		
	
	pclose(gnuplot);
}