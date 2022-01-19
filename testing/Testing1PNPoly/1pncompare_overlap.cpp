// *************************************************************************************
//					POST-NEWTONIAN POLYTROPE OVERLAP TEST
// compare_overlap.cpp
//		The 1PN Lane-Emden equation solves for the variable θ, related to ρ,P
//		The Newtonian and GR Lane-Emden equations also solve for the variable θ, related to ρ,P
//		The solution θ_{1pn} can be compared to the Newtonian and GR solutions
//			by defining an overlap based on an inner product
//		This overlap should scale with the redshift z, as follows
//			newt & 1pn		~	z
//			newt & GR 		~	z
//			1PN  & GR 		~	z^2
// Uses spline fitting to match the raddi between the 1PN and GR solutions
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"
#include "../../src/STARS/PNPolytrope.cpp"
#include "../../src/STARS/GRPolytrope.cpp"
#include "../../lib/Splinor.cpp" //it radii between 1pn (harmonic) and GR (schwarzschild)


double innerprod(int kmax, double* r, double* f1, double* f2);

int main(){
	//data for stars
	double zz[] = {
		1.e-1, 8.e-2, 6.e-2, 4.e-2, 2.e-2,
		1.e-2, 8.e-3, 6.e-3, 4.e-3, 2.e-3, 
		1.e-3, 8.e-4, 6.e-4, 4.e-4, 2.e-4, 
		1.e-4, 8.e-5, 6.e-5, 4.e-5, 2.e-5, 
		1.e-5, 1.e-6};
	int len = 1000;
	
	//initialize a series of polytropes
	Polytrope *star;
	PNPolytrope *starPN;
	GRPolytrope *starGR;
	
	char filename[255];
	
	double xcl[len], ycl[len];
	double xpn[len], ypn[len];
	double xgr[len], ygr[len];
	double dxpn, dxgr;
	double prod1, prod2, prod3, Ncl, Npn, Ngr;

	double nn[] = {0.0,1.0,1.5,3.0};
	for(double n : nn) {
		printf("Index n = %lf\n", n);
		sprintf(filename, "./1pnlimittest/compare_overlap_%0.1lf.txt", n);
		FILE *fp;
		if(!(fp=fopen(filename, "w"))){
			system("mkdir ./1pnlimittest");
			fp=fopen(filename, "w");
		}
		//begin with the classical star
		star = new Polytrope(n, len);
		double dx, XX = star->getX(len-1);
		//store the classical radius and solution
		for(int X=0; X<len; X++){
			xcl[X] = star->getX(X);
			ycl[X] = star->getY(X);
		}
		delete star;
		printf("\tclassical star completed\n");
		//now calculate the 1PN and GR polytropes
		for(double zsurf :zz){
			printf("\tsigma = %lf\n", zsurf);
			starPN = new PNPolytrope(n, zsurf, len);
			starGR = new GRPolytrope(n, zsurf, len);
			//store the solutions y and the proper radii
			dxpn = starPN->getX(2)-starPN->getX(1);
			dxgr = starGR->getX(2)-starGR->getX(1);
			//
			double int1pn = 1.0, int2pn = 1.0;
			double int1gr = 1.0, int2gr = 1.0;
			//
			xpn[0] = xgr[0] = 0.0;
			ypn[0] = ygr[0] = 1.0;
			//transform to Schwarzschild coordinates
			for(int X=1; X<len; X++){
				xpn[X] = starPN->getX(X)*sqrt(1.-2.*starPN->Phi(X));
				xgr[X] = starGR->getX(X);
 				ypn[X] = starPN->getY(X);
 				ygr[X] = starGR->getY(X);
			}
			delete starPN;
			delete starGR;
			printf("\t\tpost-newtonian star completed\n");
			printf("\t\trelativistic star completed\n");
									
			//now fit splines based on proper radii
			printf("\t\tcreating splines..."); fflush(stdout);
			Splinor *ypn_spline = new Splinor(xpn, ypn, len);
			Splinor *ygr_spline = new Splinor(xgr, ygr, len);
			printf("done\n");
						
			//evaluate based on the classical radius
			printf("\t\tinterpolating..."); fflush(stdout);
			for(int X=0; X<len; X++){
				ypn[X] = ypn_spline->interp(xcl[X]);
				ygr[X] = ygr_spline->interp(xcl[X]);
			}
			printf("done\n");
			
			//
			printf("\tprinting results for z=%lf\n", zsurf);
			Ncl = innerprod(len, xcl, ycl, ycl);
			Npn = innerprod(len, xcl, ypn, ypn);
			Ngr = innerprod(len, xcl, ygr, ygr);
			prod1 = (sqrt(Ncl*Npn)-innerprod(len, xcl, ycl, ypn))/sqrt(Ncl*Npn);
			prod2 = (sqrt(Ncl*Ngr)-innerprod(len, xcl, ycl, ygr))/sqrt(Ncl*Ngr);
			prod3 = (sqrt(Ngr*Npn)-innerprod(len, xcl, ypn, ygr))/sqrt(Npn*Ngr);
			//prod1 = 1.-innerprod(len, xcl, ycl, ypn)/sqrt(Ncl*Npn);
			//prod2 = 1.-innerprod(len, xcl, ycl, ygr)/sqrt(Ncl*Ngr);
			//prod3 = 1.-innerprod(len, xcl, ypn, ygr)/sqrt(Npn*Ngr);
			//if(prod1<0) prod1=-prod1;//1.e-16;
			//if(prod2<0) prod2=-prod2;//1.e-16;
			//if(prod3<0) prod3=-prod3;//1.e-16;
			fprintf(fp, "%0.12le\t%0.12le\t%0.12le\t%0.12le\n", zsurf, prod1, prod2, prod3);
			delete ypn_spline;
			delete ygr_spline;
		}
		fclose(fp);
		//plot everything in single graph, for simplicity
		char outname[91];
		sprintf(outname, "./1pnlimittest/compare_overlap_%0.1lf.png", n);
		FILE *gnuplot = popen("gnuplot -persist", "w");
		fprintf(gnuplot, "reset\n");
		fprintf(gnuplot, "set term png size 800,800\n");
		fprintf(gnuplot, "set output '%s'\n", outname);
		fprintf(gnuplot, "set encoding utf8\n");
		sprintf(filename, "./1pnlimittest/compare_overlap_%0.1lf.txt", n);
		fprintf(gnuplot, "set title 'Scaling of Overlap between Polytropes in Classical, GR, and 1PN, ");
		fprintf(gnuplot, "n=%0.1lf' enhanced\n", n);
		
		fprintf(gnuplot, "set logscale x\n");
		fprintf(gnuplot, "set logscale y\n");
		fprintf(gnuplot, "set format x '10^{%%L}'\n");
		fprintf(gnuplot, "set format y '10^{%%L}'\n");
		fprintf(gnuplot, "set xrange [1e-1:1e-6]\n");
		fprintf(gnuplot, "set ytics 0.1\n");
		
		double start=3.e-3, end = 7.e-6, range = log10(start/end), y0=1e-4, z0 = 1e-5;
		fprintf(gnuplot, "set arrow 1 from %le,%le to %le,%le lc rgb 'red' nohead\n", start, y0, end, y0*pow(1e-2,range));
		fprintf(gnuplot, "set arrow 2 from %le,%le to %le,%le lc rgb 'red' nohead\n", start, z0, end, z0*pow(1e-4,range));
		
		fprintf(gnuplot, "set xlabel 'z_{surf}'\n");
		fprintf(gnuplot, "set ylabel 'overlap'\n");
		fprintf(gnuplot, "plot '%s' u 1:2 w l t 'CL-PN',", filename);
		fprintf(gnuplot, "'%s' u 1:3 w l t 'CL-GR',", filename);
		fprintf(gnuplot, "'%s' u 1:4 w l t 'PN-GR'\n", filename);

		
		
		pclose(gnuplot);
		char command[255];
		//system(command);
	}
	
	return 0;
}

//trapezoidal integration, I guess?
double innerprod(int kmax, double* r, double* f1, double* f2){
	double sum=0.0, int1 = 0.0, int2 = f1[0]*f2[0];
	double dx = 0.0;

	for(int k=1; k<kmax-2; k++){
		dx = r[k]-r[k-1];	
		int1 = int2;
		int2 = f1[k]*f2[k];
		sum += 0.5*(int1+int2)*dx;
	}
	return sum;
}



