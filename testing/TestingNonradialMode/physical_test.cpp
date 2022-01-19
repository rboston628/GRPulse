// *************************************************************************************
//					NEWTONIAN MODE ERROR SCALE TEST
// physical_test.cpp
//		This Dziembowski equations solve for variables y_1,y_2,y_3,y_4
//		These are related to Drho, DP, etc.
//		This code takes the solution in terms of y_i, finds physical variabes, 
//			and puts them into the LAWE to compute a residual
//		This same process is performed by the SSR() method of each Mode,
//			However, this program calculates and graphs the residual at each grid point
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"
#include "../../src/MODES/NonradialModeDriver.cpp"
#include "../../src/MODES/Mode.cpp"

//data for stars
double n[10] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};

int main(){
	//ask user for tolerance of calculation
	int length =2048*10;
	
	//initialize a series of polytropes
	Polytrope *star;
	//Isopycnic *star;
	NonradialModeDriver *driv;
	const unsigned int num_var = NonradialModeDriver::num_var;
	Mode<num_var> *mode;
	
	double Gamma1 = 5./3.;
	
	//for each star, check ten modes in three L
	double checkCont[10][3][10];
	double checkNewt[10][3][10];
	double checkPois[10][3][10];
	
	FILE *fp;
	
	for(int i=1; i<10; i++){
		//fp = fopen("./exacttest/testof.txt", "w");
		//calculate polytrope
		star = new Polytrope(n[i], length);
		driv = new NonradialModeDriver(star, Gamma1);
		int Len = driv->length();
		//star = new Isopycnic(length);
		char mkdir[100];
		sprintf(mkdir, "rm -r ./physicaltest/polytrope_%1.1f/", n[i]);
		system(mkdir);
		sprintf(mkdir, "mkdir -p ./physicaltest/polytrope_%1.1f/", n[i]);
		system(mkdir);
		
		double sig2omeg = pow(star->Radius(),3)/(star->Mass());
			
		int l0=2;
		for(int l=3; l<4; l++){
			for(int k=1; k<8; k++){
				printf("-----------------------------------------------------\n");
				mode = new Mode<num_var>(k,l,l, driv);
				int K = mode->modeOrder();
				printf("N=%1.1f\tL=%d\tK=%d\n", n[i], l, K);
	
				char filename[100];
				sprintf(filename, "./physicaltest/polytrope_%1.1f/RMSR_%1.1f_L%d_K%d.txt", n[i],n[i],l,K);
				fp = fopen(filename,"w");
				
				double R = star->Radius();
				
				checkCont[i][l-l0][k] = 0.0;
				checkPois[i][l-l0][k] = 0.0;
				checkNewt[i][l-l0][k] = 0.0;
				double e1,e2,e3;
				double n1,n2,n3;
				double d2Phi, difxi, dDP;
				for(int X=4; X<Len-4; X++){
					int XX = 2*X;
					//stellar variables, to simplify equations
					double rho    = star->rho(XX);
					double P      = star->P(XX);
					double drhodr = star->drhodr(XX);
					double dPdr   = star->dPdr(XX);
					double g      = star->dPhidr(XX);
					double r      = star->rad(XX);
					double G1     = 5.0/3.0;
					//mode variables
					double sigma2 = mode->getOmega2()/sig2omeg;
					double L2     = double(l*l+l);
					double xi   = r  *mode->getY(0,X);	//r  *y1
					double chi  = r*g*mode->getY(1,X);	//r*g*y2
					double DPhi = r*g*mode->getY(2,X);	//r*g*y3
					double dDPhi=   g*mode->getY(3,X);	//  g*y4
					double Drho = rho*rho/(G1*P)*(chi - DPhi - g*xi) - xi*drhodr;
										
					//calculate numerical derivatives
					double b3,b2,b1,a1,a2,a3, h1;
					h1= star->rad(XX+2)-star->rad(XX );
					// d2DPhi/dr2 = (d/dr)(dDphi/dr), dDPhi/dr = g y4
					b3=star->dPhidr(XX-6)*mode->getY(3,X-3);
					b2=star->dPhidr(XX-4)*mode->getY(3,X-2);
					b1=star->dPhidr(XX-2)*mode->getY(3,X-1);
					a1=star->dPhidr(XX+2)*mode->getY(3,X+1);
					a2=star->dPhidr(XX+4)*mode->getY(3,X+2);
					a3=star->dPhidr(XX+6)*mode->getY(3,X+3);
					//d2Phi = (-a2 + 8.*a1 - 8.*b1 + b2)/(12.*h1);
					d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//next dxi/dr, xi = r y1
					b3=star->rad(XX-6)*mode->getY(0,X-3);
					b2=star->rad(XX-4)*mode->getY(0,X-2);
					b1=star->rad(XX-2)*mode->getY(0,X-1);
					a1=star->rad(XX+2)*mode->getY(0,X+1);
					a2=star->rad(XX+4)*mode->getY(0,X+2);
					a3=star->rad(XX+6)*mode->getY(0,X+3);
					//difxi = (-a2 + 8.*a1 - 8.*b1 + b2)/(12.*h1);
					difxi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//now dDP/dr, DP = rho*(chi-DPhi), chi=rg*y2, DPhi=rg*y3
					b3=star->rho(XX-6)*star->dPhidr(XX-6)*star->rad(XX-6)*(mode->getY(1,X-3)-mode->getY(2,X-3));
					b2=star->rho(XX-4)*star->dPhidr(XX-4)*star->rad(XX-4)*(mode->getY(1,X-2)-mode->getY(2,X-2));
					b1=star->rho(XX-2)*star->dPhidr(XX-2)*star->rad(XX-2)*(mode->getY(1,X-1)-mode->getY(2,X-1));
					a1=star->rho(XX+2)*star->dPhidr(XX+2)*star->rad(XX+2)*(mode->getY(1,X+1)-mode->getY(2,X+1));
					a2=star->rho(XX+4)*star->dPhidr(XX+4)*star->rad(XX+4)*(mode->getY(1,X+2)-mode->getY(2,X+2));
					a3=star->rho(XX+6)*star->dPhidr(XX+6)*star->rad(XX+6)*(mode->getY(1,X+3)-mode->getY(2,X+3));
					//dDP = (-a2 + 8.*a1 - 8.*b1 + b2)/(12.*h1);
					dDP = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//Now calculate residuals
					//Perturbed Poisson equation
					e1 = fabs(
							4.0*m_pi*Drho
							+ L2*pow(r,-2)*DPhi 
							- 2.0*dDPhi/r  - d2Phi
					);
					n1 = fabs(4.0*m_pi*Drho )
							+ fabs( L2*pow(r,-2)*DPhi )
							+ fabs( 2.0*dDPhi/r ) + fabs( d2Phi );
					//continuity equation
					e2 = fabs(
							Drho + xi*drhodr + rho*( 2.*xi/r + difxi - L2*pow(r,-2)/sigma2*chi )
					);
					n2 = fabs(Drho) + fabs(xi*drhodr) + fabs(2.*rho/r*xi)
							+ fabs(rho*difxi)
							+ fabs(rho*L2*chi*pow(r,-2)/sigma2);
					//newton's equation -- the r component. The theta component defines chi = r sigma2 xiH
					e3 = fabs( rho*sigma2*xi - g*Drho - rho*dDPhi - dDP );
					n3 = fabs( rho*sigma2*xi ) + fabs( g*Drho ) + fabs( rho*dDPhi ) + fabs( dDP );
					//normalize residuals
					e1 = e1/n1;
					e2 = e2/n2;
					e3 = e3/n3;
					//collect residuals
					checkPois[i][l-l0][k] += e1;
					checkCont[i][l-l0][k] += e2;
					checkNewt[i][l-l0][k] += e3;
					fprintf(fp, "%le\t%le\t%le\t%le\t%le\t%le\t%le\n", mode->getRad(X), e1, e2, e3, xi, chi, DPhi);
				}
				
				printf("Sum of Errors in Continuity:\t%le\n", checkCont[i][l-l0][k]);
				printf("Sum of Errors in Newton's Law:\t%le\n", checkNewt[i][l-l0][k]);
				printf("Sum of Errors in Poisson:\t%le\n", checkPois[i][l-l0][k]);
				fclose(fp);
				delete mode;
				
				int base = 10;
				char outname[100];
				sprintf(outname, "./physicaltest/polytrope_%1.1f/RMSR_%1.1f_L%d_K%d.png", n[i],n[i], l, K);
				FILE *gnuplot = popen("gnuplot -persist", "w");
				fprintf(gnuplot, "reset\n");
				fprintf(gnuplot, "set term png size 1600,800\n");
				fprintf(gnuplot, "set samples %d\n", Len);
				fprintf(gnuplot, "set output '%s'\n", outname);
				fprintf(gnuplot, "set xlabel 'r/R'\n");
				fprintf(gnuplot, "set ylabel 'residual'\n");
				fprintf(gnuplot, "set logscale y %d\n", base);
				fprintf(gnuplot, "set format y '%d^{%%L}'\n", base);
				fprintf(gnuplot, "set ytics %d\n", base);
				fprintf(gnuplot, "set title 'Normalized Residuals in Perturbation Equations, Mode L=%d, K=%d, Background polytrope n=%1.1f, Grid = %d'\n", l,K,n[i],length);
				fprintf(gnuplot, "plot '%s' u 1:3 w l t 'Continuity resid'", filename);
				fprintf(gnuplot, ",    '%s' u 1:4 w l t 'Newtons Law resid'", filename);
				fprintf(gnuplot, ",    '%s' u 1:2 w l t 'Poisson resid'", filename);
				fprintf(gnuplot, "\n");
				pclose(gnuplot);//*/			
			}
		}		
		delete star;
	}
	return 0;
}