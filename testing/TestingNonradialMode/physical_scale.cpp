// *************************************************************************************
//					NEWTONIAN MODE RMSR SCALE TEST
// physical_scale.cpp
//		This Dziembowski equations solve for variables y_1,y_2,y_3,y_4
//		These are related to Drho, DP, etc.
//		This code takes the solution in terms of y_i, finds physical variabes, 
//			and puts them into the LAWE to compute a residual
//		We expect this RMSR to scale like N^{-4}
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/Polytrope.cpp"
#include "../../src/MODES/NonradialModeDriver.cpp"
#include "../../src/MODES/Mode.cpp"

//data for stars
double n[10] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};


int main(){
	//ask user for tolerance of calculation
	int fac = 2;
	int length1 = 2048;
	int length2 = fac*length1-(fac-1);
	
	//initialize a series of polytropes
	//Polytrope *star;
	Polytrope *star1, *star2;
	NonradialModeDriver  *driv1, *driv2;
	const unsigned int num_var = NonradialModeDriver::num_var;
	Mode<num_var> *mode1, *mode2;
	
	//adiabatic exponent
	double Gamma1 = 5./3.;
		
	FILE *fp;
	for(int i=0; i<10; i++){
		//initialize the stars
		star1 = new Polytrope(n[i],length1);
		star2 = new Polytrope(n[i],length2);
		driv1 = new NonradialModeDriver(star1, Gamma1);
		driv2 = new NonradialModeDriver(star2, Gamma1);
		
		char mkdir[100];
		sprintf(mkdir, "rm -r ./physscale/poly%1.1f/", n[i]);
		system(mkdir);
		sprintf(mkdir, "mkdir -p ./physscale/poly%1.1f/", n[i]);
		system(mkdir);
	
		double sig2omeg = pow(star1->Radius(),3)/(star1->Mass());
		
		int l0=2;
		for(int l=1; l<4; l++){
			for(int k=1; k<8; k++){
				printf("-----------------------------------------------------\n");
				mode1 = new Mode<num_var>(2*(k+1),l,l, driv1);
				mode2 = new Mode<num_var>(2*(k+1),l,l, driv2);
				int Len = driv1->length();
				int K1 = mode1->modeOrder(), K2 = mode2->modeOrder();
				printf("N=%1.1f\tL=%d\tK=%d\tK=%d\n", n[i], l, K1,K2);
			
				if(K1!=K2) continue;
			
				char filename[100];
				sprintf(filename, "./physscale/poly%1.1f/test_N%1.1f_L%d_K%d.txt", n[i],n[i], l,K1);
				fp = fopen(filename,"w");
			
				double R = star1->Radius();
			
				double e1,e2,e3;
				double n1,n2,n3;
				double d2Phi, difxi, dDP;
				for(int X=4; X<Len-4; X++){
					int XX = 2*X;
					//stellar variables, to simplify equations
					double rho    = star1->rho(XX);
					double P      = star1->P(XX);
					double drhodr = star1->drhodr(XX);
					double dPdr   = star1->dPdr(XX);
					double g      = star1->dPhidr(XX);
					double r      = star1->rad(XX);
					double G1     = 5.0/3.0;
					//mode variables
					double sigma2 = mode1->getOmega2()/sig2omeg;
					double L2     = double(l*l+l);
					double xi   = r  *mode1->getY(0,X);	//r  *y1
					double chi  = r*g*mode1->getY(1,X);	//r*g*y2
					double DPhi = r*g*mode1->getY(2,X);	//r*g*y3
					double dDPhi=   g*mode1->getY(3,X);	//  g*y4
					double Drho = rho*rho/(G1*P)*(chi - DPhi - g*xi) - xi*drhodr;
									
					//calculate numerical derivatives
					double b3,b2,b1,a1,a2,a3, h1;
					h1= star1->rad(XX+2)-star1->rad(XX );
					// d2DPhi/dr2 = (d/dr)(dDphi/dr), dDPhi/dr = g y4
					b3=star1->dPhidr(XX-6)*mode1->getY(3,X-3);
					b2=star1->dPhidr(XX-4)*mode1->getY(3,X-2);
					b1=star1->dPhidr(XX-2)*mode1->getY(3,X-1);
					a1=star1->dPhidr(XX+2)*mode1->getY(3,X+1);
					a2=star1->dPhidr(XX+4)*mode1->getY(3,X+2);
					a3=star1->dPhidr(XX+6)*mode1->getY(3,X+3);
					//d2Phi = (-a2 + 8.*a1 - 8.*b1 + b2)/(12.*h1);
					d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//next dxi/dr, xi = r y1
					b3=star1->rad(XX-6)*mode1->getY(0,X-3);
					b2=star1->rad(XX-4)*mode1->getY(0,X-2);
					b1=star1->rad(XX-2)*mode1->getY(0,X-1);
					a1=star1->rad(XX+2)*mode1->getY(0,X+1);
					a2=star1->rad(XX+4)*mode1->getY(0,X+2);
					a3=star1->rad(XX+6)*mode1->getY(0,X+3);
					//difxi = (-a2 + 8.*a1 - 8.*b1 + b2)/(12.*h1);
					difxi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//now dDP/dr, DP = rho*(chi-DPhi), chi=rg*y2, DPhi=rg*y3
					b3=star1->rho(XX-6)*star1->dPhidr(XX-6)*star1->rad(XX-6)*(mode1->getY(1,X-3)-mode1->getY(2,X-3));
					b2=star1->rho(XX-4)*star1->dPhidr(XX-4)*star1->rad(XX-4)*(mode1->getY(1,X-2)-mode1->getY(2,X-2));
					b1=star1->rho(XX-2)*star1->dPhidr(XX-2)*star1->rad(XX-2)*(mode1->getY(1,X-1)-mode1->getY(2,X-1));
					a1=star1->rho(XX+2)*star1->dPhidr(XX+2)*star1->rad(XX+2)*(mode1->getY(1,X+1)-mode1->getY(2,X+1));
					a2=star1->rho(XX+4)*star1->dPhidr(XX+4)*star1->rad(XX+4)*(mode1->getY(1,X+2)-mode1->getY(2,X+2));
					a3=star1->rho(XX+6)*star1->dPhidr(XX+6)*star1->rad(XX+6)*(mode1->getY(1,X+3)-mode1->getY(2,X+3));
					//dDP = (-a2 + 8.*a1 - 8.*b1 + b2)/(12.*h1);
					dDP = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//Now calculate residuals
					//Perturbed Poisson equation
					e1 = fabs(
							4.0*m_pi*Drho
							+ L2*pow(r,-2)*DPhi 
							- 2.0*dDPhi/r  - d2Phi
					);
					//continuity equation
					e2 = fabs(
							Drho + xi*drhodr + rho*( 2.*xi/r + difxi - L2*pow(r,-2)/sigma2*chi )
					);
					//newton's equation -- the r component. The theta component defines chi = r sigma2 xiH
					e3 = fabs( rho*sigma2*xi - g*Drho - rho*dDPhi - dDP );
				
					//now go to finer mesh star
					int fXX = fac*XX;
					int fX  = fac*X;
					//stellar variables, to simplify equations
					rho    = star2->rho(fXX);
					P      = star2->P(fXX);
					drhodr = star2->drhodr(fXX);
					dPdr   = star2->dPdr(fXX);
					g      = star2->dPhidr(fXX);
					r      = star2->rad(fXX);
					G1     = 5.0/3.0;
					//mode variables
					sigma2 = mode2->getOmega2()/sig2omeg;
					L2     = double(l*l+l);
					xi   = r  *mode2->getY(0,fX);	//r  *y1
					chi  = r*g*mode2->getY(1,fX);	//r*g*y2
					DPhi = r*g*mode2->getY(2,fX);	//r*g*y3
					dDPhi=   g*mode2->getY(3,fX);	//  g*y4
					Drho = rho*rho/(G1*P)*(chi - DPhi - g*xi) - xi*drhodr;
									
					//calculate numerical derivatives
					h1= star2->rad(fXX+fac*2)-star2->rad(fXX );
					// d2DPhi/dr2 = (d/dr)(dDphi/dr), dDPhi/dr = g y4
					b3=star2->dPhidr(fXX-fac*6)*mode2->getY(3,fX-fac*3);
					b2=star2->dPhidr(fXX-fac*4)*mode2->getY(3,fX-fac*2);
					b1=star2->dPhidr(fXX-fac*2)*mode2->getY(3,fX-fac*1);
					a1=star2->dPhidr(fXX+fac*2)*mode2->getY(3,fX+fac*1);
					a2=star2->dPhidr(fXX+fac*4)*mode2->getY(3,fX+fac*2);
					a3=star2->dPhidr(fXX+fac*6)*mode2->getY(3,fX+fac*3);
					//d2Phi = (-a2 + 8.*a1 - 8.*b1 + b2)/(12.*h1);
					d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//next dxi/dr, xi = r y1
					b3=star2->rad(fXX-fac*6)*mode2->getY(0,fX-fac*3);
					b2=star2->rad(fXX-fac*4)*mode2->getY(0,fX-fac*2);
					b1=star2->rad(fXX-fac*2)*mode2->getY(0,fX-fac*1);
					a1=star2->rad(fXX+fac*2)*mode2->getY(0,fX+fac*1);
					a2=star2->rad(fXX+fac*4)*mode2->getY(0,fX+fac*2);
					a3=star2->rad(fXX+fac*6)*mode2->getY(0,fX+fac*3);
					//difxi = (-a2 + 8.*a1 - 8.*b1 + b2)/(12.*h1);
					difxi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//now dDP/dr, DP = rho*(chi-DPhi), chi=rg*y2, DPhi=rg*y3
					b3=star2->rho(fXX-fac*6)*star2->dPhidr(fXX-fac*6)*star2->rad(fXX-fac*6)*(mode2->getY(1,fX-fac*3)-mode2->getY(2,fX-fac*3));
					b2=star2->rho(fXX-fac*4)*star2->dPhidr(fXX-fac*4)*star2->rad(fXX-fac*4)*(mode2->getY(1,fX-fac*2)-mode2->getY(2,fX-fac*2));
					b1=star2->rho(fXX-fac*2)*star2->dPhidr(fXX-fac*2)*star2->rad(fXX-fac*2)*(mode2->getY(1,fX-fac*1)-mode2->getY(2,fX-fac*1));
					a1=star2->rho(fXX+fac*2)*star2->dPhidr(fXX+fac*2)*star2->rad(fXX+fac*2)*(mode2->getY(1,fX+fac*1)-mode2->getY(2,fX+fac*1));
					a2=star2->rho(fXX+fac*4)*star2->dPhidr(fXX+fac*4)*star2->rad(fXX+fac*4)*(mode2->getY(1,fX+fac*2)-mode2->getY(2,fX+fac*2));
					a3=star2->rho(fXX+fac*6)*star2->dPhidr(fXX+fac*6)*star2->rad(fXX+fac*6)*(mode2->getY(1,fX+fac*3)-mode2->getY(2,fX+fac*3));
					//dDP = (-a2 + 8.*a1 - 8.*b1 + b2)/(12.*h1);
					dDP = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
			
					//Poisson equation
					n1 = fabs(
							4.0*m_pi*Drho
							+ L2*pow(r,-2)*DPhi 
							- 2.0*dDPhi/r  - d2Phi
					);
					//continuity equation
					n2 = fabs(
							Drho + xi*drhodr + rho*( 2.*xi/r + difxi - L2*pow(r,-2)/sigma2*chi )
					);
					//newton's equation -- the r component. The theta component defines chi = r sigma2 xiH
					n3 = fabs( rho*sigma2*xi - g*Drho - rho*dDPhi - dDP );
				
					e1 = e1/n1;
					e2 = e2/n2;
					e3 = e3/n3;
				
					//print residuals
					fprintf(fp, "%le\t%le\t%le\t%le\t%le\t%le\t%le\n", mode1->getRad(X), e1, e2, e3, n1, n2, n3);
				}
				fclose(fp);
				delete mode1;
				delete mode2;
			
				char outname[100];
				sprintf(outname, "./physscale/poly%1.1f/test_N%1.1f_L%d_K%d.png", n[i],n[i], l,K1);
				FILE *gnuplot = popen("gnuplot -persist", "w");
				fprintf(gnuplot, "reset\n");
				fprintf(gnuplot, "set term png size 1600,800\n");
				fprintf(gnuplot, "set samples %d\n", Len);
				fprintf(gnuplot, "set output '%s'\n", outname);
				fprintf(gnuplot, "set xlabel 'r/R'\n");
				fprintf(gnuplot, "set ylabel 'ε_{%dK}/ε_{%dK}'\n", fac*2, 2 );
				fprintf(gnuplot, "set logscale y %d\n", fac);
				fprintf(gnuplot, "set format y '%d^{%%L}'\n", fac);
				fprintf(gnuplot, "set ytics %d\n", fac);
				fprintf(gnuplot, "set title 'Scaling of Perturbation Residuals, Mode L=%d, K=%d, polytrope star n=%1.1f, scale factor = %d'\n", l,K1,n[i], fac);
				fprintf(gnuplot, "set arrow 1 from 0.0,%le to 1.01, %le lw 4 lc rgb 'red' nohead\n", pow(fac,4), pow(fac,4));
				fprintf(gnuplot, "plot '%s' u 1:3 w l t 'Continuity resid'", filename);
				fprintf(gnuplot, ",    '%s' u 1:4 w l t 'Newtons Law resid'", filename);
				fprintf(gnuplot, ",    '%s' u 1:2 w l t 'Poisson resid'", filename);
				fprintf(gnuplot, "\n");
				pclose(gnuplot);			
			}
		}		
		delete star1;
		delete star2;
		delete driv1;
		delete driv2;
	}


	return 0;
}