// *************************************************************************************
//					POST-NEWTONIAN MODE ERROR SCALE TEST
// physical_test.cpp
//		This Dziembowski equations solve for variables y_1,y_2,y_3,y_4
//		These are related to Drho, DP, etc.
//		This code takes the solution in terms of y_i, finds physical variabes, 
//			and puts them into the LAWE to compute a residual
//		This same process is performed by the SSR() method of each Mode,
//			However, this program calculates and graphs the residual at each grid point
// *************************************************************************************

#include "../../src/STARS/Star.cpp"
#include "../../src/STARS/PNPolytrope.cpp"
#include "../../src/MODES/PNNonradialModeDriver.cpp"
#include "../../src/MODES/Mode.cpp"

//data for stars
double n[10] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};


int main(){
	typedef PNPolytrope STAR ;
	typedef PNNonradialModeDriver DRIVER;
	const unsigned int num_var = DRIVER::num_var;
	typedef Mode<num_var> MODE;

	//this will be the relativity parameter used in calculation
	double zsurf = 1.0e-4;

	//the grid size of the polytrope
	int length =2048;
	
	//initialize a series of polytropes
	STAR *star;
	DRIVER *driv;
	MODE *mode;
	
	double Gamma1 = 5./3.;
	
	//for each star, check ten modes in three L
	double checkCont[10][3][10];
	double checkNewt[10][3][10];
	double checkPois[10][3][10];
	
	FILE *fp;
	
	for(int i=0; i<10; i++){
		//calculate polytrope
		star = new STAR(n[i], zsurf, length);
		driv = new DRIVER(star, Gamma1);
		printf("g");
		star->writeStar();
		int Len = driv->length();
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
				mode = new MODE(k,l,l, driv);
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
				double cee2=1;
				for(int X=4; X<Len-4; X++){
					int XX = 2*X;
					//stellar variables, to simplify equations
					double rho    = star->rho(XX);
					double P      = star->P(XX);
					double drhodr = star->drhodr(XX);
					double dPdr   = star->dPdr(XX);
					double g      = star->dPhidr(XX);
					double gPN    = star->dPsidr(XX);
					double q      = g + gPN/cee2;
					double Phi    = star->Phi(XX);
					double r      = star->rad(XX);
					double A      = star->Schwarzschild_A(XX,Gamma1);
					//mode variables
					double freq2 = mode->getOmega2()*star->getC(XX)*q/r;
					double freq  = sqrt(freq2);
					double L2    = double(l*l+l);
					double xi   = r  *mode->getY(0,X);//r  *y1
					double chi  = r*q*mode->getY(1,X);//r*g*y2
					double DPhi = r*q*mode->getY(2,X);//r*g*y3
					double dDPhi=   q*mode->getY(3,X);//  g*y4
					double DPsi = r*q*mode->getY(4,X)*star->Phi(Len-1);
					double dDPsi=   q*mode->getY(5,X)*star->Phi(Len-1);
					double DWR  = q*star->Phi(Len-1)/(4.*freq)*mode->getY(6,X);//divide by i
					double DWH  = q*star->Phi(Len-1)/(4.*freq)*mode->getY(7,X);//divide by i
					double Drho = (
							-2.*P*DPhi*rho - 4.*r*freq*DWH*rho*rho - DPsi*rho*rho
							+2.*P*rho*chi - 4.*rho*rho*Phi*chi + P*xi*dPdr
						)/(cee2*P*Gamma1) + (
							rho*rho*(chi-DPhi) + xi*rho*dPdr - P*Gamma1*xi*drhodr
						)/(P*Gamma1);
					double DP = Gamma1*P/(rho+P/cee2)*Drho + Gamma1*P*xi*A;
		
					//calculate numerical derivatives
					double b3,b2,b1,a1,a2,a3, h1;
					h1= star->rad(XX+2)-star->rad(XX );
					// d2DPhi/dr2 = (d/dr)(dDphi/dr), dDPhi/dr = g y4
					b3=(star->dPhidr(XX-6)+star->dPsidr(XX-6)/cee2)*mode->getY(3,X-3);
					b2=(star->dPhidr(XX-4)+star->dPsidr(XX-4)/cee2)*mode->getY(3,X-2);
					b1=(star->dPhidr(XX-2)+star->dPsidr(XX-2)/cee2)*mode->getY(3,X-1);
					a1=(star->dPhidr(XX+2)+star->dPsidr(XX+2)/cee2)*mode->getY(3,X+1);
					a2=(star->dPhidr(XX+4)+star->dPsidr(XX+4)/cee2)*mode->getY(3,X+2);
					a3=(star->dPhidr(XX+6)+star->dPsidr(XX+6)/cee2)*mode->getY(3,X+3);
					double d2DPhi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//next dxi/dr, xi = r y1
					b3=star->rad(XX-6)*mode->getY(0,X-3);
					b2=star->rad(XX-4)*mode->getY(0,X-2);
					b1=star->rad(XX-2)*mode->getY(0,X-1);
					a1=star->rad(XX+2)*mode->getY(0,X+1);
					a2=star->rad(XX+4)*mode->getY(0,X+2);
					a3=star->rad(XX+6)*mode->getY(0,X+3);
					difxi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//now dDP/dr, DP = rho*(chi-DPhi), chi=rg*y2, DPhi=rg*y3
					b3 = driv->DPfunc(X-3, freq, mode);
					b2 = driv->DPfunc(X-2, freq, mode);
					b1 = driv->DPfunc(X-1, freq, mode);
					a1 = driv->DPfunc(X+1, freq, mode);
					a2 = driv->DPfunc(X+2, freq, mode);
					a3 = driv->DPfunc(X+3, freq, mode);
					dDP = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//Now calculate residuals
					//Perturbed Poisson equation
					e1 = fabs(
							4.0*m_pi*Drho
							+ L2*pow(r,-2)*DPhi 
							- 2.0*dDPhi/r  - d2DPhi
					);
					n1 = fabs(4.0*m_pi*Drho )
							+ fabs( L2*pow(r,-2)*DPhi )
							+ fabs( 2.0*dDPhi/r ) + fabs( d2DPhi )
					;
					//continuity equation
					e2 = fabs(
							Drho*(1.-2.*Phi/cee2) - 3.*rho*DPhi/cee2
							+ xi*(drhodr+dPdr/cee2-2.*drhodr*Phi/cee2-2.*rho*g/cee2)
							+ (rho+P/cee2-2.*rho*Phi/cee2)*( 2.*xi/r + difxi - L2*pow(r,-2)/freq2*chi )		
					);
					n2 = fabs(Drho*(1.-2.*Phi/cee2)) + fabs(3.*rho*DPhi/cee2)
							+ fabs(xi*drhodr)+fabs(xi*dPdr/cee2)+2.*fabs(xi*drhodr*Phi/cee2) + 2.*fabs(xi*rho*g/cee2)
							+ ( fabs(rho+P/cee2)+fabs(2.*rho*Phi)/cee2 )*( fabs(2.*xi/r)+fabs(difxi)+fabs(L2*pow(r,-2)/freq2*chi) )
					;
					//newton's equation -- the r component. The theta component defines chi
					e3 = fabs( -(rho+P-2.*rho*Phi)*freq2*xi + dDP + rho*dDPhi + g*Drho
							+ rho*dDPsi + Drho*gPN + 2.*Phi*dDP + 2.*dPdr*DPhi
							+ 4.*rho*freq*DWR + (2.*rho*Phi+P)*dDPhi + (2.*rho*DPhi+2.*Drho*Phi+DP)*g);
					n3 = ( fabs(rho+P)+fabs(2.*rho*Phi) )*fabs(freq2*xi) + fabs( dDP ) + fabs( g*Drho ) + fabs( rho*dDPhi ) 
							+ fabs(rho*dDPsi) + fabs(Drho*gPN) + fabs(2.*Phi*dDP)
							+ fabs(4.*rho*freq*DWR) + fabs(2.*rho*Phi+P)*fabs(dDPhi) + fabs(2.*Drho*Phi+DP)*fabs(g);
					//normalize residuals
					e1 = e1/n1;
					e2 = e2/n2;
					e3 = e3/n3;
					//collect residuals
					checkPois[i][l-l0][k] += e1*e1;
					checkCont[i][l-l0][k] += e2*e2;
					checkNewt[i][l-l0][k] += e3*e3;
					double R = star->Radius();
					fprintf(fp, "%le\t%le\t%le\t%le\t%le\t%le\t%le\n", mode->getRad(X), e1, e2, e3, xi, chi, DPhi);
				}
				checkPois[i][l-l0][k] = sqrt(checkPois[i][l-l0][k]/(Len-8));
				checkCont[i][l-l0][k] = sqrt(checkCont[i][l-l0][k]/(Len-8));
				checkNewt[i][l-l0][k] = sqrt(checkNewt[i][l-l0][k]/(Len-8));
				
				printf("Root-Mean-Square Fractional Error in Continuity   :\t%le\n", checkCont[i][l-l0][k]);
				printf("Root-Mean-Square Fractional Error in Newton's  Law:\t%le\n", checkNewt[i][l-l0][k]);
				printf("Root-Mean-Square Fractional Error in Poisson's Law:\t%le\n", checkPois[i][l-l0][k]);
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
		delete driv;
	}
	return 0;
}