// *************************************************************************************
//					POST-NEWTONIAN MODE RMSR SCALE TEST
// physical_scale.cpp
//		This Dziembowski equations solve for variables y_1,y_2,y_3,y_4
//		These are related to Drho, DP, etc.
//		This code takes the solution in terms of y_i, finds physical variabes, 
//			and puts them into the LAWE to compute a residual
//		We expect this RMSR to scale like N^{-4}
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

	//ask user for tolerance of calculation
	int fac = 2;
	int length1 = 2048;
	int length2 = fac*length1-(fac-1);
	
	//initialize a series of polytropes
	STAR *star1, *star2;
	DRIVER  *driv1, *driv2;
	MODE *mode1, *mode2;
	
	//adiabatic exponent
	double Gamma1 = 5./3.;
	
	//redshfit
	double zsurf = 1.e-4;
		
	FILE *fp;
	for(int i=0; i<10; i++){
		//initialize the stars
		star1 = new STAR(n[i],zsurf, length1);
		star2 = new STAR(n[i],zsurf, length2);
		driv1 = new DRIVER(star1, Gamma1);
		driv2 = new DRIVER(star2, Gamma1);
		
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
				mode1 = new MODE(2*(k+1),l,l, driv1);
				mode2 = new MODE(2*(k+1),l,l, driv2);
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
				double cee2=1;
				for(int X=4; X<Len-4; X++){
					int XX = 2*X;
					//stellar variables, to simplify equations
					double rho    = star1->rho(XX);
					double P      = star1->P(XX);
					double drhodr = star1->drhodr(XX);
					double dPdr   = star1->dPdr(XX);
					double g      = star1->dPhidr(XX);
					double gPN    = star1->dPsidr(XX);
					double q      = g + gPN/cee2;
					double Phi    = star1->Phi(XX);
					double r      = star1->rad(XX);
					double A      = star1->Schwarzschild_A(XX,Gamma1);
					//mode variables
					double freq2 = mode1->getOmega2()*star1->getC(XX)*q/r;
					double freq  = sqrt(freq2);
					double L2    = double(l*l+l);
					double xi   = r  *mode1->getY(0,X);//r  *y1
					double chi  = r*q*mode1->getY(1,X);//r*g*y2
					double DPhi = r*q*mode1->getY(2,X);//r*g*y3
					double dDPhi=   q*mode1->getY(3,X);//  g*y4
					double DPsi = r*q*mode1->getY(4,X)*star1->Phi(Len-1);
					double dDPsi=   q*mode1->getY(5,X)*star1->Phi(Len-1);
					double DWR  = q*star1->Phi(Len-1)/(4.*freq)*mode1->getY(6,X);//divide by i
					double DWH  = q*star1->Phi(Len-1)/(4.*freq)*mode1->getY(7,X);//divide by i
					double Drho = (
							-2.*P*DPhi*rho - 4.*r*freq*DWH*rho*rho - DPsi*rho*rho
							+2.*P*rho*chi - 4.*rho*rho*Phi*chi + P*xi*dPdr
						)/(cee2*P*Gamma1) + (
							rho*rho*(chi-DPhi) + xi*rho*dPdr - P*Gamma1*xi*drhodr
						)/(P*Gamma1);
					double DP = Gamma1*P/(rho+P/cee2)*Drho + Gamma1*P*xi*A;//*/
		
					//calculate numerical derivatives
					double b3,b2,b1,a1,a2,a3, h1;
					h1= star1->rad(XX+2)-star1->rad(XX );
					// d2DPhi/dr2 = (d/dr)(dDphi/dr), dDPhi/dr = g y4
					b3=(star1->dPhidr(XX-6)+star1->dPsidr(XX-6)/cee2)*mode1->getY(3,X-3);
					b2=(star1->dPhidr(XX-4)+star1->dPsidr(XX-4)/cee2)*mode1->getY(3,X-2);
					b1=(star1->dPhidr(XX-2)+star1->dPsidr(XX-2)/cee2)*mode1->getY(3,X-1);
					a1=(star1->dPhidr(XX+2)+star1->dPsidr(XX+2)/cee2)*mode1->getY(3,X+1);
					a2=(star1->dPhidr(XX+4)+star1->dPsidr(XX+4)/cee2)*mode1->getY(3,X+2);
					a3=(star1->dPhidr(XX+6)+star1->dPsidr(XX+6)/cee2)*mode1->getY(3,X+3);
					double d2DPhi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//next dxi/dr, xi = r y1
					b3=star1->rad(XX-6)*mode1->getY(0,X-3);
					b2=star1->rad(XX-4)*mode1->getY(0,X-2);
					b1=star1->rad(XX-2)*mode1->getY(0,X-1);
					a1=star1->rad(XX+2)*mode1->getY(0,X+1);
					a2=star1->rad(XX+4)*mode1->getY(0,X+2);
					a3=star1->rad(XX+6)*mode1->getY(0,X+3);
					difxi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//now dDP/dr, DP = rho*(chi-DPhi), chi=rg*y2, DPhi=rg*y3
					b3 = driv1->DPfunc(X-3, freq, mode1);
					b2 = driv1->DPfunc(X-2, freq, mode1);
					b1 = driv1->DPfunc(X-1, freq, mode1);
					a1 = driv1->DPfunc(X+1, freq, mode1);
					a2 = driv1->DPfunc(X+2, freq, mode1);
					a3 = driv1->DPfunc(X+3, freq, mode1);
					dDP = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//Now calculate residuals
					//Perturbed Poisson equation
					e1 = fabs(
							4.0*m_pi*Drho
							+ L2*pow(r,-2)*DPhi 
							- 2.0*dDPhi/r  - d2DPhi
					);
					//continuity equation
					e2 = fabs(
							Drho*(1.-2.*Phi/cee2) - 3.*rho*DPhi/cee2
							+ xi*(drhodr+dPdr/cee2-2.*drhodr*Phi/cee2-2.*rho*g/cee2)
							+ (rho+P/cee2-2.*rho*Phi/cee2)*( 2.*xi/r + difxi - L2*pow(r,-2)/freq2*chi )		
					);
					//newton's equation -- the r component. The theta component defines chi
					e3 = fabs( -(rho+P-2.*rho*Phi)*freq2*xi + dDP + rho*dDPhi + g*Drho
							+ rho*dDPsi + Drho*gPN + 2.*Phi*dDP + 2.*dPdr*DPhi
							+ 4.*rho*freq*DWR + (2.*rho*Phi+P)*dDPhi + (2.*rho*DPhi+2.*Drho*Phi+DP)*g);	
					//now for second star
					int fXX = fac*XX;
					int fX  = fac*X;
					rho    = star2->rho(fXX);
					P      = star2->P(fXX);
					drhodr = star2->drhodr(fXX);
					dPdr   = star2->dPdr(fXX);
					g      = star2->dPhidr(fXX);
					gPN    = star2->dPsidr(fXX);
					q      = g + gPN/cee2;
					Phi    = star2->Phi(fXX);
					r      = star2->rad(fXX);
					A      = star2->Schwarzschild_A(fXX,Gamma1);
					//mode variables
					xi   = r  *mode2->getY(0,fX);//r  *y1
					chi  = r*q*mode2->getY(1,fX);//r*g*y2
					DPhi = r*q*mode2->getY(2,fX);//r*g*y3
					dDPhi=   q*mode2->getY(3,fX);//  g*y4
					DPsi = r*q*mode2->getY(4,X)*star2->Phi(fac*Len-fac);
					dDPsi=   q*mode2->getY(5,X)*star2->Phi(fac*Len-fac);
					DWR  = q*star1->Phi(fac*Len-fac)/(4.*freq)*mode1->getY(6,X);//divide by i
					DWH  = q*star1->Phi(fac*Len-fac)/(4.*freq)*mode1->getY(7,X);//divide by i
					Drho = (
							-2.*P*DPhi*rho - 4.*r*freq*DWH*rho*rho - DPsi*rho*rho
							+2.*P*rho*chi - 4.*rho*rho*Phi*chi + P*xi*dPdr
						)/(cee2*P*Gamma1) + (
							rho*rho*(chi-DPhi) + xi*rho*dPdr - P*Gamma1*xi*drhodr
						)/(P*Gamma1);
					DP = Gamma1*P/(rho+P/cee2)*Drho + Gamma1*P*xi*A;//*/
					//calculate numerical derivatives
					h1= star1->rad(fXX+fac*2)-star1->rad(fXX );
					// d2DPhi/dr2 = (d/dr)(dDphi/dr), dDPhi/dr = g y4
					b3=(star2->dPhidr(fXX-fac*6)+star2->dPsidr(fXX-fac*6)/cee2)*mode2->getY(3,fX-fac*3);
					b2=(star2->dPhidr(fXX-fac*4)+star2->dPsidr(fXX-fac*4)/cee2)*mode2->getY(3,fX-fac*2);
					b1=(star2->dPhidr(fXX-fac*2)+star2->dPsidr(fXX-fac*2)/cee2)*mode2->getY(3,fX-fac*1);
					a1=(star2->dPhidr(fXX+fac*2)+star2->dPsidr(fXX+fac*2)/cee2)*mode2->getY(3,fX+fac*1);
					a2=(star2->dPhidr(fXX+fac*4)+star2->dPsidr(fXX+fac*4)/cee2)*mode2->getY(3,fX+fac*2);
					a3=(star2->dPhidr(fXX+fac*6)+star2->dPsidr(fXX+fac*6)/cee2)*mode2->getY(3,fX+fac*3);
					d2DPhi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//next dxi/dr, xi = r y1
					b3=star2->rad(fXX-fac*6)*mode2->getY(0,fX-fac*3);
					b2=star2->rad(fXX-fac*4)*mode2->getY(0,fX-fac*2);
					b1=star2->rad(fXX-fac*2)*mode2->getY(0,fX-fac*1);
					a1=star2->rad(fXX+fac*2)*mode2->getY(0,fX+fac*1);
					a2=star2->rad(fXX+fac*4)*mode2->getY(0,fX+fac*2);
					a3=star2->rad(fXX+fac*6)*mode2->getY(0,fX+fac*3);
					difxi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//now dDP/dr, DP = rho*(chi-DPhi), chi=rg*y2, DPhi=rg*y3
					b3 = driv2->DPfunc(fX-fac*3, freq, mode2);
					b2 = driv2->DPfunc(fX-fac*2, freq, mode2);
					b1 = driv2->DPfunc(fX-fac*1, freq, mode2);
					a1 = driv2->DPfunc(fX+fac*1, freq, mode2);
					a2 = driv2->DPfunc(fX+fac*2, freq, mode2);
					a3 = driv2->DPfunc(fX+fac*3, freq, mode2);
					dDP = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h1);
					//Now calculate residuals
					//Perturbed Poisson equation
					n1 = fabs(
							4.0*m_pi*Drho
							+ L2*pow(r,-2)*DPhi 
							- 2.0*dDPhi/r  - d2DPhi
					);
					//continuity equation
					n2 = fabs(
							Drho*(1.-2.*Phi/cee2) - 3.*rho*DPhi/cee2
							+ xi*(drhodr+dPdr/cee2-2.*drhodr*Phi/cee2-2.*rho*g/cee2)
							+ (rho+P/cee2-2.*rho*Phi/cee2)*( 2.*xi/r + difxi - L2*pow(r,-2)/freq2*chi )		
					);
					//newton's equation -- the r component. The theta component defines chi
					n3 = fabs( -(rho+P-2.*rho*Phi)*freq2*xi + dDP + rho*dDPhi + g*Drho
							+ rho*dDPsi + Drho*gPN + 2.*Phi*dDP + 2.*dPdr*DPhi
							+ 4.*rho*freq*DWR + (2.*rho*Phi+P)*dDPhi + (2.*rho*DPhi+2.*Drho*Phi+DP)*g);				
					//normalize residuals
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