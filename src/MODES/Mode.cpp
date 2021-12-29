//**************************************************************************************
//							The MODE Object
//	Mode.cpp
//		Equation agnostic 
//			-- all information specific to physics supplied by the driver
// 		Capable of nonradial, Cowling, 1PN modes by supplying different drivers
//  NOTE: If you make ANY changes to this file, you MUST make clean and recompile
//**************************************************************************************

#ifndef MODECLASS
#define MODECLASS

#include "Mode.h"

//This is the basic setup routine common to all constructors
// preparing all arrays for integration and matching
//The individual constructors differ only in how they initialize frequency
template <size_t  numvar>
void Mode<numvar>::basic_setup(){
	cee2 = star->light_speed2();
	Gee = star->Gee();

	//begin by setting all values to 1 -- will be changed later in code
	for(int i=0;i<numvar;i++){
		yCenter[i]  = 1.0;
		ySurface[i] = 1.0;
	}
	//we set up stellar grid as twice the mode grid to avoid interpolation at half points
	len_star = star->length();
	len      = driver->length();
	Gamma1   = driver->Gamma1();
	
	//determine the index where fit for inward and outward integrations will occur
	//some experimentation suggests this should be closer to surface for large l
	int imax = len;
	R = star->rad(len_star-1);
	double R80 = 0.8*R;	//at 80% of stellar radius
	for(int i=0;i<len-1;i++) if(star->rad(2*i) > R80) {imax=i;break;}
	double b = double(l)/double(l+1);
	//based on l, set the fit infex at an intermediary point between stellar index and R80
	xfit = star->indexFit;//int( b*imax + (1.0-b)*star->indexFit );
	//set the arrays that define the perturbation
	rad = new double[len];
	y = new double*[numvar];
	for(int i=0;i<numvar;i++) y[i] = new double[len];
	//define values of radius for the mode
	for(int X=0; X<len; X++){
		rad[X] = driver->rad(X);
	}
	//define the conversion factor from frequency in rad/s to dimensionless; w=f*sig2omeg
	sig2omeg = pow(star->Radius(),3)/(Gee*star->Mass());
	
	//set up boundary matrix and index array
	double **yy = new double*[numvar];
	int *ind = new int[numvar];
	for(int i=0;i<numvar;i++) yy[i] = new double[numvar];
	driver->getBoundaryMatrix(num_var, yCenter, ySurface, yy, ind);
	for(int i=0;i<numvar;i++) for(int j=0;j<numvar;j++) boundaryMatrix[i][j]=yy[i][j];
	for(int i=0;i<numvar;i++) indexOrder[i] = ind[i];
}

//This constructor is given an initial value of frequency to use
//This initial guess is used as a starting point in search for frequency
template <size_t numvar> 
Mode<numvar>::Mode(double omg2, int l, int m, ModeDriver *drv)
//	 : ModeBase(drv) 
	: l(l), m(m), omega2(omg2), driver(drv), star(drv->star)
{	
	basic_setup();
	//now we  converge to solution
	converge();
	this->k = verifyMode();
}


//This constructor guesses initial frequency based on a desired k,l mode
//The initial guess is taken from the analytic solutions for n=0 polytrope
template <size_t numvar> 
Mode<numvar>::Mode(int K, int L, int M, ModeDriver *drv)
//	 : ModeBase(drv)
	: l(L), m(M), k(K), driver(drv), star(drv->star)
{	
	basic_setup();
	
	double ell = double(l);
	omega2 = 0.0;
	//for p-modes
	if(k>=0){
		//for initial guess, may as well use homogeneous frequency
		//see Cox 17.7, homogeneous model
		//homogeneous model frequencies usually larger, so reduce size of n used
		// there is no magic to this guess -- it just seems to work okay
		double ks = double(k - (2*l-4));		
		// 2Dn = -4 + Gamma1*[ n*(2l+2n+5) + 2l + 3]
		double Dn = Gamma1*( ks*(ks+ell+0.5) ) - 2.0;
		if(Gamma1==0.0) Dn = star->Gamma1(0)*( ks*(ks+ell+0.5) ) - 2.0;
		//omega^2 = Dn + sqrt( Dn^2 + l(l+1) )
		omega2 = Dn + sqrt( Dn*Dn + ell*(ell+1.0) );
		//if(omega2 < 0) omega2 = Dn - sqrt( Dn*Dn + ell*(ell+1.0) );
		//from dimensional considerations, sigma^2 = (GM/R^3) *omega^2
	}
	//for g-modes
	else {
		//find base sigma2 when n=-1 (use above formula, simplify)
		double sig1 = -2.0 + sqrt(4. + ell*(ell+1.));
		//convert to freq
		sig1 = sqrt( sig1*nug ); //fq = fg*sigma
		//empirical formula freq_n = freq_-1 * exp(-0.4 sqrt(-n-1)) for n<0
		//determined for polytropes -- works less well with realistic WD stars
		sig1 = sig1*exp(-0.4*sqrt(-k));
		//now convert freq back to sigma2
		omega2 = (sig1*sig1)/nug; //sigma^2
		//omega2 = omega2/sig2omeg;
	}

	//now we  converge to solution
	//by matching outward and inward solutions in the center
	converge();
	this->k = verifyMode();
}

//This constructor is given a range within which the desired frequency is known to exist
//The range must bound exactly one eigenfrequency
template <size_t numvar> 
Mode<numvar>::Mode(double omeg2lo, double omeg2hi, int l, int m, ModeDriver *drv)
//	 : ModeBase(drv)
	: l(l), m(m), driver(drv), star(drv->star)
{	
	basic_setup();

	//BEGIN BISECTION HUNT
	converged = false;
	int a=53, b=122, r=17737, ok = 39;//for pseudorandom bracket division
	//if omegas are in wrong order, swap 'em!
	if (omeg2lo > omeg2hi){
		double tmp = omeg2hi;
		omeg2hi = omeg2lo;
		omeg2lo = tmp;
	}
	else if (omeg2lo==omeg2hi){	//this should not occur
		omeg2lo *= 0.9;
	}
	
	//the brackets given are often themselves zeros -- inch them closer to avoid this
	double dw = (omeg2hi-omeg2lo)*1.0e-3;
	double wmin=omeg2lo+dw, wmax=omeg2hi-dw;
	//now find the initial bracketing values of the Wronskian
	double Wsmin = RK4center(omeg2lo, yCenter, ySurface);
	double Wsmax = RK4center(omeg2hi, yCenter, ySurface);
	
	//if the bracketing values do not bound a zero, just use the mid-point as a starting value
	omega2 =  0.5*(omeg2lo + omeg2hi);
	if(Wsmin*Wsmax>0.0) {
		converge();
		this->k = verifyMode();
		return;
	}
	
	//if the brackets are fine, we now begin a biection search
	int stop=0;
	double w1 = 0.5*(omeg2lo + omeg2hi), w2 = omeg2hi;
	double W1 = RK4center(w1, yCenter, ySurface), W2;
	while( fabs(w2-w1) > 0.0 ){
		w2 = w1;
		W2 = W1;
		//test each bracket
		if( W1*Wsmax > 0.0 ){
			if( w1 < omeg2hi ){
				omeg2hi = w1;
				Wsmax = W1;
			}
		}
		else if( W1*Wsmin > 0.0 ){
			if( w1 > omeg2lo ){
				omeg2lo = w1;
				Wsmin = W1;
			}
		}
		w1 = 0.5*(omeg2hi+omeg2lo);
		W1 = RK4center(w1, yCenter, ySurface);
		//if the new Wronskian is same as old, pick a random value inside brackets
		if(W2==W1){
			ok = (a*ok+b)%r; //generates a psuedo-random integer in (0,r)
			w1 = omeg2lo + (double(ok)/double(r))*fabs(omeg2hi-omeg2lo);
			W1 = RK4center(w1, yCenter, ySurface);
		}
		if(++stop>20){
			converge();
			this->k = verifyMode();
			return;
		}
	}
	//we now know our value of omega2
	omega2=w1;
	linearMatch(omega2, yCenter, ySurface);
	//if the Wronskian still isn't good, just find values as normal
	if(W1 > 1e-10) converge();
	else converged = true;
	this->k = verifyMode();
}

//destructor
template <size_t numvar> 
Mode<numvar>::~Mode(){
		//free space used by y1,...,y4 from stack
		delete[] rad;
		for(int i=0; i<numvar; i++)
			delete[] y[i];
}

template <size_t numvar>
void Mode<numvar>::converge(){
	converged = false;
	convergeBisect(0.0);
	linearMatch(omega2, yCenter,ySurface);
	converged = true;
}

//Find frequency using a bisection search, based on Wronskian, up to tolerance tol
template <size_t numvar>
void Mode<numvar>::convergeBisect(double tol){
	//the apparently magical numbers are arbitrary
	double w1  = omega2, w2, W1, W2;
	//brackets
	double wmax, wmin, wdx;
	double Wmax, Wmin, Wdx;
	int stop = 0, bnd = 200;
	int a=53, b=122, r=17737, ok = 39;//for pseudorandom positioning
	
	W1 = W2 = RK4center(w1, yCenter,ySurface);	
	//while the two Ws are on same side of axis, keep reposition until zero is bound
	//we use Newton's method to the nearest zero
	w2=w1;
	while(W1*W2 >= 0.0){
		//compute numerical derivative
		wdx = 1.01*w2;
		Wdx = RK4center(wdx, yCenter,ySurface);\
		wdx = w2 - W2*(wdx-w2)/(Wdx-W2);
		//Newton's method can fail at this if it only approaches the zero from one direction
		//In such a case, slightly broaden the bracket to get to other side of zero
		if( wdx == w2 ) {
			if(wdx>w1) wdx *= 1.01;
			if(wdx<w1) wdx *= 0.99;
		}
		//limit amount of change allowed in single step
		if(wdx > 1.1*w2) wdx = 1.1*w2;	//if we increased, don't increase too much
		if(wdx < 0.9*w2) wdx = 0.9*w2;	//if we decreased, don't decrease too much
		Wdx = W2;
		W2 = RK4center(wdx, yCenter, ySurface);
		w2 = wdx;
	}
	//sometimes the brackets will be on either end of hyperbolic divergence
	//check that solution changes continuously between the two brackets
	double w3 = 0.5*(w1+w2), W3 = RK4center(w3, ySurface, yCenter);
	double scale = 1.01;
	while( (W3*W1>=0.0 & fabs(W3)>fabs(W1)) | (W3*W2>=0.0 & fabs(W3)>fabs(W2))){
		//this can sometimes be solved by taking broader steps in secant method
		scale *= 1.1;
		w2=w1;
		W2=W1;
		while(W1*W2 >= 0.0){
			//compute numerical derivative
			wdx = scale*w2;
			Wdx = RK4center(wdx, yCenter,ySurface);
			wdx = w2 - W2*(wdx-w2)/(Wdx-W2);
			//Newton's method can fail at this if it only approaches the zero from one direction
			//In such a case, slightly broaden the bracket to get to other side of zero
			if( wdx == w2 ) {
				if(wdx>w1) wdx *= 1.01;
				if(wdx<w1) wdx *= 0.99;
			}
			//limit amount of change allowed in single step
			//if(wdx > 1.1*w2) wdx = 1.1*w2;	//if we increased, don't increase too much
			//if(wdx < 0.9*w2) wdx = 0.9*w2;	//if we decreased, don't decrease too much
			Wdx = W2;
			W2 = RK4center(wdx, yCenter, ySurface);
			w2 = wdx;
		}
		w3 = 0.5*(w1+w2);
		W3 = RK4center(w3, yCenter, ySurface);
	}
	
	if(w2 > w1){
		wmax = w2; Wmax = W2;
		wmin = w1; Wmin = W1; 
	}
	else {
		wmax = w1; Wmax = W1;
		wmin = w2; Wmin = W2;
	}
	//swap if they are backwards
	if(wmin > wmax){
		double temp = wmax;
		wmax = wmin; wmin = temp;
		temp = Wmax;
		Wmax = Wmin; Wmin = temp;
	}
	
	//now bisect search inside the brackets
	w1 = 0.5*(wmin+wmax);
	W2 = RK4center(w1, yCenter,ySurface);
	stop=0;
	while( fabs(w2-w1) > tol ){
		w2 = w1;
		W1 = W2;
		//test s
		if( W1*Wmax > 0.0 ){
			if( w1 < wmax ){
				wmax = w1;
				Wmax = W1;
			}
		}
		else if( W1*Wmin > 0.0 ){
			if( w1 > wmin ){
				wmin = w1;
				Wmin = W1;
			}
		}
		w1 = 0.5*(wmin+wmax);
		W2 = RK4center(w1, yCenter, ySurface);
		if(W2==W1){//if the problem becomes stuck in a loop
			ok = (a*ok+b)%r; //generates a psuedo-random integer in (0,r)
			w1 = wmin + (double(ok)/double(r))*fabs(wmax-wmin);
			W2 = RK4center(w1, yCenter,ySurface);
		}
	}
	omega2 = w1;
}

//This is deprecated
//Newton convergence on omega2, up to tolerance tol, max number of steps term
//if term = 0, then no maximum number of steps (until integer overflow, anyway)
template <size_t numvar> 
void Mode<numvar>::convergeNewton(double tol, int term){
	int stop=1, fix=0, kick=3;
	double w1  = 0.0, w2 = omega2;
	double W1 = RK4center(omega2, yCenter,ySurface), W2;
	while( fabs(w1-w2) > tol){
		w1=w2;		
		w2 = 1.01*w1;
		W2 = RK4center(w2, yCenter,ySurface);
		if(W2==W1){
			printf("stux\n");
			omega2 = w1;
			return;
		}
		w2 = w1 - W1*(w2-w1)/(W2-W1);
		W1 = RK4center(w2, yCenter,ySurface);
		
		if(isnan(W1)) printf("NaN in W1 at step %d\n", stop);
		if(isnan(w2)) printf("NaN in w2 at step %d\n", stop);

		if( stop++ == term) break;
	}
	omega2 = w2;
	printf("W2 = %le", W1);
}

//linearly match inward and outward solutions
template <size_t numvar> 
void Mode<numvar>::linearMatch(double w2, double y0[numvar], double ys[numvar]){
	double DY[numvar][numvar];
	for(int i=0;i<numvar/2;i++){
		RK4out(xfit, w2, boundaryMatrix[i]  );
		for(int j=0;j<numvar;j++) DY[i][j] = y[j][xfit];
	}
	for(int i=numvar/2; i<numvar; i++){
		RK4in( xfit, w2, boundaryMatrix[i]);
		for(int j=0;j<numvar;j++) DY[i][j] = y[j][xfit];
	}
		
	//ALROGITHM TO FIND COEFFICIENTS
	double A[numvar][numvar];
	for(int i=0;i<numvar;i++){
		for(int j=0; j<numvar/2; j++)      A[i][j] = DY[j][i]; //outward solutions +
		for(int j=numvar/2; j<numvar; j++) A[i][j] =-DY[j][i]; // inward solutions -
	}
	
	double aa[numvar] = {0.0}; //for(int i=0; i<num_var; i++) aa[i] = 0.0;
	double bb[numvar] = {0.0}; //for(int i=0; i<num_var; i++) bb[i] = 0.0;
	invertMatrix(A, bb, aa);
	
	//for(int i=0; i<numvar; i++) printf("%lf ", aa[i]);
	//printf("\n");
	//for the basis BCs we chose, this will be the properly scaled physical solution
	//if we change the BCs, we must change these results to match
	for(int a=0; a<numvar/2; a++){
		y0[indexOrder[a]] *= aa[a];
	}
	for(int a=numvar/2; a<numvar-1; a++){
		ys[indexOrder[a]] *= aa[a];
	}
	ys[indexOrder[0]] = 1.0;
	
	
	RK4out(xfit, w2, y0);
	RK4in( xfit, w2, ys);
}

//integrates outward from interior to xmax
template <size_t numvar> 
void Mode<numvar>::RK4out(int xmax, double w2, double y0[numvar]){
	//intermediate values used in equations
	double XC=0.0, YC[num_var]={0.0};
	//contain shifts dy1 = dx*y1', dy2 = dx*y2', etc.
	double K[4][num_var];
	//Butcher tableau for RK4
	static const double b[4] = {0.5, 0.5, 1.0, 0.0};
	static const double c[4] = {1.0, 2.0, 2.0, 1.0};
	static const int    d[4] = {0, 1, 1, 2};
	
	//step size, based on star's grid
	double dx = rad[1];
	//coefficients for derivatives
	double coeff[num_var][num_var];

	//set initial values
	int start = driver->CentralBC(y, y0, w2, l);
	
	//now begin RK4
	for(int x = start; x<xmax; x++){
		//begin using grid values for first correction
		XC = rad[x];
		for(int i=0;i<num_var;i++) YC[i] = y[i][x];
		//dx = log(rad[x+1]/rad[x]);
		dx = rad[x+1]-rad[x];
		for(int a = 0; a<4; a++){
			driver->getCoeff( &coeff[0][0],x,d[a],w2,l);
			//calculate shift for next step
			for(int i=0;i<num_var;i++){
				K[a][i] = 0.0;
				for(int j=0;j<num_var;j++) K[a][i] += dx*coeff[i][j]*YC[j]/XC;
			}
			//evaluate next correction from shift
			XC = rad[x] + b[a]*dx;
			for(int i=0;i<num_var;i++) YC[i] = y[i][x] + b[a]*K[a][i];
		}
		//calculate value of solution at next step by averaging corrections
		for(int i=0;i<num_var;i++){
			y[i][x+1] = y[i][x] + (K[0][i]/6.+K[1][i]/3.+K[2][i]/3.+K[3][i]/6.);
		}
	}
}

//integrates inward from surface toward xmin
template <size_t numvar> 
void Mode<numvar>::RK4in( int xmin, double w2, double ys[numvar]){
	//the XC, Y1C-Y4C are intermediate values used in equations
	double YC[numvar], XC;
	//contain shifts dy1 = dx*y1', dy2 = dx*y2', etc.
	double K[4][numvar];
	//Butcher tableau for RK4
	static const double b[4] = {0.5, 0.5, 1.0, 0.0};
	static const double c[4] = {1.0, 2.0, 2.0, 1.0};
	//step backward, not forward
	static const int d[4] = {0,-1,-1,-2};

	//step size, based on star's grid
	double dx = rad[len-2] - rad[len-1];
	//coefficients for derivatives
	double coeff[numvar][numvar];

	//set initial values
	int start = driver->SurfaceBC(y, ys, w2, l);
	//now begin RK4
	for(int x = start; x>xmin; x--){
		//begin using grid values for first correction
		XC = rad[x];
		for(int i=0;i<num_var;i++) YC[i] = y[i][x];
		//dx = log(rad[x-1]/rad[x]);
		dx = rad[x-1] - rad[x];
		for(int a = 0; a<4; a++){
			driver->getCoeff(&coeff[0][0],x,d[a],w2,l);
			//calculate shift for next step
			for(int i=0;i<num_var;i++){
				K[a][i] = 0.0;
				for(int j=0;j<num_var;j++) K[a][i] += dx*coeff[i][j]*YC[j]/XC;
			}
			//evaluate next correction from shift
			XC = rad[x] + b[a]*dx;
			for(int i=0;i<num_var;i++) YC[i] = y[i][x] + b[a]*K[a][i];
		}
		//calculate value of solution at next step by averaging shifts
		for(int i=0;i<num_var;i++){
			y[i][x-1] = y[i][x] + (K[0][i]/6.+K[1][i]/3.+K[2][i]/3.+K[3][i]/6.);
		}
	}
}

//integrates from both edges toward center and returns Wronskian
template <size_t numvar> 
double Mode<numvar>::RK4center(double w2, double y0[numvar], double ys[numvar]){		
	//prepare a matrix
	double DY[numvar][numvar];
	//store values of outward solution at fitting point in rows
	for(int i=0;i<numvar/2;i++){
		RK4out(xfit, w2, boundaryMatrix[i]);		
		for(int j=0;j<numvar;j++) DY[i][j] = y[j][xfit];
	}
	//store values of inward solution at fitting point in rows
	for(int i=numvar/2; i<numvar; i++){
		RK4in( xfit, w2, boundaryMatrix[i]);		
		for(int j=0;j<numvar;j++) DY[i][j] = y[j][xfit];
	}
	
	//the Wronskian is the determinant of this matrix
	return determinant(DY);
}

//this method serves to verify that the n is indeed the desired mode number
// according to Scuflaire and Osaki classification scheme, counting zeros of xi
template <size_t numvar> 
int Mode<numvar>::verifyMode(){
	//trace through solution in (xi, chi) plane, parameterized by index
	int quad=0, quadP, N=0;
	for(int x=1; x<len-2; x++){
		quadP = quad;
		//determine quadrant of (xi, chi)
		quad = (y[0][x]>=0 ? (y[1][x]>0? 1 : 2 ) : (y[1][x]>=0? 4 : 3));
		if(quadP==quad) continue;
		//if solution rotates clockwise, count as negative modes (g modes)
		if(quadP == 1 & quad == 2) N--;
		if(quadP == 3 & quad == 4) N--;
		//if solution rotates counter-clockwise, positive modes (p modes)
		if(quadP == 2 & quad == 1) N++;
		if(quadP == 4 & quad == 3) N++;	
	}
	return N;
}


//print out the mode information and plot it on gnuplot
template <size_t numvar> 
void Mode<numvar>::writeMode(char *c){
	//create names for files to be opened
	char filename[256];
	char rootname[256];
	char txtname[256];
	char outname[256];
	if(c==NULL)	sprintf(filename, "./out/%s/mode", star->name);
	else{
		sprintf(filename, "./%s/modes", c);
	}
	sprintf(rootname, "%s/mode_%d.%d", filename, l,k);
	if(!converged) sprintf(rootname, "%sXX", rootname); 
	//save data to folder to avoid clutter - make sure folder exists
	sprintf(txtname, "%s.txt", rootname);
	sprintf(outname, "%s.png", rootname);
	FILE *fp;
	if(!(fp = fopen(txtname, "w")) ){
		char command[256];
		sprintf(command, "mkdir -p %s", filename);
		system(command);
		fp = fopen(txtname, "w");
	}
	double R = rad[len-1];
	double M = star->Mass();
	for(int x=0; x<len; x++){
		fprintf(fp, "%0.16le", rad[x]/R);
		//fprintf(fp, "%0.16le", log(1.-star->mr(x)/M));
		for(int a=0; a<num_var; a++) fprintf(fp, "\t%0.16le", y[a][x]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	//plot file in png in gnuplot, and open png
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set output '%s'\n", outname);
	char title[100]; star->graph_title(title);
	fprintf(gnuplot, "set title 'full mode %d,%d in %s, period=%0.5lf s'\n", l,k, title, getPeriod());
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'log|y|'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	//fprintf(gnuplot, "set xrange [-8:0]\n");
	std::string varname[numvar]; driver->varnames(varname);
	//fprintf(gnuplot, "set arrow 1 from 5e3,1e-6 to 2e5,15.6e-12 lc rgb 'red' nohead\n");
	fprintf(gnuplot, "set arrow 1 from %le, graph 0 to %le, graph 1 lc rgb 'red' nohead\n", rad[xfit]/R, rad[xfit]/R);
	fprintf(gnuplot, "plot ");
	for(int a=0; a<numvar; a++){
		fprintf(gnuplot, "%c '%s' u 1:(abs($%d)) w l t '%s'", (a==0? ' ':','),txtname, a+2, varname[a].c_str());
	}
	fprintf(gnuplot, "\n");

	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
}

//write the mode, and then open the png to screen for easy viewing
template <size_t numvar> 
void Mode<numvar>::printMode(char *c){
	writeMode(c);
	char outname[258];
	if(c==NULL) sprintf(outname, "./out/%s/mode/mode_%d.%d", star->name, l,k);
	else {
		sprintf(outname, "./%s/modes/mode_%d.%d", c, l,k);
	}
	if(!converged) sprintf(outname, "%sXX", outname);
	sprintf(outname, "%s.png", outname);
	char openmyplot[248];
	sprintf(openmyplot, "open %s", outname);
	system(openmyplot);
}

//ways to access the frequency
template <size_t numvar> 
double Mode<numvar>::getOmega2(){
	return omega2;
}
template <size_t numvar> 
double Mode<numvar>::getFreq(){
	return sqrt(omega2/sig2omeg);
}
template <size_t numvar> 
double Mode<numvar>::getPeriod(){
	return 2.*m_pi/getFreq();
}


//we want to be able to call SSR() on each Mode object
//however, SSR() requires equations from ModeDriver
//this is the best compromise
template <size_t numvar> 
double Mode<numvar>::SSR(){
	return driver->SSR(omega2, l, this);
}

template <size_t numvar>
double Mode<numvar>::tidal_overlap(){
	return driver->tidal_overlap(this);
}



#endif