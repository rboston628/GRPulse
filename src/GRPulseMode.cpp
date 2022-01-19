//**************************************************************************************
//							GR PULSE MODE
// GRPulseMode.cpp
//		Handles all functionality with finding modes and calculating their properies
//**************************************************************************************


#include "GRPulseMain.h"
#include "GRPulseIO.h"

double compare_JCD(double, int, int, double);
double compare_Pekeris(double, int, int, double);

// the following arrays come from tables in JCD-DJM paper:
// 		Christensen-Dalsgaard & Mullan (1994), MNRAS 270, 921-935
// these were very useful for testing in early stages of the code

//columns l=2 and l=3 from Table 1 in JCD-DJM paper
const double JCD1_5[3][35] = {
	{256.8042, 424.9943, 582.9618, 736.2824, 886.9146, 1035.8195, 1183.5437, 1330.4248, 1476.6835, 1622.4704,
		1767.8912, 1913.0229, 2057.9220, 2202.6318, 2347.1851, 2491.6079, 2635.9207, 2780.1399, 2924.2788,
		3068.3483, 3212.3576, 3356.3141, 3500.2242, 3644.0932, 3787.9257, 3931.7257, 4075.4965, 4219.2411,
		4362.9620, 4506.6614, 4650.3412, 4794.0033, 4937.6492, 5081.2801, 5224.8974},
	{320.2718, 484.2210, 641.7601, 795.5461, 946.8945, 1096.5604, 1245.0155, 1392.5718, 1539.4439, 1685.7838,
		1831.7026, 1977.2825, 2122.5862, 2267.6620, 2412.5474, 2557.2726, 2701.8615, 2846.3336, 2990.7050,
		3131.9888, 3279.1960, 3423.3360, 3567.4165, 3711.4442, 3855.4249, 3999.3634, 4143.2640, 4287.1305,
		4430.9660, 4574.7734, 4718.5552, 4862.3137, 5006.0508, 5149.7683, 5293.4677},
	{369.0420, 534.0912, 693.1482, 848.4539, 1001.2131, 1152.1531, 1301.7468, 1450.3169, 1598.0915, 1745.2365,
		1891.8754, 2038.1018, 2183.9880, 2329.5905, 2474.9541, 2620.1149, 2765.1019, 2909.9392, 3054.6465,
		3199.2403, 3343.7343, 3488.1403, 3632.4682, 3776.7265, 3920.9226, 4065.0627, 4209.1525, 4353.1996,
		4497.1996, 4641.1649, 4785.0958, 4928.9953, 5072.8660, 5216.7101, 5360.5299}
};

//columns l=2 and l=3 from Table 2 in JCD-DJM paper
const double JCD3_0[3][35] = {
	{337.2152, 463.5718, 590.0694, 716.6289, 843.1066,  969.4331, 1095.5832, 1221.5530, 1347.3484, 1472.9799,
		1598.4594, 1723.7990, 1849.0105, 1974.1051, 2099.0928, 2223.9831, 2348.7844, 2473.5043, 2598.1500,
		2722.7276, 2847.2430, 2971.7011, 3096.1067, 3220.4639, 3344.7767, 3469.0483, 3593.2820, 3717.4807,
		3841.6468, 3965.7827, 4089.8906, 4213.9725, 4338.0301, 4462.0650, 4586.0789},
	{390.1223, 516.1992, 643.0677, 769.7802, 896.2910, 1022.6135, 1148.7622, 1274.7491, 1400.5849, 1526.2795,
		1651.8424, 1777.2825, 1902.6082, 2027.8277, 2152.9485, 2277.9776, 2402.9216, 2527.7867, 2652.5785,
		2777.3023, 2901.9629, 3026.5648, 3151.1121, 3275.6085, 3400.0577, 3524.4627, 3648.8266, 3773.1522,
		3897.4418, 4021.6979, 4145.9226, 4270.1179, 4394.2856, 4518.4275 ,4642.5450},
	{428.8391, 558.2981, 686.8732, 814.7027, 942.0267, 1068.9866, 1195.6639, 1322.1090, 1448.3552, 1574.4265,
		1700.3411, 1826.1136, 1951.7561, 2077.2793, 2202.6923, 2328.0032, 2453.2197, 2578.3482, 2703.3949,
		2828.3655, 2953.2651, 3078.0983, 3202.8695, 3327.5828, 3452.2418, 3576.8500, 3701.4106, 3825.9264,
		3950.4002, 4074.8346, 4199.2318, 4323.5941, 4447.9235, 4572.2218, 4696.4908}
};

//columns l=2 and l=3 from Table 3 in JCD-DJM paper
const double JCD4_0[3][35] = {
	{507.0621, 571.8190, 625.9711, 736.0301, 848.1925,  960.7919, 1073.6477, 1186.7036, 1299.9278, 1413.2959,
		1526.7873, 1640.3842, 1754.0710, 1867.8342, 1981.6625, 2095.5458, 2209.4757, 2323.4448, 2437.4470,
		2551.4770, 2665.5302, 2779.6026, 2893.6909, 3007.7923, 3121.9042, 3236.0246, 3350.1515, 3464.2834,
		3578.4189, 3692.5569, 3806.6962, 3920.8360, 4034.9756, 4149.1143, 4263.2516},
	{648.2038, 711.8195, 791.8287, 874.4382, 933.9077, 1024.6134, 1131.6998, 1242.3945, 1354.4090, 1467.0790,
		1580.1333, 1693.4374, 1806.9155, 1920.5210, 2034.2233, 2148.0012, 2261.8392, 2375.7260, 2489.6527,
		2603.6125, 2717.5999, 2831.6102, 2945.6399, 3059.6858, 3173.7454, 3287.8163, 3401.8968, 3515.9851,
		3630.0799, 3744.1800, 3858.2842, 3972.3917, 4086.5016, 4200.6132, 4314.7259},
	{713.1345, 782.3745, 833.5703, 936.7684, 996.2500, 1059.2311, 1169.2703, 1282.1336, 1395.5967, 1509.3076,
		1623.1501, 1737.0730, 1851.0500, 1965.0657, 2079.1108, 2193.1790, 2307.2658, 2421.3680, 2535.4830,
		2649.6088, 2763.7438, 2877.8864, 2992.0355, 3106.1901, 3220.3492, 3334.5120, 3448.6777, 3562.8457,
		3677.0153, 3791.1861, 3905.3575, 4019.5291, 4133.7005, 4247.8714, 4362.0413}
};

//this function will find eigenmode and frequency for the desired k,l in the output data
// requires a mode driver type (the <class>)
//  defined below
template <class MODE> int mode_finder(CalculationOutputData&);

//this will create the correct mode driver and pass the correct driver type to the mode_finder
int create_modes(CalculationOutputData &data_out){	
	switch(data_out.regime){
		case regime::PN0:
			switch(data_out.modetype){
				case modetype::cowling:
					data_out.modetype = modetype::cowling;
					data_out.driver = new CowlingModeDriver(data_out.star, data_out.adiabatic_index);
					mode_finder<CowlingModeDriver>(data_out);
					break;
				case modetype::radial:
				case modetype::quasinormal:
				case modetype::nonradial:
					data_out.modetype = modetype::nonradial;
					data_out.driver = new NonradialModeDriver(data_out.star, data_out.adiabatic_index);			
					mode_finder<NonradialModeDriver>(data_out);
					break;
			}
			break;
		case regime::PN1:
			//there is only one type of mode implemented in the 1pn case
			data_out.modetype = modetype::nonradial;
			data_out.driver = new PNNonradialModeDriver(reinterpret_cast<PNStar*>(data_out.star), data_out.adiabatic_index);
			mode_finder<PNNonradialModeDriver>(data_out);
			break;
		case regime::GR:
			//there is only one type of mode implemented in the GR case
			data_out.modetype = modetype::cowling;
			data_out.driver = new GRCowlingModeDriver(reinterpret_cast<GRStar*>(data_out.star), data_out.adiabatic_index);
			mode_finder<GRCowlingModeDriver>(data_out);
			break;
	}
		
	return 0;
}

//**************************************************************************************
// THE MODE FINDER ALGORITHM	
// This function is the true workhorse of GRPulse.
// This will organize the specified mode numbers, then begin looking through the list
// If it does not first find the desired mode, it will use a bisection search to find it
// If an accidentally discovered mode is in a list, it will save it
// Mode output is printed at the end of each run of constant L
//**************************************************************************************
template <class MODEDRIVER> int mode_finder(CalculationOutputData &data){
	//ensure that the passed class is a daughter of Mode
	static_assert((std::is_base_of<ModeDriver,MODEDRIVER>::value), "failure");
	
	typedef Mode<MODEDRIVER::num_var> MODE;
	static_assert((std::is_base_of<ModeBase,MODE>::value), "failure");
	
	printf("num_var = %d\n", MODE::num_var);
	
	//information on the nodes that need to be found
	int num = data.mode_num;
	int l_list[num];
	//for psuedorandom points inside bracket intervals
	int a=53, b=122, r=17737, ok = 39; //these are arbitrary but not magical -- cycle ~ r-1
	//good ol' 2*pi
	double twopi = 6.283185307179586;

	//flags indicating whether error columns should be calculated for test cases
	double index=data.input_params[0];
	bool isIsopycnic = (data.model==model::polytrope) & (index==0.0);
	bool isJCD = (data.model==model::polytrope) &
		(data.regime==regime::PN0) &
		(index==1.5 | index==3.0 | index==4.0) &
		(fabs(data.adiabatic_index - 5./3.)<1.e-5);
	bool comp1PN = (data.model==model::polytrope) & (data.regime==regime::PN1 | data.regime==regime::GR);

	
	//produce a list of the different L asked for, each represented once
	int nextl=0;
	l_list[0] = data.l[0];
	for(int j=1; j<num; j++){
		//if the last mode had the same L, keep going
		if(data.l[j] == l_list[nextl]) continue;
		//otherwise look through rest of list to see if same L occured earlier
		else{
			bool already = false;
			for(int i=0; (i<=nextl) & !already; i++){
				if(data.l[j] == l_list[i]) already = true;
			}
			//if it occured earlier keep going, otherwise save it
			if(already) continue;
			else {
				nextl++;
				l_list[nextl] = data.l[j];
			}
		}
	}
	//now for each L, produce a sublist of desired modes
	int nextmode = 0, klo=100000, khi=-100000;
	for(int j=0; j<=nextl; j++){
		int nextk=0, kl[num]; //kl is the sublist of desired modes
		for(int i=0; i<num; i++){	//save in list kl all K-values with this L
			if(data.l[i] == l_list[j]){
				kl[nextk] = data.k[i];
				nextk++;			//the number of k to calculate
			}
		}
		//also save values for the min,max values of K in list
		for(int i=0; i<nextk; i++){
			if(kl[i]>khi) khi = kl[i];
			if(kl[i]<klo) klo = kl[i];
		}

		//****************************************************************************
		// Here the real work begins
		// This will run through each k sublist, attempting new modes
		//   The calculation does not always return the desired mode
		//   However, the returned mode is often in the list later on
		//   This adjusts the search for correct k, while also saving those found by
		//     accident to reduce the amount of recalculation
		//   Search uses a bracketed bisection search on omega2.
		//****************************************************************************
		printf("L=%d,K\n", l_list[j]); fflush(stdout);
		//now we need to find the (L,K) modes given by (l_list[j], kl[i])
		int ktry[nextk];	 //the list of our calculated K
		double w2try[nextk]; //the corresponding frequencies
		MODE *modetry[nextk];
		for(int i=0; i<nextk; i++) ktry[i] = 0;//zero all elements for easy identification
		for(int i=0; i<nextk; i++){
			//STEP 1:  check if we have already been filled from earlier
			if(ktry[i] != 0 && w2try[i] > 0.0) {
				//if we have then leave
				printf("\tK=%d\tk=%d\talready found!\n", kl[i], ktry[i]);
				continue;
			}
			//STEP 2:  if not filled from earlier, try creating a mode
			//make a mode and find its mode-order k
			printf("\tK=%d\t", kl[i]); fflush(stdout);
			modetry[i] = new MODE(kl[i], l_list[j], 0, data.driver);
			ktry[i]  = modetry[i]->modeOrder();
			w2try[i] = modetry[i]->getOmega2();
			//is this the one we want?
			if(ktry[i] == kl[i]){
				printf("\t%lf\n", w2try[i]);
				continue;
			}
			//STEP 3:  if our trial mode is not the one we want, keep looking
			else {
				//STEP 3a:  check if the mode we found appears later in the list
				printf("%d\tlooking in list. ", ktry[i]); fflush(stdout);
				bool inlist = false;
				for(int ii=i+1; (ii<nextk) & !inlist; ii++){
					//if it is in the list, save it
					if(ktry[i] == kl[ii]) { 
						printf("found in list. "); fflush(stdout);
						inlist = true;
						modetry[ii] = modetry[i];
						modetry[i] = NULL; //necessary for proper memory management
						ktry[ii] = ktry[i];
						w2try[ii] = w2try[i];
					}
				}
				//STEP 3b: create brackets to be used in our bisection search
				//  we may be able to use a previously-found omega2 as one bracket
				printf("finding brackets.\n"); fflush(stdout);
				double w2min=0.0, w2max=0.0, dw2=0.0, w2in=0.0, w2out=0.0;
				int kk = ktry[i], kmax=khi+1, kmin=klo-1;
				//scan list of modes to find brackets
				for(int j=0; j<nextk; j++){
					//if the list of previously found modes (ktry) contains a lower mode
					//then use that mode as a minimum bracket
					if(ktry[j]!=0 & ktry[j] < kl[i] & ktry[j]>kmin){
						kmin  = ktry[j];
						w2min = w2try[j];
					}
					//if the lsit of previous found modes contains a higher mode
					//then use that mode as a maximum bracket
					if(ktry[j]!=0 & ktry[j] > kl[i] & ktry[j]<kmax){
						kmax  = ktry[j];
						w2max = w2try[j];
					}
				}
				//bracket search when discovered mode is HIGHER than desired mode
				if(kk > kl[i]){		
					//if the max bracket was not found in list
					//then use the current mode as a max bracket (since it is high)
					if(kmax==khi+1){
						kmax  = kk;
						w2max = w2try[i];
					}
					//if the moin bracket was not found in list
					//then use an absolute minimum as a min bracket
					// (note: this puts limits on allowed gmodes at -999999999)
					if(kmin==klo-1){
						w2min = 0.0;
						kmin = -1000000000;
					}	
				}
				//bracket search when discovered mode is LOWER than desired mode
				else if(kk < kl[i]){
					//printf("low bracket: k=%d, w2=%le\n", ktry[i], w2try[i]);
					//if the min bracket was not found in list
					//then use the current mode as a min bracket
					if(kmin==klo-1){
						kmin  = kk;
						w2min = w2try[i];
					}
					dw2 = w2min;
					//if the max bracket was not found in list, search for it
					if(kmax==khi+1){
						//start at current mode and increase
						kmax  = (kk>kmin? kk : kmin);
						w2max = (kk>kmin? w2try[i]+dw2 : w2min+dw2);
					}
					//increase and search until we find a max bracket
					double incr = 2.0;
					while(kmax < kl[i]){
						delete modetry[i];
						kmin = kmax;
						w2min = w2max;
						w2max = incr*w2max;
						modetry[i] = new MODE(w2max, l_list[j],0,data.driver);
						kmax = modetry[i]->modeOrder();
						printf("\t\t(%d,%d) in (%f,%f)\n",kmin, kmax, w2min, w2max);
						if(isnan(w2min) | isnan(w2max)) return 1;
						//if we found it, quit
						inlist = false;
						if(kmax == kl[i]){
							kk = kmax;
							w2try[i] = modetry[i]->getOmega2();
							break;
						}
						//check if our trial mode is anywhere else in the list
						else for(int ii=i+1;(ii<nextk) & !inlist; ii++){
							if((kmax==kl[ii]) & (ktry[ii]==0)) {
								printf("\t\t\t%d in list\n", kmax);
								inlist = true;
								modetry[ii] = modetry[i];
								modetry[i] = NULL;
								ktry[ii] = kmax;
								w2try[ii] = modetry[ii]->getOmega2();
							}
						}
						if(w2max > 1e6) {
							printf("TOO LARGE\n");
							incr = 0.5;
						}
						if(w2max < 1e-6){
							printf("No Dice\n");
							break;
						}
					}
				}
				//swap brackets if backward
				if(kmin >= kmax){
					int tk = kmax;
					kmax = kmin;
					kmin = tk;
					double tw = w2max;
					w2max = w2min;
					w2min = tw;
				}
				
				//STEP 3c:  now we have brackets -- these SHOULD put bounds in frequency
				//  for w2 in (w2min, w2max), will produce k in (kmin, kmax)
				//  reiterate bisection search until desired mode is found
				printf("\t\tbracketed (%d,%d) in (%f,%f)\n", kmin, kmax, w2min, w2max);
				double prevmin=w2min, prevmax=w2max; //value of previous frequency
				int stop=0; //integer to limit number of iterations
				while(kk != kl[i]){
					//STEP 3c (i): bisect the brackets 
					w2in = 0.5*(w2min+w2max); //bisect the brackets
					//STEP 3c (ii): create a trial mode
					delete modetry[i];
					modetry[i] = new MODE(w2in, l_list[j],0,data.driver);
					kk = modetry[i]->modeOrder();
					w2out = modetry[i]->getOmega2();
					printf("%d\t\t(%d,%d) in (%f,%f)\t %d\t %d %lf-->%lf\n",
							l_list[j],kmin, kmax, w2min, w2max, kl[i], kk, w2in, w2out);
					//STEP 3c (iii): compare the trial mode to desired mode
					//if we found it, move on to next
					if(kk == kl[i]) {
						break;
					}
					//if we didn't find it, see if either bracket can be moved
					if(     kk > kl[i] & kk <= kmax & w2max>w2out & w2out>0.0){
						kmax = kk;
						w2max = w2in;
					}
					else if(kk < kl[i] & kk >= kmin & w2min<w2out & w2out>0.0){
						kmin = kk;
						w2min = w2in;
					}
					
					//STEP 3c (iv): check if the trial mode is anywhere in the list
					inlist = false;
					for(int ii=i+1;(ii<nextk)&!inlist; ii++){
						if((kk==kl[ii]) & (ktry[ii]==0) & (kl[ii]!=kl[i]) & (w2out>0.0) & (kl[ii]!=0)) {
							printf("\t\t\t%d in list\n", kk);
							inlist = true;
							modetry[ii] = modetry[i];
							modetry[i ] = NULL;
							ktry[ii] = kk;
							w2try[ii] = modetry[ii]->getOmega2();
							stop=0;
						}
					}	
					
					//if the sought mode is bracketed between two known frequencies
					// then we can use special Mode constructor
					if((kmin == kl[i]-1) & (kmax == kl[i]+1)){
						delete modetry[i];
						modetry[i] = new MODE(w2min, w2max, l_list[j],0,data.driver);
						kk = modetry[i]->modeOrder();
						w2out = modetry[i]->getOmega2();
					}
					
					//STEP 3c (v):  check if the brackets have moved since last iteration
					//if not, pick a random location within brackets to test
					//repeat until one of the brackets can move
					//this prevents us from getting stuck
					double w2maxT = w2max, w2minT = w2min;
					int enough=0;
					while(kk!=kl[i] && (w2min==prevmin)&&(w2max==prevmax)){
						//pick a pseudo-random place in brackets
						ok = (a*ok+b)%r; //generates a psuedo-random integer in (0,r)
						w2in = w2minT + (double(ok)/double(r))*(w2maxT-w2minT);
						delete modetry[i];
						modetry[i] = new MODE(w2in, l_list[j],0,data.driver);
						kk = modetry[i]->modeOrder();
						w2out = modetry[i]->getOmega2();
						printf("\tR\t(%d,%d) in (%f,%f)\t %d\t %d %lf-->%lf\n",
								kmin, kmax, w2minT, w2maxT,kl[i], kk, w2in, w2out);
						//check if we found mode, or if we can move brackets
						if(kk==kl[i]) break;
						if(     kk > kl[i] & kk <= kmax & w2out>0.0 & w2out>w2min){
							kmax = kk;
							w2maxT = w2out;
						}
						else if(kk < kl[i] & kk >= kmin & w2out>0.0 & w2out<w2max){
							kmin = kk;
							w2minT = w2out;
						}
						//accounts for fact multiple w2in lead to same k
						if(w2in < w2maxT & kk == kmax & w2in>0.0){
							w2maxT = w2in;
						}
						else if(w2in > w2minT & kk == kmin & w2in>0.0){
							w2minT = w2in;
						}
						//sometimes zeros are inaccessible
						//move brackets and try again
						if(fabs(w2minT-w2maxT)<1e-2*w2minT) {
							//kk = kl[i];
							w2min = w2minT;
							w2max = w2maxT;
							//w2out = (w2min+w2max)*0.5;
							break;
						}
						//cancel if we have tried more than 5 random spots
						if(enough++ > 5) {
							w2min=w2minT;
							w2max=w2maxT;
							break;
						}
					}
					//update past values
					prevmin = w2min, prevmax=w2max;
										
					//if we are just unable to find the mode, say so
					if(kk!=kl[i] & fabs(w2max-w2min) < 1e-2*w2min){
						modetry[i] = NULL;
						w2try[i] = 0.0;
						ktry[i] = 0;
						printf("too close\t%le\n", fabs((w2max-w2min)/w2max) );
						break;
					}
					
					if(stop++ > 10) {
						break;
					}
				}
				//if we were unable to find the mode, the w2try[i] is set to 0
				//  do not save the values in that case
				if(w2try[i]==0.0) continue;
				//now save the mode numbers, frequency
				ktry[i] = kk;
				w2try[i] = modetry[i]->getOmega2();
			}
			printf("\t\tk=%d\n", ktry[i]);		
		}
		
		//STEP 4:  return through list to pick up any missed modes
		//  this calculation follows same steps as STEP 3 but because more
		//   modes have been filled, we are more likely to have good brackets
		printf("%d\tpicking up stragglers...\n", l_list[j]);
		for(int i=0; i<nextk; i++){
			//STEP 4a:  check if we have already been filled from earlier
			if(ktry[i] == kl[i] && w2try[i] > 0.0) {
				continue;
			}
			//STEP 4b: if we have not, fill
			else {
				//STEP 4b (i): find brackets to use in bisection by scanning list of modes
				//  this pass through, we are more likely to have bracketing modes
				printf("\t\t%d\t", kl[i]); fflush(stdout);
				double w2min=0.0, w2max=0.0, dw2=0.0, w2in=0.0, w2out=0.0;
				int kk = ktry[i], kmax=khi+1, kmin=klo-1;
				for(int j=0; j<nextk; j++){
					//if the list of previously found modes (ktry) contains a lower mode
					//then use that mode as a minimum bracket
					if(ktry[j]!=0 & ktry[j] < kl[i] & ktry[j]>kmin){
						kmin  = ktry[j];
						w2min = w2try[j];
					}
					//if the list of previous found modes contains a higher mode
					//then use that mode as a maximum bracket
					if(ktry[j]!=0 & ktry[j] > kl[i] & ktry[j]<kmax){
						kmax  = ktry[j];
						w2max = w2try[j];
					}
				}
				//if we don't have a lower bound from the list, use absolute minimum
				if(kmin==klo-1){ 
						w2min = 0.0;
						kmin = -1000000000;
				}
				//if we don't have a higher bound from calculated list, just leave
				if(kmax==khi+1) { printf("no brackets\n"); continue;}
				printf("(%d,%d)\t", kmin,kmax);
				
				//STEP 4b (ii): now bisect the brackets until correct mode is found
				double prevmin=w2min, prevmax=w2max; //value of previous frequency
				int stop=0; //integer to limit number of iterations
				while(kk != kl[i]){
					w2in = 0.5*(w2min+w2max); //bisect the brackets
					//create a trial mode
					delete modetry[i];
					modetry[i] = new MODE(w2in, l_list[j],0,data.driver);
					kk = modetry[i]->modeOrder();
					w2out = modetry[i]->getOmega2();
					//if we found it, then great.  move on to next
					if(kk == kl[i]) {
						break;
					}
					//if we didn't find it, see if either bracket can be moved
					if(     kk > kl[i] & kk <= kmax & w2max >w2out & w2out>0.0){
						kmax = kk;
						w2max = w2out;
					}
					else if(kk < kl[i] & kk >= kmin & w2min<w2out & w2out>0.0){
						kmin = kk;
						w2min = w2out;
					}
					//if the sought mode is bracketed between two known frequencies
					// then we can use special Mode constructor
					if((kmin == kl[i]-1) & (kmax == kl[i]+1)){
						delete modetry[i];
						modetry[i] = new MODE(w2min, w2max, l_list[j],0,data.driver);
						kk = modetry[i]->modeOrder();
						w2out = modetry[i]->getOmega2();
					}
					
					//check if the brackets have moved since last time
					//if not, pick a random location within brackets to test
					//repeat until one of the brackets can move
					double w2maxT = w2max, w2minT = w2min;
					int enough=0;
					while(kk!=kl[i] && (w2min==prevmin)&&(w2max==prevmax)){
						//if scanning didn't work, pick a pseudo-random place in brackets
						ok = (a*ok+b)%r; //generates a psuedo-random integer in (0,r)
						w2in = w2minT + (double(ok)/double(r))*(w2maxT-w2minT);
						delete modetry[i];
						modetry[i] = new MODE(w2in, l_list[j],0,data.driver);
						kk = modetry[i]->modeOrder();
						w2out = modetry[i]->getOmega2();
						//check if we found mode, or if we can move brackets
						if(kk==kl[i]) break;
						if(kk > kl[i] & kk <= kmax & w2out>0.0){
							kmax = kk;
							w2maxT = w2out;
						}
						else if(kk < kl[i] & kk >= kmin & w2out>0.0){
							kmin = kk;
							w2minT = w2out;
						}
						//accounts for fact multiple w2in lead to same k
						if(w2in < w2maxT & kk == kmax & w2in>0.0){
							w2maxT = w2in;
						}
						else if(w2in > w2minT & kk == kmin & w2in>0.0){
							w2minT = w2in;
						}
						//sometimes zeros are inaccessible
						//move brackets and try again
						if(fabs(w2minT-w2maxT)<1e-2*w2minT) {
							w2min = w2minT;
							w2max = w2maxT;
							break;
						}
						//cancel if we have tried more than 20 random spots
						if(enough++ > 5) {
							w2min=w2minT;
							w2max=w2maxT;
							break;
						}
					}		
						
					//update past values
					prevmin = w2min, prevmax=w2max;
					if(stop++ > 5) {
						break;
					}
				}
				//now save the mode numbers, frequency
				ktry[i] = kk;
				w2try[i] = modetry[i]->getOmega2();
				if(kk==kl[i]) printf("found\n");
				else printf("not found\n");
			}
		}
		printf("\tdone\n");
		
		//STEP 5: organize the data lists
		int enext = data.mode_done;
		for(int i=0; i<nextk; i++){
			//STEP 5a: if a mode was not found, say so and save nothing to it
			if(ktry[i]!=kl[i] | w2try[i]==0.0){
				printf("unable to find mode %d,%d\n", l_list[j], kl[i]);
				data.l[nextmode] = l_list[j];
				data.k[nextmode] = kl[i];
				data.w[nextmode] = 0.0;
				data.mode[nextmode] = NULL;
				data.f[nextmode] = 0.0;
				data.period[nextmode] = 0.0;
				data.mode_SSR[nextmode] = 0.0;
				nextmode++;
				data.mode_done++;
				continue;
			}
			//STEP 5b: otherwise save all data
			data.l[nextmode] = l_list[j];
			data.k[nextmode] = kl[i];
			data.w[nextmode] = sqrt(w2try[i]);
			if(sqrt(w2try[i])<0.0) data.w[nextmode] = -sqrt(-w2try[i]);
			data.mode[nextmode] = modetry[i];
			//freq0 converts to rad/s, divide by 2pi to get frequency in Hertz
			data.f[nextmode] = data.freq0*data.w[nextmode]/(twopi);
			//period in seconds is 1/f
			data.period[nextmode] = 1./data.f[nextmode];
			data.mode_SSR[nextmode] = modetry[i]->SSR();			
			nextmode++;
			data.mode_done++;
		}
	
		//STEP 6: calculate desired errors
		int e=0;
		//STEP 6a: for n=0 polytropes, we can compare frequencies to the exact Pekeris formula
		if(isIsopycnic){
			for(int i=enext;i<data.mode_done; i++){
				if(data.k[i] >=0)
					data.err[e][i] = compare_Pekeris(data.w[i], data.l[i], data.k[i], data.adiabatic_index);
				else data.err[e][i] = nan("");
			}
			e++;
		}
		//STEP 6b: for certain polytrope frequencies, we can compare to tables in JCD-DJM paper
		if(isJCD) {
			for(int i=enext;i<data.mode_done; i++){
				//if((data.l[i]==1 | data.l[i]==2 | data.l[i]==3) & (data.k[i]>0 & data.k[i]<36))
					data.err[e][i] = compare_JCD(index, data.l[i], data.k[i], data.w[i]);
				//else data.err[e][i] = nan("");
			}
			e++;
		}
		//STEP 6c: for 1PN polytropes, we can compare against frequencies from Newtonian polytropes
		//  will calculate a Newtonian polytrope, then try to match the corresponding mode
		//  this will only try to match it once, using the 1PN frequency to start
		if(comp1PN){
			Polytrope *star0PN = new Polytrope(1.0,1.0, index, data.Ngrid);
			NonradialModeDriver *drv0PN = new NonradialModeDriver(star0PN, data.adiabatic_index);
			const int nn = NonradialModeDriver::num_var;
			Mode<nn> *testmode;
			for(int i=enext; i<data.mode_done; i++){
				testmode = new Mode<nn>(data.w[i]*data.w[i], data.l[i], 0, drv0PN);
				double w0pn = sqrt(testmode->getOmega2());
				if(testmode->modeOrder() == data.k[i]){
					data.err[e][i] = 2.*fabs(w0pn-data.w[i])/(w0pn+data.w[i]);
				}
				else {
					data.err[e][i] = nan("");
				}
				delete testmode;
			}
			delete star0PN;
			delete drv0PN;
			e++;
		}
		if(e!=data.i_err) printf("Error in mode error listing...\n");
		
		//STEP 7: at the end of each L, print all data to the output file 
		write_mode_output(data);
	}
		
	return 0;
}


double compare_JCD(double n, int l, int k, double w){
	//convert dimensionless freuqueny to same scale used in JCD-DJM paper
	w = round(w*nug*10000.0)/10000.0;
	double fJCD = 0.0;
	//for polytropes n=1.5, n=3, n=4, compare to tables,l=1,2,3, with 1<=k<=35
	if((l==1 | l==2 | l==3) & (k>0 & k<36)){
		switch(int(n*10)){
			case 15:
				fJCD = JCD1_5[l-1][k-1];
				break;
			case 30:
				fJCD = JCD3_0[l-1][k-1];
				break;
			case 40:
				fJCD = JCD4_0[l-1][k-1];
				break;
		}
		return abs(w-fJCD);
	}
	else return nan("");
}

//a single function to check against the Pekeris formula
double compare_Pekeris(double w, int l, int k, double Gam1){
	double wPek2 = 0.0, dnl = Gam1*double(k)*(double(k+l)+0.5)-2.;
	wPek2 = dnl + sqrt(dnl*dnl + double(l*l+l));
	//there is also a formula for f-modes due to Chandrasekhar (1964), c.f. Cox (1980) eq 17.80
	if(k==0) wPek2 = double(2.*l*(l-1))/double(2*l+1);
	return fabs(w*w - wPek2)/wPek2;
}