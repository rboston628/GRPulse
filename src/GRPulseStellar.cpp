//**************************************************************************************
//							GR PULSE STELLAR
// GRPulseStellar.cpp
//		Handles all functionality with finding stars and calculating their properies
//**************************************************************************************


#include "GRPulseMain.h"

//functions to create the different stellar models
int create_classical_polytrope(CalculationOutputData &);
int create_classical_CHWD(CalculationOutputData &);
int create_classical_MESA(CalculationOutputData &);
int create_1pn_polytrope(CalculationOutputData &);
int create_1pn_CHWD(CalculationOutputData &);
int create_gr_polytrope(CalculationOutputData &);

//will select the correct creation fuction from above based on user input
int create_star(CalculationOutputData &data_out){
	printf("Creating star...\n");
	switch(data_out.regime){
		case regime::PN0:
			switch(data_out.model) {
				case model::polytrope:
					create_classical_polytrope(data_out);
					break;
				case model::CHWD:
					create_classical_CHWD(data_out);
					break;
				case model::MESA:
					create_classical_MESA(data_out);
					break;
			}
			break;
		case regime::PN1:
			switch(data_out.model) {
				case model::polytrope:
					create_1pn_polytrope(data_out);
					break;
				case model::CHWD:
					create_1pn_CHWD(data_out);
					break;
				default:
					create_1pn_polytrope(data_out);
					break;
			}
			break;
		case regime::GR:
			switch(data_out.model){
				case model::polytrope:
					create_gr_polytrope(data_out);
					break;
				default:
					printf("ERROR: This stellar model is not implemented for GR\n");
					return 1;
			}
			break;
			
	}
	return 0;
}

//create a classical polytrope
int create_classical_polytrope(CalculationOutputData& data){
	double mass_CGS = data.mass*data.unitset.base_mass;
	double radius_CGS = data.radius*data.unitset.base_length;
	//create star	
	data.star = new Polytrope(mass_CGS, radius_CGS, data.input_params[0], data.Ngrid);

	//calculate error estimate for this model
	data.star_SSR = data.star->SSR();
	printf("SSR = %le\n", data.star_SSR);
	
	return 0;
}

//create a post-newtonian polytrope
int create_1pn_polytrope(CalculationOutputData& data){
	//the user-supplied radius should be taken as the Schwarzschild radius
	//the PNPolytrope code uses isotropic coordinates
	//to rectify, convert the radius being used in the frequency scale
	double Rschwarz = data.radius;
	double zsurf = data.zsurf;
	double Riso = Rschwarz*(1.+zsurf);
	double mass_CGS = data.mass*data.unitset.base_mass;
	double radius_CGS = Riso*data.unitset.base_length;
	data.freq0 = sqrt(G_CGS*mass_CGS*pow(radius_CGS,-3));
	printf("\tfrequency scale = %lf\n", data.freq0);
	//now create the star	
	data.star = new PNPolytrope(data.input_params[0], data.zsurf, data.Ngrid);
	
	//calculate error estimate for this model
	data.star_SSR = data.star->SSR();
	printf("SSR = %le\n", data.star_SSR);
	
	return 0;
}

//create a general relativistic polytrope
int create_gr_polytrope(CalculationOutputData& data){
	data.star = new GRPolytrope(data.input_params[0], data.zsurf, data.Ngrid);

	//calculate error estimate for this model
	data.star_SSR = data.star->SSR();
	printf("SSR = %le\n", data.star_SSR);
	
	return 0;
}


//create a classical WD following Chandrasekhar
int create_classical_CHWD(CalculationOutputData& data){
	data.star = new ChandrasekharWD(data.input_params[0], data.input_params[1], data.Ngrid);

	//adjust inputs to match the actual values
	data.units = units::CGS;
	data.mass = data.star->Mass();
	data.radius = data.star->Radius();
	format_units(data);

	//calculate error estimate for this model
	data.star_SSR = data.star->SSR();
	printf("SSR = %le\n", data.star_SSR);
	
	return 0;
}

//create a post-newtonian WD following Chandrasekhar
int create_1pn_CHWD(CalculationOutputData& data){
	//the user-supplied radius should be taken as the Schwarzschild radius
	//the PNPolytrope code uses isotropic coordinates
	//to rectify, convert the radius being used in the frequency scale
	double Rschwarz = data.radius;
	double zsurf = data.zsurf;
	double Riso = Rschwarz*(1.+zsurf);
	double mass_CGS = data.mass*data.unitset.base_mass;
	double radius_CGS = Riso*data.unitset.base_length;
	data.freq0 = sqrt(G_CGS*mass_CGS*pow(radius_CGS,-3));
	//now create the star	
	data.star = new PNChandrasekharWD(data.input_params[0], data.zsurf, data.Ngrid);
	
	//calculate error estimate for this model
	data.star_SSR = data.star->SSR();
	printf("SSR = %le\n", data.star_SSR);
	
	return 0;
}

//create a classical WD found in MESA
int create_classical_MESA(CalculationOutputData& data){
	char inputname[255];
	sprintf(inputname, "%s.dat", data.str_input_param.c_str() );
	printf("source=%s.dat\n", data.str_input_param.c_str());
	data.star = new MESA(inputname, data.Ngrid);
	
	//adjust the inputs around the fact this is a MESA object
	data.units = units::CGS;			//MESA uses CGS units
	data.mass = data.star->Mass();		//the mass is determined by model
	data.radius = data.star->Radius();	//the radius is determined by model
	data.Ngrid = data.star->length();
	data.params = units::pmass|units::pradius;
	format_units(data);
	
	//calculate error estimate for this MESA model
	data.star_SSR = data.star->SSR();
	printf("SSR = %le\n", data.star_SSR);
	
	return 0;
}



