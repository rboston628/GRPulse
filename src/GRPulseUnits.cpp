//**************************************************************************************
//GRPulseUnits.cpp
//	This file contains all fuctions for handling unit conversions for the program
//**************************************************************************************

#include "GRPulseMain.h"

#ifndef GRPULSEUNITSH
#define GRPULSEUNITSH


int format_units(CalculationOutputData& data){	
	//fix the unitset based on which united were specified
	switch(data.units){
		case units::astro:
			data.unitset.G = G_astro;
			data.unitset.C = C_astro;
			data.unitset.base_mass = MSOLAR;
			data.unitset.base_length = 100000.0; // 1 km in cm
			data.unitset.base_time = 1.0; //use base time of 1 second
			break;
		case units::geo:
			data.unitset.G = 1.0;
			data.unitset.C = 1.0;
			data.unitset.base_time = 1.0; //use base time of 1 second
			data.unitset.base_length = 1.0/C_CGS;// 1 lightsecond, in cm
			data.unitset.base_mass = 1.0/C_CGS;// mass equivalent to 1 lightsecond, in cm
			break;
		case units::SI:
			data.unitset.G = G_CGS*1.0e-3;
			data.unitset.C = C_CGS*1.0e-2;
			data.unitset.base_time = 1.0; //use base time of 1 second;
			data.unitset.base_length = 100.0; //use base length of 1m, in cm
			data.unitset.base_mass = 1000.0; //use base mass of 1kg, in g
			break;
		case units::CGS:	
			data.unitset.G = G_CGS;
			data.unitset.C = C_CGS;
			data.unitset.base_time = 1.0; //use base time of 1 second;
			data.unitset.base_length = 1.0; //use base length of 1cm
			data.unitset.base_mass = 1.0; //use base mass of 1g
			break;
	}
	
	//the user can specify certain parameters of the star for the model; this calculates the others
	switch(data.params){
		case (units::pmass|units::pradius):
			data.logg = log10(G_CGS*data.mass*data.unitset.base_mass*pow(data.unitset.base_length*data.radius,-2));
			data.zsurf = 1./sqrt( 1. - 2.*data.unitset.G*data.mass/(data.radius*pow(data.unitset.C,2)) ) - 1.;
			break;
		case (units::pmass|units::pzsurf):
			data.radius = 2.*data.unitset.G*data.mass*pow(data.unitset.C,-2)/(1.-pow(1.+data.zsurf,-2));
			data.logg = log10(G_CGS*data.mass*data.unitset.base_mass*pow(data.unitset.base_length*data.radius,-2));
			break;
		case (units::pmass|units::plogg):
			data.radius = sqrt(data.unitset.G*data.mass*pow(10.0,-data.logg));
			data.zsurf = 1./sqrt( 1. - 2.*pow(10.,data.logg)*data.radius*pow(data.unitset.C,-2) ) - 1.;
			break;
		case (units::pradius|units::pzsurf):
			data.mass = 0.5*pow(data.unitset.C,2)*data.radius/data.unitset.G*(1.-pow(1.+data.zsurf,-2));
			data.logg = log10(G_CGS*data.mass*data.unitset.base_mass*pow(data.unitset.base_length*data.radius,-2));
			break;
		case (units::pradius|units::plogg):
			data.mass = pow(data.radius,2)*pow(10.0,data.logg)/data.unitset.G;
			data.zsurf = 1./sqrt( 1. - 2.*pow(10.,data.logg)*data.radius*pow(data.unitset.C,-2) ) - 1.;
			break;
		case (units::pzsurf|units::plogg):
			//I don't know why this isn't implemented yet -- I don't feel like doing it now
			data.mass   = 0.0;
			data.radius = 0.0;
			break;
	}
	
	printf("\tmass=%lg\n", data.mass);
	printf("\tradius=%lg\n", data.radius);
	printf("\tzsurf=%lg\n\tlogg=%lg\n", data.zsurf, data.logg);
	double mass_CGS = data.mass*data.unitset.base_mass;
	double radius_CGS = data.radius*data.unitset.base_length;
	data.freq0 = sqrt(G_CGS*mass_CGS*pow(radius_CGS,-3));
	printf("\tfrequency scale = %lf\n", data.freq0);
	
	return 0;
}



#endif