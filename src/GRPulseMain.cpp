//**************************************************************************************
//
//							GR PULSE
//	Main program.  All other functionality is called from here.
//	This program is not necessary to use the stellar models or the mode functions
//		defined in the src files, but provides a single and simple access to them
//
//**************************************************************************************


#include "distribution.h"
#include "GRPulseMain.h"
#include "GRPulseIO.h"

int main(int argc, char* argv[]){
	//Handle command-line information to extract the input file
	FILE* input_file;
	CalculationInputData calcdataIn;
	CalculationOutputData calcdataOut;
	char input_file_name[128];
	if(argc == 2) {
		//the calculation filename is sent as command-line argument
		sprintf(input_file_name, "%s", argv[1]);
		//call the read_input routine on the given filename (see GRPulseIO.h)
		if(!read_input(input_file_name, calcdataIn)) printf("file read\n");
		else return 1;
	}
	else {
		printf("Usage is GRPulse.out <input file>\n");
		return 1;
	}
	
	//remove any information from past calculations, and create directory for the calculation
	char command[256];
	sprintf(command, "rm -r ./output/%s", calcdataIn.calcname.c_str());
	system(command);
	sprintf(command, "mkdir ./output/%s", calcdataIn.calcname.c_str());
	system(command);
	sprintf(command, "mkdir ./output/%s/star", calcdataIn.calcname.c_str());
	system(command);
	sprintf(command, "mkdir ./output/%s/modes", calcdataIn.calcname.c_str());
	system(command);
	//print a copy of the input file for future reference
	echo_input(calcdataIn);
	//prepare the output, based on the input
	setup_output(calcdataIn, calcdataOut);

	//make the background stellar model
	create_star(calcdataOut);
	//write output on the bakground stellar model
	write_stellar_output(calcdataOut);
	
	create_modes(calcdataOut);
	write_mode_output(calcdataOut);
	if(calcdataOut.regime==regime::PN0 && calcdataOut.modetype==modetype::nonradial)
		write_tidal_overlap(calcdataOut);
	return 0;
}