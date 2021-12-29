// **************************************************************************************
// GRPulseIO.h
//	This defines the protocols for handling use input and program output
//	Will take a filename form the commandline as an input file
//**************************************************************************************

#ifndef GRPULSEIOH
#define GRPULSEIOH

#include "GRPulseMain.h"

//will read user input from the specified file and organize the input data
int read_input(char input_file_name[128], CalculationInputData&);

//will re-print the user input in order to re-run the same calculation
int echo_input(CalculationInputData &);

//prepare to write output
int setup_output(CalculationInputData&, CalculationOutputData&);

//just prints a boc around certain parts of the output file
void print_splash(FILE *fp, char*, int WIDTH);

//write output information about the star itself
int write_stellar_output(CalculationOutputData&);

//write output data about the calculated modes
int write_mode_output(CalculationOutputData&);

//write both stellar and mode output
int write_output(CalculationInputData&);

//write tidal overlap coefficients in separate file
int write_tidal_overlap(CalculationOutputData&);

#endif