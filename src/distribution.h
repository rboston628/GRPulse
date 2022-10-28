// **************************************************************************************
// constants.h
//		This file describes dependencies and includes many broadly-used 
//      physical constants, and certain broadly-used for the GRPulse code.  
// **************************************************************************************

//constants and dependencies
#ifndef DISTRIBUTION
#define DISTRIBUTION
#include "STARS/Star.h"
#include "MODES/Mode.h"
//specific Newtonian stars
#include "STARS/Polytrope.h"
#include "STARS/ChandrasekharWD++.h"
#include "STARS/MESA.h"
//for 1PN models of stars
#include "STARS/PNStar.h"
#include "STARS/PNPolytrope.h"
#include "STARS/PNChandrasekharWD++.h"
//for GR models of stars
#include "STARS/GRStar.h"
#include "STARS/GRPolytrope.h"
//mode drivers
#include "MODES/ModeDriver.h"
#include "MODES/CowlingModeDriver.h"
#include "MODES/NonradialModeDriver.h"
//for 1PN mode calculations
#include "MODES/PNNonradialModeDriver.h"
//for GR mode calculations
#include "MODES/GRCowlingModeDriver.h"
#endif



