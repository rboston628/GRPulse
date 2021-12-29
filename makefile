#IDIR = include
IDIR=src
LDIR=lib
SDIR=src
ODIR=obj

.SUFFIXES:
.SUFFIXES: .cpp .o

CC=g++ -std=c++11
CFLAGS=-I$(IDIR) -Wuninitialized -Weffc++

## files needed to compile stellar models
#  dependencies
_STARDEPS = constants.h\
 	 STARS/Star.h \
		STARS/Polytrope.h \
		STARS/ChandrasekharWD.h \
		STARS/MESA.h\
 	 STARS/PNStar.h \
		STARS/PNPolytrope.h \
		STARS/PNChandrasekharWD.h \
	 STARS/GRStar.h \
		STARS/GRPolytrope.h
STARDEPS = $(patsubst %, $(IDIR)/%, $(_STARDEPS))
#  source
_STARSRC = STARS/Star.cpp \
		STARS/Polytrope.cpp STARS/PNPolytrope.cpp STARS/GRPolytrope.cpp \
		STARS/ChandrasekharWD.cpp STARS/PNChandrasekharWD.cpp \
		STARS/MESA.cpp
STARSRC  = $(patsubst %, $(SDIR)/%, $(_STARSRC))

## files needed to compile mode drivers
#  dependencies
_DRVDEPS = constants.h\
	STARS/Star.h STARS/PNStar.h\
	MODES/ModeDriver.h \
		MODES/NonradialModeDriver.h \
		MODES/CowlingModeDriver.h\
		MODES/PNNonradialModeDriver.h\
		MODES/GRCOwlingModeDriver.h 
DRVDEPS = $(patsubst %, $(IDIR)/%, $(_DRVDEPS))
#  source
_DRVSRC  = MODES/NonradialModeDriver.cpp \
		MODES/PNNonradialModeDriver.cpp \
		MODES/CowlingModeDriver.cpp \
		MODES/GRCowlingModeDriver.cpp
DRVSRC   = $(patsubst %, $(SDIR)/%, $(_DRVSRC))


## files needed to compile the mode object
#  dependencies
_MODEDEPS = constants.h\
	STARS/Star.h MODES/ModeDriver.h MODES/Mode.h
MODEDEPS = $(patsubst %, $(IDIR)/%, $(_MODEDEPS))
#  source
_MODESRC = MODES/Mode.cpp
MODESRC  = $(patsubst %, $(SDIR)/%, $(_MODESRC))


## files needed to compile main program
#  dependencies
_MAINDEPS = constants.h GRPulseMain.h GRPulseIO.h\
	STARS/Star.h \
		STARS/Polytrope.h \
		STARS/ChandrasekharWD.h \
		STARS/MESA.h \
  	STARS/PNStar.h \
		STARS/PNPolytrope.h \
		STARS/PNChandrasekharWD.h \
	STARS/GRStar.h \
		STARS/GRPolytrope.h \
  	MODES/Mode.h\
  	MODES/ModeDriver.h \
		MODES/CowlingModeDriver.h \
		MODES/NonradialModeDriver.h \
		MODES/PNNonradialModeDriver.h \
		MODES/GRCowlingModeDriver.h 
MAINDEPS = $(patsubst %, $(IDIR)/%, $(_MAINDEPS))
#  soure
_MAINSRC = GRPulseMain.cpp GRPulseIO.cpp GRPulseStellar.cpp GRPulseMode.cpp GRPulseUnits.cpp
MAINSRC  = $(patsubst %, $(SDIR)/%, $(_MAINSRC))


## prepare object names
STAROBJ = $(patsubst %.cpp,$(ODIR)/%.o, $(_STARSRC))
DRVOBJ  = $(patsubst %.cpp,$(ODIR)/%.o, $(_DRVSRC))
MODEOBJ = $(patsubst %.cpp,$(ODIR)/%.o, $(_MODESRC))
MAINOBJ = $(patsubst %.cpp,$(ODIR)/%.o, $(_MAINSRC))


## Main rule for GRPulse program
GRPulse:  $(MAINOBJ) $(MODEOBJ) $(STAROBJ) $(DRVOBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LDIR)/mylib.a

## Rules for each subsection -- only update if their dependencies change
$(STAROBJ): $(ODIR)/%.o: $(SDIR)/%.cpp $(STARDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(DRVOBJ): $(ODIR)/%.o: $(SDIR)/%.cpp $(DRVDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(MODEOBJ): $(ODIR)/%.o: $(SDIR)/%.cpp $(MODEDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(MAINOBJ): $(ODIR)/%.o: $(SDIR)/%.cpp  $(MAINDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)



.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(ODIR)/STARS/*.o $(ODIR)/MODES/*.o

