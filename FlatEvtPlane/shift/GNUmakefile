OBJS = miniDst.o Param.o CentralityMaker.o StRefMultCorr.o produceShiftFactor.o
EXE =  produceShiftFactor

#OBJS = MuDst.o createphiwgt.o
#EXE = createphiwgt

#OBJS = MuDst.o StPhysicalHelixD.o StHelixD.o StTofGeometry.o createphiwgt.o
#EXE = createphiwgt

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++
FLAGS = -m32 -std=c++11 -fno-inline -Wall -g $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -m32 -std=c++11 -fno-inline -o $(EXE) $(OBJS) -lSpectrum $(ROOTFLAGS) $(ROOTLIBS)
%.o: %.cxx
	$(COMPILE)  $< 

