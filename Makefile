# anaroot directory
TARTSYS=/opt/anaroot

all: makeraw

ROOTCFLAGS  = $(shell $(CFG_ROOT_BIN)root-config --cflags)
ROOTLIBS    = $(shell $(CFG_ROOT_BIN)root-config --libs)

CFLAGS = -I$(TARTSYS)/include
LFLAGS = -L$(TARTSYS)/lib

GXX = g++

.SUFFIXES: .cpp .o

.cpp.o:
	$(GXX) $(CFLAGS) $(ROOTCFLAGS) -c $<

makeraw: and.h makeraw.cpp ./src/ribf123_rawvar.h
	$(GXX) $(CFLAGS) $(LFLAGS) $(ROOTCFLAGS) $(ROOTLIBS) -lanacore -lXMLParser -o $@ makeraw.cpp

clean:
	rm -f *.cpp~ *.o makeraw
