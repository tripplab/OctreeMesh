include ../BuildFlags.mk

SRCDIR:=../../../source
EXE:=../../../gid/problemtypes/ElectricPotential.gid/ElectricPotential

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/METIS -I$(SRCDIR)/GiDProblemTypes/ElectricPotential
LDFLAGS:=-L.. -lFEMT -lMETIS

CXX:=g++

all: $(EXE)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE)

######################################################################
# ElectricPotential ##################################################
######################################################################
$(EXE): main.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

main.o: $(SRCDIR)/GiDProblemTypes/ElectricPotential/main.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

