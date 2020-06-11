include ../BuildFlags.mk

SRCDIR:=../../../source
EXE:=../../../gid/problemtypes/Solid.gid/Solid.Schur

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/METIS -I$(SRCDIR)/GiDProblemTypes/Solid
LDFLAGS:=-L.. -lFEMT -lMETIS

CXX:=mpic++

all: $(EXE)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE)

######################################################################
# Solid ##############################################################
######################################################################
$(EXE): main.schur.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

main.schur.o: $(SRCDIR)/GiDProblemTypes/Solid/main.schur.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

