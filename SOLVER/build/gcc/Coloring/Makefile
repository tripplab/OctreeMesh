include ../BuildFlags.mk

SRCDIR:=../../../source
EXE:=../../../gid/problemtypes/Coloring.gid/Coloring

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/METIS -I$(SRCDIR)/GiDProblemTypes/Coloring
LDFLAGS:=-L.. -lFEMT -lMETIS

CXX:=g++

all: $(EXE)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE)

######################################################################
# Heat ##############################################################
######################################################################
$(EXE): main.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

main.o: $(SRCDIR)/GiDProblemTypes/Coloring/main.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

