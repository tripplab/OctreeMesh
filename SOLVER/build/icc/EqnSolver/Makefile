include ../BuildFlags.mk

SRCDIR:=../../../source
EXE:=../../../tools/EqnSolver

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/METIS -I$(SRCDIR)/Tools/EqnSolver
LDFLAGS:=-L.. -lFEMT -lMETIS

CXX:=icpc

all: $(EXE)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE)

######################################################################
# Solid ##############################################################
######################################################################
$(EXE): main.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

main.o: $(SRCDIR)/Tools/EqnSolver/main.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

