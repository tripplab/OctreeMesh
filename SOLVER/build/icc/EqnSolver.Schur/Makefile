include ../BuildFlags.mk

SRCDIR:=../../../source
EXE:=../../../tools/EqnSolver.Schur

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/METIS -I$(SRCDIR)/Tools/EqnSolver
LDFLAGS:=-L.. -lFEMT -lMETIS

CXX:=export OMPI_CC:=icpc; mpic++

all: $(EXE)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE)

######################################################################
# EqnSolver ##########################################################
######################################################################
$(EXE): main.schur.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

main.schur.o: $(SRCDIR)/Tools/EqnSolver/main.schur.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

