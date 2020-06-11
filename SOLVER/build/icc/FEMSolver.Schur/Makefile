include ../BuildFlags.mk

SRCDIR:=../../../source
EXE:=../../../tools/FEMSolver.Schur

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/METIS -I$(SRCDIR)/Tools/FEMSolver
LDFLAGS:=-L.. -lFEMT -lMETIS

CXX:=export OMPI_CC:=icpc; mpic++

all: $(EXE)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE)

######################################################################
# FEMSolver ##########################################################
######################################################################
$(EXE): main.schur.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

main.schur.o: $(SRCDIR)/Tools/FEMSolver/main.schur.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

