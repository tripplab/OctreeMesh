include ../BuildFlags.mk

SRCDIR:=../../../source
EXE_CPP:=../../../tools/FEMSolverExample_cpp
EXE_FOR:=../../../tools/FEMSolverExample_for

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/Examples/FEMSolverExample
LDFLAGS:=

all: $(EXE_CPP) $(EXE_FOR)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE_CPP) $(EXE_FOR)

######################################################################
# FEMSolverExample ###################################################
######################################################################
$(EXE_CPP): FEMSolverExample_cpp.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

FEMSolverExample_cpp.o: $(SRCDIR)/Examples/FEMSolverExample/FEMSolverExample.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

$(EXE_FOR): FEMSolverExample_for.o
	echo LD $@
	$(FC) $(CPPFLAGS) $(FFLAGS) -o $@ $+ $(LDFLAGS)

FEMSolverExample_for.o: $(SRCDIR)/Examples/FEMSolverExample/FEMSolverExample.for
	echo FC $<
	$(FC) $(FFLAGS) -o $@ -c $<
