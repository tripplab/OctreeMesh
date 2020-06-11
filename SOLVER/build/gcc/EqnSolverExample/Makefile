include ../BuildFlags.mk

SRCDIR:=../../../source
EXE_CPP:=../../../tools/EqnSolverExample_cpp
EXE_FOR:=../../../tools/EqnSolverExample_for

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/Examples/EqnSolverExample
LDFLAGS:=

all: $(EXE_CPP) $(EXE_FOR)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE_CPP) $(EXE_FOR)

######################################################################
# EqnSolverExample ###################################################
######################################################################
$(EXE_CPP): EqnSolverExample_cpp.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

EqnSolverExample_cpp.o: $(SRCDIR)/Examples/EqnSolverExample/EqnSolverExample.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

$(EXE_FOR): EqnSolverExample_for.o
	echo LD $@
	$(FC) $(CPPFLAGS) $(FFLAGS) -o $@ $+ $(LDFLAGS)

EqnSolverExample_for.o: $(SRCDIR)/Examples/EqnSolverExample/EqnSolverExample.for
	echo FC $<
	$(FC) $(FFLAGS) -o $@ -c $<
