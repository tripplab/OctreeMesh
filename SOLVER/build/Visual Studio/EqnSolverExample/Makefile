!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

EXE_PREFIX=..\..\..\tools\EqnSolverExample
EXE=$(EXE_PREFIX).exe


all: $(EXE)

clean:
	del /Q *.obj $(EXE) $(EXE).manifest $(EXE_PREFIX).pdb

######################################################################
# EqnSolverExample ###################################################
######################################################################
$(EXE): EqnSolverExample.obj
	$(LINK) /OUT:$@ $** $(LDFLAGS)

EqnSolverExample.obj: $(SRCDIR)\Examples\EqnSolverExample\EqnSolverExample.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?
