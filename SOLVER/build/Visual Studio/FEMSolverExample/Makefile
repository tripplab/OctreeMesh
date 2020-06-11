!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

EXE_PREFIX=..\..\..\tools\FEMSolverExample
EXE=$(EXE_PREFIX).exe


all: $(EXE)

clean:
	del /Q *.obj $(EXE) $(EXE).manifest $(EXE_PREFIX).pdb

######################################################################
# FEMSolverExample ###################################################
######################################################################
$(EXE): FEMSolverExample.obj
	$(LINK) /OUT:$@ $** $(LDFLAGS)

FEMSolverExample.obj: $(SRCDIR)\Examples\FEMSolverExample\FEMSolverExample.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?
