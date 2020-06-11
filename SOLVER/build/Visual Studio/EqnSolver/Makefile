!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

CPPFLAGS=$(CPPFLAGS) /I$(SRCDIR) /I$(SRCDIR)\METIS /I$(SRCDIR)\Tools\EqnSolver
LDFLAGS=$(LDFLAGS) ..\FEMT.lib ..\METIS.lib

EXE_PREFIX=..\..\..\tools\EqnSolver
EXE=$(EXE_PREFIX).exe


all: $(EXE)

clean:
	del /Q *.obj $(EXE) $(EXE).manifest $(EXE_PREFIX).pdb

######################################################################
# EqnSolver ##########################################################
######################################################################
$(EXE): main.obj
	$(LINK) /OUT:$@ $** $(LDFLAGS)

main.obj: $(SRCDIR)\Tools\EqnSolver\main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

