!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

CPPFLAGS=$(CPPFLAGS) /I$(SRCDIR) /I$(SRCDIR)\METIS /I$(SRCDIR)\Tools\EqnSolver $(MPI_INCLUDE)
LDFLAGS=$(LDFLAGS) ..\FEMT.lib ..\METIS.lib $(MPI_LIBS)

EXE_PREFIX=..\..\..\tools\EqnSolver.Schur
EXE=$(EXE_PREFIX).exe


all: $(EXE)

clean:
	del /Q *.obj $(EXE) $(EXE).manifest $(EXE_PREFIX).pdb

######################################################################
# EqnSolver ##########################################################
######################################################################
$(EXE): main.schur.obj
	$(LINK) /OUT:$@ $** $(LDFLAGS)

main.schur.obj: $(SRCDIR)\Tools\EqnSolver\main.schur.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?
