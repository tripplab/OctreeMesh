!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

CPPFLAGS=$(CPPFLAGS) /I$(SRCDIR) /I$(SRCDIR)\METIS /I$(SRCDIR)\GiDProblemTypes\Solid $(MPI_INCLUDE)
LDFLAGS=$(LDFLAGS) ..\FEMT.lib ..\METIS.lib $(MPI_LIBS)

EXE_PREFIX=..\..\..\gid\problemtypes\Solid.gid\Solid.Schur
EXE=$(EXE_PREFIX).exe


all: $(EXE)

clean:
	del /Q *.obj $(EXE) $(EXE).manifest $(EXE_PREFIX).pdb

######################################################################
# Solid ##############################################################
######################################################################
$(EXE): main.schur.obj
	$(LINK) /OUT:$@ $** $(LDFLAGS)

main.schur.obj: $(SRCDIR)\GiDProblemTypes\Solid\main.schur.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?
