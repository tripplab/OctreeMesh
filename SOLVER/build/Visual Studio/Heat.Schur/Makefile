!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

CPPFLAGS=$(CPPFLAGS) /I$(SRCDIR) /I$(SRCDIR)\METIS /I$(SRCDIR)\GiDProblemTypes\Heat $(MPI_INCLUDE)
LDFLAGS=$(LDFLAGS) ..\FEMT.lib ..\METIS.lib $(MPI_LIBS)

EXE_PREFIX=..\..\..\gid\problemtypes\Heat.gid\Heat.Schur
EXE=$(EXE_PREFIX).exe


all: $(EXE)

clean:
	del /Q *.obj $(EXE) $(EXE).manifest $(EXE_PREFIX).pdb

######################################################################
# Heat ###############################################################
######################################################################
$(EXE): main.schur.obj
	$(LINK) /OUT:$@ $** $(LDFLAGS)

main.schur.obj: $(SRCDIR)\GiDProblemTypes\Heat\main.schur.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?
