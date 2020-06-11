!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

CPPFLAGS=$(CPPFLAGS) /I$(SRCDIR) /I$(SRCDIR)\METIS /I$(SRCDIR)\GiDProblemTypes\Partitioning
LDFLAGS=$(LDFLAGS) ..\FEMT.lib ..\METIS.lib

EXE_PREFIX=..\..\..\gid\problemtypes\Partitioning.gid\Partitioning
EXE=$(EXE_PREFIX).exe


all: $(EXE)

clean:
	del /Q *.obj $(EXE) $(EXE).manifest $(EXE_PREFIX).pdb

######################################################################
# Partitioning #######################################################
######################################################################
$(EXE): main.obj
	$(LINK) /OUT:$@ $** $(LDFLAGS)

main.obj: $(SRCDIR)\GiDProblemTypes\Partitioning\main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?
