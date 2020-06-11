!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

CPPFLAGS=$(CPPFLAGS) /I$(SRCDIR) /I$(SRCDIR)\METIS /I$(SRCDIR)\GiDProblemTypes\Coloring
LDFLAGS=$(LDFLAGS) ..\FEMT.lib ..\METIS.lib

EXE_PREFIX=..\..\..\gid\problemtypes\Coloring.gid\Coloring
EXE=$(EXE_PREFIX).exe


all: $(EXE)

clean:
	del /Q *.obj $(EXE) $(EXE).manifest $(EXE_PREFIX).pdb

######################################################################
# Coloring ###########################################################
######################################################################
$(EXE): main.obj
	$(LINK) /OUT:$@ $** $(LDFLAGS)

main.obj: $(SRCDIR)\GiDProblemTypes\Coloring\main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?
