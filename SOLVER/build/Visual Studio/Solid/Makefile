!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

CPPFLAGS=$(CPPFLAGS) /I$(SRCDIR) /I$(SRCDIR)\METIS /I$(SRCDIR)\GiDProblemTypes\Solid
LDFLAGS=$(LDFLAGS) ..\FEMT.lib ..\METIS.lib

EXE_PREFIX=..\..\..\gid\problemtypes\Solid.gid\Solid
EXE=$(EXE_PREFIX).exe


all: $(EXE)

clean:
	del /Q *.obj $(EXE) $(EXE).manifest $(EXE_PREFIX).pdb

######################################################################
# Solid ##############################################################
######################################################################
$(EXE): main.obj
	$(LINK) /OUT:$@ $** $(LDFLAGS)

main.obj: $(SRCDIR)\GiDProblemTypes\Solid\main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?
