!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

CPPFLAGS=$(CPPFLAGS) /I$(SRCDIR) /I$(SRCDIR)\Tools\Diff
LDFLAGS=$(LDFLAGS) ..\FEMT.lib

EXE_PREFIX=..\..\..\tools\Diff
EXE=$(EXE_PREFIX).exe


all: $(EXE)

clean:
	del /Q *.obj $(EXE) $(EXE).manifest $(EXE_PREFIX).pdb

######################################################################
# Diff ###############################################################
######################################################################
$(EXE): main.obj
	$(LINK) /OUT:$@ $** $(LDFLAGS)

main.obj: $(SRCDIR)\Tools\Diff\main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?
