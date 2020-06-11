include ../BuildFlags.mk

SRCDIR:=../../../source
EXE:=../../../gid/problemtypes/Heat.gid/Heat

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/METIS -I$(SRCDIR)/GiDProblemTypes/Heat
LDFLAGS:=-L.. -lFEMT -lMETIS

CXX:=icpc

all: $(EXE)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE)

######################################################################
# Heat ###############################################################
######################################################################
$(EXE): main.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

main.o: $(SRCDIR)/GiDProblemTypes/Heat/main.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

