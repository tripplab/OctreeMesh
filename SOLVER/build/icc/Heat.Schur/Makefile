include ../BuildFlags.mk

SRCDIR:=../../../source
EXE:=../../../gid/problemtypes/Heat.gid/Heat.Schur

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/METIS -I$(SRCDIR)/GiDProblemTypes/Heat
LDFLAGS:=-L.. -lFEMT -lMETIS

CXX:=export OMPI_CC:=icpc; mpic++

all: $(EXE)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE)

######################################################################
# Heat ###############################################################
######################################################################
$(EXE): main.schur.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

main.schur.o: $(SRCDIR)/GiDProblemTypes/Heat/main.schur.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

