include ../BuildFlags.mk

SRCDIR:=../../../source
EXE:=../../../tools/Diff

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR) -I$(SRCDIR)/Tools/Diff
LDFLAGS:=-L.. -lFEMT

CXX:=icpc

all: $(EXE)

clean:
	echo CLEAN
	rm -f *.o
	rm -f $(EXE)

######################################################################
# Diff ###############################################################
######################################################################
$(EXE): main.o
	echo LD $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

main.o: $(SRCDIR)/Tools/Diff/main.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

