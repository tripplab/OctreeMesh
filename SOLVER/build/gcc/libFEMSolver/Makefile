include ../BuildFlags.mk

SRCDIR:=../../../source

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR)/ -I$(SRCDIR)/METIS
LDFLAGS:=

CXX:=g++

all: ../libFEMSolver.a

clean:
	echo CLEAN
	rm -f *.o ../libFEMSolver.a

######################################################################
# libFEMSolver #######################################################
######################################################################
Basic:=Debug.o Float.o Format.o Integer.o Log.o Memory.o Time.o
FiniteElement:=Assembler.o Mesh.o Shape.o
Library:=libFEMSolver.o
../libFEMSolver.a: $(Basic) $(FiniteElement) $(Library)
	echo AR $@
	$(AR) rcs $@ $+

# Basic ##############################################################
Debug.o: $(SRCDIR)/Basic/Debug.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Float.o: $(SRCDIR)/Basic/Float.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Format.o: $(SRCDIR)/Basic/Format.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Integer.o: $(SRCDIR)/Basic/Integer.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Log.o: $(SRCDIR)/Basic/Log.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Memory.o: $(SRCDIR)/Basic/Memory.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Time.o: $(SRCDIR)/Basic/Time.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

# FiniteElement ######################################################
Assembler.o: $(SRCDIR)/FiniteElement/Assembler.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Mesh.o: $(SRCDIR)/FiniteElement/Mesh.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Shape.o: $(SRCDIR)/FiniteElement/Shape.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

# Library ############################################################
libFEMSolver.o: $(SRCDIR)/libFEMSolver/libFEMSolver.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<
