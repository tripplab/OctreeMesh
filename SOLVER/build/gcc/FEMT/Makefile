include ../BuildFlags.mk

SRCDIR:=../../../source

CPPFLAGS:=$(CPPFLAGS) -I$(SRCDIR)/ -I$(SRCDIR)/METIS
LDFLAGS:=

ifdef USE_MPI
CXX:=mpic++
else
CXX:=g++
endif

all: ../libFEMT.a

clean:
	echo CLEAN
	rm -f *.o ../libFEMT.a

######################################################################
# FEMT ###############################################################
######################################################################
Basic:=Debug.o Float.o Format.o Integer.o Log.o Memory.o Random.o String.o Time.o
ifdef USE_MPI
Communication:=MPI.o
endif
File:=File.o MatFile.o
FiniteElement:=Assembler.o Mesh.o Partition.o Shape.o ShapeFunctions.o
Math:=Formula.o
../libFEMT.a: $(Basic) $(Communication) $(File) $(FiniteElement) $(Math)
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

Random.o: $(SRCDIR)/Basic/Random.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

String.o: $(SRCDIR)/Basic/String.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Time.o: $(SRCDIR)/Basic/Time.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

# Communication ######################################################
MPI.o: $(SRCDIR)/Communication/MPI.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

# File ###############################################################
File.o: $(SRCDIR)/File/File.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

MatFile.o: $(SRCDIR)/File/MatFile.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

# FiniteElement ######################################################
Assembler.o: $(SRCDIR)/FiniteElement/Assembler.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Mesh.o: $(SRCDIR)/FiniteElement/Mesh.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Partition.o: $(SRCDIR)/FiniteElement/Partition.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

Shape.o: $(SRCDIR)/FiniteElement/Shape.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

ShapeFunctions.o: $(SRCDIR)/FiniteElement/ShapeFunctions.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

# Math ###############################################################
Formula.o: $(SRCDIR)/Math/Formula.cpp
	echo CC $<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<
