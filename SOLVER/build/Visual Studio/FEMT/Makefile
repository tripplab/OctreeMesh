!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

!IF $(USE_MPI)==1
CPPFLAGS=$(CPPFLAGS) /I$(SRCDIR) /I$(SRCDIR)\METIS $(MPI_INCLUDE)
!ELSE
CPPFLAGS=$(CPPFLAGS) /I$(SRCDIR) /I$(SRCDIR)\METIS
!ENDIF

all: ..\FEMT.lib

clean:
	del /Q *.obj *.pdb ..\FEMT.lib

######################################################################
# FEMT ###############################################################
######################################################################
Basic=Debug.obj Float.obj Format.obj Integer.obj Log.obj Memory.obj Random.obj String.obj Time.obj
!IF $(USE_MPI)==1
Communication=MPI.obj
!ENDIF
File=File.obj MatFile.obj
FiniteElement=Assembler.obj Mesh.obj Partition.obj Shape.obj ShapeFunctions.obj
Math=Formula.obj
..\FEMT.lib: $(Basic) $(Communication) $(File) $(FiniteElement) $(Math)
	$(AR) /NOLOGO /OUT:$@ $**

# Basic ##############################################################
Debug.obj: $(SRCDIR)\Basic\Debug.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

Float.obj: $(SRCDIR)\Basic\Float.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

Format.obj: $(SRCDIR)\Basic\Format.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

Integer.obj: $(SRCDIR)\Basic\Integer.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

Log.obj: $(SRCDIR)\Basic\Log.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

Memory.obj: $(SRCDIR)\Basic\Memory.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

Random.obj: $(SRCDIR)\Basic\Random.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

String.obj: $(SRCDIR)\Basic\String.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

Time.obj: $(SRCDIR)\Basic\Time.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

# Communication ######################################################
MPI.obj: $(SRCDIR)\Communication\MPI.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

# File ###############################################################
File.obj: $(SRCDIR)\File\File.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

MatFile.obj: $(SRCDIR)\File\MatFile.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

# FiniteElement ######################################################
Assembler.obj: $(SRCDIR)\FiniteElement\Assembler.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

Mesh.obj: $(SRCDIR)\FiniteElement\Mesh.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

Partition.obj: $(SRCDIR)\FiniteElement\Partition.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

Shape.obj: $(SRCDIR)\FiniteElement\Shape.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

ShapeFunctions.obj: $(SRCDIR)\FiniteElement\ShapeFunctions.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?

# Math ###############################################################
Formula.obj: $(SRCDIR)\Math\Formula.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo$@ /c $?
