.SILENT:
.PHONY: all clean

VERSION:=beta36

#USE_DEBUG:=1
#USE_MPI:=1

ifdef USE_DEBUG
CPPFLAGS:=-Wall -Wextra -D$(VERSION) -DVERSION=$(VERSION) -DMEMORY_USAGE
CXXFLAGS:=-openmp -g
FFLAGS:=-g -free -cpp
else
CPPFLAGS:=-Wall -Wextra -D$(VERSION) -DVERSION=$(VERSION) -DNDEBUG
CXXFLAGS:=-openmp -O3
FFLAGS:=-O3 -free -cpp
endif

FC:=ifort
CXX:=icpc
