CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic

SRC = src/main.cc src/mesh_reader.cc src/topology.cc src/stl_writer.cc
OBJ = $(SRC:.cc=.o)

all: mesh2stl

mesh2stl: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) mesh2stl
