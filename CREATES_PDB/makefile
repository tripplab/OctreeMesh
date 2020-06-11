CCOMP=g++
CFLAGS=-std=c++11 -fopenmp -Wall
LFLAGS=-fopenmp -lm
EJECUTABLE=./output/pdb

all: folders $(EJECUTABLE)
folders: ./input ./obj ./obj/boundary ./obj/capside ./obj/file ./obj/geometry ./obj/interpolator ./obj/mesh ./obj/octree ./obj/results ./output
./input:
	mkdir ./input
./obj:
	mkdir ./obj
./obj/boundary:
	mkdir ./obj/boundary
./obj/capside:
	mkdir ./obj/capside
./obj/file:
	mkdir ./obj/file
./obj/geometry:
	mkdir ./obj/geometry
./obj/interpolator:
	mkdir ./obj/interpolator
./obj/mesh:
	mkdir ./obj/mesh
./obj/octree:
	mkdir ./obj/octree
./obj/results:
	mkdir ./obj/results
./output:
	mkdir ./output
$(EJECUTABLE): ./obj/results/results.o ./obj/file/file.o ./obj/geometry/geometry.o ./obj/boundary/boundary.o ./obj/mesh/mesh.o ./obj/capside/capside.o ./obj/octree/octree_driver.o ./obj/interpolator/interpolator.o ./obj/main.o
	$(CCOMP) $(LFLAGS) ./obj/results/results.o ./obj/file/file.o ./obj/geometry/geometry.o ./obj/boundary/boundary.o ./obj/mesh/mesh.o ./obj/capside/capside.o ./obj/octree/octree_driver.o ./obj/interpolator/interpolator.o ./obj/main.o -o $(EJECUTABLE)
./obj/results/results.o: ./src/results/results.cc
	$(CCOMP) $(CFLAGS) -c ./src/results/results.cc -o ./obj/results/results.o
./obj/file/file.o: ./src/file/file.cc
	$(CCOMP) $(CFLAGS) -c ./src/file/file.cc -o ./obj/file/file.o
./obj/geometry/geometry.o: ./src/geometry/geometry.cc
	$(CCOMP) $(CFLAGS) -c ./src/geometry/geometry.cc -o ./obj/geometry/geometry.o
./obj/boundary/boundary.o: ./src/boundary/boundary.cc
	$(CCOMP) $(CFLAGS) -c ./src/boundary/boundary.cc -o ./obj/boundary/boundary.o
./obj/mesh/mesh.o: ./src/mesh/mesh.cc
	$(CCOMP) $(CFLAGS) -c ./src/mesh/mesh.cc -o ./obj/mesh/mesh.o
./obj/capside/capside.o: ./src/capside/capside.cc
	$(CCOMP) $(CFLAGS) -c ./src/capside/capside.cc -o ./obj/capside/capside.o
./obj/octree/octree_driver.o: ./src/octree/octree_driver.cc
	$(CCOMP) $(CFLAGS) -c ./src/octree/octree_driver.cc -o ./obj/octree/octree_driver.o
./obj/interpolator/interpolator.o: ./src/interpolator/interpolator.cc
	$(CCOMP) $(CFLAGS) -c ./src/interpolator/interpolator.cc -o ./obj/interpolator/interpolator.o
./obj/main.o: ./src/main.cc
	$(CCOMP) $(CFLAGS) -c ./src/main.cc -o ./obj/main.o
clean:
	rm -f ./*.txt
	rm -f ./*.log
	rm -f ./*.pdb
	rm -f ./*.msh
	rm -f ./*.res
	rm -f ./output/pdb
	rm -f ./obj/main.o
	rm -f ./obj/results/results.o
	rm -f ./obj/file/file.o
	rm -f ./obj/geometry/geometry.o
	rm -f ./obj/boundary/boundary.o
	rm -f ./obj/mesh/mesh.o
	rm -f ./obj/octree/octree_driver.o
	rm -f ./obj/interpolator/interpolator.o
	rm -f ./obj/capside/capside.o
run:
	./output/pdb ./input/T3_1cwp_full.vdb ./2_T3_1CWP.post.msh ./2_T3_1CWP.post.res ./2_T3_1CWP.pdb 5
