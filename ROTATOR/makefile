CCOMP=g++
CFLAGS=-std=c++11 -Wall
LFLAGS=-lm
EJECUTABLE=./rotator

all: $(EJECUTABLE)
$(EJECUTABLE): ./main.o
	$(CCOMP) $(LFLAGS) ./main.o -o $(EJECUTABLE)
./obj/main.o: ./main.cc
	$(CCOMP) $(CFLAGS) -c ./main.cc -o ./main.o
clean:
	rm -f ./rotator
	rm -f ./main.o
run:
	./rotator ./2_T3_1CWP.geometry.dat ./2_T3_1CWP.geometry.dat
	#./rotator ../MESHER/2_T3_1CWP.geometry.dat ../MESHER/2_T3_1CWP.geometry.dat

