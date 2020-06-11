include BuildFlags.mk

ifdef USE_MPI
ALL:=libs tools tools_mpi examples problems problems_mpi
else
ALL:=libs tools examples problems
endif
CLEAN:=clean_libs clean_tools clean_tools_mpi clean_examples clean_problems clean_problems_mpi

PARALLEL=4

all: $(ALL)

libs:
	make -j $(PARALLEL) -C METIS
	make -j $(PARALLEL) -C FEMT
	make -j $(PARALLEL) -C libFEMSolver
tools:
	make -j $(PARALLEL) -C FEMSolver
	make -j $(PARALLEL) -C EqnSolver
	make -j $(PARALLEL) -C MatSolver
tools_mpi:
	make -j $(PARALLEL) -C EqnSolver.Schur
	make -j $(PARALLEL) -C FEMSolver.Schur
examples:
	make -j $(PARALLEL) -C EqnSolverExample
	make -j $(PARALLEL) -C FEMSolverExample
problems:
	make -j $(PARALLEL) -C Coloring
	make -j $(PARALLEL) -C Diff
	make -j $(PARALLEL) -C ElectricPotential
	make -j $(PARALLEL) -C Heat
	make -j $(PARALLEL) -C Partitioning
	make -j $(PARALLEL) -C Solid
problems_mpi:
	make -j $(PARALLEL) -C Heat.Schur
	make -j $(PARALLEL) -C Solid.Schur

clean: $(CLEAN)

clean_libs:
	make -C METIS clean
	make -C FEMT clean
	make -C libFEMSolver clean
clean_tools:
	make -C FEMSolver clean
	make -C EqnSolver clean
	make -C MatSolver clean
clean_tools_mpi:
	make -C EqnSolver.Schur clean
	make -C FEMSolver.Schur clean
clean_examples:
	make -C EqnSolverExample clean
	make -C FEMSolverExample clean
clean_problems:
	make -C Coloring clean
	make -C Diff clean
	make -C ElectricPotential clean
	make -C Heat clean
	make -C Partitioning clean
	make -C Solid clean
clean_problems_mpi:
	make -C Heat.Schur clean
	make -C Solid.Schur clean
