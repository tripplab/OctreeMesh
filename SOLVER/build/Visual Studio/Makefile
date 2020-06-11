!INCLUDE BuildFlags.mk

!IF $(USE_MPI)==1
ALL=libs tools tools_mpi examples problems problems_mpi
!ELSE
ALL=libs tools examples problems
!ENDIF

CLEAN=clean_libs clean_tools clean_tools_mpi clean_examples clean_problems clean_problems_mpi

all: $(ALL)

libs:
	cd METIS
	nmake /nologo
	cd ..\FEMT
	nmake /nologo
	cd ..
tools:
	cd FEMSolver
	nmake /nologo
	cd ..\EqnSolver
	nmake /nologo
	cd ..\MatSolver
	nmake /nologo
	cd ..
tools_mpi:
	cd EqnSolver.Schur
	nmake /nologo
	cd ..\FEMSolver.Schur
	nmake /nologo
	cd ..
examples:
	cd EqnSolverExample
	nmake /nologo
	cd ..\FEMSolverExample
	nmake /nologo
	cd ..
problems:
	cd Coloring
	nmake /nologo
	cd ..\Diff
	nmake /nologo
	cd ..\ElectricPotential
	nmake /nologo
	cd ..\Heat
	nmake /nologo
	cd ..\Heat.Multigrid
	nmake /nologo
	cd ..\Partitioning
	nmake /nologo
	cd ..\Solid
	nmake /nologo
	cd ..
problems_mpi:
	cd Solid.Schur
	nmake /nologo
	cd ..\Heat.Schur
	nmake /nologo
	cd ..

clean: $(CLEAN)

clean_libs:
	cd METIS
	nmake /nologo clean
	cd ..\FEMT
	nmake /nologo clean
	cd ..
clean_tools:
	cd FEMSolver
	nmake /nologo clean
	cd ..\EqnSolver
	nmake /nologo clean
	cd ..\MatSolver
	nmake /nologo clean
	cd ..
clean_tools_mpi:
	cd EqnSolver.Schur
	nmake /nologo clean
	cd ..\FEMSolver.Schur
	nmake /nologo clean
	cd ..
clean_examples:
	cd EqnSolverExample
	nmake /nologo clean
	cd ..\FEMSolverExample
	nmake /nologo clean
	cd ..
clean_problems:
	cd Coloring
	nmake /nologo clean
	cd ..\Diff
	nmake /nologo clean
	cd ..\ElectricPotential
	nmake /nologo clean
	cd ..\Heat
	nmake /nologo clean
	cd ..\Heat.Multigrid
	nmake /nologo clean
	cd ..\Partitioning
	nmake /nologo clean
	cd ..\Solid
	nmake /nologo clean
	cd ..
clean_problems_mpi:
	cd Solid.Schur
	nmake /nologo clean
	cd ..\Heat.Schur
	nmake /nologo clean
	cd ..
