@echo off

echo "Clean GCC"

del /Q *.a

pushd METIS
del /Q *.o
popd

pushd FEMT
del /Q *.o
popd

pushd FEMSolver
del /Q *.o
popd

pushd EqnSolver
del /Q *.o
popd

pushd MatSolver
del /Q *.o
popd

pushd FEMSolver.Schur
del /Q *.o
popd

pushd EqnSolver.Schur
del /Q *.o
popd

pushd EqnSolverExample
del /Q *.o
popd

pushd FEMSolverExample
del /Q *.o
popd

pushd Coloring
del /Q *.o
popd

pushd Diff
del /Q *.o
popd

pushd ElectricPotential
del /Q *.o
popd

pushd Heat
del /Q *.o
popd

pushd Heat.Multigrid
del /Q *.o
popd

pushd Partitioning
del /Q *.o
popd

pushd Solid
del /Q *.o
popd

pushd Heat.Schur
del /Q *.o
popd

pushd Solid.Schur
del /Q *.o
popd
