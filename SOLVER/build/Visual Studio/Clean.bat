@echo off

echo "Clean Visual Studio"

del /Q *.lib

pushd METIS
del /Q *.obj
popd

pushd FEMT
del /Q *.obj
popd

pushd Solid
del /Q *.obj
popd

pushd Solid.Schur
del /Q *.obj
popd

pushd Heat
del /Q *.obj
popd

pushd Heat.Multigrid
del /Q *.obj
popd

pushd Heat.Schur
del /Q *.obj
popd

pushd Heat.DCG
del /Q *.obj
popd

pushd StructureOptimization
del /Q *.obj
popd

pushd Coloring
del /Q *.o
popd

pushd Diff
del /Q *.obj
popd

pushd FEMSolver
del /Q *.obj
popd

pushd FEMSolver.Schur
del /Q *.obj
popd

pushd Partitioning
del /Q *.obj
popd

pushd EqnSolver
del /Q *.obj
popd

pushd EqnSolverExample
del /Q *.obj
popd

pushd FEMSolverExample
del /Q *.obj
popd
