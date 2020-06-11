@echo off

rem LOG_LEVEL
rem   0  Only fatal error messages
rem   1  Description of the current process
rem   2  Solvers iterations are shown (default)
set LOG_LEVEL=2

rem SOLVER_THREADS
rem Sets the maximum number of threads used by solvers
set SOLVER_THREADS=1

rem SUBSTRUCTURING_THREADS
rem Number of threads used in the master to assemble matrices
set SUBSTRUCTURING_THREADS=1

rem SUBSTRUCTURING_TOLERANCE
rem Tolerance for convergence of Schur matrix
set SUBSTRUCTURING_TOLERANCE=1e-5

set named_pipe_data="FEMData"
set named_pipe_result="FEMResult"

FEMSolver.Schur.exe %named_pipe_data% %named_pipe_result%
