#!/bin/sh

# LOG_LEVEL
#   0  Only fatal error messages
#   1  Description of the current process
#   2  Solvers iterations are shown (default)
export LOG_LEVEL=2

# SOLVER_THREADS
# Sets the maximum number of threads used by solvers
export SOLVER_THREADS=1

# SUBSTRUCTURING_THREADS
# Number of threads used in the master to assemble matrices
export SUBSTRUCTURING_THREADS=1

# SUBSTRUCTURING_TOLERANCE
# Tolerance for convergence of Schur matrix
export SUBSTRUCTURING_TOLERANCE=1e-5

named_pipe_data="FEMData"
named_pipe_result="FEMResult"

./FEMSolver.Schur $named_pipe_data $named_pipe_result
