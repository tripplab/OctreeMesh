#!/bin/sh

# LOG_LEVEL
#   0  Only fatal error messages
#   1  Description of the current process
#   2  Solvers iterations are shown (default)
export LOG_LEVEL=2

# SOLVER_TYPE
# Chooses solver for all session
#   1  Conjugate gradient (default)
#   2  Cholesky decomposition
#   3  Cholesky2 decomposition
#   4  Biconjugate gradient
#   5  LU decomposition
export SOLVER_TYPE=1

# SOLVER_THREADS
# Sets the maximum number of threads used by solvers
export SOLVER_THREADS=1

# SOLVER_TOLERANCE
# For iterative solvers, this defines the convergence criteria
export SOLVER_TOLERANCE=1e-5

# SOLVER_MAX_STEPS
# For iterative solvers, it is the maximum number of iterations allowed
export SOLVER_MAX_STEPS=10000

# PRECONDITIONER_TYPE
# For iterative solvers.
#   0=None
#   1=Jacobi (default)
#   2=Incomplete_Cholesky (Only for Conjugate gradient)
#   3=Incomplete_Cholesky2 (Only for Conjugate gradient)
#   4=Incomplete_LU (Only for Biconjugate gradient)
#   5=Sparse_Approximate_Inverse
export PRECONDITIONER_TYPE=1

# PRECONDITIONER_LEVEL
# Set the preconditioner level for Incomplete_Cholesky, Incomplete_Cholesky2, Incomplete_LU and Sparse_Approximate_Inverse
export PRECONDITIONER_LEVEL=1

# PRECONDITIONER_THRESHOLD
# Set the preconditioner threshold for Sparse_Approximate_Inverse
export PRECONDITIONER_THRESHOLD=0.0

named_pipe_data="EqnData"
named_pipe_result="EqnResult"

./EqnSolver $named_pipe_data $named_pipe_result
