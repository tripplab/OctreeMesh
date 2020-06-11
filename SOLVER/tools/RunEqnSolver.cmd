@echo off

rem LOG_LEVEL
rem   0  Only fatal error messages
rem   1  Description of the current process
rem   2  Solvers iterations are shown (default)
set LOG_LEVEL=2

rem SOLVER_TYPE
rem Chooses solver for all session
rem   1  Conjugate gradient (default)
rem   2  Cholesky decomposition
rem   3  Cholesky2 decomposition
rem   4  Biconjugate gradient
rem   5  LU decomposition
set SOLVER_TYPE=1

rem SOLVER_THREADS
rem Sets the maximum number of threads used by solvers
set SOLVER_THREADS=1

rem SOLVER_TOLERANCE
rem For iterative solvers, this defines the convergence criteria
set SOLVER_TOLERANCE=1e-5

rem SOLVER_MAX_STEPS
rem For iterative solvers, it is the maximum number of iterations allowed
set SOLVER_MAX_STEPS=10000

rem PRECONDITIONER_TYPE
rem For iterative solvers.
rem   0=None
rem   1=Jacobi (default)
rem   2=Incomplete_Cholesky (Only for Conjugate gradient)
rem   3=Incomplete_Cholesky2 (Only for Conjugate gradient)
rem   4=Incomplete_LU (Only for Biconjugate gradient)
rem   5=Sparse_Approximate_Inverse
set PRECONDITIONER_TYPE=1

rem PRECONDITIONER_LEVEL
rem Set the preconditioner level for Incomplete_Cholesky, Incomplete_Cholesky2, Incomplete_LU and Sparse_Approximate_Inverse
set PRECONDITIONER_LEVEL=1

rem PRECONDITIONER_THRESHOLD
rem Set the preconditioner threshold for Sparse_Approximate_Inverse
set PRECONDITIONER_THRESHOLD=0.0

set named_pipe_data="EqnData"
set named_pipe_result="EqnResult"

EqnSolver %named_pipe_data% %named_pipe_result%
