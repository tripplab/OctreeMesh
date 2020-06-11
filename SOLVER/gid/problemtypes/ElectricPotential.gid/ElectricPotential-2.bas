*intformat "%i "
*realformat "%.7g "
*if(strcmp(gendata(Solver_type),"Conjugate_gradient")==0)
*set var solvertype=1
*elseif(strcmp(gendata(Solver_type),"Cholesky_decomposition")==0)
*set var solvertype=2
*elseif(strcmp(gendata(Solver_type),"Cholesky2_decomposition")==0)
*set var solvertype=3
*elseif(strcmp(gendata(Solver_type),"Biconjugate_gradient")==0)
*set var solvertype=4
*elseif(strcmp(gendata(Solver_type),"LU_decomposition")==0)
*set var solvertype=5
*endif
*if(strcmp(gendata(Preconditioner),"None")==0)
*set var preconditioner=0
*elseif(strcmp(gendata(Preconditioner),"Jacobi")==0)
*set var preconditioner=1
*elseif(strcmp(gendata(Preconditioner),"Incomplete_Cholesky")==0)
*set var preconditioner=2
*elseif(strcmp(gendata(Preconditioner),"Incomplete_Cholesky2")==0)
*set var preconditioner=3
*elseif(strcmp(gendata(Preconditioner),"Incomplete_LU")==0)
*set var preconditioner=4
*elseif(strcmp(gendata(Preconditioner),"Sparse_approximate_inverse")==0)
*set var preconditioner=5
*endif
*if((solvertype==4)&&((preconditioner==2)||(preconditioner==3)))
*messagebox Error: Allowed Biconjugate gradient preconditioners are Jacobi, Incomplete LU and Sparse_Approximate_Inverse.
*endif
; Solver file

{Solver}
*solvertype; Type (1=Conjugate_gradient, 2=Cholesky_decomposition, 3=Cholesky2_decomposition, 4=Biconjugate_gradient, 5=LU_decomposition)
*gendata(Threads,int); Threads
; Parameters for iterative solvers
*gendata(Tolerance,real); Tolerance
*gendata(Max_steps,int); Max_steps
*preconditioner; Preconditioner (0=None, 1=Jacobi, 2=Incomplete_Cholesky, 3=Incomplete_Cholesky2, 4=Incomplete_LU, 5=Sparse_Approximate_Inverse)
*gendata(Preconditioner_level,int); Preconditioner_level
*gendata(Preconditioner_threshold,real); Preconditioner_threshold

{Substructuring}
*gendata(Substructuring_threads,int); Substructuring_threads
*gendata(Substructuring_tolerance,real); Substructuring_tolerance
*gendata(Multiple_results,int); Multiple_results (0=No, 1=Yes)
