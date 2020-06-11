! test.for
! Copyright (C) 2014 Miguel Vargas (miguel.vargas@gmail.com)
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Library General Public License for more details.
!
! You should have received a copy of the GNU Library General Public
! License along with this library; if not, write to the Free
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

!
! Compilation
! GCC:
!     gfortran -ffree-line-length-none test.f90 -L. -lFEMSolver -lMETIS -lstdc++ -fopenmp
! Intel:
!     ifort test.f90 -L. -lFEMSolver -lMETIS -lstdc++ -openmp
!

PROGRAM test

	IMPLICIT NONE

! Solver ID
	INTEGER*4 ID
	
! Solver parameters
	INTEGER*4 solver_type
	INTEGER*4 solver_threads
	REAL*8    solver_tolerance
	INTEGER*4 solver_max_steps
	INTEGER*4 preconditioner_type
	INTEGER*4 preconditioner_level
	REAL*8    preconditioner_threshold
	INTEGER*4 message_level

! Mesh data
	INTEGER*4 number_of_nodes
	INTEGER*4 number_of_elements
	INTEGER*4 element_type
	INTEGER*4 nodes_per_element
	INTEGER*4 degrees_of_freedom
	INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: connectivity

! Elemental matrix
   REAL*8, DIMENSION(3*3) :: Ke

! Vectors
	INTEGER*1, DIMENSION(:), ALLOCATABLE :: fixed
	REAL*8, DIMENSION(:), ALLOCATABLE :: fixed_values
	REAL*8, DIMENSION(:), ALLOCATABLE :: x
	REAL*8, DIMENSION(:), ALLOCATABLE :: b ! Vector of independent terms

! Valid solution
	LOGICAL*1 valid_solution

! Selecto solver to use
	ID = 1
	
! Initialize solver
	solver_type = 1 ! 1=Conjugate_gradient, 2=Cholesky_decomposition, 3=Cholesky2_decomposition, 4=Biconjugate_gradient, 5=LU_decomposition
	solver_threads = 1
	solver_tolerance = 1.0e-5
	solver_max_steps = 10000
	preconditioner_type = 1 ! 0=None, 1=Jacobi, 2=Incomplete_Cholesky, 3=Incomplete_Cholesky2, 4=Incomplete_LU, 5=Sparse_Approximate_Inverse
	preconditioner_level = 1
	preconditioner_threshold = 0.0
	message_level = 2 ! 0=None, 1=Information_messages, 2=All_messages
	CALL FEMSolverInit(ID, solver_type, solver_threads, solver_tolerance, solver_max_steps, preconditioner_type, preconditioner_level, preconditioner_threshold, message_level)

! Set connectivity
	number_of_nodes = 17
	number_of_elements = 21
	element_type = 2 ! 2=Triangle, 3=Quadrilateral, 4=Tetrahedra, 5=Hexahedra
	nodes_per_element = 3
	degrees_of_freedom = 1
	ALLOCATE(connectivity(number_of_elements, nodes_per_element))
	connectivity(1,:) = (/3, 8, 5/)
	connectivity(2,:) = (/3, 5, 1/)
	connectivity(3,:) = (/5, 8, 12/)
	connectivity(4,:) = (/7, 2, 4/)
	connectivity(5,:) = (/7, 4, 11/)
	connectivity(6,:) = (/4, 2, 1/)
	connectivity(7,:) = (/17, 16, 13/)
	connectivity(8,:) = (/17, 13, 15/)
	connectivity(9,:) = (/13, 16, 14/)
	connectivity(10,:) = (/12, 15, 10/)
	connectivity(11,:) = (/10, 15, 13/)
	connectivity(12,:) = (/12, 10, 5/)
	connectivity(13,:) = (/14, 11, 9/)
	connectivity(14,:) = (/9, 11, 4/)
	connectivity(15,:) = (/14, 9, 13/)
	connectivity(16,:) = (/10, 13, 6/)
	connectivity(17,:) = (/6, 13, 9/)
	connectivity(18,:) = (/10, 6, 5/)
	connectivity(19,:) = (/6, 9, 4/)
	connectivity(20,:) = (/5, 6, 4/)
	connectivity(21,:) = (/1, 5, 4/)
	CALL FEMSolverSetConnectivity(ID, number_of_nodes, number_of_elements, element_type, nodes_per_element, degrees_of_freedom, connectivity)
	DEALLOCATE(connectivity)

! Initialize A
	CALL FEMSolverFillA(ID, 0.0)
	
! Set elemental matrices
	Ke = (/6.4699850468326987e-005, -2.6213892626938322e-005, -3.8485957841388651e-005, -2.6213892626938322e-005, 6.4699852119895029e-005, -3.8485959492956721e-005, -3.8485957841388651e-005, -3.8485959492956721e-005, 7.6971917334345372e-005/)
	CALL FEMSolverSetAe(ID, 1, Ke)
	Ke = (/0.0001194163, -6.4851481229273816e-005, -5.4564862420669250e-005, -6.4851481229273816e-005, 6.4518950865214794e-005, 3.3253036405902742e-007, -5.4564862420669250e-005, 3.3253036405902742e-007, 5.4232332056610221e-005/)
	CALL FEMSolverSetAe(ID, 2, Ke)
	Ke = (/6.4518952570947138e-005, -6.4851496361322073e-005, 3.3254379037494387e-007, -6.4851496361322073e-005, 0.0001194164, -5.4564874551608474e-005, 3.3254379037494387e-007, -5.4564874551608474e-005, 5.4232330761233532e-005/)
	CALL FEMSolverSetAe(ID, 3, Ke)
	Ke = (/6.4699857166259221e-005, -2.6213888552978070e-005, -3.8485968613281151e-005, -2.6213888552978070e-005, 6.4699842120739703e-005, -3.8485953567761646e-005, -3.8485968613281151e-005, -3.8485953567761646e-005, 7.6971922181042797e-005/)
	CALL FEMSolverSetAe(ID, 4, Ke)
	Ke = (/0.0001194164, -6.4851506330956674e-005, -5.4564867583856409e-005, -6.4851506330956674e-005, 6.4518961777500545e-005, 3.3254455345613586e-007, -5.4564867583856409e-005, 3.3254455345613586e-007, 5.4232323030400286e-005/)
	CALL FEMSolverSetAe(ID, 5, Ke)
	Ke = (/6.4518956933392303e-005, -6.4851494166733798e-005, 3.3253723334149865e-007, -6.4851494166733798e-005, 0.0001194164, -5.4564864260067519e-005, 3.3253723334149865e-007, -5.4564864260067519e-005, 5.4232327026726017e-005/)
	CALL FEMSolverSetAe(ID, 6, Ke)
	Ke = (/6.4699848987102081e-005, -2.6213890805761957e-005, -3.8485958181340134e-005, -2.6213890805761957e-005, 6.4699852125378693e-005, -3.8485961319616740e-005, -3.8485958181340134e-005, -3.8485961319616740e-005, 7.6971919500956867e-005/)
	CALL FEMSolverSetAe(ID, 7, Ke)
	Ke = (/0.0001194164, -6.4851485993385854e-005, -5.4564874361501674e-005, -6.4851485993385854e-005, 6.4518947014268201e-005, 3.3253897911763883e-007, -5.4564874361501674e-005, 3.3253897911763883e-007, 5.4232335382384034e-005/)
	CALL FEMSolverSetAe(ID, 8, Ke)
	Ke = (/6.4518961173541726e-005, -6.4851491940559977e-005, 3.3253076701823620e-007, -6.4851491940559977e-005, 0.0001194163, -5.4564854162971068e-005, 3.3253076701823620e-007, -5.4564854162971068e-005, 5.4232323395952824e-005/)
	CALL FEMSolverSetAe(ID, 9, Ke)
	Ke = (/6.4699842120739703e-005, -2.6213888552978070e-005, -3.8485953567761646e-005, -2.6213888552978070e-005, 6.4699857166259221e-005, -3.8485968613281151e-005, -3.8485953567761646e-005, -3.8485968613281151e-005, 7.6971922181042797e-005/)
	CALL FEMSolverSetAe(ID, 10, Ke)
	Ke = (/0.000118198, -6.4851513039450068e-005, -5.3346449325098674e-005, -6.4851513039450068e-005, 6.5184045310075873e-005, -3.3253227062582225e-007, -5.3346449325098674e-005, -3.3253227062582225e-007, 5.3678981595724479e-005/)
	CALL FEMSolverSetAe(ID, 11, Ke)
	Ke = (/6.5184037957109425e-005, -6.4851504353269536e-005, -3.3253360383990270e-007, -6.4851504353269536e-005, 0.000118198, -5.3346454060646762e-005, -3.3253360383990270e-007, -5.3346454060646762e-005, 5.3678987664486678e-005/)
	CALL FEMSolverSetAe(ID, 12, Ke)
	Ke = (/6.4699852125378693e-005, -2.6213890805761964e-005, -3.8485961319616733e-005, -2.6213890805761964e-005, 6.4699848987102095e-005, -3.8485958181340127e-005, -3.8485961319616733e-005, -3.8485958181340127e-005, 7.6971919500956867e-005/)
	CALL FEMSolverSetAe(ID, 13, Ke)
	Ke = (/0.0001181979, -6.4851479284895551e-005, -5.3346442685054684e-005, -6.4851479284895551e-005, 6.5184030546843529e-005, -3.3255126194796938e-007, -5.3346442685054684e-005, -3.3255126194796938e-007, 5.3678993947002662e-005/)
	CALL FEMSolverSetAe(ID, 14, Ke)
	Ke = (/6.5184022707578200e-005, -6.4851491940559963e-005, -3.3253076701823625e-007, -6.4851491940559963e-005, 0.000118198, -5.3346469426500777e-005, -3.3253076701823625e-007, -5.3346469426500777e-005, 5.3679000193519014e-005/)
	CALL FEMSolverSetAe(ID, 15, Ke)
	Ke = (/6.4434406098309372e-005, -2.5551726981809021e-005, -3.8882679116500358e-005, -2.5551726981809021e-005, 6.4434416237837695e-005, -3.8882689256028667e-005, -3.8882679116500358e-005, -3.8882689256028667e-005, 7.7765368372529039e-005/)
	CALL FEMSolverSetAe(ID, 16, Ke)
	Ke = (/7.6740651804875536e-005, -3.9583147051491118e-005, -3.7157504753384418e-005, -3.9583147051491118e-005, 6.6011012444926317e-005, -2.6427865393435195e-005, -3.7157504753384418e-005, -2.6427865393435195e-005, 6.3585370146819603e-005/)
	CALL FEMSolverSetAe(ID, 17, Ke)
	Ke = (/6.6010999096737498e-005, -3.9583145859760614e-005, -2.6427853236976891e-005, -3.9583145859760614e-005, 7.6740665893494642e-005, -3.7157520033734028e-005, -2.6427853236976891e-005, -3.7157520033734028e-005, 6.3585373270710926e-005/)
	CALL FEMSolverSetAe(ID, 18, Ke)
	Ke = (/7.4932933309631666e-005, -3.8450383464282361e-005, -3.6482549845349306e-005, -3.8450383464282361e-005, 6.6423858769906367e-005, -2.7973475305623986e-005, -3.6482549845349306e-005, -2.7973475305623986e-005, 6.4456025150973285e-005/)
	CALL FEMSolverSetAe(ID, 19, Ke)
	Ke = (/6.0475784682483771e-005, -6.7551660037995739e-005, 7.0758753555119602e-006, -6.7551660037995739e-005, 0.0001333117, -6.5760034998968276e-005, 7.0758753555119602e-006, -6.5760034998968276e-005, 5.8684159643456314e-005/)
	CALL FEMSolverSetAe(ID, 20, Ke)
	Ke = (/7.1839453851444221e-005, -3.5919728741874401e-005, -3.5919725109569793e-005, -3.5919728741874401e-005, 6.6664336196327976e-005, -3.0744607454453576e-005, -3.5919725109569793e-005, -3.0744607454453576e-005, 6.6664332564023382e-005/)
	CALL FEMSolverSetAe(ID, 21, Ke)

! Allocate vectors
	ALLOCATE(fixed(number_of_nodes))
	ALLOCATE(fixed_values(number_of_nodes))
	ALLOCATE(b(number_of_nodes))
	ALLOCATE(x(number_of_nodes))

! Set fixed conditions
	fixed = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0/)
	CALL FEMSolverSetFixed(ID, fixed)

! Fixed values
	fixed_values = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 273, 0, 273, 0/)

! Right side vector
	b = (/0.000281733, 0.000281733, 0.0, 0.0, 0.0, 0.0, 0.000281733, 0.000281733, 0.0, 0.0, 0.000281733, 0.000281733, 0.0, 0.0, 0.000281733, 0.0, 0.000281733/)

! Compensate fixed values and set B
	CALL FEMSolverCompensateFixed(ID, fixed_values, b)

! Run solver
	CALL FEMSolverRun(ID, x, valid_solution)

! Set right side vector
	b = (/0.000563466, 0.000563466, 0.0, 0.0, 0.0, 0.0, 0.000563466, 0.000563466, 0.0, 0.0, 0.000563466, 0.000563466, 0.0, 0.0, 0.000563466, 0.0, 0.000563466/)
	CALL FEMSolverSetB(ID, b)

! Run solver again
	CALL FEMSolverRun(ID, x, valid_solution)
	
! Deallocate vectors
	DEALLOCATE(x)
	DEALLOCATE(b)
	DEALLOCATE(fixed_values)
	DEALLOCATE(fixed)

! End solver
	CALL FEMSolverEnd(ID)

END
