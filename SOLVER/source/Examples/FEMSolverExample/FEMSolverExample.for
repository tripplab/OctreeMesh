C FEMSolverExample.for
C Copyright (C) 2012 Miguel Vargas (miguel.vargas@gmail.com)
C
C This library is free software; you can redistribute it and/or
C modify it under the terms of the GNU Library General Public
C License as published by the Free Software Foundation; either
C version 2 of the License, or (at your option) any later version.
C
C This library is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Library General Public License for more details.
C
C You should have received a copy of the GNU Library General Public
C License along with this library; if not, write to the Free
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

C
C Compilation
C GCC:
C     gfortran -ffixed-line-length-none -cpp FEMSolverExample.for
C Intel:
C     ifort -free -cpp FEMSolverExample.for
C

      PROGRAM FEMSolverExample
      IMPLICIT NONE
C Commands
      INTEGER command_end, command_set_connectivity, command_fill_A, command_set_Ae, command_set_all_Ae, command_set_x, command_set_b, command_set_fixed, command_solver_init, command_solver_run
      PARAMETER (
     . command_end              = 0,
     . command_set_connectivity = 1,
     . command_fill_A           = 2,
     . command_set_Ae           = 3,
     . command_set_all_Ae       = 4,
     . command_set_x            = 5,
     . command_set_b            = 6,
     . command_set_fixed        = 7,
     . command_solver_init      = 8,
     . command_solver_run       = 9)
C Number of elements
      INTEGER*4 E /21/
C Number of nodes
      INTEGER*4 M /17/
C Element type (2=Triangle, 3=Quadrilateral, 4=Tetrahedra, 5=Hexahedra)
      INTEGER*4 T /2/
C Nodes per element
      INTEGER*4 N /3/
C Degrees of freedom
      INTEGER*4 D /1/
C Connectivity
      INTEGER*4 connectivity(3*21) /
     . 3, 8, 5,
     . 3, 5, 1,
     . 5, 8, 12,
     . 7, 2, 4,
     . 7, 4, 11,
     . 4, 2, 1,
     . 17, 16, 13,
     . 17, 13, 15,
     . 13, 16, 14,
     . 12, 15, 10,
     . 10, 15, 13,
     . 12, 10, 5,
     . 14, 11, 9,
     . 9, 11, 4,
     . 14, 9, 13,
     . 10, 13, 6,
     . 6, 13, 9,
     . 10, 6, 5,
     . 6, 9, 4,
     . 5, 6, 4,
     . 1, 5, 4/
C Elemental matrices
      REAL*8 Ke(3*3,21) /
     . 6.4699850468326987e-005, -2.6213892626938322e-005, -3.8485957841388651e-005, -2.6213892626938322e-005, 6.4699852119895029e-005, -3.8485959492956721e-005, -3.8485957841388651e-005, -3.8485959492956721e-005, 7.6971917334345372e-005,
     . 0.0001194163, -6.4851481229273816e-005, -5.4564862420669250e-005, -6.4851481229273816e-005, 6.4518950865214794e-005, 3.3253036405902742e-007, -5.4564862420669250e-005, 3.3253036405902742e-007, 5.4232332056610221e-005,
     . 6.4518952570947138e-005, -6.4851496361322073e-005, 3.3254379037494387e-007, -6.4851496361322073e-005, 0.0001194164, -5.4564874551608474e-005, 3.3254379037494387e-007, -5.4564874551608474e-005, 5.4232330761233532e-005,
     . 6.4699857166259221e-005, -2.6213888552978070e-005, -3.8485968613281151e-005, -2.6213888552978070e-005, 6.4699842120739703e-005, -3.8485953567761646e-005, -3.8485968613281151e-005, -3.8485953567761646e-005, 7.6971922181042797e-005,
     . 0.0001194164, -6.4851506330956674e-005, -5.4564867583856409e-005, -6.4851506330956674e-005, 6.4518961777500545e-005, 3.3254455345613586e-007, -5.4564867583856409e-005, 3.3254455345613586e-007, 5.4232323030400286e-005,
     . 6.4518956933392303e-005, -6.4851494166733798e-005, 3.3253723334149865e-007, -6.4851494166733798e-005, 0.0001194164, -5.4564864260067519e-005, 3.3253723334149865e-007, -5.4564864260067519e-005, 5.4232327026726017e-005,
     . 6.4699848987102081e-005, -2.6213890805761957e-005, -3.8485958181340134e-005, -2.6213890805761957e-005, 6.4699852125378693e-005, -3.8485961319616740e-005, -3.8485958181340134e-005, -3.8485961319616740e-005, 7.6971919500956867e-005,
     . 0.0001194164, -6.4851485993385854e-005, -5.4564874361501674e-005, -6.4851485993385854e-005, 6.4518947014268201e-005, 3.3253897911763883e-007, -5.4564874361501674e-005, 3.3253897911763883e-007, 5.4232335382384034e-005,
     . 6.4518961173541726e-005, -6.4851491940559977e-005, 3.3253076701823620e-007, -6.4851491940559977e-005, 0.0001194163, -5.4564854162971068e-005, 3.3253076701823620e-007, -5.4564854162971068e-005, 5.4232323395952824e-005,
     . 6.4699842120739703e-005, -2.6213888552978070e-005, -3.8485953567761646e-005, -2.6213888552978070e-005, 6.4699857166259221e-005, -3.8485968613281151e-005, -3.8485953567761646e-005, -3.8485968613281151e-005, 7.6971922181042797e-005,
     . 0.000118198, -6.4851513039450068e-005, -5.3346449325098674e-005, -6.4851513039450068e-005, 6.5184045310075873e-005, -3.3253227062582225e-007, -5.3346449325098674e-005, -3.3253227062582225e-007, 5.3678981595724479e-005,
     . 6.5184037957109425e-005, -6.4851504353269536e-005, -3.3253360383990270e-007, -6.4851504353269536e-005, 0.000118198, -5.3346454060646762e-005, -3.3253360383990270e-007, -5.3346454060646762e-005, 5.3678987664486678e-005,
     . 6.4699852125378693e-005, -2.6213890805761964e-005, -3.8485961319616733e-005, -2.6213890805761964e-005, 6.4699848987102095e-005, -3.8485958181340127e-005, -3.8485961319616733e-005, -3.8485958181340127e-005, 7.6971919500956867e-005,
     . 0.0001181979, -6.4851479284895551e-005, -5.3346442685054684e-005, -6.4851479284895551e-005, 6.5184030546843529e-005, -3.3255126194796938e-007, -5.3346442685054684e-005, -3.3255126194796938e-007, 5.3678993947002662e-005,
     . 6.5184022707578200e-005, -6.4851491940559963e-005, -3.3253076701823625e-007, -6.4851491940559963e-005, 0.000118198, -5.3346469426500777e-005, -3.3253076701823625e-007, -5.3346469426500777e-005, 5.3679000193519014e-005,
     . 6.4434406098309372e-005, -2.5551726981809021e-005, -3.8882679116500358e-005, -2.5551726981809021e-005, 6.4434416237837695e-005, -3.8882689256028667e-005, -3.8882679116500358e-005, -3.8882689256028667e-005, 7.7765368372529039e-005,
     . 7.6740651804875536e-005, -3.9583147051491118e-005, -3.7157504753384418e-005, -3.9583147051491118e-005, 6.6011012444926317e-005, -2.6427865393435195e-005, -3.7157504753384418e-005, -2.6427865393435195e-005, 6.3585370146819603e-005,
     . 6.6010999096737498e-005, -3.9583145859760614e-005, -2.6427853236976891e-005, -3.9583145859760614e-005, 7.6740665893494642e-005, -3.7157520033734028e-005, -2.6427853236976891e-005, -3.7157520033734028e-005, 6.3585373270710926e-005,
     . 7.4932933309631666e-005, -3.8450383464282361e-005, -3.6482549845349306e-005, -3.8450383464282361e-005, 6.6423858769906367e-005, -2.7973475305623986e-005, -3.6482549845349306e-005, -2.7973475305623986e-005, 6.4456025150973285e-005,
     . 6.0475784682483771e-005, -6.7551660037995739e-005, 7.0758753555119602e-006, -6.7551660037995739e-005, 0.0001333117, -6.5760034998968276e-005, 7.0758753555119602e-006, -6.5760034998968276e-005, 5.8684159643456314e-005,
     . 7.1839453851444221e-005, -3.5919728741874401e-005, -3.5919725109569793e-005, -3.5919728741874401e-005, 6.6664336196327976e-005, -3.0744607454453576e-005, -3.5919725109569793e-005, -3.0744607454453576e-005, 6.6664332564023382e-005/
C Vector with constrain values
      REAL*8 x(17) /0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 273, 0, 273, 0/
C Vector of independent terms
      REAL*8 b(17) /0.000281733, 0.000281733, 0, 0, 0, 0, 0.000281733, 0.000281733, 0, 0, 0.000281733, 0.000281733, 0, 0, 0.000281733, 0, 0.000281733/
C Vector that indicates where the constrains are
      INTEGER*1 fixed(17) /0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0/
C Result vector
      REAL*8 r(17)
C Variable to indicate command
      INTEGER*1 command
C Data pipe
      INTEGER*4 FEMData /100/
C Result pipe
      INTEGER*4 FEMResult /101/
C Indexes
      INTEGER*4 i, j
C Paths for pipes
#ifdef WIN32
      CHARACTER*128 dataname /'\\.\pipe\FEMData'/
      CHARACTER*128 resultname /'\\.\pipe\FEMResult'/
#else
      CHARACTER*128 dataname /'/tmp/FEMData'/
      CHARACTER*128 resultname /'/tmp/FEMResult'/
#endif

C Open data and result pipes
      OPEN(FEMData, FILE=dataname, ACCESS='STREAM')
      OPEN(FEMResult, FILE=resultname, ACCESS='STREAM')

C Send mesh data
      command = command_set_connectivity
      WRITE(FEMData) command, M, E, T, N, D, (connectivity(i), i = 1, E*N)

C Send elemental matrices
      command = command_set_Ae;
      DO i = 1, E
        WRITE(FEMData) command, i, (Ke(j, i), j = 1, N*N)
      ENDDO

C Send vector with constrain values
      command = command_set_x
      WRITE(FEMData) command, (x(j), j = 1, M)

C Send vector of independent terms
      command = command_set_b
      WRITE(FEMData) command, (b(j), j = 1, M)

C Send vector that indicates where the constrains are
      command = command_set_fixed
      WRITE(FEMData) command, (fixed(j), j = 1, M)

C Initialize solver
      command = command_solver_init
      WRITE(FEMData) command

C Run solver and read result
      command = command_solver_run
      WRITE(FEMData) command
      flush(FEMData)
      READ(FEMResult) (r(j), j = 1, M)

C Display result
      DO i = 1, M
        WRITE(6, *) r(i)
      ENDDO

C Send end session command
      command = command_end
      WRITE(FEMData) command

C Close pipes
      CLOSE(FEMResult)
      CLOSE(FEMData)

      END
