C EqnSolverExample.for
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
C     gfortran -ffixed-line-length-none -cpp EqnSolverExample.for
C Intel:
C     ifort -free -cpp EqnSolverExample.for
C

      PROGRAM EqnSolverExample
      IMPLICIT NONE
C Commands
      INTEGER command_end, command_set_size, command_set_row, command_set_x, command_set_b, command_set_fixed, command_solver_init, command_solver_run
      PARAMETER (
     . command_end         = 0,
     . command_set_size    = 1,
     . command_set_row     = 2,
     . command_set_x       = 3,
     . command_set_b       = 4,
     . command_set_fixed   = 5,
     . command_solver_init = 6,
     . command_solver_run  = 7)
C Number of equations
      INTEGER*4 M /22/
C Matrix values
      REAL*8 values(120) /
     . 0.000123262649322419, -3.60513546547159e-005, -3.39967111532538e-005, -5.32145835144495e-005,
     . -3.60513546547159e-005, 0.000170510996033461, -4.34952789122033e-005, -6.08778004323386e-005, -3.0086562034203e-005,
     . -3.39967111532538e-005, 0.000170318953739018, -5.57481717486078e-005, -4.56652509159795e-005, -3.49088199211772e-005,
     . -5.32145835144495e-005, -4.34952789122033e-005, -5.57481717486078e-005, 0.000341662307627607, -7.11430412224589e-005, -6.86630215458964e-005, -4.93982106839912e-005,
     . -6.08778004323386e-005, -7.11430412224589e-005, 0.000342667616164504, -5.86038603040437e-005, -4.42802365962771e-005, -4.31018496289288e-005, -6.46608279804565e-005,
     . -4.56652509159795e-005, -6.86630215458964e-005, 0.000341547174453471, -5.62295144334364e-005, -4.42702807061953e-005, -5.89712307716659e-005, -6.77478760802972e-005,
     . -3.0086562034203e-005, -5.86038603040437e-005, 0.000122792254967864, -3.41018326296169e-005,
     . -3.49088199211772e-005, -5.62295144334364e-005, 0.000122915210160465, -3.17768758058511e-005,
     . -4.93982106839912e-005, -4.42802365962771e-005, -4.42702807061953e-005, 0.000329610569383882, -4.72254168572741e-005, -5.27761219876592e-005, -5.23122253623367e-005, -3.9348077190148e-005,
     . -5.89712307716659e-005, -3.17768758058511e-005, 0.000170291900663079, -4.46349962484123e-005, -3.49087978371502e-005,
     . -4.31018496289288e-005, -3.41018326296169e-005, 0.000170883134605931, -6.64646115422989e-005, -2.72148408050865e-005,
     . -6.46608279804565e-005, -4.72254168572741e-005, -6.64646115422989e-005, 0.000342335610854092, -5.09176371517581e-005, -6.2797286976477e-005, -5.02698303458277e-005,
     . -6.77478760802972e-005, -5.27761219876592e-005, -4.46349962484123e-005, 0.000341041263674863, -5.81729055099591e-005, -5.62295230426859e-005, -6.14798408058491e-005,
     . -5.23122253623367e-005, -5.09176371517581e-005, 0.000342524050812863, -7.14963888156276e-005, -5.57119746311947e-005, -5.11077281633579e-005, -6.09780966885886e-005,
     . -3.9348077190148e-005, -5.81729055099591e-005, -7.14963888156276e-005, 0.000344682167058638, -4.53552672668097e-005, -6.39237466853512e-005, -6.63857815907427e-005,
     . -3.49087978371502e-005, -5.62295230426859e-005, 0.000122915205551077, -3.17768846712407e-005,
     . -2.72148408050865e-005, -6.2797286976477e-005, 0.000122870955324627, -3.28588275430632e-005,
     . -6.14798408058491e-005, -4.53552672668097e-005, -3.17768846712407e-005, 0.000170574963019778, -3.19629702758784e-005,
     . -5.02698303458277e-005, -5.57119746311947e-005, -3.28588275430632e-005, 0.000170881035319409, -3.20404027993234e-005,
     . -5.11077281633579e-005, -6.39237466853512e-005, 0.000169496564745875, -2.48620285243623e-005, -2.96030613728036e-005,
     . -6.63857815907427e-005, -3.19629702758784e-005, -2.48620285243623e-005, 0.000123210780390983,
     . -6.09780966885886e-005, -3.20404027993234e-005, -2.96030613728036e-005, 0.000122621560860716/
C Matrix indexes
      INTEGER*4 indexes(120) /
     . 1, 2, 3, 4,
     . 1, 2, 4, 5, 7,
     . 1, 3, 4, 6, 8,
     . 1, 2, 3, 4, 5, 6, 9,
     . 2, 4, 5, 7, 9, 11, 12,
     . 3, 4, 6, 8, 9, 10, 13,
     . 2, 5, 7, 11,
     . 3, 6, 8, 10,
     . 4, 5, 6, 9, 12, 13, 14, 15,
     . 6, 8, 10, 13, 16,
     . 5, 7, 11, 12, 17,
     . 5, 9, 11, 12, 14, 17, 19,
     . 6, 9, 10, 13, 15, 16, 18,
     . 9, 12, 14, 15, 19, 20, 22,
     . 9, 13, 14, 15, 18, 20, 21,
     . 10, 13, 16, 18,
     . 11, 12, 17, 19,
     . 13, 15, 16, 18, 21,
     . 12, 14, 17, 19, 22,
     . 14, 15, 20, 21, 22,
     . 15, 18, 20, 21,
     . 14, 19, 20, 22/
C Row sizes
      INTEGER*4 count(22) /4, 5, 5, 7, 7, 7, 4, 4, 8, 5, 5, 7, 7, 7, 7, 4, 4, 5, 5, 5, 4, 4/
C Vector with constrain values
      REAL*8 x(22) /273, 273, 0, 0, 0, 0, 273, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 293, 0/
C Vector of independent terms
      REAL*8 b(22) /0, 0, 0, 0, 0, 0, 0, -0.000433883773336144, 0, -0.000867767546672287, 0, 0, 0, 0, 0, -0.000433883773336144, -0.000216941883549719, 0, -0.000433883718353043, 0, 0, -0.000216941834803323/
C Vector that indicates where the constrains are
      INTEGER*1 fixed(22) /1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0/
C Result vector
      REAL*8 r(22)
C Variable to indicate command
      INTEGER*1 command
C Data pipe
      INTEGER*4 EqnData /100/
C Result pipe
      INTEGER*4 EqnResult /101/
C Indexes
      INTEGER*4 i, j, j1, jn
C Paths for pipes
#ifdef WIN32
      CHARACTER*128 dataname /'\\.\pipe\EqnData'/
      CHARACTER*128 resultname /'\\.\pipe\EqnResult'/
#else
      CHARACTER*128 dataname /'/tmp/EqnData'/
      CHARACTER*128 resultname /'/tmp/EqnResult'/
#endif

C Open data and result pipes
      OPEN(EqnData, FILE=dataname, ACCESS='STREAM')
      OPEN(EqnResult, FILE=resultname, ACCESS='STREAM')

C Send number of equations
      command = command_set_size
      WRITE(EqnData) command, M

C Send rows
      command = command_set_row
      j1 = 1
      DO i = 1, M
        jn = j1 + count(i) - 1
        WRITE(EqnData) command, i, count(i), (indexes(j), j = j1, jn), (values(j), j = j1, jn)
        j1 = jn + 1
      ENDDO

C Send vector with constrain values
      command = command_set_x
      WRITE(EqnData) command, (x(j), j = 1, M)

C Send vector of independent terms
      command = command_set_b
      WRITE(EqnData) command, (b(j), j = 1, M)

C Send vector that indicates where the constrains are
      command = command_set_fixed
      WRITE(EqnData) command, (fixed(j), j = 1, M)

C Initialize solver
      command = command_solver_init
      WRITE(EqnData) command

C Run solver and read result
      command = command_solver_run
      WRITE(EqnData) command
      flush(EqnData)
      READ(EqnResult) (r(j), j = 1, M)

C Display result
      DO i = 1, M
        WRITE(6, *) r(i)
      ENDDO

C Send end session command
      command = command_end
      WRITE(EqnData) command

C Close pipes
      CLOSE(EqnResult)
      CLOSE(EqnData)

      END
