# EqnSolverExample.py
# Copyright (C) 2012 Miguel Vargas (miguel.vargas@gmail.com)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, write to the Free
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import sys
import array

# Commands
command_end         = array.array('b', [0])
command_set_size    = array.array('b', [1])
command_set_row     = array.array('b', [2])
command_set_x       = array.array('b', [3])
command_set_b       = array.array('b', [4])
command_set_fixed   = array.array('b', [5])
command_solver_init = array.array('b', [6])
command_solver_run  = array.array('b', [7])
# Number of equations
M = 22
# Matrix values
values = [ \
    array.array('d', [0.000123262649322419, -3.60513546547159e-005, -3.39967111532538e-005, -5.32145835144495e-005]), \
    array.array('d', [-3.60513546547159e-005, 0.000170510996033461, -4.34952789122033e-005, -6.08778004323386e-005, -3.0086562034203e-005]), \
    array.array('d', [-3.39967111532538e-005, 0.000170318953739018, -5.57481717486078e-005, -4.56652509159795e-005, -3.49088199211772e-005]), \
    array.array('d', [-5.32145835144495e-005, -4.34952789122033e-005, -5.57481717486078e-005, 0.000341662307627607, -7.11430412224589e-005, -6.86630215458964e-005, -4.93982106839912e-005]), \
    array.array('d', [-6.08778004323386e-005, -7.11430412224589e-005, 0.000342667616164504, -5.86038603040437e-005, -4.42802365962771e-005, -4.31018496289288e-005, -6.46608279804565e-005]), \
    array.array('d', [-4.56652509159795e-005, -6.86630215458964e-005, 0.000341547174453471, -5.62295144334364e-005, -4.42702807061953e-005, -5.89712307716659e-005, -6.77478760802972e-005]), \
    array.array('d', [-3.0086562034203e-005, -5.86038603040437e-005, 0.000122792254967864, -3.41018326296169e-005]), \
    array.array('d', [-3.49088199211772e-005, -5.62295144334364e-005, 0.000122915210160465, -3.17768758058511e-005]), \
    array.array('d', [-4.93982106839912e-005, -4.42802365962771e-005, -4.42702807061953e-005, 0.000329610569383882, -4.72254168572741e-005, -5.27761219876592e-005, -5.23122253623367e-005, -3.9348077190148e-005]), \
    array.array('d', [-5.89712307716659e-005, -3.17768758058511e-005, 0.000170291900663079, -4.46349962484123e-005, -3.49087978371502e-005]), \
    array.array('d', [-4.31018496289288e-005, -3.41018326296169e-005, 0.000170883134605931, -6.64646115422989e-005, -2.72148408050865e-005]), \
    array.array('d', [-6.46608279804565e-005, -4.72254168572741e-005, -6.64646115422989e-005, 0.000342335610854092, -5.09176371517581e-005, -6.2797286976477e-005, -5.02698303458277e-005]), \
    array.array('d', [-6.77478760802972e-005, -5.27761219876592e-005, -4.46349962484123e-005, 0.000341041263674863, -5.81729055099591e-005, -5.62295230426859e-005, -6.14798408058491e-005]), \
    array.array('d', [-5.23122253623367e-005, -5.09176371517581e-005, 0.000342524050812863, -7.14963888156276e-005, -5.57119746311947e-005, -5.11077281633579e-005, -6.09780966885886e-005]), \
    array.array('d', [-3.9348077190148e-005, -5.81729055099591e-005, -7.14963888156276e-005, 0.000344682167058638, -4.53552672668097e-005, -6.39237466853512e-005, -6.63857815907427e-005]), \
    array.array('d', [-3.49087978371502e-005, -5.62295230426859e-005, 0.000122915205551077, -3.17768846712407e-005]), \
    array.array('d', [-2.72148408050865e-005, -6.2797286976477e-005, 0.000122870955324627, -3.28588275430632e-005]), \
    array.array('d', [-6.14798408058491e-005, -4.53552672668097e-005, -3.17768846712407e-005, 0.000170574963019778, -3.19629702758784e-005]), \
    array.array('d', [-5.02698303458277e-005, -5.57119746311947e-005, -3.28588275430632e-005, 0.000170881035319409, -3.20404027993234e-005]), \
    array.array('d', [-5.11077281633579e-005, -6.39237466853512e-005, 0.000169496564745875, -2.48620285243623e-005, -2.96030613728036e-005]), \
    array.array('d', [-6.63857815907427e-005, -3.19629702758784e-005, -2.48620285243623e-005, 0.000123210780390983]), \
    array.array('d', [-6.09780966885886e-005, -3.20404027993234e-005, -2.96030613728036e-005, 0.000122621560860716])]
# Matrix indexes
indexes = [ \
    array.array('i', [1, 2, 3, 4]), \
    array.array('i', [1, 2, 4, 5, 7]), \
    array.array('i', [1, 3, 4, 6, 8]), \
    array.array('i', [1, 2, 3, 4, 5, 6, 9]), \
    array.array('i', [2, 4, 5, 7, 9, 11, 12]), \
    array.array('i', [3, 4, 6, 8, 9, 10, 13]), \
    array.array('i', [2, 5, 7, 11]), \
    array.array('i', [3, 6, 8, 10]), \
    array.array('i', [4, 5, 6, 9, 12, 13, 14, 15]), \
    array.array('i', [6, 8, 10, 13, 16]), \
    array.array('i', [5, 7, 11, 12, 17]), \
    array.array('i', [5, 9, 11, 12, 14, 17, 19]), \
    array.array('i', [6, 9, 10, 13, 15, 16, 18]), \
    array.array('i', [9, 12, 14, 15, 19, 20, 22]), \
    array.array('i', [9, 13, 14, 15, 18, 20, 21]), \
    array.array('i', [10, 13, 16, 18]), \
    array.array('i', [11, 12, 17, 19]), \
    array.array('i', [13, 15, 16, 18, 21]), \
    array.array('i', [12, 14, 17, 19, 22]), \
    array.array('i', [14, 15, 20, 21, 22]), \
    array.array('i', [15, 18, 20, 21]), \
    array.array('i', [14, 19, 20, 22])]
# Row sizes
count = array.array('i', [4, 5, 5, 7, 7, 7, 4, 4, 8, 5, 5, 7, 7, 7, 7, 4, 4, 5, 5, 5, 4, 4])
# Vector with constrain values
x = array.array('d', [273, 273, 0, 0, 0, 0, 273, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 293, 0])
# Vector of independent terms
b = array.array('d', [0, 0, 0, 0, 0, 0, 0, -0.000433883773336144, 0, -0.000867767546672287, 0, 0, 0, 0, 0, -0.000433883773336144, -0.000216941883549719, 0, -0.000433883718353043, 0, 0, -0.000216941834803323])
# Vector that indicates where the constrains are
fixed = array.array('b', [1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0])
# Result vector
r = array.array('d')

# Paths for pipes
if sys.platform == "win32":
    dataname = r'\\.\pipe\EqnData'
    resultname = r'\\.\pipe\EqnResult'
else:
    dataname = r'/tmp/EqnData'
    resultname = r'/tmp/EqnResult'

# Open data and result pipes
EqnData = open(dataname, "wb")
EqnResult = open(resultname, "rb")

# Send number of equations
command_set_size.tofile(EqnData)
array.array('i', [M]).tofile(EqnData)

# Send rows
for i in range(0, M):
    command_set_row.tofile(EqnData)
    array.array('i', [i + 1]).tofile(EqnData)
    array.array('i', [count[i]]).tofile(EqnData)
    indexes[i].tofile(EqnData)
    values[i].tofile(EqnData)

# Send vector with constrain values
command_set_x.tofile(EqnData)
x.tofile(EqnData)

# Send vector of independent terms
command_set_b.tofile(EqnData)
b.tofile(EqnData)

# Send vector that indicates where the constrains are
command_set_fixed.tofile(EqnData)
fixed.tofile(EqnData)

# Initialize solver
command_solver_init.tofile(EqnData)

# Run solver and read result
command_solver_run.tofile(EqnData)
EqnData.flush()
r.fromfile(EqnResult, M)

# Display result
for v in r[:]:
    print(v)

# Send end session command
command_end.tofile(EqnData)

# Close pipes
EqnResult.close()
EqnData.close()
