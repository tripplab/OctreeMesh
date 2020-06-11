// EqnSolverExample.cpp
// Copyright (C) 2011 Miguel Vargas (miguel.vargas@gmail.com)
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

//
// Compilation
// GCC:
//     g++ EqnSolverExample.cpp
// Intel:
//     icpc EqnSolverExample.cpp
//

#include <stdio.h>

int main()
{

// Commands
	enum Command {
		command_end         = 0,
		command_set_size    = 1,
		command_set_row     = 2,
		command_set_x       = 3,
		command_set_b       = 4,
		command_set_fixed   = 5,
		command_solver_init = 6,
		command_solver_run  = 7};
// Number of equations
	int M = 22;
// Matrix values
	double values[120] = {
		0.000123262649322419, -3.60513546547159e-005, -3.39967111532538e-005, -5.32145835144495e-005,
		-3.60513546547159e-005, 0.000170510996033461, -4.34952789122033e-005, -6.08778004323386e-005, -3.0086562034203e-005,
		-3.39967111532538e-005, 0.000170318953739018, -5.57481717486078e-005, -4.56652509159795e-005, -3.49088199211772e-005,
		-5.32145835144495e-005, -4.34952789122033e-005, -5.57481717486078e-005, 0.000341662307627607, -7.11430412224589e-005, -6.86630215458964e-005, -4.93982106839912e-005,
		-6.08778004323386e-005, -7.11430412224589e-005, 0.000342667616164504, -5.86038603040437e-005, -4.42802365962771e-005, -4.31018496289288e-005, -6.46608279804565e-005,
		-4.56652509159795e-005, -6.86630215458964e-005, 0.000341547174453471, -5.62295144334364e-005, -4.42702807061953e-005, -5.89712307716659e-005, -6.77478760802972e-005,
		-3.0086562034203e-005, -5.86038603040437e-005, 0.000122792254967864, -3.41018326296169e-005,
		-3.49088199211772e-005, -5.62295144334364e-005, 0.000122915210160465, -3.17768758058511e-005,
		-4.93982106839912e-005, -4.42802365962771e-005, -4.42702807061953e-005, 0.000329610569383882, -4.72254168572741e-005, -5.27761219876592e-005, -5.23122253623367e-005, -3.9348077190148e-005,
		-5.89712307716659e-005, -3.17768758058511e-005, 0.000170291900663079, -4.46349962484123e-005, -3.49087978371502e-005,
		-4.31018496289288e-005, -3.41018326296169e-005, 0.000170883134605931, -6.64646115422989e-005, -2.72148408050865e-005,
		-6.46608279804565e-005, -4.72254168572741e-005, -6.64646115422989e-005, 0.000342335610854092, -5.09176371517581e-005, -6.2797286976477e-005, -5.02698303458277e-005,
		-6.77478760802972e-005, -5.27761219876592e-005, -4.46349962484123e-005, 0.000341041263674863, -5.81729055099591e-005, -5.62295230426859e-005, -6.14798408058491e-005,
		-5.23122253623367e-005, -5.09176371517581e-005, 0.000342524050812863, -7.14963888156276e-005, -5.57119746311947e-005, -5.11077281633579e-005, -6.09780966885886e-005,
		-3.9348077190148e-005, -5.81729055099591e-005, -7.14963888156276e-005, 0.000344682167058638, -4.53552672668097e-005, -6.39237466853512e-005, -6.63857815907427e-005,
		-3.49087978371502e-005, -5.62295230426859e-005, 0.000122915205551077, -3.17768846712407e-005,
		-2.72148408050865e-005, -6.2797286976477e-005, 0.000122870955324627, -3.28588275430632e-005,
		-6.14798408058491e-005, -4.53552672668097e-005, -3.17768846712407e-005, 0.000170574963019778, -3.19629702758784e-005,
		-5.02698303458277e-005, -5.57119746311947e-005, -3.28588275430632e-005, 0.000170881035319409, -3.20404027993234e-005,
		-5.11077281633579e-005, -6.39237466853512e-005, 0.000169496564745875, -2.48620285243623e-005, -2.96030613728036e-005,
		-6.63857815907427e-005, -3.19629702758784e-005, -2.48620285243623e-005, 0.000123210780390983,
		-6.09780966885886e-005, -3.20404027993234e-005, -2.96030613728036e-005, 0.000122621560860716};
// Matrix indexes
	int indexes[120]= {
		1, 2, 3, 4,
		1, 2, 4, 5, 7,
		1, 3, 4, 6, 8,
		1, 2, 3, 4, 5, 6, 9,
		2, 4, 5, 7, 9, 11, 12,
		3, 4, 6, 8, 9, 10, 13,
		2, 5, 7, 11,
		3, 6, 8, 10,
		4, 5, 6, 9, 12, 13, 14, 15,
		6, 8, 10, 13, 16,
		5, 7, 11, 12, 17,
		5, 9, 11, 12, 14, 17, 19,
		6, 9, 10, 13, 15, 16, 18,
		9, 12, 14, 15, 19, 20, 22,
		9, 13, 14, 15, 18, 20, 21,
		10, 13, 16, 18,
		11, 12, 17, 19,
		13, 15, 16, 18, 21,
		12, 14, 17, 19, 22,
		14, 15, 20, 21, 22,
		15, 18, 20, 21,
		14, 19, 20, 22};
// Row sizes
	int count[22] = {4, 5, 5, 7, 7, 7, 4, 4, 8, 5, 5, 7, 7, 7, 7, 4, 4, 5, 5, 5, 4, 4};
// Vector with constrain values
	double x[22] = {273, 273, 0, 0, 0, 0, 273, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 293, 0};
// Vector of independent terms
	double b[22] = {0, 0, 0, 0, 0, 0, 0, -0.000433883773336144, 0, -0.000867767546672287, 0, 0, 0, 0, 0, -0.000433883773336144, -0.000216941883549719, 0, -0.000433883718353043, 0, 0, -0.000216941834803323};
// Vector that indicates where the constrains are
	bool fixed[22] = {1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0};
// Result vector
	double r[22];
// Variable to indicate command
	char command;
// Data pipe
	FILE* EqnData;
// Result pipe
	FILE* EqnResult;
// Indexes
	int i, j1;
// Names for the pipes
#ifdef WIN32
	const char* dataname = "\\\\.\\pipe\\EqnData";
	const char* resultname = "\\\\.\\pipe\\EqnResult";
#else
	const char* dataname = "/tmp/EqnData";
	const char* resultname = "/tmp/EqnResult";
#endif

// Open data and result pipes
	EqnData = fopen(dataname, "wb");
	EqnResult = fopen(resultname, "rb");

// Send number of equations
	command = command_set_size;
	fwrite(&command, 1, 1, EqnData);
	fwrite(&M, sizeof(int), 1, EqnData);

// Send rows
	command = command_set_row;
	j1 = 0;
	for (i = 0; i < M; ++i)
	{
		int row = i + 1;
		fwrite(&command, 1, 1, EqnData);
		fwrite(&row, sizeof(int), 1, EqnData);
		fwrite(&count[i], sizeof(int), 1, EqnData);
		fwrite(indexes + j1, sizeof(int), count[i], EqnData);
		fwrite(values + j1, sizeof(double), count[i], EqnData);
		j1 += count[i];
	}

	// Send vector with constrain values
	command = command_set_x;
	fwrite(&command, 1, 1, EqnData);
	fwrite(x, sizeof(double), M, EqnData);

	// Send vector of independent terms
	command = command_set_b;
	fwrite(&command, 1, 1, EqnData);
	fwrite(b, sizeof(double), M, EqnData);

	// Send vector that indicates where the constrains are
	command = command_set_fixed;
	fwrite(&command, 1, 1, EqnData);
	fwrite(fixed, sizeof(bool), M, EqnData);

	// Initialize solver
	command = command_solver_init;
	fwrite(&command, 1, 1, EqnData);

// Run solver and read result
	command = command_solver_run;
	fwrite(&command, 1, 1, EqnData);
	fflush(EqnData);
	fread(r, sizeof(double), M, EqnResult);

// Display result
	for (i = 0; i < M; ++i)
	{
		printf("%f\n", r[i]);
	}

// Send end session command
	command = command_end;
	fwrite(&command, 1, 1, EqnData);
	fflush(EqnData);

// Close pipes
	fclose(EqnResult);
	fclose(EqnData);

 	return 0;
}
