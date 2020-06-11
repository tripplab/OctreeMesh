// main.cpp
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

#include <Basic/Log.h>
#include <Basic/Macros.h>
#include <File/File.h>
#include <FiniteElement/Assembler.h>
#include <FiniteElement/Mesh.h>
#include <Solver/SolverParameters.h>


// Default environment variables values
#define LOG_LEVEL                2
#define SOLVER_TYPE              1 // 1=Conjugate_gradient, 2=Cholesky_decomposition, 3=Cholesky2_decomposition, 4=Biconjugate_gradient, 5=LU_decomposition
#define SOLVER_THREADS           1
#define SOLVER_TOLERANCE         1e-5
#define SOLVER_MAX_STEPS         10000
#define PRECONDITIONER_TYPE      1 // 0=None, 1=Jacobi, 2=Incomplete_Cholesky, 3=Incomplete_Cholesky2, 4=Incomplete_LU, 5=Sparse_Approximate_Inverse
#define PRECONDITIONER_LEVEL     1
#define PRECONDITIONER_THRESHOLD 0.0


enum Command
{
	command_end         = 0,
	command_set_size    = 1,
	command_set_row     = 2,
	command_set_x       = 3,
	command_set_b       = 4,
	command_set_fixed   = 5,
	command_solver_init = 6,
	command_solver_run  = 7
};


int main(int argc, char** argv)
{
	if (argc > 4)
	{
		fprintf(stderr, "Invalid number of arguments. Use:\n  %s [parameters_file] [<data_pipe_name> <results_pipe_name>]\n", argv[0]);
		fputs("In Windows systems, pipe names must have the form: \\\\.\\pipe\\[pipe_name]\n", stderr);
		return 1;
	}

	// Set log level from environment
	log_level = GetEnvInteger(LOG_LEVEL);

	try
	{
		Log(1, "------------------------------ EqnSolver -----------------------------");
		Log(1, "Version:   %s", MacroValueToString(VERSION));

		// Load solver parameters
		SolverParameters<double>* solver_parameters = (SolverParameters<double>*)0;
		if ((argc == 2) || (argc == 4))
		{
			solver_parameters = new SolverParameters<double>(argv[1]);
		}
		else
		{
			SolverType solver_type = (SolverType)GetEnvInteger(SOLVER_TYPE);
			int solver_threads = GetEnvInteger(SOLVER_THREADS);
			double solver_tolerance = GetEnvFloat(SOLVER_TOLERANCE);
			int solver_max_steps = GetEnvInteger(SOLVER_MAX_STEPS);
			PreconditionerType preconditioner_type = (PreconditionerType)GetEnvInteger(PRECONDITIONER_TYPE);
			int preconditioner_level = GetEnvInteger(PRECONDITIONER_LEVEL);
			double preconditioner_threshold = GetEnvFloat(PRECONDITIONER_THRESHOLD);

			solver_parameters = new SolverParameters<double>(solver_type, solver_threads, solver_tolerance, solver_max_steps, preconditioner_type, preconditioner_level, preconditioner_threshold);
		}
		solver_parameters->PrintInfo();

		// Create and open pipes
		Log(1, "Create pipes ---------------------------------------------------------");
		const char* data_path = (char*)0;
		const char* result_path = (char*)0;
		File data_pipe;
		File result_pipe;
		if (argc == 3)
		{
			data_path = argv[1];
			result_path = argv[2];
		}
		else if (argc == 4)
		{
			data_path = argv[2];
			result_path = argv[3];
		}
		else
		{
			data_path = "EqnData";
			result_path = "EqnResult";
		}
		#ifdef WIN32
			Log(1, "Data pipe: \\\\.\\pipe\\%s", data_path);
			Log(1, "Result pipe: \\\\.\\pipe\\%s", result_path);
		#else
			Log(1, "Data pipe: /tmp/%s", data_path);
			Log(1, "Result pipe: /tmp/%s", result_path);
		#endif
		#pragma omp parallel sections num_threads(2)
		{
			#pragma omp section
			{
				data_pipe.NamedPipe(data_path, named_pipe_read_only);
			}
			#pragma omp section
			{
				result_pipe.NamedPipe(result_path, named_pipe_write_only);
			}
		}

		// Solver data
		CSRMatrix<double> A;
		Vector<double> x;
		Vector<double> b;
		Vector<bool> fixed;

		// For systems that require reordering
		Vector<int> index;
		Vector<int> inverse_index;
		CSRMatrix<double> A_r;
		Vector<double> x_r;
		Vector<double> b_r;
		Vector<bool> fixed_r;

		Solver<double>* solver = (Solver<double>*)0;
		for (bool go = true; go; )
		{
			char command;
			data_pipe.Read(command);
			switch ((Command)command)
			{
				case command_end:
				{
					result_pipe.Close();
					data_pipe.Close();
					go = false;
					continue;
				}
				case command_set_size:
				{
					int size;

					data_pipe.Read(&size, 1);
					A.Resize(size, size);
					x.Resize(size);
					b.Resize(size);
					fixed.Resize(size);
					if (solver_parameters->reorder)
					{
						A_r.Resize(size, size);
						x_r.Resize(size);
						b_r.Resize(size);
						fixed_r.Resize(size);
					}
					continue;
				}
				case command_set_row:
				{
					int row;
					int count;
					data_pipe.Read(&row, 1);
					data_pipe.Read(&count, 1);
					A.AllocateRow(row, count);
					data_pipe.Read(&A.index[row][1], count);
					data_pipe.Read(&A.entry[row][1], count);
					continue;
				}
				case command_set_x:
				{
					data_pipe.Read(x.data, x.size);
					continue;
				}
				case command_set_b:
				{
					data_pipe.Read(b.data, b.size);
					continue;
				}
				case command_set_fixed:
				{
					data_pipe.Read(fixed.data, fixed.size);
					continue;
				}
				case command_solver_init:
				{
					delete solver;
					if (solver_parameters->reorder)
					{
						FindingAnOrdering(A, index, inverse_index);
						Reorder(index, inverse_index, A, A_r, solver_parameters->threads);
						Reorder(index, x, x_r, solver_parameters->threads);
						Reorder(index, fixed, fixed_r, solver_parameters->threads);
						solver = solver_parameters->Instanciate(A_r, x_r, b_r, fixed_r, true);
					}
					else
					{
						solver = solver_parameters->Instanciate(A, x, b, fixed, true);
					}
					continue;
				}
				case command_solver_run:
				{
					Assert(solver);

					if (solver_parameters->reorder)
					{
						Reorder(index, b, b_r, solver_parameters->threads);
					}
					solver->CompensateFixed();

					bool valid_solution = solver->Calculate();
					Log(1, "Solution: %s", valid_solution ? "valid" : "INVALID");

					if (solver_parameters->reorder)
					{
						Reorder(inverse_index, x_r, x, solver_parameters->threads);
					}
					result_pipe.Write(x.data, x.size);
					result_pipe.Flush();
					continue;
				}
				default:
				{
					Log(1, "Invalid command %i", command);
					result_pipe.Close();
					data_pipe.Close();
					go = false;
					continue;
				}
			}
		}

		delete solver;
		delete solver_parameters;

		if (Memory::memory_usage)
		{
			Log(1, "Peak allocated memory: %lu bytes", (unsigned long)Memory::peak_usage);
		}
	}
	catch (Exception&)
	{
		DebugPosition("Catch fatal exception");
	}

	if (Memory::current_usage != 0)
	{
		fprintf(stderr, "[Error] Memory leak: %lu bytes\n", (unsigned long)Memory::current_usage);
	}

 	return 0;
}
