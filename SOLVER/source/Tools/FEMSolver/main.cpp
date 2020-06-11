// main.cpp
// Copyright (C) 2012 Miguel Vargas (miguel.vargas@gmail.com)
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

#include <Basic/Debug.h>
#include <Basic/Integer.h>
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
	command_end              = 0,
	command_set_connectivity = 1,
	command_fill_A           = 2,
	command_set_Ae           = 3,
	command_set_all_Ae       = 4,
	command_set_x            = 5,
	command_set_b            = 6,
	command_set_fixed        = 7,
	command_solver_init      = 8,
	command_solver_run       = 9
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
		Log(1, "FEMSolver ------------------------------------------------------------");
		Log(1, "-Version: " MacroValueToString(VERSION));

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
			data_path = "FEMData";
			result_path = "FEMResult";
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

		// Geometry data
		Mesh mesh;
		Vector<int> element_index;
		Vector<int> node_index;
		int degrees_of_freedom = 0;
		Assembler* assembler = (Assembler*)0;

		// Solver data
		CSRMatrix<double> A;
		Matrix<double> Ae;
		Vector<double> x;
		Vector<double> b;
		Vector<bool> fixed;

		// Independent terms
		Vector<double> global_b(b.size);

		// Solution
		Vector<double> global_x;

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
				case command_set_connectivity:
				{
					int nodes_count;
					int elements_count;
					int element_type;
					int nodes_per_element;

					data_pipe.Read(&nodes_count, 1);
					data_pipe.Read(&elements_count, 1);
					data_pipe.Read(&element_type, 1);
					data_pipe.Read(&nodes_per_element, 1);
					data_pipe.Read(&degrees_of_freedom, 1);

					mesh.Resize((ShapeType)element_type, nodes_per_element, elements_count, nodes_count);
					data_pipe.Read(mesh.connectivity.data, elements_count*nodes_per_element);

					// Check geometry
					{
						Sequence<int, 32> error_nodes;
						if (!mesh.CheckNodes(error_nodes))
						{
							Log(1, "ERROR in mesh, %i unconnected nodes, check 'nodes.txt'", error_nodes.size);
							FormatInteger format(true, false, 1, FormatInteger::decimal);
							File nodes_file;
							nodes_file.Create("nodes.txt");
							for (SequenceItem<int>* item = error_nodes.first; item; item = item->next)
							{
								nodes_file.Write(item->value, format);
							}
							nodes_file.Close();
						}

						Vector<Sequence<int, 32> > elements_group;
						mesh.CheckElements(elements_group);
						if (elements_group.size > 1)
						{
							Log(1, "ERROR in mesh, %i unconnected groups of elements, check 'elements.txt'", elements_group.size);
							FormatInteger format(true, false, 1, FormatInteger::decimal);
							File elements_file;
							elements_file.Create("elements.txt");
							elements_file.Write("group element_1 element_2 ...\n");
							for (int g = 1; g <= elements_group.size; ++g)
							{
								elements_file.Write(g, format);
								for (SequenceItem<int>* item = elements_group.entry[g].first; item; item = item->next)
								{
									elements_file.Write(item->value, format);
								}
								elements_file.Write("\n");
							}
							elements_file.Close();
						}
					}

					element_index.Resize(mesh.elements_count);
					element_index.FillSeries(1, 1); // We will use all elements (no partitioning)

					node_index.Resize(mesh.nodes_count);
					mesh.GetElementsNodes(element_index, solver_parameters->reorder, node_index); // Get nodes for elements (reorder if necesary)

					delete assembler;
					assembler = new Assembler(mesh, element_index, node_index, degrees_of_freedom); // Assembler class
					if (!assembler)
					{
						Throw(Memory::exception);
					}

					int M = node_index.size*degrees_of_freedom;
					A.Resize(M, M);
					assembler->AllocateMatrix(A); // Allocate memory for stiffness matrix (reordered if necesary)
					A.Fill(0);
					x.Resize(M);
					b.Resize(M);
					fixed.Resize(M);
					fixed.Fill(false);
					global_x.Resize(M);
					global_b.Resize(M);

					int N = mesh.nodes_per_element*degrees_of_freedom;
					Ae.Resize(N, N);

					continue;
				}
				case command_fill_A:
				{
					double value;
					data_pipe.Read(&value, 1);
					A.Fill(value);
					continue;
				}
				case command_set_Ae:
				{
					int element_id;
					data_pipe.Read(&element_id, 1);
					data_pipe.Read(Ae.data, Ae.rows*Ae.columns);

					Assert(assembler);
					Assert((element_id >= 1) && (element_id <= mesh.elements_count));
					assembler->AssembleAe(element_id, Ae, A);
					continue;
				}
				case command_set_all_Ae:
				{
					Assert(assembler);

					for (int element_id = 1; element_id <= element_index.size; ++element_id)
					{
						data_pipe.Read(Ae.data, Ae.rows*Ae.columns);
						assembler->AssembleAe(element_id, Ae, A);
					}
					continue;
				}
				case command_set_x:
				{
					Assert(assembler);

					data_pipe.Read(global_x.data, global_x.size);
					assembler->AssembleV(global_x, x);
					continue;
				}
				case command_set_b:
				{
					Assert(assembler);

					data_pipe.Read(global_b.data, global_b.size);
					assembler->AssembleV(global_b, b);
					continue;
				}
				case command_set_fixed:
				{
					Assert(assembler);

					Vector<bool> global_fixed(fixed.size);
					data_pipe.Read(global_fixed.data, global_fixed.size);
					assembler->AssembleV(global_fixed, fixed);
					continue;
				}
				case command_solver_init:
				{
					delete solver;
					solver = solver_parameters->Instanciate(A, x, b, fixed, false);
					continue;
				}
				case command_solver_run:
				{
					Assert(solver);

					solver->CompensateFixed();
					bool valid_solution = solver->Calculate();
					Log(1, "Solution: %s", valid_solution ? "valid" : "INVALID");

					// Store solution
					for (int n = 1; n <= node_index.size; ++n)
					{
						int i = (n - 1)*degrees_of_freedom;
						int gi = (node_index.entry[n] - 1)*degrees_of_freedom;
						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							global_x.entry[gi + d] = x.entry[i + d]; // Fill solution (inverse order if necesary)
						}
					}
					result_pipe.Write(global_x.data, global_x.size);
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
		delete assembler;
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
