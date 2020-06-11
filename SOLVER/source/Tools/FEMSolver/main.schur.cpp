// main.schur.cpp
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

#include <Basic/Integer.h>
#include <Basic/Log.h>
#include <Basic/Macros.h>
#include <Communication/MPI.h>
#include <File/File.h>
#include <FiniteElement/Assembler.h>
#include <FiniteElement/Mesh.h>
#include <FiniteElement/Partition.h>
#include <Solver/SolverParameters.h>
#include <Solver/SchurComplement.h>


// Default environment variables values
#define LOG_LEVEL                2
#define SOLVER_TYPE              1 // 1=Conjugate_gradient, 2=Cholesky_decomposition, 3=Cholesky2_decomposition, 4=Biconjugate_gradient, 5=LU_decomposition
#define SOLVER_THREADS           1
#define SOLVER_TOLERANCE         1e-5
#define SOLVER_MAX_STEPS         10000
#define PRECONDITIONER_TYPE      1 // 0=None, 1=Jacobi, 2=Incomplete_Cholesky, 3=Incomplete_Cholesky2, 4=Incomplete_LU, 5=Sparse_Approximate_Inverse
#define PRECONDITIONER_LEVEL     1
#define PRECONDITIONER_THRESHOLD 0.0
#define SUBSTRUCTURING_THREADS   1
#define SUBSTRUCTURING_TOLERANCE 1e-5


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
		return 1;
	}

	// Set log level from environment
	log_level = GetEnvInteger(LOG_LEVEL);

	try
	{
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
			int substructuring_threads = GetEnvInteger(SUBSTRUCTURING_THREADS);
			double substructuring_tolerance = GetEnvFloat(SUBSTRUCTURING_TOLERANCE);

			solver_parameters = new SolverParameters<double>(solver_type, solver_threads, solver_tolerance, solver_max_steps, preconditioner_type, preconditioner_level, preconditioner_threshold, substructuring_threads, substructuring_tolerance, false);
		}

		// MPI object
		MPI mpi(argc, argv);
		if (mpi.size < 3)
		{
			if (mpi.rank == 0)
			{
				fputs("[Error] Invalid number processes to start, minimum 3 processes should be started.\n", stderr);
			}
			return 1;
		}

		if (mpi.rank == 0) // Master
		{
			int partitions_count = mpi.size - 1;

			Log(1, "FEMSolver.Schur ------------------------------------------------------");
			Log(1, "-Version:      " MacroValueToString(VERSION));
			Log(1, "-Partitions:    %i", partitions_count);

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

			int M = 0; // Matrix size
			int N = 0; // Elemental matrix size

			// Geometry data
			Mesh mesh;
			Vector<int> element_index;
			Vector<int> node_index;
			int degrees_of_freedom = 0;

			// Elemental matrices
			Vector<Matrix<double> > elemental_matrices;

			// Global vectors
			Vector<double> global_b;
			Vector<double> global_x;
			Vector<bool> global_fixed;

			// Substructuration data
			Vector<Substructure> substructures;
			Boundary boundary;

			// Main loop
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

						M = node_index.size*degrees_of_freedom;
						global_x.Resize(M);
						global_b.Resize(M);
						global_fixed.Resize(M);
						global_fixed.Fill(false);

						N = mesh.nodes_per_element*degrees_of_freedom;
						elemental_matrices.Resize(elements_count);
						for (int element_id = 1; element_id <= element_index.size; ++element_id)
						{
							elemental_matrices.entry[element_id].Resize(N, N);
						}
						continue;
					}
					case command_fill_A:
					{
						continue;
					}
					case command_set_Ae:
					{
						int element_id;
						data_pipe.Read(&element_id, 1);
						data_pipe.Read(elemental_matrices.entry[element_id].data, N*N);
						continue;
					}
					case command_set_all_Ae:
					{
						for (int element_id = 1; element_id <= element_index.size; ++element_id)
						{
							data_pipe.Read(elemental_matrices.entry[element_id].data, N*N);
						}
						continue;
					}
					case command_set_x:
					{
						data_pipe.Read(global_x.data, global_x.size);
						continue;
					}
					case command_set_b:
					{
						data_pipe.Read(global_b.data, global_b.size);
						continue;
					}
					case command_set_fixed:
					{
						data_pipe.Read(global_fixed.data, global_fixed.size);
						continue;
					}
					case command_solver_init:
					{
						StructurePartitioning(mesh, partitions_count, substructures, boundary);
						Log(1, "Partitioning done");
						continue;
					}
					case command_solver_run:
					{
						// Start Schur complement routines
						SchurComplementMaster<double> schur_complement_master(mpi, substructures, boundary, degrees_of_freedom);

						Log(1, "Systems of equations -------------------------------------------------");
						Log(1, "-Local solver threads: %i", solver_parameters->threads);

						// Systems of equations
						#pragma omp parallel for default(shared) schedule(dynamic,1) num_threads(solver_parameters->substructuring_threads)
						for (int p = 1; p <= partitions_count; ++p)
						{
							const Vector<int>& element_index = substructures.entry[p].element_index;
							const Vector<int>& node_index = substructures.entry[p].node_index;

							// Assembler class
							Assembler assembler(mesh, element_index, node_index, degrees_of_freedom);

							// Matrices size
							int Mi = node_index.size*degrees_of_freedom; // Matrix size
							int Bi = substructures.entry[p].boundary_count*degrees_of_freedom;
							int Ii = Mi - Bi;

							// Assemble K
							CSRMatrix<double> K(Mi, Mi);
							assembler.AllocateMatrix(K);
							K.Fill(0);
							for (int i = 1; i <= element_index.size; ++i)
							{
								int element_id = element_index.entry[i];
								assembler.AssembleAe(element_id, elemental_matrices.entry[element_id], K);
							}

							// Assemble x
							Vector<double> x(Mi);
							assembler.AssembleV(global_x, x);

							// Assemble b
							Vector<double> b(Mi);
							assembler.AssembleV(global_b, b);

							// Assemble fixed
							Vector<bool> fixed(Mi);
							assembler.AssembleV(global_fixed, fixed);

							// Initialize Schur method for this partition
							schur_complement_master.SendSystemOfEquations(p, K, x, b, fixed);

							Log(1, "Partition %i:", p);
							Log(1, "-Degrees of freedom: %i (%i + %i)", Mi, Ii, Bi);
							Log(1, "-nnz(K):             %i", K.NonZero());
						}

						// Allocate ABB, xB, bB and fixedB
						int B = boundary.node_index.size*degrees_of_freedom;
						CSRMatrix<double> ABB(B, B);
						Vector<double> xB(B);
						Vector<double> bB(B);
						Vector<bool> fixedB(B);

						// Assembler class
						Assembler assemblerB(mesh, boundary.element_index, boundary.node_index, degrees_of_freedom);

						// Assemble ABB
						assemblerB.AllocateMatrix(ABB);
						ABB.Fill(0.0);
						for (int i = 1; i <= boundary.element_index.size; ++i)
						{
							int element_id = boundary.element_index.entry[i];
							assemblerB.AssembleAe(element_id, elemental_matrices.entry[element_id], ABB);
						}

						// Assemble xB
						assemblerB.AssembleV(global_x, xB);

						// Assemble bB
						assemblerB.AssembleV(global_b, bB);

						// Assemble fixedB
						assemblerB.AssembleV(global_fixed, fixedB);

						Log(1, "Boundaries -----------------------------------------------------------");
						Log(1, "-Degrees of freedom: %i", B);
						Log(1, "-nnz(ABB):           %i", ABB.NonZero());

						Log(1, "Substructuring -------------------------------------------------------");
						Log(1, "-Threads:          %i", solver_parameters->substructuring_threads);
						Log(1, "-Tolerance:        %.5g", solver_parameters->substructuring_tolerance);

						schur_complement_master.InitializeSolver(ABB, xB, bB, fixedB);
						schur_complement_master.ComplesateFixed(xB, bB, fixedB, solver_parameters->substructuring_threads);
						schur_complement_master.Solver(ABB, xB, bB, solver_parameters->substructuring_tolerance, solver_parameters->max_steps, solver_parameters->substructuring_threads);
						schur_complement_master.CompleteSolution(xB, global_x);

						result_pipe.Write(global_x.data, global_x.size);
						result_pipe.Flush();

						if (Memory::memory_usage)
						{
							Log(1, "Peak allocated memory ------------------------------------------------");
							Log(1, "-Master:        %12lu bytes", (unsigned long)Memory::peak_usage);
							Vector<size_t> partition_memory_peak_usage(partitions_count);
							schur_complement_master.GetMemoryPeakUsage(partition_memory_peak_usage);
							size_t total_memory_peak_usage = Memory::peak_usage;
							for (int p = 1; p <= partitions_count; ++p)
							{
								total_memory_peak_usage += partition_memory_peak_usage.entry[p];
								Log(1, "-Partition %3i: %12lu bytes", p, (unsigned long)partition_memory_peak_usage.entry[p]);
							}
							Log(1, "-Total:         %12lu bytes", (unsigned long)total_memory_peak_usage);
						}

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
		}
		else // rank != 0
		{
			log_level = 0; // Prevents slaves logs to appear

			SchurComplementSlave<double> schur_complement_slave(mpi, solver_parameters->threads);
			schur_complement_slave.ReceiveSystemOfEquations();
			schur_complement_slave.Solver();
			schur_complement_slave.CompleteSolution();
			if (Memory::memory_usage)
			{
				schur_complement_slave.ReportMemoryPeakUsage();
			}
		}

		delete solver_parameters;
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
