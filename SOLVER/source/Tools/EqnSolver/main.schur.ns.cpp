// main.schur.ns.cpp
// Copyright (C) 2013 Miguel Vargas (miguel.vargas@gmail.com)
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
#include <Communication/MPI.h>
#include <File/File.h>
#include <FiniteElement/Assembler.h>
#include <FiniteElement/Mesh.h>
#include <FiniteElement/Partition.h>
#include <Solver/SolverParameters.h>
#include <Solver/SchurNS.h>


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

			Log(1, "EqnSolver.Schur ------------------------------------------------------");
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
			Vector<double> global_x;
			Vector<double> global_b;
			Vector<bool> global_fixed;

			// Substructuration data
			Vector<Substructure> substructures;
			Boundary boundary;

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
						global_x.Resize(size);
						global_b.Resize(size);
						global_fixed.Resize(size);
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
						StructurePartitioning(A, partitions_count, substructures, boundary);
						Log(1, "Partitioning done");
						continue;
					}
					case command_solver_run:
					{
						// Start Schur complement routines
						SchurNSMaster<double> schur_ns_master(mpi, substructures, boundary, 1);

						Log(1, "Systems of equations -------------------------------------------------");
						Log(1, "-Local solver threads: %i", solver_parameters->threads);

						// Systems of equations
						#pragma omp parallel for default(shared) schedule(dynamic,1) num_threads(solver_parameters->substructuring_threads)
						for (int p = 1; p <= partitions_count; ++p)
						{
							const Vector<int>& node_index = substructures.entry[p].node_index;

							Vector<int> global_to_local(A.rows);
							global_to_local.Fill(0);
							for (int i = 1; i <= node_index.size; ++i)
							{
								int n = node_index.entry[i];
								global_to_local.entry[n] = i;
							}

							// Sizes
							int Mi = node_index.size; // Matrix size
							int Bi = substructures.entry[p].boundary_count;
							int Ii = Mi - Bi;

							// Subsystem of equations
							CSRMatrix<double> K(Mi, Mi);
							Vector<double> x(Mi);
							Vector<double> b(Mi);
							Vector<bool> fixed(Mi);
							Vector<int> K_index_i(Mi);
							Vector<double> K_entry_i(Mi);
							for (int i = 1; i <= node_index.size; ++i)
							{
								int n = node_index.entry[i];

								int K_count_i = 0;
								int A_count_n = A.count[n];
								for (int k = 1; k <= A_count_n; ++k)
								{
									int gj = A.index[n][k];
									int lj = global_to_local.entry[gj];
									if (lj != 0)
									{
										++K_count_i;
										K_index_i.entry[K_count_i] = lj;
										K_entry_i.entry[K_count_i] = A.entry[n][k];
									}
								}
								K.AllocateRow(i, K_count_i);
								for (int k = 1; k <= K_count_i; ++k)
								{
									K.index[i][k] = K_index_i.entry[k];
									K.entry[i][k] = K_entry_i.entry[k];
								}

								x.entry[i] = global_x.entry[n];
								b.entry[i] = global_b.entry[n];
								fixed.entry[i] = global_fixed.entry[n];
							}

							// Initialize Schur method for this partition
							schur_ns_master.SendSystemOfEquations(p, K, x, b, fixed);

							Log(1, "Partition %i:", p);
							Log(1, "-Degrees of freedom: %i (%i + %i)", Mi, Ii, Bi);
							Log(1, "-nnz(K):             %i", K.NonZero());
						}

						// Allocate ABB, xB, bB and fixedB
						int B = boundary.node_index.size;
						CSRMatrix<double> ABB(B, B);
						Vector<double> xB(B);
						Vector<double> bB(B);
						Vector<bool> fixedB(B);

						{
							Vector<int> global_to_local(A.rows);
							global_to_local.Fill(0);
							for (int i = 1; i <= B; ++i)
							{
								int n = boundary.node_index.entry[i];
								global_to_local.entry[n] = i;
							}

							Vector<int> ABB_index_i(B);
							Vector<double> ABB_entry_i(B);
							for (int i = 1; i <= B; ++i)
							{
								int n = boundary.node_index.entry[i];

								int ABB_count_i = 0;
								int A_count_n = A.count[n];
								for (int k = 1; k <= A_count_n; ++k)
								{
									int gj = A.index[n][k];
									int lj = global_to_local.entry[gj];
									if (lj != 0)
									{
										++ABB_count_i;
										ABB_index_i.entry[ABB_count_i] = lj;
										ABB_entry_i.entry[ABB_count_i] = A.entry[n][k];
									}
								}
								ABB.AllocateRow(i, ABB_count_i);
								for (int k = 1; k <= ABB_count_i; ++k)
								{
									ABB.index[i][k] = ABB_index_i.entry[k];
									ABB.entry[i][k] = ABB_entry_i.entry[k];
								}

								xB.entry[i] = global_x.entry[n];
								bB.entry[i] = global_b.entry[n];
								fixedB.entry[i] = global_fixed.entry[n];
							}
						}

						Log(1, "Boundaries -----------------------------------------------------------");
						Log(1, "-Degrees of freedom: %i", B);
						Log(1, "-nnz(ABB):           %i", ABB.NonZero());

						Log(1, "Substructuring -------------------------------------------------------");
						Log(1, "-Threads:          %i", solver_parameters->substructuring_threads);
						Log(1, "-Tolerance:        %.5g", solver_parameters->substructuring_tolerance);

						schur_ns_master.InitializeSolver(ABB, xB, bB, fixedB);
						schur_ns_master.ComplesateFixed(xB, bB, fixedB, solver_parameters->substructuring_threads);

						CSRMatrix<double> ABBt(B, B);
						Transpose(ABB, ABBt);
						schur_ns_master.Solver(ABB, ABBt, xB, bB, solver_parameters->substructuring_tolerance, solver_parameters->max_steps, solver_parameters->substructuring_threads);
						schur_ns_master.CompleteSolution(xB, global_x);

						result_pipe.Write(global_x.data, global_x.size);
						result_pipe.Flush();

						if (Memory::memory_usage)
						{
							Log(1, "Peak allocated memory ------------------------------------------------");
							Log(1, "-Master:        %12lu bytes", (unsigned long)Memory::peak_usage);
							Vector<size_t> partition_memory_peak_usage(partitions_count);
							schur_ns_master.GetMemoryPeakUsage(partition_memory_peak_usage);
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

			SchurNSSlave<double> schur_ns_slave(mpi, solver_parameters->threads);
			schur_ns_slave.ReceiveSystemOfEquations();
			schur_ns_slave.Solver();
			schur_ns_slave.CompleteSolution();
			if (Memory::memory_usage)
			{
				schur_ns_slave.ReportMemoryPeakUsage();
			}
		}

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
