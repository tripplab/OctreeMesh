// main.schur.cpp
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

#include "Solid.h"

#include <Basic/Log.h>
#include <Basic/Macros.h>
#include <Container/CSRMatrix.h>
#include <Container/Vector.h>
#include <Communication/MPI.h>
#include <FiniteElement/Assembler.h>
#include <FiniteElement/Geometry.h>
#include <FiniteElement/Partition.h>
#include <FiniteElement/ShapeFunctions.h>
#include <Solver/SchurComplement.h>
#include <Solver/SolverParameters.h>

#include <omp.h>


#define PREFIX_MAXIMUM_LENGTH 900
#define PATH_MAXIMUM_LENGTH 1000


int main(int argc, char** argv)
{
	if ((argc < 2) || (argc > 4))
	{
		fprintf(stderr, "Invalid number of arguments. Use:\n  %s <problem_files_prefix> [log_level] [log_file]\n", argv[0]);
		fprintf(stderr, " Example:\n  %s gid_examples/problem.gid/problem\n\n", argv[0]);
		fputs(" Required problem files are:\n", stderr);
		fputs("  <problem_file_prefix>.solver.dat\n", stderr);
		fputs("  <problem_file_prefix>.geometry.dat\n", stderr);
		fputs("  <problem_file_prefix>.problem.dat\n\n", stderr);
		fputs(" By default, log level is set to 2:\n", stderr);
		return 1;
	}

	char solver_file_name[PATH_MAXIMUM_LENGTH];
	char geometry_file_name[PATH_MAXIMUM_LENGTH];
	char problem_file_name[PATH_MAXIMUM_LENGTH];
	char results_file_name[PATH_MAXIMUM_LENGTH];

	sprintf(solver_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.solver.dat", argv[1]);
	sprintf(geometry_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.geometry.dat", argv[1]);
	sprintf(problem_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.problem.dat", argv[1]);

	log_level = (argc > 2) ? atoi(argv[2]) : 2;
	if (argc == 4)
	{
		log_stream = fopen(argv[3], "wb");
		if (!log_stream)
		{
			log_stream = stdout;
			fprintf(stderr, "[Warning] Failed to create log file: %s", argv[3]);
		}
	}

	try
	{
		// Load solver parameters
		SolverParameters<double> solver_parameters(solver_file_name);

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

			Log(1, "Solid.Schur ----------------------------------------------------------");
			Log(1, "-Version:       " MacroValueToString(VERSION));
			Log(1, "-Solver file:   %s", solver_file_name);
			Log(1, "-Geometry file: %s", geometry_file_name);
			Log(1, "-Problem file:  %s", problem_file_name);
			Log(1, "-Partitions:    %i", partitions_count);

			// Load geometry
			Geometry<double> geometry(geometry_file_name);
			geometry.PrintInfo();

			// Partition mesh
			Vector<Substructure> substructures;
			Boundary boundary;
			StructurePartitioning(geometry.mesh, partitions_count, substructures, boundary);
			Log(1, "Partitioning done");

			// Set shape functions
			ShapeFunctions<double> shape_functions(geometry.mesh, geometry.nodes);

			// Load problem parameters
			Solid<double> solid(problem_file_name, shape_functions, geometry.material_index);

			// Assemble and send systems of equations and links
			Log(1, "Systems of equations -------------------------------------------------");
			Log(1, "-Local solver threads: %i", solver_parameters.threads);

			// Assemble global vectors (displacements, forces and fixed)
			int degrees_of_freedom = geometry.nodes.dimension;
			int M = geometry.mesh.nodes_count*degrees_of_freedom;
			Vector<double> global_u(M);
			Vector<double> global_f(M);
			Vector<bool> global_fixed(M);

			// Fill global vectors
			solid.FillU(global_u);
			solid.FillFixed(global_fixed);
			solid.FillF(global_f);

			// Elemental matriz size
			int N = geometry.mesh.nodes_per_element*degrees_of_freedom;

			// Start Schur complement routines
			SchurComplementMaster<double> schur_complement_master(mpi, substructures, boundary, degrees_of_freedom);

			// Systems of equations
			#pragma omp parallel for default(shared) schedule(dynamic,1) num_threads(solver_parameters.substructuring_threads)
			for (int p = 1; p <= partitions_count; ++p)
			{
				const Vector<int>& element_index = substructures.entry[p].element_index;
				const Vector<int>& node_index = substructures.entry[p].node_index;

				// Assembler class
				Assembler assembler(geometry.mesh, element_index, node_index, degrees_of_freedom);

				// Matrices size
				int Mi = node_index.size*degrees_of_freedom; // Matrix size
				int Bi = substructures.entry[p].boundary_count*degrees_of_freedom;
				int Ii = Mi - Bi;

				// Assemble K
				CSRMatrix<double> K(Mi, Mi);
				assembler.AllocateMatrix(K);
				K.Fill(0);
				Matrix<double> Ke(N, N);
				for (int i = 1; i <= element_index.size; ++i)
				{
					int element_id = element_index.entry[i];

					solid.FillKe(element_id, Ke);
					assembler.AssembleAe(element_id, Ke, K);
				}

				// Assemble u
				Vector<double> u(Mi);
				assembler.AssembleV(global_u, u);

				// Assemble f
				Vector<double> f(Mi);
				assembler.AssembleV(global_f, f);

				// Assemble fixed
				Vector<bool> fixed(Mi);
				assembler.AssembleV(global_fixed, fixed);

				// Initialize Schur method for this partition
				schur_complement_master.SendSystemOfEquations(p, K, u, f, fixed);

				Log(1, "Partition %i:", p);
				Log(1, "-Degrees of freedom: %i (%i + %i)", Mi, Ii, Bi);
				Log(1, "-nnz(K):             %i", K.NonZero());
			}

			// Allocate KBB, uB, fB and fixedB
			int B = boundary.node_index.size*degrees_of_freedom;
			CSRMatrix<double> KBB(B, B);
			Vector<double> uB(B);
			Vector<double> fB(B);
			Vector<bool> fixedB(B);

			// Assembler class
			Assembler assemblerB(geometry.mesh, boundary.element_index, boundary.node_index, degrees_of_freedom);

			// Assemble KBB
			assemblerB.AllocateMatrix(KBB);
			KBB.Fill(0.0);
			Matrix<double> Ke(N, N);
			for (int i = 1; i <= boundary.element_index.size; ++i)
			{
				int element_id = boundary.element_index.entry[i];

				solid.FillKe(element_id, Ke);
				assemblerB.AssembleAe(element_id, Ke, KBB);
			}

			// Assemble uB
			assemblerB.AssembleV(global_u, uB);

			// Assemble fB
			assemblerB.AssembleV(global_f, fB);

			// Assemble fixedB
			assemblerB.AssembleV(global_fixed, fixedB);

			Log(1, "Boundaries -----------------------------------------------------------");
			Log(1, "-Degrees of freedom: %i", B);
			Log(1, "-nnz(KBB):           %i", KBB.NonZero());

			Log(1, "Substructuring -------------------------------------------------------");
			Log(1, "-Threads:          %i", solver_parameters.substructuring_threads);
			Log(1, "-Tolerance:        %.5g", solver_parameters.substructuring_tolerance);
			Log(1, "-Multiple results: %s", solver_parameters.multiple_results ? "yes" : "no");

			schur_complement_master.InitializeSolver(KBB, uB, fB, fixedB);
			schur_complement_master.ComplesateFixed(uB, fB, fixedB, solver_parameters.substructuring_threads);
			schur_complement_master.Solver(KBB, uB, fB, solver_parameters.substructuring_tolerance, solver_parameters.max_steps, solver_parameters.substructuring_threads);
			schur_complement_master.CompleteSolution(uB, global_u);

			// Save results
			Log(1, "Save results ---------------------------------------------------------");
			if (solver_parameters.multiple_results)
			{
				char mesh_file_name[PATH_MAXIMUM_LENGTH];
				char mesh_name[31];

				for (int p = 1; p <= partitions_count; ++p)
				{
					Vector<int>& node_index = substructures.entry[p].node_index;
					Vector<int>& element_index = substructures.entry[p].element_index;

					if (solid.save_mesh)
					{
						sprintf(mesh_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.%03i.post.msh", argv[1], p);
						sprintf(mesh_name, "Solid_%03i", p);
						geometry.SaveMesh(mesh_file_name, mesh_name, element_index, node_index, p);
						Log(1, "Mesh saved: %s", mesh_file_name);
					}

					Vector<double> u(node_index.size*degrees_of_freedom);
					for (int i = 1; i <= node_index.size; ++i)
					{
						int l = (i -1)*degrees_of_freedom;
						int g = (node_index.entry[i] - 1)*degrees_of_freedom;
						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							u.entry[l + d] = global_u.entry[g + d];
						}
					}

					// Save results
					double displacement_max = 0;
					double von_mises_max = 0;

					sprintf(results_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.%03i.post.res", argv[1], p);
					File results_file;
					results_file.Create(results_file_name);
					results_file.Write("GiD Post Result File 1.0\n\n");
					solid.WriteStrainStressGaussPoints(results_file);
					if (solid.calculate_displacement)
					{
						solid.WriteDisplacements(results_file, node_index, u, 1, displacement_max);
						Log(1, "Maximum displacement: %0.3g m", displacement_max);
					}
					if (solid.calculate_strain)
					{
						solid.WriteStrain(results_file, element_index, node_index, u, 1);
					}
					if (solid.calculate_stress)
					{
						solid.WriteStress(results_file, element_index, node_index, u, 1);
					}
					if (solid.calculate_von_mises)
					{
						solid.WriteVonMises(results_file, element_index, node_index, u, 1, von_mises_max);
						Log(1, "Maximum Von Mises: %0.3f", von_mises_max);
					}
					results_file.Close();
					Log(1, "Results saved: %s", results_file_name);
				}
			}
			else // Save single file results
			{
				Vector<int> global_element_index(geometry.mesh.elements_count);
				global_element_index.FillSeries(1, 1);

				Vector<int> global_node_index(geometry.mesh.nodes_count);
				global_node_index.FillSeries(1, 1);

				// Save mesh
				if (solid.save_mesh)
				{
					char mesh_file_name[PATH_MAXIMUM_LENGTH];

					sprintf(mesh_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.post.msh", argv[1]);
					geometry.SaveMesh(mesh_file_name, "Solid", global_element_index, global_node_index, 1);
					Log(1, "Mesh saved: %s", mesh_file_name);
				}

				// Save results
				double displacement_max = 0;
				double von_mises_max = 0;

				sprintf(results_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.post.res", argv[1]);
				File results_file;
				results_file.Create(results_file_name);
				results_file.Write("GiD Post Result File 1.0\n\n");
				solid.WriteStrainStressGaussPoints(results_file);
				if (solid.calculate_displacement)
				{
					solid.WriteDisplacements(results_file, global_node_index, global_u, 1, displacement_max);
					Log(1, "Maximum displacement: %0.3g m", displacement_max);
				}
				if (solid.calculate_strain)
				{
					solid.WriteStrain(results_file, global_element_index, global_node_index, global_u, 1);
				}
				if (solid.calculate_stress)
				{
					solid.WriteStress(results_file, global_element_index, global_node_index, global_u, 1);
				}
				if (solid.calculate_von_mises)
				{
					solid.WriteVonMises(results_file, global_element_index, global_node_index, global_u, 1, von_mises_max);
					Log(1, "Maximum Von Mises: %0.3f", von_mises_max);
				}
				results_file.Close();
				Log(1, "Results saved: %s", results_file_name);
			}

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
		}
		else // rank != 0
		{
			log_level = 0; // Prevents slaves logs to appear

			SchurComplementSlave<double> schur_complement_slave(mpi, solver_parameters.threads);
			schur_complement_slave.ReceiveSystemOfEquations();
			schur_complement_slave.Solver();
			schur_complement_slave.CompleteSolution();
			if (Memory::memory_usage)
			{
				schur_complement_slave.ReportMemoryPeakUsage();
			}
		}
	}
	catch (Exception&)
	{
		DebugPosition("Catch fatal exception");
	}

	if (Memory::current_usage != 0)
	{
		fprintf(stderr, "[Error] Memory leak: %lu bytes", (unsigned long)Memory::current_usage);
	}

	if (log_stream != stdout)
	{
		fclose(log_stream);
	}

	return 0;
}
