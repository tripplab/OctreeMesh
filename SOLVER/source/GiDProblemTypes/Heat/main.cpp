// main.cpp
// Copyright (C) 2010 Miguel Vargas (miguel.vargas@gmail.com)
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

#include "Heat.h"

#include <Basic/Log.h>
#include <Basic/Macros.h>
#include <Container/CSRMatrix.h>
#include <Container/Vector.h>
#include <File/MatFile.h>
#include <FiniteElement/Assembler.h>
#include <FiniteElement/Geometry.h>
#include <FiniteElement/ShapeFunctions.h>
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

	char geometry_file_name[PATH_MAXIMUM_LENGTH];
	char problem_file_name[PATH_MAXIMUM_LENGTH];
	char solver_file_name[PATH_MAXIMUM_LENGTH];
	char results_file_name[PATH_MAXIMUM_LENGTH];

	sprintf(solver_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.solver.dat", argv[1]);
	sprintf(geometry_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.geometry.dat", argv[1]);
	sprintf(problem_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.problem.dat", argv[1]);
	sprintf(results_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.post.res", argv[1]);

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
		Log(1, "Heat ----------------------------------------------------------------");
		Log(1, "-Version:       " MacroValueToString(VERSION));
		Log(1, "-Solver file:   %s", solver_file_name);
		Log(1, "-Geometry file: %s", geometry_file_name);
		Log(1, "-Problem file:  %s", problem_file_name);

		// Load solver parameters
		SolverParameters<double> solver_parameters(solver_file_name);
		solver_parameters.PrintInfo();

		// Load geometry
		Geometry<double> geometry(geometry_file_name);
		geometry.PrintInfo();

		// Check geometry
		{
			Sequence<int, 32> error_nodes;
			if (!geometry.mesh.CheckNodes(error_nodes))
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

			{
				Vector<Sequence<int, 32> > elements_group;
				geometry.mesh.CheckElements(elements_group);
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
		}

		// We will use all elements (no partitioning)
		Vector<int> element_index(geometry.mesh.elements_count);
		element_index.FillSeries(1, 1);

		// Get nodes for elements (reorder if necesary)
		Vector<int> node_index(geometry.mesh.nodes_count);
		geometry.mesh.GetElementsNodes(element_index, solver_parameters.reorder, node_index);

		// Set shape functions
		ShapeFunctions<double> shape_functions(geometry.mesh, geometry.nodes);

		// Load problem parameters
		Heat<double> heat(problem_file_name, shape_functions, geometry.material_index);

		// Assembler class
		Assembler assembler(geometry.mesh, element_index, node_index, 1);

		// Matrices size
		int M = node_index.size; // Globlal matrix size
		int N = geometry.mesh.nodes_per_element; // Elemental matriz size

		if (heat.problem_type == time_independent)
		{
			Log(1, "System of equations --------------------------------------------------");

			// Assemble system of equations
			CSRMatrix<double> K(M, M);
			assembler.AllocateMatrix(K);
			K.Fill(0);
			Matrix<double> Ke(N, N);
			for (int i = 1; i <= element_index.size; ++i)
			{
				int element_id = element_index.entry[i];
		
				heat.FillKe(element_id, Ke);
				assembler.AssembleAe(element_id, Ke, K);
			}

			// Assemble displacements and fixed
			Vector<double> u(M);
			Vector<bool> fixed(M);
			{
				Vector<double> global_u(M);
				Vector<bool> global_fixed(M);
				heat.FillU(global_u);
				heat.FillFixed(global_fixed);
				assembler.AssembleV(global_u, u);
				assembler.AssembleV(global_fixed, fixed);
			}

			// Assemble forces
			Vector<double> f(M);
			{
				Vector<double> global_f(M);
				heat.FillF(global_f);
				assembler.AssembleV(global_f, f);
			}

			Log(1, "-Degrees of freedom: %i", M);
			Log(1, "-nnz(K):             %i", K.NonZero());

			// Initialize solver
			Log(1, "Running solver -------------------------------------------------------");
			Solver<double>* solver = solver_parameters.Instanciate(K, u, f, fixed, false);
	
			// Solve system of equations
			solver->CompensateFixed();
			bool valid_solution = solver->Calculate();
			Log(1, "Solution: %s", valid_solution ? "valid" : "INVALID");
	
			// Save system of equations
			if (heat.save_system_of_equations)
			{
				// Generate permutation matrix
				CSRMatrix<int> P(M, M);
				for (int n = 1; n <= node_index.size; ++n)
				{
					P.AllocateRow(n, 1);
					P.entry[n][1] = 1;
					P.index[n][1] = node_index.entry[n];
				}
	
				char mat4_file_name[PATH_MAXIMUM_LENGTH];
				sprintf(mat4_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.post.mat", argv[1]);

				MatFile mat_file;
				mat_file.Create(mat4_file_name);
				mat_file.Store("A", K);
				mat_file.Store("x", u);
				mat_file.Store("b", f);
				mat_file.Store("P", P);
				mat_file.Close();
				Log(1, "System of equations saved: %s", mat4_file_name);
			}
	
			// Save results
			File results_file;
			results_file.Create(results_file_name);
			results_file.Write("GiD Post Result File 1.0\n\n");
			heat.WriteFluxGaussPoints(results_file);
			if (heat.calculate_temperature)
			{
				heat.WriteTemperatures(results_file, node_index, u, 1);
			}
			if (heat.calculate_flux)
			{
				heat.WriteFlux(results_file, element_index, node_index, u, 1);
			}
			results_file.Close();
			Log(1, "Results saved: %s", results_file_name);

			// Free solver resources
			delete solver;
		}
		else // time_dependent
		{
			// R. W. Lewis, P. Nithiarasu, K. N. Seetharamu
			// Fundamentals of the Finite Element Method for Heat and Fluid Flow
			// Wiley, 2004, pp 156-160

			Log(1, "System of equations --------------------------------------------------");

			// Assemble matrices
			CSRMatrix<double> CpK(M, M);
			CSRMatrix<double> CmK(M, M);
			assembler.AllocateMatrix(CpK);
			CmK.CopyStructure(CpK);
			CpK.Fill(0);
			CmK.Fill(0);

			Matrix<double> Ke(N, N);
			Matrix<double> Ce(N, N);
			for (int i = 1; i <= element_index.size; ++i)
			{
				int element_id = element_index.entry[i];
		
				heat.FillKeCe(element_id, Ke, Ce);
				assembler.AssembleAe(element_id, Ke, CpK);
				assembler.AssembleAe(element_id, Ce, CmK);
			}

			for (int i = 1; i <= M; ++i)
			{
				int k_max = CpK.Count(i);
				for (int k = 1; k <= k_max; ++k)
				{
					double K = CpK.entry[i][k];
					double C = CmK.entry[i][k];
					CpK.entry[i][k] = C + heat.time_scheme_factor*heat.time_per_step*K;
					CmK.entry[i][k] = C - (1 - heat.time_scheme_factor)*heat.time_per_step*K;
				}
			}

			// Assemble fixed
			Vector<bool> fixed(M);
			{
				Vector<bool> global_fixed(M);
				heat.FillFixed(global_fixed);
				assembler.AssembleV(global_fixed, fixed);
			}

			Log(1, "-Degrees of freedom: %i", M);
			Log(1, "-nnz(K):             %i", CpK.NonZero());

			Log(1, "Running solver -------------------------------------------------------");

			// Initialize solver
			Vector<double> global_u(M);   // Unordered u
			Vector<double> u(M);          // Reordered u
			Vector<double> global_f(M);   // Unordered f
			Vector<double> f(M);          // Reordered f
			Vector<double> f_previous(M); // Previous f
			Vector<double> b(M);
			Solver<double>* solver = solver_parameters.Instanciate(CpK, u, b, fixed, true);

			// Calculate initial u and f
			heat.time = 0;
			heat.FillU(global_u);
			assembler.AssembleV(global_u, u);
			heat.FillF(global_f);
			assembler.AssembleV(global_f, f_previous);

			// Iterate and save results
			File results_file;
			results_file.Create(results_file_name);
			results_file.Write("GiD Post Result File 1.0\n\n");
			heat.WriteFluxGaussPoints(results_file);
			for (int step = 1; step <= heat.steps; ++step)
			{
				// Set time
				heat.time = heat.time_per_step*step;

				// Assemble f
				heat.FillF(global_f);
				assembler.AssembleV(global_f, f);

				// Calculate b
				for (int i = 1; i <= M; ++i)
				{
					double sum = 0;
					int k_max = CmK.Count(i);
					for (int k = 1; k <= k_max; ++k)
					{
						sum += CmK.entry[i][k]*u.entry[CmK.index[i][k]];
					}
					b.entry[i] = sum + heat.time_per_step*(heat.time_scheme_factor*f.entry[i] + (1 - heat.time_scheme_factor)*f_previous.entry[i]);
					f_previous.entry[i] = f.entry[i];
				}

				// Assemble u
				heat.FillU(global_u);
				assembler.AssembleV(global_u, u);
				solver->CompensateFixed();

				// Solve system of equations
				if (!solver->Calculate())
				{
					Log(1, "Solution: INVALID");
					break;
				}

				if ((step - 1) % heat.result_every_steps == 0)
				{
					if (heat.calculate_temperature)
					{
						heat.WriteTemperatures(results_file, node_index, u, step);
					}
					if (heat.calculate_flux)
					{
						heat.WriteFlux(results_file, element_index, node_index, u, step);
					}
				}
				Log(2, "Step %i", step);
			}
			results_file.Close();
			Log(1, "Results saved: %s", results_file_name);

			// Free solver resources
			delete solver;
		}

		// Save mesh
		if (heat.save_mesh)
		{
			char mesh_file_name[PATH_MAXIMUM_LENGTH];

			sprintf(mesh_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.post.msh", argv[1]);
			geometry.SaveMesh(mesh_file_name, "Heat", element_index, node_index, 1);
			Log(1, "Mesh saved: %s", mesh_file_name);
		}

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
		fprintf(stderr, "[Error] Memory leak: %lu bytes", (unsigned long)Memory::current_usage);
	}

	if (log_stream != stdout)
	{
		fclose(log_stream);
	}

	return 0;
}
