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

#include "Solid.h"

#include <Basic/Log.h>
#include <Basic/Macros.h>
#include <Container/CSRMatrix.h>
#include <Container/Vector.h>
#include <File/MatFile.h>
#include <FiniteElement/Assembler.h>
#include <FiniteElement/Geometry.h>
#include <FiniteElement/ShapeFunctions.h>
#include <Solver/SolverParameters.h>

#include <stdlib.h>


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
		Log(1, "Solid ----------------------------------------------------------------");
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

		// We will use all elements (no partitioning)
		Vector<int> element_index(geometry.mesh.elements_count);
		element_index.FillSeries(1, 1);

		// Get nodes for elements (reorder if necesary)
		Vector<int> node_index(geometry.mesh.nodes_count);
		geometry.mesh.GetElementsNodes(element_index, solver_parameters.reorder, node_index);

		// Set shape functions
		ShapeFunctions<double> shape_functions(geometry.mesh, geometry.nodes);

		// Load problem parameters
		Solid<double> solid(problem_file_name, shape_functions, geometry.material_index);

		// Assembler class
		Assembler assembler(geometry.mesh, element_index, node_index, geometry.nodes.dimension);

		// Matrices size
		int M = node_index.size*geometry.nodes.dimension; // Globlal matrix size
		int N = geometry.mesh.nodes_per_element*geometry.nodes.dimension; // Elemental matriz size

		if (solid.problem_type == stationary_problem)
		{
			Log(1, "System of equations --------------------------------------------------");

			// Assemble stiffness matrix
			CSRMatrix<double> K(M, M);
			assembler.AllocateMatrix(K);
			K.Fill(0);
			Matrix<double> Ke(N, N);
			for (int i = 1; i <= element_index.size; ++i)
			{
				int element_id = element_index.entry[i];
		
				solid.FillKe(element_id, Ke);
				assembler.AssembleAe(element_id, Ke, K);
			}

			// Assemble displacements and fixed
			Vector<double> u(M);
			Vector<bool> fixed(M);
			{
				Vector<double> global_u(M);
				Vector<bool> global_fixed(M);
				solid.FillU(global_u);
				solid.FillFixed(global_fixed);
				assembler.AssembleV(global_u, u);
				assembler.AssembleV(global_fixed, fixed);
			}

			// Assemble forces
			Vector<double> f(M);
			{
				Vector<double> global_f(M);
				solid.FillF(global_f);
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
			if (solid.save_system_of_equations)
			{
				// Generate permutation matrix
				CSRMatrix<int> P(M, M);
				for (int n = 1, i = 1; n <= node_index.size; ++n)
				{
					for (int d = 1; d <= geometry.nodes.dimension; ++d, ++i)
					{
						P.AllocateRow(i, 1);
						P.entry[i][1] = 1;
						P.index[i][1] = (node_index.entry[n] - 1)*geometry.nodes.dimension + d;
					}
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
			double displacement_max = 0;
			double von_mises_max = 0;

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
			Log(1, "Total weight: %0.3f kg", solid.TotalWeight());
			Log(1, "Results saved: %s", results_file_name);

			// Free solver resources
			delete solver;
		}
		else // dynamic_problem
		{
			Log(1, "System of equations --------------------------------------------------");

			double alpha = solid.time_scheme_factor;
			double beta = (1 - alpha)*(1 - alpha)/4;
			double gamma = (1 - 2*alpha)/2;
			double h = solid.time_per_step;

			// Assemble matrices
			CSRMatrix<double> P(M, M);
			CSRMatrix<double> Q(M, M);
			CSRMatrix<double> R(M, M);
			CSRMatrix<double> S(M, M);
			assembler.AllocateMatrix(P);
			Q.CopyStructure(P);
			R.CopyStructure(P);
			S.CopyStructure(P);
			P.Fill(0);
			Q.Fill(0);
			R.Fill(0);

			Matrix<double> Ke(N, N);
			Matrix<double> Ce(N, N);
			Matrix<double> Me(N, N);
			for (int i = 1; i <= element_index.size; ++i)
			{
				int element_id = element_index.entry[i];
		
				solid.FillKeCeMe(element_id, Ke, Ce, Me);
				assembler.AssembleAe(element_id, Ke, P);
				assembler.AssembleAe(element_id, Ce, Q);
				assembler.AssembleAe(element_id, Me, R);
			}

			for (int i = 1; i <= M; ++i)
			{
				int k_max = P.Count(i);
				for (int k = 1; k <= k_max; ++k)
				{
					double K = P.entry[i][k];
					double C = Q.entry[i][k];
					double M = R.entry[i][k];
					P.entry[i][k] = M/(h*h*beta) + (1 + alpha)*C*gamma/(h*beta) + (1 + alpha)*K; // *dn
					Q.entry[i][k] = M*(1/(2*beta) - 1) - (1 + alpha)*C*h*(1 - gamma/(2*beta)); // *a
					R.entry[i][k] = M/(h*beta) - C*((1 + alpha)*(1 - gamma/beta) + alpha); // *v
					S.entry[i][k] = M/(h*h*beta) + (1 + alpha)*C*gamma/(h*beta) + alpha*K; // *d
				}
			}

			// Assemble fixed
			Vector<bool> fixed(M);
			{
				Vector<bool> global_fixed(M);
				solid.FillFixed(global_fixed);
				assembler.AssembleV(global_fixed, fixed);
			}

			Log(1, "-Degrees of freedom: %i", M);
			Log(1, "-nnz(K):             %i", P.NonZero());

			Log(1, "Running solver -------------------------------------------------------");

			// Initialize solver
			Vector<double> global_u(M); // Unordered u
			Vector<double> u(M);        // Reordered u
			Vector<double> global_f(M); // Unordered f
			Vector<double> f(M);        // Reordered f
			Vector<double> u_prev(M);   // Previous u
			Vector<double> b(M);
			Solver<double>* solver = solver_parameters.Instanciate(P, u, b, fixed, true);

			// Calculate initial u
			solid.time = 0;
			solid.FillU(global_u);
			assembler.AssembleV(global_u, u);

			// Velocity and acceleration
			Vector<double> a(M);
			Vector<double> v(M);
			a.Fill(0);
			v.Fill(0);

			// Iterate and save results
			File results_file;
			results_file.Create(results_file_name);
			results_file.Write("GiD Post Result File 1.0\n\n");
			solid.WriteStrainStressGaussPoints(results_file);
			for (int step = 1; step <= solid.steps; ++step)
			{
				// Assemble f
				solid.time = solid.time_per_step*step + (1 + alpha)*solid.time_per_step;
				solid.FillF(global_f);
				assembler.AssembleV(global_f, f);

				// Calculate b
				for (int i = 1; i <= M; ++i)
				{
					double sum = f.entry[i];
					int k_max = S.Count(i);
					for (int k = 1; k <= k_max; ++k)
					{
						int j = S.index[i][k];

						sum += Q.entry[i][k]*a.entry[j];
						sum += R.entry[i][k]*v.entry[j];
						sum += S.entry[i][k]*u.entry[j];
					}
					b.entry[i] = sum;
					u_prev.entry[i] = u.entry[i];
				}

				// Assemble u
				solid.time = solid.time_per_step*step;
				solid.FillU(global_u);
				assembler.AssembleV(global_u, u);
				solver->CompensateFixed();

				// Solve system of equations
				if (!solver->Calculate())
				{
					Log(1, "Solution: INVALID");
					break;
				}

				// Update a, v
				for (int i = 1; i <= M; ++i)
				{
					double V = v.entry[i];
					double A = a.entry[i];
					double UN = u.entry[i];
					double U = u_prev.entry[i];
					a.entry[i] = UN/(h*h*beta) - U/(h*h*beta) - V/(h*beta) - A*(1/(2*beta) - 1);
					v.entry[i] = V + h*(1 - gamma)*A + h*gamma*a.entry[i];
				}

				Log(2, "Step %i", step);
				if ((step - 1) % solid.result_every_steps == 0)
				{
					if (solid.calculate_displacement)
					{
						double maximum;
						solid.WriteDisplacements(results_file, node_index, u, step, maximum);
						Log(1, "Maximum displacement: %0.2f", maximum);
					}
					if (solid.calculate_strain)
					{
						solid.WriteStrain(results_file, element_index, node_index, u, step);
					}
					if (solid.calculate_stress)
					{
						solid.WriteStress(results_file, element_index, node_index, u, step);
					}
					if (solid.calculate_von_mises)
					{
						double maximum;
						solid.WriteVonMises(results_file, element_index, node_index, u, step, maximum);
						Log(1, "Maximum Von Mises: %0.2f", maximum);
					}
				}
			}
			results_file.Close();
			Log(1, "Results saved: %s", results_file_name);

			// Free solver resources
			delete solver;
		}

		// Save mesh
		if (solid.save_mesh)
		{
			char mesh_file_name[PATH_MAXIMUM_LENGTH];

			sprintf(mesh_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.post.msh", argv[1]);
			geometry.SaveMesh(mesh_file_name, "Solid", element_index, node_index, 1);
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
