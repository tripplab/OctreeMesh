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

#include "ElectricPotential.h"

#include <Basic/Log.h>
#include <Basic/Macros.h>
#include <Container/CSRMatrix.h>
#include <Container/Vector.h>
#include <File/MatFile.h>
#include <FiniteElement/Assembler.h>
#include <FiniteElement/Geometry.h>
#include <FiniteElement/ShapeFunctions.h>
#include <Solver/SolverParameters.h>


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
		fputs("  <problem_file_prefix>.electric_potential.dat\n\n", stderr);
		fputs(" By default, log level is set to 2:\n", stderr);
		return 1;
	}

	char geometry_file_name[PATH_MAXIMUM_LENGTH];
	char problem_file_name[PATH_MAXIMUM_LENGTH];
	char solver_file_name[PATH_MAXIMUM_LENGTH];
	char results_file_name[PATH_MAXIMUM_LENGTH];
	char capacitances_file_name[PATH_MAXIMUM_LENGTH];

	sprintf(solver_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.solver.dat", argv[1]);
	sprintf(geometry_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.geometry.dat", argv[1]);
	sprintf(problem_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.problem.dat", argv[1]);
	sprintf(results_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.post.res", argv[1]);
	sprintf(capacitances_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.cap.res", argv[1]);

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
		Log(1, "ElectricPotential ----------------------------------------------------------------");
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
		ElectricPotential<double> electric_potential(problem_file_name, shape_functions, geometry.material_index);

		// Assembler class
		Assembler assembler(geometry.mesh, element_index, node_index, 1);

		// Matrices size
		int M = node_index.size; // Globlal matrix size
		int N = geometry.mesh.nodes_per_element; // Elemental matriz size

		switch (electric_potential.problem_type)
		{
			case problem_simple:
			{
				Log(1, "System of equations --------------------------------------------------");

				// Allocate K
				CSRMatrix<double> K(M, M);
				assembler.AllocateMatrix(K);
				K.Fill(0);

				// Assemble elemental matrices
				Matrix<double> Ke(N, N);
				for (int e = 1; e <= element_index.size; ++e)
				{
					int element_id = element_index.entry[e];

					electric_potential.FillKe(element_id, Ke);
					assembler.AssembleAe(element_id, Ke, K);
				}

				// Assemble u
				Vector<double> global_u(M);   // Unordered u
				Vector<double> u(M);          // Reordered u
				electric_potential.FillU(global_u);
				assembler.AssembleV(global_u, u);

				// Assemble f
				Vector<double> global_f(M);
				Vector<double> f(M);
				electric_potential.FillF(global_f);
				assembler.AssembleV(global_f, f);

				// Assemble fixed
				Vector<bool> global_fixed(M);
				Vector<bool> fixed(M);
				electric_potential.FillFixed(global_fixed);
				assembler.AssembleV(global_fixed, fixed);

				Log(1, "-Degrees of freedom: %i", M);
				Log(1, "-nnz(K):             %i", K.NonZero());

				// Run solver
				Log(1, "Running solver -------------------------------------------------------");
				Solver<double>* solver = solver_parameters.Instanciate(K, u, f, fixed, false);

				// Solve system of equations
				solver->CompensateFixed();
				bool valid_solution = solver->Calculate();
				Log(1, "Solution: %s", valid_solution ? "valid" : "INVALID");

				// Save system of equations
				if (electric_potential.save_system_of_equations)
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
					mat_file.Store("fixed", fixed);
					mat_file.Close();
					Log(1, "System of equations saved: %s", mat4_file_name);
				}

				// Save results
				File results_file;
				results_file.Create(results_file_name);
				results_file.Write("GiD Post Result File 1.0\n\n");
				electric_potential.WriteFluxGaussPoints(results_file);
				if (electric_potential.calculate_potential)
				{
					electric_potential.WritePotentials(results_file, node_index, u, 1);
				}
				if (electric_potential.calculate_electric_field)
				{
					electric_potential.WriteElectricField(results_file, element_index, node_index, u, 1);
				}
				results_file.Close();
				Log(1, "Results saved: %s", results_file_name);

				// Free solver resources
				delete solver;
				break;
			}
			case problem_capacitance_matrix:
			{
				Log(1, "System of equations --------------------------------------------------");

				// Allocate K
				CSRMatrix<double> K(M, M);
				assembler.AllocateMatrix(K);
				K.Fill(0);

				// Assemble elemental matrices
				Matrix<double> Ke(N, N);
				for (int e = 1; e <= element_index.size; ++e)
				{
					int element_id = element_index.entry[e];

					electric_potential.FillKe(element_id, Ke);
					assembler.AssembleAe(element_id, Ke, K);
				}

				// Allocate global_u and u
				Vector<double> global_u(M);   // Unordered u
				Vector<double> u(M);          // Reordered u

				// Fill global_f
				Vector<double> global_f(M);
				Vector<double> f(M);
				electric_potential.FillF(global_f);

				// Assemble fixed
				Vector<bool> global_fixed(M);
				Vector<bool> fixed(M);
				electric_potential.FillFixed(global_fixed);
				assembler.AssembleV(global_fixed, fixed);

				Log(1, "-Degrees of freedom: %i", M);
				Log(1, "-nnz(K):             %i", K.NonZero());

				// Keep rows of K to calculate charges
				CSRMatrix<double> Kc(K.rows, K.columns);
				for (int i = 1; i <= K.rows; ++i)
				{
					if (fixed.entry[i])
					{
						Kc.AllocateRow(i, K.Count(i));
						Kc.CopyRowIndexes(i, K);
						Kc.CopyRowValues(i, K);
					}
				}

				// Initialize solver
				Solver<double>* solver = solver_parameters.Instanciate(K, u, f, fixed, true);

				// Inverse node index
				Vector<int> node_index_inverse(node_index.size);
				for (int n = 1; n <= node_index.size; ++n)
				{
					node_index_inverse.entry[node_index.entry[n]] = n;
				}

				// Iterate and save results
				Log(1, "Calculating capacitances ---------------------------------------------");
				File capacitances_file;
				capacitances_file.Create(capacitances_file_name);

				File results_file;
				results_file.Create(results_file_name);
				results_file.Write("GiD Post Result File 1.0\n\n");
				electric_potential.WriteFluxGaussPoints(results_file);
				for (int ea = 1; ea <= electric_potential.electrodes_count; ea += electric_potential.segments_per_step)
				{
					// Assemble u
					electric_potential.electrode_active = ea;
					electric_potential.FillU(global_u);
					assembler.AssembleV(global_u, u);

					// Assemble f
					assembler.AssembleV(global_f, f);

					// Run solver
					solver->CompensateFixed();
					if (!solver->Calculate())
					{
						Log(1, "Solution: INVALID");
						break;
					}

					if (electric_potential.calculate_potential)
					{
						electric_potential.WritePotentials(results_file, node_index, u, ea);
					}
					electric_potential.WriteCapacitance(capacitances_file, node_index_inverse, Kc, u);
					Log(2, "Active electrode: %i", ea);
				}
				results_file.Close();
				Log(1, "Results saved: %s", results_file_name);
				capacitances_file.Close();
				Log(1, "Capacitance matrix saved: %s", capacitances_file_name);

				// Free solver resources
				delete solver;
				break;
			}
			case problem_sensitivity_analysis:
			{
				// Keep old log_level
				int old_log_level = log_level;

				// Fill global_f
				Vector<double> global_f(M);
				electric_potential.FillF(global_f);

				// Assemble fixed
				Vector<bool> global_fixed(M);
				Vector<bool> fixed(M);
				electric_potential.FillFixed(global_fixed);
				assembler.AssembleV(global_fixed, fixed);

				Log(1, "Calculating sensitivity analysis -------------------------------------");

				// Inverse node index
				Vector<int> node_index_inverse(node_index.size);
				for (int n = 1; n <= node_index.size; ++n)
				{
					node_index_inverse.entry[node_index.entry[n]] = n;
				}

				// Sensitivity analysis material
				int sensitivity_analysis_material_id = 0;
				for (int m = 1; m <= electric_potential.materials.size; ++m)
				{
					if (electric_potential.materials.entry[m].sensitivity_use)
					{
						sensitivity_analysis_material_id = m;
						break;
					}
				}

				// Sensitivity analisys pertmitivitties
				int high_permittivity_material_id = electric_potential.materials.size - 1;
				int low_permittivity_material_id = electric_potential.materials.size;
				double high_permittivity = electric_potential.materials.entry[high_permittivity_material_id].permittivity;
				double low_permittivity = electric_potential.materials.entry[low_permittivity_material_id].permittivity;

				Log(1, "Calculating capacitance_high");
				Matrix<double> capacitance_high(electric_potential.electrodes_count, electric_potential.electrodes_count);
				{
					log_level = 0;

					// Allocate K
					CSRMatrix<double> K(M, M);
					assembler.AllocateMatrix(K);
					K.Fill(0);

					// Assemble elemental matrices
					electric_potential.materials.entry[sensitivity_analysis_material_id].permittivity = high_permittivity;
					Matrix<double> Ke(N, N);
					for (int e = 1; e <= element_index.size; ++e)
					{
						int element_id = element_index.entry[e];

						electric_potential.FillKe(element_id, Ke);
						assembler.AssembleAe(element_id, Ke, K);
					}

					// Keep rows of K to calculate charges
					CSRMatrix<double> Kc(K.rows, K.columns);
					for (int i = 1; i <= K.rows; ++i)
					{
						if (fixed.entry[i])
						{
							Kc.AllocateRow(i, K.Count(i));
							Kc.CopyRowIndexes(i, K);
							Kc.CopyRowValues(i, K);
						}
					}

					// Allocate global_u, u and f
					Vector<double> global_u(M);
					Vector<double> u(M);
					Vector<double> f(M);

					Solver<double>* solver = solver_parameters.Instanciate(K, u, f, fixed, true);

					Vector<double> capacitance(electric_potential.electrodes_count);
					for (int ea = 1; ea <= electric_potential.electrodes_count; ++ea)
					{
						// Assemble u
						electric_potential.electrode_active = ea;
						electric_potential.FillU(global_u);
						assembler.AssembleV(global_u, u);

						// Assemble f
						assembler.AssembleV(global_f, f);

						// Run solver
						solver->CompensateFixed();
						if (!solver->Calculate())
						{
							log_level = old_log_level;
							Log(1, "Solution: INVALID");
							break;
						}

						// Fill capacitance with high permitivitty
						electric_potential.GetCapacitance(node_index_inverse, Kc, u, capacitance);
						for (int et = 1; et <= electric_potential.electrodes_count; ++et)
						{
							capacitance_high.entry[ea][et] = capacitance.entry[et];
						}
					}

					log_level = old_log_level;

					// Free solver resources
					delete solver;
				}

				Log(1, "Calculating capacitance_low");
				Matrix<double> capacitance_low(electric_potential.electrodes_count, electric_potential.electrodes_count);
				{
					log_level = 0;

					// Allocate K
					CSRMatrix<double> K(M, M);
					assembler.AllocateMatrix(K);
					K.Fill(0);

					// Assemble elemental matrices
					electric_potential.materials.entry[sensitivity_analysis_material_id].permittivity = low_permittivity;
					Matrix<double> Ke(N, N);
					for (int e = 1; e <= element_index.size; ++e)
					{
						int element_id = element_index.entry[e];

						electric_potential.FillKe(element_id, Ke);
						assembler.AssembleAe(element_id, Ke, K);
					}

					// Keep rows of K to calculate charges
					CSRMatrix<double> Kc(K.rows, K.columns);
					for (int i = 1; i <= K.rows; ++i)
					{
						if (fixed.entry[i])
						{
							Kc.AllocateRow(i, K.Count(i));
							Kc.CopyRowIndexes(i, K);
							Kc.CopyRowValues(i, K);
						}
					}

					// Allocate global_u, u and f
					Vector<double> global_u(M);
					Vector<double> u(M);
					Vector<double> f(M);

					Solver<double>* solver = solver_parameters.Instanciate(K, u, f, fixed, true);

					Vector<double> capacitance(electric_potential.electrodes_count);
					for (int ea = 1; ea <= electric_potential.electrodes_count; ++ea)
					{
						// Assemble u
						electric_potential.electrode_active = ea;
						electric_potential.FillU(global_u);
						assembler.AssembleV(global_u, u);

						// Assemble f
						assembler.AssembleV(global_f, f);

						// Run solver
						solver->CompensateFixed();
						if (!solver->Calculate())
						{
							log_level = old_log_level;
							Log(1, "Solution: INVALID");
							break;
						}

						// Fill capacitance with low permitivitty
						electric_potential.GetCapacitance(node_index_inverse, Kc, u, capacitance);
						for (int et = 1; et <= electric_potential.electrodes_count; ++et)
						{
							capacitance_low.entry[ea][et] = capacitance.entry[et];
						}
					}

					log_level = old_log_level;

					// Free solver resources
					delete solver;
				}

				Log(1, "Calculating sensitivity");
				Vector<Matrix<double> > sensitivity(geometry.mesh.elements_count);
				{
					// Set solver_threads to 1 and parallelize sensistivity calculation
					int threads = solver_parameters.threads;
					solver_parameters.threads = 1;
					omp_set_num_threads(1); // Set number of threads for matrix operations
					log_level = 0;

					int tests_count = 0;

					// Calculate element_size and element_size_max
					Vector<double> element_size(element_index.size);
					double element_size_max = 0;
					for (int e = 1; e <= element_index.size; ++e)
					{
						int element_id = element_index.entry[e];
						if (electric_potential.material_index.entry[element_id] == sensitivity_analysis_material_id)
						{
							double sum = 0;
							for (int q = 1; q <= electric_potential.element_integration_rules.size; ++q)
							{
								double det_J = shape_functions.ElementDetJ(element_id, electric_potential.element_integration_rules.entry[q].point);

								sum += electric_potential.element_integration_rules.entry[q].weight*det_J;
							}
							element_size.entry[element_id] = sum;
							if (element_size_max < sum)
							{
								element_size_max = sum;
							}
						}
						else
						{
							element_size.entry[element_id] = 0;
						}
					}

					#pragma omp parallel for schedule(dynamic,1) firstprivate(electric_potential) num_threads(threads)
					for (int e = 1; e <= element_index.size; ++e)
					{
						int element_id = element_index.entry[e];
						if (electric_potential.material_index.entry[element_id] == sensitivity_analysis_material_id)
						{
							// Allocate K
							CSRMatrix<double> K(M, M);
							assembler.AllocateMatrix(K);
							K.Fill(0);

							// Assemble elemental matrices
							electric_potential.material_index.entry[element_id] = high_permittivity_material_id;
							Matrix<double> Ke(N, N);
							for (int e = 1; e <= element_index.size; ++e)
							{
								int element_id = element_index.entry[e];

								electric_potential.FillKe(element_id, Ke);
								assembler.AssembleAe(element_id, Ke, K);
							}
							electric_potential.material_index.entry[element_id] = sensitivity_analysis_material_id;

							// Keep rows of K to calculate charges
							CSRMatrix<double> Kc(K.rows, K.columns);
							for (int i = 1; i <= K.rows; ++i)
							{
								if (fixed.entry[i])
								{
									Kc.AllocateRow(i, K.Count(i));
									Kc.CopyRowIndexes(i, K);
									Kc.CopyRowValues(i, K);
								}
							}

							// Allocate global_u, u and f
							Vector<double> global_u(M);
							Vector<double> u(M);
							Vector<double> f(M);

							Solver<double>* solver = solver_parameters.Instanciate(K, u, f, fixed, true);

							Vector<double> capacitance(electric_potential.electrodes_count);
							sensitivity.entry[element_id].Resize(electric_potential.electrodes_count, electric_potential.electrodes_count);
							for (int ea = 1; ea <= electric_potential.electrodes_count; ++ea)
							{
								// Assemble u
								electric_potential.electrode_active = ea;
								electric_potential.FillU(global_u);
								assembler.AssembleV(global_u, u);

								// Assemble f
								assembler.AssembleV(global_f, f);

								// Run solver
								solver->CompensateFixed();
								if (!solver->Calculate())
								{
									Log(0, "Solution: INVALID");
									break;
								}

								// Fill sensitivity
								electric_potential.GetCapacitance(node_index_inverse, Kc, u, capacitance);
								for (int et = 1; et <= electric_potential.electrodes_count; ++et)
								{
									sensitivity.entry[element_id].entry[ea][et] = (capacitance.entry[et] - capacitance_low.entry[ea][et])/(capacitance_high.entry[ea][et] - capacitance_low.entry[ea][et])/(high_permittivity - low_permittivity)*element_size_max/element_size.entry[element_id];
								}
							}

							// Free solver resources
							delete solver;

							#pragma omp critical
							{
								++tests_count;
								if (old_log_level >= 2)
								{
									if (!(tests_count % 100))
									{
										log_level = old_log_level;
										Log(2, "-Tests: %i", tests_count);
										log_level = 0;
									}
								}
							}
						}
					}

					log_level = old_log_level;
				}

				File results_file;
				results_file.Create(results_file_name);
				results_file.Write("GiD Post Result File 1.0\n\n");

				// Write sensitivity on elements
				electric_potential.WriteSensitivityGaussPoints(results_file);
				for (int ea = 1; ea <= electric_potential.electrodes_count; ++ea)
				{
					for (int et = 1; et <= electric_potential.electrodes_count; ++et)
					{
						electric_potential.WriteSensitivity(results_file, element_index, sensitivity, ea, et);
					}
				}

				// Write sensitivity on nodes
				if (electric_potential.result_on_nodes)
				{
					Mesh& mesh = geometry.mesh;
					Vector<int> count(node_index.size);
					count.Fill(0);
					for (int e = 1; e <= element_index.size; ++e)
					{
						int element_id = element_index.entry[e];
						if (electric_potential.material_index.entry[element_id] == sensitivity_analysis_material_id)
						{
							for (int c = 1; c <= mesh.nodes_per_element; ++c)
							{
								int node_id = mesh.connectivity.entry[element_id][c];
								++count.entry[node_id];
							}
						}
					}
					Vector<double> value(node_index.size);
					for (int ea = 1; ea <= electric_potential.electrodes_count; ++ea)
					{
						for (int et = 1; et <= electric_potential.electrodes_count; ++et)
						{
							value.Fill(0);
							for (int e = 1; e <= element_index.size; ++e)
							{
								int element_id = element_index.entry[e];
								if (electric_potential.material_index.entry[element_id] == sensitivity_analysis_material_id)
								{
									for (int c = 1; c <= mesh.nodes_per_element; ++c)
									{
										int node_id = mesh.connectivity.entry[element_id][c];
										value.entry[node_id] += sensitivity.entry[element_id].entry[ea][et];
									}
								}
							}
							electric_potential.WriteSensitivityOnNodes(results_file, value, count, ea, et);
						}
					}
				}
				results_file.Close();
				Log(1, "Results saved: %s", results_file_name);

				break;
			}
		}

		// Save mesh
		if (electric_potential.save_mesh)
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
