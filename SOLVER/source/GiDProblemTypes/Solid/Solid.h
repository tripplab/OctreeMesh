// Solid.h
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

#ifndef _Solid_h_
#define _Solid_h_

#include <Basic/Float.h>
#include <Basic/Integer.h>
#include <Basic/Log.h>
#include <Basic/System.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>
#include <FiniteElement/Shape.h>
#include <FiniteElement/ShapeFunctions.h>
#include <FiniteElement/ShapeIntegrationRule.h>
#include <File/File.h>
#include <Math/Formula.h>

#include <string.h>


#define SOLID_SECTION_NAME_MAX_SIZE 256
#define SOLID_USER_FUNCTIONS_COUNT 20
#define SOLID_USER_FUNCTION_MAX_SIZE 65536


enum PlaneProblem
{
	plane_stress = 1,
	plane_strain = 2
};


enum ProblemType
{
	stationary_problem = 1,
	dynamic_problem = 2
};


template <typename T>
class Solid
{
	public:

		struct Material
		{
			T poisson_ratio;
			T young_modulus;
			T density;
			T thickness;
		};

		struct Displacement
		{
			int node_id;
			int user_function_x_id;
			T displacement_x;
			bool fixed_x;
			int user_function_y_id;
			T displacement_y;
			bool fixed_y;
			int user_function_z_id;
			T displacement_z;
			bool fixed_z;
		};

		struct NodalForce
		{
			int node_id;
			int user_function_x_id;
			T force_x;
			int user_function_y_id;
			T force_y;
			int user_function_z_id;
			T force_z;
		};

		struct NormalForce
		{
			int element_id;
			int user_function_id;
			T force;
			Vector<int> facet_nodes;
		};

		// General
		PlaneProblem plane_problem;
		bool calculate_displacement;
		bool calculate_strain;
		bool calculate_stress;
		bool calculate_von_mises;
		ProblemType problem_type;
		bool save_mesh;
		bool save_system_of_equations;

		// Dynamic
		T time_per_step;
		int steps;
		int result_every_steps;
		T time_scheme_factor;
		T rayleigh_damping_a;
		T rayleigh_damping_b;

		// Gravity
		bool use_mass_forces;
		T gravity;

		// User functions
		T x;
		T y;
		T z;
		T time;
		Vector<Formula<T> > user_function;

		Vector<Material> materials;
		Vector<Displacement> displacements;
		Vector<NodalForce> nodal_forces;
		Vector<NormalForce> normal_forces;
		Vector<IntegrationRule<T> > element_integration_rules;
		Vector<IntegrationRule<T> > facet_integration_rules;

		const ShapeFunctions<T>& shape_functions;
		const Vector<int>& material_index;


		Solid(const char* file_name, const ShapeFunctions<T>& shape_functions, const Vector<int>& material_index) throw(Memory::Exception, File::Exception)
		:	plane_problem(),
			calculate_displacement(),
			calculate_strain(),
			calculate_stress(),
			calculate_von_mises(),
			problem_type(),
			save_mesh(),
			save_system_of_equations(),
			time_per_step(),
			steps(),
			result_every_steps(),
			time_scheme_factor(),
			rayleigh_damping_a(),
			rayleigh_damping_b(),
			use_mass_forces(),
			gravity(),
			x(),
			y(),
			z(),
			time(0),
			user_function(SOLID_USER_FUNCTIONS_COUNT),
			materials(),
			displacements(),
			nodal_forces(),
			normal_forces(),
			element_integration_rules(),
			facet_integration_rules(),
			shape_functions(shape_functions),
			material_index(material_index)
		{
			try
			{
				char section_name[SOLID_SECTION_NAME_MAX_SIZE];
				int tmp_integer;

				// Load problem
				File file;
				file.Open(file_name);
				file.SkipComments(';');

				// Load general parameters
				file.Read(section_name, SOLID_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{General}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(tmp_integer);
				plane_problem = (PlaneProblem)tmp_integer;
				file.SkipComments(';');
				file.Read(calculate_displacement);
				file.SkipComments(';');
				file.Read(calculate_strain);
				file.SkipComments(';');
				file.Read(calculate_stress);
				file.SkipComments(';');
				file.Read(calculate_von_mises);
				file.SkipComments(';');
				file.Read(tmp_integer);
				problem_type = (ProblemType)tmp_integer;
				file.SkipComments(';');
				file.Read(save_mesh);
				file.SkipComments(';');
				file.Read(save_system_of_equations);
				file.SkipComments(';');

				// Load dynamic problem parameters
				file.Read(section_name, SOLID_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{DynamicParameters}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(time_per_step);
				file.SkipComments(';');
				file.Read(steps);
				file.SkipComments(';');
				file.Read(result_every_steps);
				file.SkipComments(';');
				file.Read(time_scheme_factor);
				file.SkipComments(';');
				file.Read(rayleigh_damping_a);
				file.SkipComments(';');
				file.Read(rayleigh_damping_b);
				file.SkipComments(';');

				// Load gravity parameters
				file.Read(section_name, SOLID_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Gravity}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(use_mass_forces);
				file.SkipComments(';');
				file.Read(gravity);
				file.SkipComments(';');

				// Load user defined functions
				file.Read(section_name, SOLID_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{UserFunctions}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				for (int i = 1; i <= SOLID_USER_FUNCTIONS_COUNT; ++i)
				{
					char function_text[SOLID_USER_FUNCTION_MAX_SIZE];
					file.ReadLine(function_text, SOLID_USER_FUNCTION_MAX_SIZE - 10);
					strncat(function_text, ";x;y;z;t", 10);
					try
					{
						user_function.entry[i].Define(function_text);
					}
					catch (Memory::Exception&)
					{
						ReThrow();
					}
					catch (Exception&)
					{
						Throw(File::exception_format);
					}
				}

				// Load materials
				int materials_count;
				file.Read(section_name, SOLID_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Materials}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(materials_count);
				file.SkipComments(';');
				materials.Resize(materials_count);
				for (int m = 1; m <= materials_count; ++m)
				{
					file.Read(materials.entry[m].poisson_ratio);
					file.Read(materials.entry[m].young_modulus);
					file.Read(materials.entry[m].density);
					file.Read(materials.entry[m].thickness);
					file.SkipComments(';');
				}

				// Load displacements conditions
				int displacement_count;
				file.Read(section_name, SOLID_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Displacement}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(displacement_count);
				file.SkipComments(';');
				displacements.Resize(displacement_count);
				for (int c = 1; c <= displacement_count; ++c)
				{
					displacements.entry[c].fixed_x = false;
					displacements.entry[c].fixed_y = false;
					displacements.entry[c].fixed_z = false;
				}
				for (int c = 1; c <= displacement_count; ++c)
				{
					int fixed_x;
					int fixed_y;
					int fixed_z;

					file.Read(displacements.entry[c].node_id);
					file.Read(displacements.entry[c].user_function_x_id);
					file.Read(displacements.entry[c].displacement_x);
					file.Read(fixed_x);
					displacements.entry[c].fixed_x |= fixed_x == 1;
					file.Read(displacements.entry[c].user_function_y_id);
					file.Read(displacements.entry[c].displacement_y);
					file.Read(fixed_y);
					displacements.entry[c].fixed_y |= fixed_y == 1;
					file.Read(displacements.entry[c].user_function_z_id);
					file.Read(displacements.entry[c].displacement_z);
					file.Read(fixed_z);
					displacements.entry[c].fixed_z |= fixed_z == 1;
					file.SkipComments(';');
				}

				// Load nodal force conditions
				int nodal_force_count;
				file.Read(section_name, SOLID_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{NodalForce}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(nodal_force_count);
				file.SkipComments(';');
				nodal_forces.Resize(nodal_force_count);
				for (int c = 1; c <= nodal_force_count; ++c)
				{
					file.Read(nodal_forces.entry[c].node_id);
					file.Read(nodal_forces.entry[c].user_function_x_id);
					file.Read(nodal_forces.entry[c].force_x);
					file.Read(nodal_forces.entry[c].user_function_y_id);
					file.Read(nodal_forces.entry[c].force_y);
					file.Read(nodal_forces.entry[c].user_function_z_id);
					file.Read(nodal_forces.entry[c].force_z);
					file.SkipComments(';');
				}

				// Load normal force conditions
				int normal_force_count;
				file.Read(section_name, SOLID_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{NormalForce}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(normal_force_count);
				file.SkipComments(';');
				normal_forces.Resize(normal_force_count);
				for (int c = 1; c <= normal_force_count; ++c)
				{
					file.Read(normal_forces.entry[c].element_id);
					file.Read(normal_forces.entry[c].user_function_id);
					file.Read(normal_forces.entry[c].force);
					normal_forces.entry[c].facet_nodes.Resize(shape_functions.mesh.nodes_per_facet);
					for (int f = 1; f <= shape_functions.mesh.nodes_per_facet; ++f)
					{
						file.Read(normal_forces.entry[c].facet_nodes.entry[f]);
					}
					file.SkipComments(';');
				}

				file.Close();

				// Fill integration rules
				int element_integration_points_count;
				int facet_integration_points_count;
				switch (shape_functions.mesh.element_type)
				{
					case shape_triangle:
					{
						element_integration_points_count = (shape_functions.mesh.nodes_per_element == 3) ? 1 : 3;
						facet_integration_points_count = (shape_functions.mesh.nodes_per_facet == 2) ? 1 : 2;
						break;
					}
					case shape_quadrilateral:
					{
						element_integration_points_count = (shape_functions.mesh.nodes_per_element == 4) ? 2 : (shape_functions.mesh.nodes_per_element == 8) ? 2 : 3;
						facet_integration_points_count = (shape_functions.mesh.nodes_per_facet == 2) ? 1 : 2;
						break;
					}
					case shape_tetrahedron:
					{
						element_integration_points_count = (shape_functions.mesh.nodes_per_element == 4) ? 1 : 4;
						facet_integration_points_count = (shape_functions.mesh.nodes_per_facet == 3) ? 1 : 3;
						break;
					}
					case shape_hexahedron:
					{
						element_integration_points_count = (shape_functions.mesh.nodes_per_element == 8) ? 2 : (shape_functions.mesh.nodes_per_element == 20) ? 3 : 3;
						facet_integration_points_count = (shape_functions.mesh.nodes_per_facet == 4) ? 2 : (shape_functions.mesh.nodes_per_facet == 8) ? 2 : 3;
						break;
					}
					default:
					{
						Throw(File::exception_format);
					}
				}
				ShapeIntegrationRule(shape_functions.mesh.element_type, element_integration_points_count, element_integration_rules);
				ShapeIntegrationRule(shape_functions.mesh.facet_type, facet_integration_points_count, facet_integration_rules);

				Log(1, "Problem --------------------------------------------------------------");
				if (shape_functions.nodes.dimension == 2)
				{
					Log(1, "-Plane problem :           %s", plane_problem == plane_stress ? "Plane stress" : "Plane strain");
				}
				Log(1, "-Save mesh:                %s", save_mesh ? "yes" : "no");
				Log(1, "-Save system of equations: %s", save_system_of_equations ? "yes" : "no");
				Log(1, "-Materials used:           %i", materials.size);
				Log(1, "-Displacement conditions:  %i", displacements.size);
				Log(1, "-Nodal force conditions:   %i", nodal_forces.size);
				Log(1, "-Normal force conditions:  %i", normal_forces.size);
				if (use_mass_forces)
				{
					Log(1, "Mass forces:");
					Log(1, "-Gravity applied:          %g", gravity);
				}
				if (problem_type == stationary_problem)
				{
					Log(1, "-Problem type:             Stationary");
				}
				else
				{
					Log(1, "-Problem type:             Dynamic");
					Log(1, "Dynamic parameters:");
					Log(1, "-Time per step:            %g", time_per_step);
					Log(1, "-Steps:                    %i", steps);
					Log(1, "-Result every steps:       %i", result_every_steps);
					Log(1, "-Time alpha factor:        %g", time_scheme_factor);
					Log(1, "-Rayleigh damping a        %g", rayleigh_damping_a);
					Log(1, "-Rayleigh damping b        %g", rayleigh_damping_b);
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void FillKe(const int element_id, Matrix<T>& Ke) const throw(Memory::Exception)
		{
			// S. W. Sloan
			// A fast stiffness formulation for finite element analysis of two-dimensional solids
			// International Journal for Numerical Methods in Engineering, Vol. 17, pp. 1313-1323, 1981

			try
			{
				int m = shape_functions.mesh.nodes_per_element*shape_functions.nodes.dimension;

				Vector<T> N(shape_functions.mesh.nodes_per_element);
				Matrix<T> dN(shape_functions.mesh.nodes_per_element, shape_functions.nodes.dimension);

				int material_id = material_index.entry[element_id];
				Material& material = materials.entry[material_id];

				if (shape_functions.nodes.dimension == 2)
				{
					Matrix<T> Bt(shape_functions.mesh.nodes_per_element*2, 3);
					Matrix<T> DB(3, shape_functions.mesh.nodes_per_element*2);
					Ke.Fill(T(0));

					T E = material.young_modulus;
					T v = material.poisson_ratio;

					if (plane_problem == plane_stress)
					{
						T a = E/(1 - v*v);
						T b = a*v;
						T c = a*(1 - v)/2;
						T t = material.thickness;

						for (int q = 1; q <= element_integration_rules.size; ++q)
						{
							T det_J;
							shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

							for (register int i = 1; i <= shape_functions.mesh.nodes_per_element; ++i)
							{
								register int i2 = i*2;
								register int i1 = i2 - 1;
								T p = dN.entry[i][1];
								T q = dN.entry[i][2];
								Bt.entry[i1][1] = p; Bt.entry[i2][1] = 0;
								Bt.entry[i1][2] = 0; Bt.entry[i2][2] = q;
								Bt.entry[i1][3] = q; Bt.entry[i2][3] = p;
								DB.entry[1][i1] = a*p; DB.entry[1][i2] = b*q;
								DB.entry[2][i1] = b*p; DB.entry[2][i2] = a*q;
								DB.entry[3][i1] = c*q; DB.entry[3][i2] = c*p;
							}

							T weight = element_integration_rules.entry[q].weight;
							for (register int i = 1; i <= m; ++i)
							{
								for (register int j = 1; j <= m; ++j)
								{
									T sum = 0;
									for (register int k = 1; k <= 3; ++k)
									{
										sum += Bt.entry[i][k]*DB.entry[k][j];
									}
									Ke.entry[i][j] += weight*t*sum*det_J;
								}
							}
						}
					}
					else // plane_problem == plane_strain
					{
						T a = E/((1 + v)*(1 - 2*v));
						T b = a*(1 - v);
						T c = a*v;
						T d = a*(1 - 2*v)/2;

						for (int q = 1; q <= element_integration_rules.size; ++q)
						{
							T det_J;
							shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

							for (register int i = 1; i <= shape_functions.mesh.nodes_per_element; ++i)
							{
								register int i2 = i*2;
								register int i1 = i2 - 1;
								T p = dN.entry[i][1];
								T q = dN.entry[i][2];
								Bt.entry[i1][1] = p; Bt.entry[i2][1] = 0;
								Bt.entry[i1][2] = 0; Bt.entry[i2][2] = q;
								Bt.entry[i1][3] = q; Bt.entry[i2][3] = p;
								DB.entry[1][i1] = b*p; DB.entry[1][i2] = c*q;
								DB.entry[2][i1] = c*p; DB.entry[2][i2] = b*q;
								DB.entry[3][i1] = d*q; DB.entry[3][i2] = d*p;
							}

							T weight = element_integration_rules.entry[q].weight;
							for (register int i = 1; i <= m; ++i)
							{
								for (register int j = 1; j <= m; ++j)
								{
									register T sum = 0;
									for (register int k = 1; k <= 3; ++k)
									{
										sum += Bt.entry[i][k]*DB.entry[k][j];
									}
									Ke.entry[i][j] += weight*sum*det_J;
								}
							}
						}
					}
				}
				else if (shape_functions.nodes.dimension == 3)
				{
					Matrix<T> Bt(shape_functions.mesh.nodes_per_element*3, 6);
					Matrix<T> DB(6, shape_functions.mesh.nodes_per_element*3);

					T z = material.young_modulus/((1 + material.poisson_ratio)*(1 - 2*material.poisson_ratio));
					T a = z*(1 - material.poisson_ratio);
					T b = z*material.poisson_ratio;
					T c = z*(1 - 2*material.poisson_ratio)/2;

					Ke.Fill(T(0));
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						T det_J;
						shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

						for (register int i = 1; i <= shape_functions.mesh.nodes_per_element; ++i)
						{
							register int i3 = i*3;
							register int i2 = i3 - 1;
							register int i1 = i3 - 2;
							T p = dN.entry[i][1];
							T q = dN.entry[i][2];
							T r = dN.entry[i][3];
							Bt.entry[i1][1] = p; Bt.entry[i2][1] = 0; Bt.entry[i3][1] = 0;
							Bt.entry[i1][2] = 0; Bt.entry[i2][2] = q; Bt.entry[i3][2] = 0;
							Bt.entry[i1][3] = 0; Bt.entry[i2][3] = 0; Bt.entry[i3][3] = r;
							Bt.entry[i1][4] = q; Bt.entry[i2][4] = p; Bt.entry[i3][4] = 0;
							Bt.entry[i1][5] = 0; Bt.entry[i2][5] = r; Bt.entry[i3][5] = q;
							Bt.entry[i1][6] = r; Bt.entry[i2][6] = 0; Bt.entry[i3][6] = p;
							DB.entry[1][i1] = a*p; DB.entry[1][i2] = b*q; DB.entry[1][i3] = b*r;
							DB.entry[2][i1] = b*p; DB.entry[2][i2] = a*q; DB.entry[2][i3] = b*r;
							DB.entry[3][i1] = b*p; DB.entry[3][i2] = b*q; DB.entry[3][i3] = a*r;
							DB.entry[4][i1] = c*q; DB.entry[4][i2] = c*p; DB.entry[4][i3] = 0;
							DB.entry[5][i1] = 0;   DB.entry[5][i2] = c*r; DB.entry[5][i3] = c*q;
							DB.entry[6][i1] = c*r; DB.entry[6][i2] = 0;   DB.entry[6][i3] = c*p;
						}

						T weight = element_integration_rules.entry[q].weight;
						for (register int i = 1; i <= m; ++i)
						{
							for (register int j = 1; j <= m; ++j)
							{
								T sum = 0;
								for (register int k = 1; k <= 6; ++k)
								{
									sum += Bt.entry[i][k]*DB.entry[k][j];
								}
								Ke.entry[i][j] += weight*sum*det_J;
							}
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void FillKeCeMe(const int element_id, Matrix<T>& Ke, Matrix<T>& Ce, Matrix<T>& Me) const throw(Memory::Exception)
		{
			// S. W. Sloan
			// A fast stiffness formulation for finite element analysis of two-dimensional solids
			// International Journal for Numerical Methods in Engineering, Vol. 17, pp. 1313-1323, 1981

			try
			{
				int m = shape_functions.mesh.nodes_per_element*shape_functions.nodes.dimension;

				Vector<T> N(shape_functions.mesh.nodes_per_element);
				Matrix<T> dN(shape_functions.mesh.nodes_per_element, shape_functions.nodes.dimension);

				int material_id = material_index.entry[element_id];
				Material& material = materials.entry[material_id];

				if (shape_functions.nodes.dimension == 2)
				{
					Matrix<T> Bt(shape_functions.mesh.nodes_per_element*2, 3);
					Matrix<T> DB(3, shape_functions.mesh.nodes_per_element*2);
					Ke.Fill(T(0));
					Ce.Fill(T(0));
					Me.Fill(T(0));

					T E = material.young_modulus;
					T v = material.poisson_ratio;

					if (plane_problem == plane_stress)
					{
						T a = E/(1 - v*v);
						T b = a*v;
						T c = a*(1 - v)/2;
						T t = material.thickness;

						for (int q = 1; q <= element_integration_rules.size; ++q)
						{
							T det_J;
							shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

							for (register int i = 1; i <= shape_functions.mesh.nodes_per_element; ++i)
							{
								register int i2 = i*2;
								register int i1 = i2 - 1;
								T p = dN.entry[i][1];
								T q = dN.entry[i][2];
								Bt.entry[i1][1] = p; Bt.entry[i2][1] = 0;
								Bt.entry[i1][2] = 0; Bt.entry[i2][2] = q;
								Bt.entry[i1][3] = q; Bt.entry[i2][3] = p;
								DB.entry[1][i1] = a*p; DB.entry[1][i2] = b*q;
								DB.entry[2][i1] = b*p; DB.entry[2][i2] = a*q;
								DB.entry[3][i1] = c*q; DB.entry[3][i2] = c*p;
							}

							T weight = element_integration_rules.entry[q].weight;
							for (register int i = 1; i <= m; ++i)
							{
								for (register int j = 1; j <= m; ++j)
								{
									T sum = 0;
									for (register int k = 1; k <= 3; ++k)
									{
										sum += Bt.entry[i][k]*DB.entry[k][j];
									}
									Ke.entry[i][j] += weight*t*sum*det_J;
									if (((i - j) % 2) == 0)
									{
										int ni = (i - 1)/2 + 1;
										int nj = (j - 1)/2 + 1;
										Me.entry[i][j] += weight*material.density*t*N.entry[ni]*N.entry[nj]*det_J;
									}
									Ce.entry[i][j] += rayleigh_damping_a*Me.entry[i][j] + rayleigh_damping_b*Ke.entry[i][j];
								}
							}
						}
					}
					else // plane_problem == plane_strain
					{
						T a = E/((1 + v)*(1 - 2*v));
						T b = a*(1 - v);
						T c = a*v;
						T d = a*(1 - 2*v)/2;

						for (int q = 1; q <= element_integration_rules.size; ++q)
						{
							T det_J;
							shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

							for (register int i = 1; i <= shape_functions.mesh.nodes_per_element; ++i)
							{
								register int i2 = i*2;
								register int i1 = i2 - 1;
								T p = dN.entry[i][1];
								T q = dN.entry[i][2];
								Bt.entry[i1][1] = p; Bt.entry[i2][1] = 0;
								Bt.entry[i1][2] = 0; Bt.entry[i2][2] = q;
								Bt.entry[i1][3] = q; Bt.entry[i2][3] = p;
								DB.entry[1][i1] = b*p; DB.entry[1][i2] = c*q;
								DB.entry[2][i1] = c*p; DB.entry[2][i2] = b*q;
								DB.entry[3][i1] = d*q; DB.entry[3][i2] = d*p;
							}

							T weight = element_integration_rules.entry[q].weight;
							for (register int i = 1; i <= m; ++i)
							{
								for (register int j = 1; j <= m; ++j)
								{
									T sum = 0;
									for (register int k = 1; k <= 3; ++k)
									{
										sum += Bt.entry[i][k]*DB.entry[k][j];
									}
									Ke.entry[i][j] += weight*sum*det_J;
									if (((i - j) % 2) == 0)
									{
										int ni = (i - 1)/2 + 1;
										int nj = (j - 1)/2 + 1;
										Me.entry[i][j] += weight*material.density*N.entry[ni]*N.entry[nj]*det_J;
									}
									Ce.entry[i][j] += rayleigh_damping_a*Me.entry[i][j] + rayleigh_damping_b*Ke.entry[i][j];
								}
							}
						}
					}
				}
				else if (shape_functions.nodes.dimension == 3)
				{
					Matrix<T> Bt(shape_functions.mesh.nodes_per_element*3, 6);
					Matrix<T> DB(6, shape_functions.mesh.nodes_per_element*3);

					T z = material.young_modulus/((1 + material.poisson_ratio)*(1 - 2*material.poisson_ratio));
					T a = z*(1 - material.poisson_ratio);
					T b = z*material.poisson_ratio;
					T c = z*(1 - 2*material.poisson_ratio)/2;

					Ke.Fill(T(0));
					Ce.Fill(T(0));
					Me.Fill(T(0));
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						T det_J;
						shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

						for (register int i = 1; i <= shape_functions.mesh.nodes_per_element; ++i)
						{
							register int i3 = i*3;
							register int i2 = i3 - 1;
							register int i1 = i3 - 2;
							T p = dN.entry[i][1];
							T q = dN.entry[i][2];
							T r = dN.entry[i][3];
							Bt.entry[i1][1] = p; Bt.entry[i2][1] = 0; Bt.entry[i3][1] = 0;
							Bt.entry[i1][2] = 0; Bt.entry[i2][2] = q; Bt.entry[i3][2] = 0;
							Bt.entry[i1][3] = 0; Bt.entry[i2][3] = 0; Bt.entry[i3][3] = r;
							Bt.entry[i1][4] = q; Bt.entry[i2][4] = p; Bt.entry[i3][4] = 0;
							Bt.entry[i1][5] = 0; Bt.entry[i2][5] = r; Bt.entry[i3][5] = q;
							Bt.entry[i1][6] = r; Bt.entry[i2][6] = 0; Bt.entry[i3][6] = p;
							DB.entry[1][i1] = a*p; DB.entry[1][i2] = b*q; DB.entry[1][i3] = b*r;
							DB.entry[2][i1] = b*p; DB.entry[2][i2] = a*q; DB.entry[2][i3] = b*r;
							DB.entry[3][i1] = b*p; DB.entry[3][i2] = b*q; DB.entry[3][i3] = a*r;
							DB.entry[4][i1] = c*q; DB.entry[4][i2] = c*p; DB.entry[4][i3] = 0;
							DB.entry[5][i1] = 0;   DB.entry[5][i2] = c*r; DB.entry[5][i3] = c*q;
							DB.entry[6][i1] = c*r; DB.entry[6][i2] = 0;   DB.entry[6][i3] = c*p;
						}

						T weight = element_integration_rules.entry[q].weight;
						for (register int i = 1; i <= m; ++i)
						{
							for (register int j = 1; j <= m; ++j)
							{
								T sum = 0;
								for (register int k = 1; k <= 6; ++k)
								{
									sum += Bt.entry[i][k]*DB.entry[k][j];
								}
								Ke.entry[i][j] += weight*sum*det_J;
								if (((i - j) % 3) == 0)
								{
									int ni = (i - 1)/3 + 1;
									int nj = (j - 1)/3 + 1;
									Me.entry[i][j] += weight*material.density*N.entry[ni]*N.entry[nj]*det_J;
								}
								Ce.entry[i][j] += rayleigh_damping_a*Me.entry[i][j] + rayleigh_damping_b*Ke.entry[i][j];
							}
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void FillU(Vector<T>& u) throw()
		{
			u.Fill(0);

			if (shape_functions.nodes.dimension == 2)
			{
				for (register int i = 1; i <= displacements.size; ++i)
				{
					Displacement& displacement = displacements.entry[i];

					register int i2 = 2*displacement.node_id;
					register int i1 = i2 - 1;
					if (displacement.user_function_x_id == 0)
					{
						u.entry[i1] = displacement.displacement_x;
					}
					else
					{
						T x = shape_functions.nodes.coordinate.entry[displacement.node_id][1];
						T y = shape_functions.nodes.coordinate.entry[displacement.node_id][2];
						T z = 0;
						u.entry[i1] = user_function.entry[displacement.user_function_x_id](x, y, z, time);
					}
					if (displacement.user_function_y_id == 0)
					{
						u.entry[i2] = displacement.displacement_y;
					}
					else
					{
						T x = shape_functions.nodes.coordinate.entry[displacement.node_id][1];
						T y = shape_functions.nodes.coordinate.entry[displacement.node_id][2];
						T z = 0;
						u.entry[i2] = user_function.entry[displacement.user_function_y_id](x, y, z, time);
					}
				}
			}
			else if (shape_functions.nodes.dimension == 3)
			{
				for (register int i = 1; i <= displacements.size; ++i)
				{
					Displacement& displacement = displacements.entry[i];

					register int i3 = 3*displacement.node_id;
					register int i2 = i3 - 1;
					register int i1 = i2 - 1;
					if (displacement.user_function_x_id == 0)
					{
						u.entry[i1] = displacement.displacement_x;
					}
					else
					{
						T x = shape_functions.nodes.coordinate.entry[displacement.node_id][1];
						T y = shape_functions.nodes.coordinate.entry[displacement.node_id][2];
						T z = shape_functions.nodes.coordinate.entry[displacement.node_id][3];
						u.entry[i1] = user_function.entry[displacement.user_function_x_id](x, y, z, time);
					}
					if (displacement.user_function_y_id == 0)
					{
						u.entry[i2] = displacement.displacement_y;
					}
					else
					{
						T x = shape_functions.nodes.coordinate.entry[displacement.node_id][1];
						T y = shape_functions.nodes.coordinate.entry[displacement.node_id][2];
						T z = shape_functions.nodes.coordinate.entry[displacement.node_id][3];
						u.entry[i2] = user_function.entry[displacement.user_function_y_id](x, y, z, time);
					}
					if (displacement.user_function_z_id == 0)
					{
						u.entry[i3] = displacement.displacement_z;
					}
					else
					{
						T x = shape_functions.nodes.coordinate.entry[displacement.node_id][1];
						T y = shape_functions.nodes.coordinate.entry[displacement.node_id][2];
						T z = shape_functions.nodes.coordinate.entry[displacement.node_id][3];
						u.entry[i3] = user_function.entry[displacement.user_function_z_id](x, y, z, time);
					}
				}
			}
		}


		void FillFixed(Vector<bool>& fixed) throw()
		{
			fixed.Fill(false);

			if (shape_functions.nodes.dimension == 2)
			{
				for (register int i = 1; i <= displacements.size; ++i)
				{
					Displacement& displacement = displacements.entry[i];

					register int i2 = 2*displacement.node_id;
					register int i1 = i2 - 1;
					fixed.entry[i1] = displacement.fixed_x;
					fixed.entry[i2] = displacement.fixed_y;
				}
			}
			else if (shape_functions.nodes.dimension == 3)
			{
				for (register int i = 1; i <= displacements.size; ++i)
				{
					Displacement& displacement = displacements.entry[i];

					register int i3 = 3*displacement.node_id;
					register int i2 = i3 - 1;
					register int i1 = i2 - 1;
					fixed.entry[i1] = displacement.fixed_x;
					fixed.entry[i2] = displacement.fixed_y;
					fixed.entry[i3] = displacement.fixed_z;
				}
			}
		}


		void FillF(Vector<T>& f) throw(Memory::Exception)
		{
			try
			{
				f.Fill(0);

				// Fill nodal forces
				if (shape_functions.nodes.dimension == 2)
				{
					for (register int i = 1; i <= nodal_forces.size; ++i)
					{
						NodalForce& node_force = nodal_forces.entry[i];

						register int i2 = 2*node_force.node_id;
						register int i1 = i2 - 1;
						if (node_force.user_function_x_id == 0)
						{
							f.entry[i1] = node_force.force_x;
						}
						else
						{
							T x = shape_functions.nodes.coordinate.entry[node_force.node_id][1];
							T y = shape_functions.nodes.coordinate.entry[node_force.node_id][2];
							T z = 0;
							f.entry[i1] = user_function.entry[node_force.user_function_x_id](x, y, z, time);
						}
						if (node_force.user_function_y_id == 0)
						{
							f.entry[i2] = node_force.force_y;
						}
						else
						{
							T x = shape_functions.nodes.coordinate.entry[node_force.node_id][1];
							T y = shape_functions.nodes.coordinate.entry[node_force.node_id][2];
							T z = 0;
							f.entry[i2] = user_function.entry[node_force.user_function_y_id](x, y, z, time);
						}
					}
				}
				else if (shape_functions.nodes.dimension == 3)
				{
					for (register int i = 1; i <= nodal_forces.size; ++i)
					{
						NodalForce& node_force = nodal_forces.entry[i];

						register int i3 = 3*node_force.node_id;
						register int i2 = i3 - 1;
						register int i1 = i2 - 1;
						if (node_force.user_function_x_id == 0)
						{
							f.entry[i1] = node_force.force_x;
						}
						else
						{
							T x = shape_functions.nodes.coordinate.entry[node_force.node_id][1];
							T y = shape_functions.nodes.coordinate.entry[node_force.node_id][2];
							T z = shape_functions.nodes.coordinate.entry[node_force.node_id][3];
							f.entry[i1] = user_function.entry[node_force.user_function_x_id](x, y, z, time);
						}
						if (node_force.user_function_y_id == 0)
						{
							f.entry[i2] = node_force.force_y;
						}
						else
						{
							T x = shape_functions.nodes.coordinate.entry[node_force.node_id][1];
							T y = shape_functions.nodes.coordinate.entry[node_force.node_id][2];
							T z = shape_functions.nodes.coordinate.entry[node_force.node_id][3];
							f.entry[i2] = user_function.entry[node_force.user_function_y_id](x, y, z, time);
						}
						if (node_force.user_function_z_id == 0)
						{
							f.entry[i3] = node_force.force_z;
						}
						else
						{
							T x = shape_functions.nodes.coordinate.entry[node_force.node_id][1];
							T y = shape_functions.nodes.coordinate.entry[node_force.node_id][2];
							T z = shape_functions.nodes.coordinate.entry[node_force.node_id][3];
							f.entry[i3] = user_function.entry[node_force.user_function_z_id](x, y, z, time);
						}
					}
				}

				// Fill normal forces
				Vector<T> facet_N(shape_functions.mesh.nodes_per_facet);

				if (shape_functions.nodes.dimension == 2)
				{
					for (register int i = 1; i <= normal_forces.size; ++i)
					{
						NormalForce& facet_force = normal_forces.entry[i];
						int material_id = material_index.entry[facet_force.element_id];
						Material& material = materials.entry[material_id];

						double t = plane_problem == plane_stress ? material.thickness : 1;

						for (int q = 1; q <= facet_integration_rules.size; ++q)
						{
							// Normal vector components
							int n1 = facet_force.facet_nodes.entry[1];
							int n2 = facet_force.facet_nodes.entry[2];
							T nfx = -(shape_functions.nodes.coordinate.entry[n2][2] - shape_functions.nodes.coordinate.entry[n1][2]);
							T nfy = shape_functions.nodes.coordinate.entry[n2][1] - shape_functions.nodes.coordinate.entry[n1][1];
							T sn = sqrt(nfx*nfx + nfy*nfy);
							nfx /= sn;
							nfy /= sn;

							T det_J;
							shape_functions.FacetShapeFunctions(facet_force.facet_nodes, facet_integration_rules.entry[q].point, facet_N, det_J);

							T w = facet_integration_rules.entry[q].weight;
							for (register int j = 1; j <= shape_functions.mesh.nodes_per_facet; ++j)
							{
								register int n = (facet_force.facet_nodes.entry[j] - 1)*2;
								if (facet_force.user_function_id == 0)
								{
									f.entry[n + 1] += w*t*(facet_force.force*nfx)*facet_N.entry[j]*det_J;
									f.entry[n + 2] += w*t*(facet_force.force*nfy)*facet_N.entry[j]*det_J;
								}
								else
								{
									T x = shape_functions.nodes.coordinate.entry[facet_force.facet_nodes.entry[j]][1];
									T y = shape_functions.nodes.coordinate.entry[facet_force.facet_nodes.entry[j]][2];
									T z = 0;
									T force = user_function.entry[facet_force.user_function_id](x, y, z, time);
									f.entry[n + 1] += w*t*(force*nfx)*facet_N.entry[j]*det_J;
									f.entry[n + 2] += w*t*(force*nfy)*facet_N.entry[j]*det_J;
								}
							}
						}
					}
				}
				else if (shape_functions.nodes.dimension == 3)
				{
					for (register int i = 1; i <= normal_forces.size; ++i)
					{
						NormalForce& facet_force = normal_forces.entry[i];

						for (int q = 1; q <= facet_integration_rules.size; ++q)
						{
							// Normal vector components
							int n1 = facet_force.facet_nodes.entry[1];
							int n2 = facet_force.facet_nodes.entry[2];
							int n3 = facet_force.facet_nodes.entry[3];
							T ax = shape_functions.nodes.coordinate.entry[n2][1] - shape_functions.nodes.coordinate.entry[n1][1];
							T ay = shape_functions.nodes.coordinate.entry[n2][2] - shape_functions.nodes.coordinate.entry[n1][2];
							T az = shape_functions.nodes.coordinate.entry[n2][3] - shape_functions.nodes.coordinate.entry[n1][3];
							T bx = shape_functions.nodes.coordinate.entry[n3][1] - shape_functions.nodes.coordinate.entry[n1][1];
							T by = shape_functions.nodes.coordinate.entry[n3][2] - shape_functions.nodes.coordinate.entry[n1][2];
							T bz = shape_functions.nodes.coordinate.entry[n3][3] - shape_functions.nodes.coordinate.entry[n1][3];
							T nfx = ay*bz - az*by;
							T nfy = az*bx - ax*bz;
							T nfz = ax*by - ay*bx;
							T sn = sqrt(nfx*nfx + nfy*nfy + nfz*nfz);
							nfx /= sn;
							nfy /= sn;
							nfz /= sn;

							T det_J;
							shape_functions.FacetShapeFunctions(facet_force.facet_nodes, facet_integration_rules.entry[q].point, facet_N, det_J);

							T w = facet_integration_rules.entry[q].weight;
							for (register int j = 1; j <= shape_functions.mesh.nodes_per_facet; ++j)
							{
								register int n = (facet_force.facet_nodes.entry[j] - 1)*3;
								if (facet_force.user_function_id == 0)
								{
									f.entry[n + 1] += w*(facet_force.force*nfx)*facet_N.entry[j]*det_J;
									f.entry[n + 2] += w*(facet_force.force*nfy)*facet_N.entry[j]*det_J;
									f.entry[n + 3] += w*(facet_force.force*nfz)*facet_N.entry[j]*det_J;
								}
								else
								{
									T x = shape_functions.nodes.coordinate.entry[facet_force.facet_nodes.entry[j]][1];
									T y = shape_functions.nodes.coordinate.entry[facet_force.facet_nodes.entry[j]][2];
									T z = shape_functions.nodes.coordinate.entry[facet_force.facet_nodes.entry[j]][3];
									T force = user_function.entry[facet_force.user_function_id](x, y, z, time);
									f.entry[n + 1] += w*(force*nfx)*facet_N.entry[j]*det_J;
									f.entry[n + 2] += w*(force*nfy)*facet_N.entry[j]*det_J;
									f.entry[n + 3] += w*(force*nfz)*facet_N.entry[j]*det_J;
								}
							}
						}
					}
				}

				// Mass forces
				if (use_mass_forces)
				{
					T g = gravity/9.80665;

					Vector<T> N(shape_functions.mesh.nodes_per_element);
					Matrix<T> dN(shape_functions.mesh.nodes_per_element, shape_functions.nodes.dimension);

					const Matrix<int>& connectivity = shape_functions.mesh.connectivity;

					if (shape_functions.nodes.dimension == 2)
					{
						for (int element_id = 1; element_id <= shape_functions.mesh.elements_count; ++element_id)
						{
							int material_id = material_index.entry[element_id];
							Material& material = materials.entry[material_id];

							double t = plane_problem == plane_stress ? material.thickness : 1;

							for (int q = 1; q <= element_integration_rules.size; ++q)
							{
								T det_J;
								shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

								T w = element_integration_rules.entry[q].weight;
								for (register int j = 1; j <= shape_functions.mesh.nodes_per_element; ++j)
								{
									register int n = connectivity.entry[element_id][j]*2; // y axis
									f.entry[n] -= w*t*material.density*g*N.entry[j]*det_J;
								}
							}
						}
					}
					else if (shape_functions.nodes.dimension == 3)
					{
						for (int element_id = 1; element_id <= shape_functions.mesh.elements_count; ++element_id)
						{
							int material_id = material_index.entry[element_id];
							Material& material = materials.entry[material_id];

							for (int q = 1; q <= element_integration_rules.size; ++q)
							{
								T det_J;
								shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

								T w = element_integration_rules.entry[q].weight;
								for (register int j = 1; j <= shape_functions.mesh.nodes_per_element; ++j)
								{
									register int n = connectivity.entry[element_id][j]*3; // z axis
									f.entry[n] -= w*material.density*g*N.entry[j]*det_J;
								}
							}
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		double TotalWeight() throw(Memory::Exception)
		{
			try
			{
				T g = gravity/9.80665;

				Vector<T> N(shape_functions.mesh.nodes_per_element);
				Matrix<T> dN(shape_functions.mesh.nodes_per_element, shape_functions.nodes.dimension);

				double total_weight = 0;

				if (shape_functions.nodes.dimension == 2)
				{
					for (int element_id = 1; element_id <= shape_functions.mesh.elements_count; ++element_id)
					{
						int material_id = material_index.entry[element_id];
						Material& material = materials.entry[material_id];
						double t = plane_problem == plane_stress ? material.thickness : 1;

						for (int q = 1; q <= element_integration_rules.size; ++q)
						{
							T det_J;
							shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

							T w = element_integration_rules.entry[q].weight;
							for (register int j = 1; j <= shape_functions.mesh.nodes_per_element; ++j)
							{
								total_weight += w*t*material.density*g*N.entry[j]*det_J;
							}
						}
					}
				}
				else if (shape_functions.nodes.dimension == 3)
				{
					for (int element_id = 1; element_id <= shape_functions.mesh.elements_count; ++element_id)
					{
						int material_id = material_index.entry[element_id];
						Material& material = materials.entry[material_id];

						for (int q = 1; q <= element_integration_rules.size; ++q)
						{
							T det_J;
							shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

							T w = element_integration_rules.entry[q].weight;
							for (register int j = 1; j <= shape_functions.mesh.nodes_per_element; ++j)
							{
								total_weight += w*material.density*g*N.entry[j]*det_J;
							}
						}
					}
				}

				return total_weight;
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		double GetElementSize(int element_id) throw(Memory::Exception)
		{
			try
			{
				Vector<T> N(shape_functions.mesh.nodes_per_element);
				Matrix<T> dN(shape_functions.mesh.nodes_per_element, shape_functions.nodes.dimension);

				if (shape_functions.nodes.dimension == 2)
				{
					double area = 0;
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						T det_J;
						shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

						T w = element_integration_rules.entry[q].weight;
						for (register int j = 1; j <= shape_functions.mesh.nodes_per_element; ++j)
						{
							area += w*N.entry[j]*det_J;
						}
					}
					return area;
				}
				else if (shape_functions.nodes.dimension == 3)
				{
					double volume = 0;
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						T det_J;
						shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

						T w = element_integration_rules.entry[q].weight;
						for (register int j = 1; j <= shape_functions.mesh.nodes_per_element; ++j)
						{
							volume += w*N.entry[j]*det_J;
						}
					}
					return volume;
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		double GetTotalSize() throw(Memory::Exception)
		{
			try
			{
				Vector<T> N(shape_functions.mesh.nodes_per_element);
				Matrix<T> dN(shape_functions.mesh.nodes_per_element, shape_functions.nodes.dimension);

				if (shape_functions.nodes.dimension == 2)
				{
					double area = 0;
					for (int element_id = 1; element_id <= shape_functions.mesh.elements_count; ++element_id)
					{
						for (int q = 1; q <= element_integration_rules.size; ++q)
						{
							T det_J;
							shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

							T w = element_integration_rules.entry[q].weight;
							for (register int j = 1; j <= shape_functions.mesh.nodes_per_element; ++j)
							{
								area += w*N.entry[j]*det_J;
							}
						}
					}
					return area;
				}
				else if (shape_functions.nodes.dimension == 3)
				{
					double volume = 0;
					for (int element_id = 1; element_id <= shape_functions.mesh.elements_count; ++element_id)
					{
						for (int q = 1; q <= element_integration_rules.size; ++q)
						{
							T det_J;
							shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

							T w = element_integration_rules.entry[q].weight;
							for (register int j = 1; j <= shape_functions.mesh.nodes_per_element; ++j)
							{
								volume += w*N.entry[j]*det_J;
							}
						}
					}
					return volume;
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void WriteDisplacements(File& file, const Vector<int>& node_index, const Vector<T>& u, const int step, T& maximum) const throw(Memory::Exception, File::Exception)
		{
			try
			{
				maximum = 0;

				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				if (shape_functions.nodes.dimension == 2)
				{
					file.Write("Result \"Displacement\" \"Solid\" ");
					file.Write(step, format_integer);
					file.Write(" Vector OnNodes\n");
					file.Write("Values\n");
					for (int n = 1, i = 1; n <= node_index.size; ++n, i += 2)
					{
						file.Write(node_index.entry[n], format_integer);
						file.Write(" ");
						file.Write(u.entry[i], format_float);
						file.Write(" ");
						file.Write(u.entry[i + 1], format_float);
						file.Write(" ");
						file.Write((T)0, format_float);
						file.Write("\n");

						T displacement = sqrt(u.entry[i]*u.entry[i] + u.entry[i + 1]*u.entry[i + 1]);
						if (maximum < displacement)
						{
							maximum = displacement;
						}
					}
					file.Write("End Values\n\n");
				}
				else if (shape_functions.nodes.dimension == 3)
				{
					file.Write("Result \"Displacement\" \"Solid\" ");
					file.Write(step, format_integer);
					file.Write(" Vector OnNodes\n");
					file.Write("Values\n");
					for (int n = 1, i = 1; n <= node_index.size; ++n, i += 3)
					{
						file.Write(node_index.entry[n], format_integer);
						file.Write(" ");
						file.Write(u.entry[i], format_float);
						file.Write(" ");
						file.Write(u.entry[i + 1], format_float);
						file.Write(" ");
						file.Write(u.entry[i + 2], format_float);
						file.Write("\n");

						T displacement = sqrt(u.entry[i]*u.entry[i] + u.entry[i + 1]*u.entry[i + 1] + u.entry[i + 2]*u.entry[i + 2]);
						if (maximum < displacement)
						{
							maximum = displacement;
						}
					}
					file.Write("End Values\n\n");
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void GetStrain(const int element_id, const Vector<T>& u, Vector<Vector<T> >& strain) throw(Memory::Exception)
		{
			Assert(strain.size == element_integration_rules.size);
			Assert(strain.entry[1].size == (shape_functions.nodes.dimension - 1)*3);

			try
			{
				Vector<T> N(shape_functions.mesh.nodes_per_element);
				Matrix<T> dN(shape_functions.mesh.nodes_per_element, shape_functions.nodes.dimension);

				if (shape_functions.nodes.dimension == 2)
				{
					Vector<T> epsilon(3);
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						T det_J;
						shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

						epsilon.Fill(0);
						for (register int n = 1; n <= shape_functions.mesh.nodes_per_element; ++n)
						{
							int node_id = shape_functions.mesh.connectivity.entry[element_id][n];
							T ui = u.entry[(node_id - 1)*2 + 1];
							T vi = u.entry[node_id*2];
							epsilon.entry[1] += dN.entry[n][1]*ui;
							epsilon.entry[2] += dN.entry[n][2]*vi;
							epsilon.entry[3] += dN.entry[n][2]*ui + dN.entry[n][1]*vi;
						}
						strain.entry[q].entry[1] = epsilon.entry[1];
						strain.entry[q].entry[2] = epsilon.entry[2];
						strain.entry[q].entry[3] = epsilon.entry[3];
					}
				}
				else // shape_functions.nodes.dimension == 3
				{
					Vector<T> epsilon(6);
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						T det_J;
						shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

						epsilon.Fill(0);
						for (register int n = 1; n <= shape_functions.mesh.nodes_per_element; ++n)
						{
							int node_id = shape_functions.mesh.connectivity.entry[element_id][n];
							T ui = u.entry[(node_id - 1)*3 + 1];
							T vi = u.entry[(node_id - 1)*3 + 2];
							T wi = u.entry[node_id*3];
							epsilon.entry[1] += dN.entry[n][1]*ui;
							epsilon.entry[2] += dN.entry[n][2]*vi;
							epsilon.entry[3] += dN.entry[n][3]*wi;
							epsilon.entry[4] += dN.entry[n][2]*ui + dN.entry[n][1]*vi;
							epsilon.entry[5] += dN.entry[n][3]*vi + dN.entry[n][2]*wi;
							epsilon.entry[6] += dN.entry[n][1]*wi + dN.entry[n][3]*ui;
						}
						strain.entry[q].entry[1] = epsilon.entry[1];
						strain.entry[q].entry[2] = epsilon.entry[2];
						strain.entry[q].entry[3] = epsilon.entry[3];
						strain.entry[q].entry[4] = epsilon.entry[4];
						strain.entry[q].entry[5] = epsilon.entry[5];
						strain.entry[q].entry[6] = epsilon.entry[6];
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void GetStress(const int element_id, const Vector<Vector<T> >& strain, Vector<Vector<T> >& stress) throw()
		{
			Assert(stress.size == element_integration_rules.size);
			Assert(stress.entry[1].size == (shape_functions.nodes.dimension == 3) ? 6 : (plane_problem == plane_stress) ? 3 : 6);

			int material_id = material_index.entry[element_id];
			Material& material = materials.entry[material_id];
			T E = material.young_modulus;
			T v = material.poisson_ratio;

			if (shape_functions.nodes.dimension == 2)
			{
				if (plane_problem == plane_stress)
				{
					T a = E/(1 - v*v);
					T b = (1 - v)/2;

					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						stress.entry[q].entry[1] = a*(strain.entry[q].entry[1] + strain.entry[q].entry[2]*v);
						stress.entry[q].entry[2] = a*(strain.entry[q].entry[1]*v + strain.entry[q].entry[2]);
						stress.entry[q].entry[3] = a*b*strain.entry[q].entry[3];
					}
				}
				else // plane_problem == plane_strain
				{
					T a = E/((1 + v)*(1 - 2*v));
					T b = a*(1 - v);
					T c = a*v;
					T d = a*(1 - 2*v)/2;

					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						stress.entry[q].entry[1] = b*strain.entry[q].entry[1] + c*strain.entry[q].entry[2];
						stress.entry[q].entry[2] = c*strain.entry[q].entry[1] + b*strain.entry[q].entry[2];
						stress.entry[q].entry[3] = v*(stress.entry[q].entry[1] + stress.entry[q].entry[2]);
						stress.entry[q].entry[4] = d*strain.entry[q].entry[3];
						stress.entry[q].entry[5] = 0;
						stress.entry[q].entry[6] = 0;
					}
				}
			}
			else if (shape_functions.nodes.dimension == 3)
			{
				T a = E*(1 - v)/((1 + v)*(1 - 2*v));
				T b = v/(1 - v);
				T c = (1 - 2*v)/(2 - 2*v);

				for (int q = 1; q <= element_integration_rules.size; ++q)
				{
					stress.entry[q].entry[1] = a*b*strain.entry[q].entry[3] + a*b*strain.entry[q].entry[2] + a*strain.entry[q].entry[1];
					stress.entry[q].entry[2] = a*b*strain.entry[q].entry[3] + a*strain.entry[q].entry[2] + a*b*strain.entry[q].entry[1];
					stress.entry[q].entry[3] = a*strain.entry[q].entry[3] + a*b*strain.entry[q].entry[2] + a*b*strain.entry[q].entry[1];
					stress.entry[q].entry[4] = a*c*strain.entry[q].entry[4];
					stress.entry[q].entry[5] = a*c*strain.entry[q].entry[5];
					stress.entry[q].entry[6] = a*c*strain.entry[q].entry[6];
				}
			}
		}


		void GetVonMises(const Vector<Vector<T> >& stress, Vector<T>& von_mises) throw()
		{
			Assert(von_mises.size == element_integration_rules.size);

			if (shape_functions.nodes.dimension == 2)
			{
				if (plane_problem == plane_stress)
				{
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						T sx = stress.entry[q].entry[1];
						T sy = stress.entry[q].entry[2];
						T txy = stress.entry[q].entry[3];
						von_mises.entry[q] = sqrt(sx*sx - sx*sy + sy*sy + 3*txy*txy);
					}
				}
				else // plane_strain
				{
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						T sx = stress.entry[q].entry[1];
						T sy = stress.entry[q].entry[2];
						T sz = stress.entry[q].entry[3];
						T txy = stress.entry[q].entry[4];
						von_mises.entry[q] = sqrt(0.5*((sx - sy)*(sx - sy) + (sy - sz)*(sy - sz) + (sz - sx)*(sz - sx) + 6*txy*txy));
					}
				}
			}
			else if (shape_functions.nodes.dimension == 3)
			{
				for (int q = 1; q <= element_integration_rules.size; ++q)
				{
					// E. Oñate. Structural Analysis with the Finite Element Method. Linear Statics. Volume 1. Basis and Solids. Springer. 2009. pp. 255-256
					T sx = stress.entry[q].entry[1];
					T sy = stress.entry[q].entry[2];
					T sz = stress.entry[q].entry[3];
					T txy = stress.entry[q].entry[4];
					T tyz = stress.entry[q].entry[5];
					T txz = stress.entry[q].entry[6];
					T J2 = ((sx - sy)*(sx - sy) + (sx - sz)*(sx - sz) + (sy - sz)*(sy - sz) + 6*(txy*txy + tyz*tyz + txz*txz))/6;
					von_mises.entry[q] = sqrt(3*J2);
				}
			}
		}


		void WriteStrainStressGaussPoints(File& file) throw(Memory::Exception, File::Exception)
		{
			static const char* element_name[] =
			{
				"(Undefined)",
				"Linear",
				"Triangle",
				"Quadrilateral",
				"Tetrahedra",
				"Hexahedra"
			};

			try
			{
				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("GaussPoints \"StrainStressGaussPoints\" ElemType ");
				file.Write(element_name[shape_functions.mesh.element_type]);
				file.Write("\n");

				file.Write("Number Of Gauss Points: ");
				file.Write(element_integration_rules.size, format_integer);
				file.Write("\n");
 
				file.Write("Natural Coordinates: Given\n");
				for (int q = 1; q <= element_integration_rules.size; ++q)
				{
					for (int d = 1; d <= shape_functions.nodes.dimension; ++d)
					{
						file.Write(" ");							
						file.Write(element_integration_rules.entry[q].point.entry[d], format_float);
					}
					file.Write("\n");
				}

				file.Write("End GaussPoints\n\n");
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void WriteStrain(File& file, const Vector<int>& element_index, const Vector<int>& node_index, const Vector<T>& u, int step) throw(Memory::Exception, File::Exception)
		{
			try
			{
				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("Result \"Strain\" \"Solid\" ");
				file.Write(step, format_integer);
				file.Write(" Matrix OnGaussPoints \"StrainStressGaussPoints\"\n");
				file.Write("Values\n");

				int N = (shape_functions.nodes.dimension == 3) ? 6 : 3;

				Vector<T> global_u(shape_functions.mesh.nodes_count*shape_functions.nodes.dimension);
				for (int n = 1; n <= node_index.size; ++n)
				{
					int register node_id = node_index.entry[n];
					for (register int d = 1; d <= shape_functions.nodes.dimension; ++d)
					{
						register int l = (n - 1)*shape_functions.nodes.dimension + d;
						register int g = (node_id - 1)*shape_functions.nodes.dimension + d;
						global_u.entry[g] = u.entry[l];
					}
				}

				Vector<Vector<T> > strain(element_integration_rules.size);
				for (int i = 1; i <= element_integration_rules.size; ++i)
				{
					strain.entry[i].Resize(N);
				}
				for (int e = 1; e <= element_index.size; ++e)
				{
					int element_id = element_index.entry[e];


					GetStrain(element_id, global_u, strain);
					file.Write(element_id, format_integer);
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						for (int d = 1; d <= N; ++d)
						{
							file.Write(" ");							
							file.Write(strain.entry[q].entry[d], format_float);
						}
						file.Write("\n");
					}
				}

				file.Write("End Values\n\n");
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void WriteStress(File& file, const Vector<int>& element_index, const Vector<int>& node_index, const Vector<T>& u, int step) throw(Memory::Exception, File::Exception)
		{
			try
			{
				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("Result \"Stress\" \"Solid\" ");
				file.Write(step, format_integer);
				file.Write(" Matrix OnGaussPoints \"StrainStressGaussPoints\"\n");
				file.Write("Values\n");

				int N = (shape_functions.nodes.dimension == 3) ? 6 : 3;
				int M = (shape_functions.nodes.dimension == 3) ? 6 : (plane_problem == plane_stress) ? 3 : 6;

				Vector<T> global_u(shape_functions.mesh.nodes_count*shape_functions.nodes.dimension);
				for (int n = 1; n <= node_index.size; ++n)
				{
					int register node_id = node_index.entry[n];
					for (register int d = 1; d <= shape_functions.nodes.dimension; ++d)
					{
						register int l = (n - 1)*shape_functions.nodes.dimension + d;
						register int g = (node_id - 1)*shape_functions.nodes.dimension + d;
						global_u.entry[g] = u.entry[l];
					}
				}

				Vector<Vector<T> > strain(element_integration_rules.size);
				Vector<Vector<T> > stress(element_integration_rules.size);
				for (int i = 1; i <= element_integration_rules.size; ++i)
				{
					strain.entry[i].Resize(N);
					stress.entry[i].Resize(M);
				}
				for (int e = 1; e <= element_index.size; ++e)
				{
					int element_id = element_index.entry[e];

					GetStrain(element_id, global_u, strain);
					GetStress(element_id, strain, stress);
					file.Write(element_id, format_integer);
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						for (int d = 1; d <= M; ++d)
						{
							file.Write(" ");							
							file.Write(stress.entry[q].entry[d], format_float);
						}
						file.Write("\n");
					}
				}

				file.Write("End Values\n\n");
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void WriteVonMises(File& file, const Vector<int>& element_index, const Vector<int>& node_index, const Vector<T>& u, int step, T& maximum) throw(Memory::Exception, File::Exception)
		{
			try
			{
				maximum = Float<T>::minimum;

				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("Result \"Von Mises\" \"Solid\" ");
				file.Write(step, format_integer);
				file.Write(" Scalar OnGaussPoints \"StrainStressGaussPoints\"\n");
				file.Write("Values\n");

				int N = (shape_functions.nodes.dimension == 3) ? 6 : 3;
				int M = (shape_functions.nodes.dimension == 3) ? 6 : (plane_problem == plane_stress) ? 3 : 6;

				Vector<T> global_u(shape_functions.mesh.nodes_count*shape_functions.nodes.dimension);
				for (int n = 1; n <= node_index.size; ++n)
				{
					int register node_id = node_index.entry[n];
					for (register int d = 1; d <= shape_functions.nodes.dimension; ++d)
					{
						register int l = (n - 1)*shape_functions.nodes.dimension + d;
						register int g = (node_id - 1)*shape_functions.nodes.dimension + d;
						global_u.entry[g] = u.entry[l];
					}
				}

				Vector<Vector<T> > strain(element_integration_rules.size);
				Vector<Vector<T> > stress(element_integration_rules.size);
				Vector<T> von_mises(element_integration_rules.size);
				for (int i = 1; i <= element_integration_rules.size; ++i)
				{
					strain.entry[i].Resize(N);
					stress.entry[i].Resize(M);
				}
				for (int e = 1; e <= element_index.size; ++e)
				{
					int element_id = element_index.entry[e];

					GetStrain(element_id, global_u, strain);
					GetStress(element_id, strain, stress);
					GetVonMises(stress, von_mises);
					file.Write(element_id, format_integer);
					file.Write(" ");
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						file.Write(von_mises.entry[q], format_float);
						file.Write("\n");

						if (maximum < von_mises.entry[q])
						{
							maximum = von_mises.entry[q];
						}
					}
				}

				file.Write("End Values\n\n");
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void WriteDensityGaussPoints(File& file) throw(Memory::Exception, File::Exception)
		{
			static const char* element_name[] =
			{
				"(Undefined)",
				"Linear",
				"Triangle",
				"Quadrilateral",
				"Tetrahedra",
				"Hexahedra"
			};

			try
			{
				file.Write("GaussPoints \"DensityGaussPoints\" ElemType ");
				file.Write(element_name[shape_functions.mesh.element_type]);
				file.Write("\n");
				file.Write("Number Of Gauss Points: 1\n");
				file.Write("Natural Coordinates: internal\n");
				file.Write("End GaussPoints\n\n");
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void WriteDensity(File& file, const Vector<int>& element_index, int step) throw(Memory::Exception, File::Exception)
		{
			try
			{
				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("Result \"Density\" \"Solid\" ");
				file.Write(step, format_integer);
				file.Write(" Scalar OnGaussPoints \"DensityGaussPoints\"\n");
				file.Write("Values\n");

				for (int e = 1; e <= element_index.size; ++e)
				{
					int element_id = element_index.entry[e];
					int material_id = material_index.entry[element_id];

					file.Write(element_id, format_integer);
					file.Write(" ");							
					file.Write(materials.entry[material_id].density, format_float);
					file.Write("\n");
				}

				file.Write("End Values\n\n");
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


	private:

		Solid& operator = (const Solid&) throw()
		{
			return *this;
		}
};

#endif
