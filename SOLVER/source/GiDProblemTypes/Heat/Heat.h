// Heat.h
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

#ifndef _Heat_h_
#define _Heat_h_

#include <Basic/Float.h>
#include <Basic/Integer.h>
#include <Basic/Log.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>
#include <FiniteElement/Shape.h>
#include <FiniteElement/ShapeFunctions.h>
#include <FiniteElement/ShapeIntegrationRule.h>
#include <File/File.h>
#include <Math/Formula.h>

#include <string.h>


#define HEAT_SECTION_NAME_MAX_SIZE 256
#define HEAT_USER_FUNCTIONS_COUNT 10
#define HEAT_USER_FUNCTION_MAX_SIZE 65536


enum ProblemType
{
	time_independent = 1,
	time_dependent = 2
};


template <typename T>
class Heat
{
	public:

		struct Material
		{
			T thermal_conductivity;
			T mass_density;
			T specific_heat_capacity;
		};

		struct Temperature
		{
			int node_id;
			int user_function_id;
			T temperature;
			bool fixed;
		};

		struct HeatFlow
		{
			int user_function_id;
			T flux;
			Vector<int> facet_nodes;
		};

		struct Source
		{
			int element_id;
			int user_function_id;
			T heat;
		};

		// General
		bool calculate_temperature;
		bool calculate_flux;
		ProblemType problem_type;
		bool save_mesh;
		bool save_system_of_equations;

		// Dynamic
		T time_per_step;
		int steps;
		int result_every_steps;
		T time_scheme_factor;

		// User functions
		T time;
		Vector<Formula<T> > user_function;

		Vector<Material> materials;
		Vector<Temperature> temperatures;
		Vector<HeatFlow> heat_flows;
		Vector<Source> sources;
		Vector<IntegrationRule<T> > element_integration_rules;
		Vector<IntegrationRule<T> > facet_integration_rules;

		const ShapeFunctions<T>& shape_functions;
		const Vector<int>& material_index;


		Heat(const char* file_name, const ShapeFunctions<T>& shape_functions, const Vector<int>& material_index) throw(Memory::Exception, File::Exception)
		:	calculate_temperature(),
			calculate_flux(),
			problem_type(),
			save_mesh(),
			save_system_of_equations(),
			time_per_step(),
			steps(),
			result_every_steps(),
			time_scheme_factor(),
			time(0),
			user_function(HEAT_USER_FUNCTIONS_COUNT),
			materials(),
			temperatures(),
			heat_flows(),
			sources(),
			element_integration_rules(),
			facet_integration_rules(),
			shape_functions(shape_functions),
			material_index(material_index)
		{
			try
			{
				char section_name[HEAT_SECTION_NAME_MAX_SIZE];
				int tmp_integer;

				// Load problem
				File file;
				file.Open(file_name);
				file.SkipComments(';');

				// Load general parameters
				file.Read(section_name, HEAT_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{General}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(calculate_temperature);

				file.SkipComments(';');
				file.Read(calculate_flux);
				file.SkipComments(';');
				file.Read(tmp_integer);
				problem_type = (ProblemType)tmp_integer;
				file.SkipComments(';');
				file.Read(save_mesh);
				file.SkipComments(';');
				file.Read(save_system_of_equations);
				file.SkipComments(';');

				// Load dynamic problem parameters
				file.Read(section_name, HEAT_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{TimeDependent}") != 0)
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

				// Load user defined functions
				file.Read(section_name, HEAT_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{UserFunctions}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				for (int i = 1; i <= HEAT_USER_FUNCTIONS_COUNT; ++i)
				{
					char function_text[HEAT_USER_FUNCTION_MAX_SIZE];
					file.ReadLine(function_text, HEAT_USER_FUNCTION_MAX_SIZE - 10);
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
				file.Read(section_name, HEAT_SECTION_NAME_MAX_SIZE);
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
					file.Read(materials.entry[m].thermal_conductivity);
					file.Read(materials.entry[m].mass_density);
					file.Read(materials.entry[m].specific_heat_capacity);
					file.SkipComments(';');
				}

				// Load temperatures conditions
				int temperature_count;
				file.Read(section_name, HEAT_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Temperature}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(temperature_count);
				file.SkipComments(';');
				temperatures.Resize(temperature_count);
				for (int c = 1; c <= temperature_count; ++c)
				{
					temperatures.entry[c].fixed = false;
				}
				for (int c = 1; c <= temperature_count; ++c)
				{
					file.Read(temperatures.entry[c].node_id);
					file.Read(temperatures.entry[c].user_function_id);
					file.Read(temperatures.entry[c].temperature);
					file.Read(temperatures.entry[c].fixed);
					file.SkipComments(';');
				}

				// Load heat flow conditions
				int heat_flow_count;
				file.Read(section_name, HEAT_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{HeatFlow}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(heat_flow_count);
				file.SkipComments(';');
				heat_flows.Resize(heat_flow_count);
				for (int c = 1; c <= heat_flow_count; ++c)
				{
					file.Read(heat_flows.entry[c].user_function_id);
					file.Read(heat_flows.entry[c].flux);
					heat_flows.entry[c].facet_nodes.Resize(shape_functions.mesh.nodes_per_facet);
					for (int f = 1; f <= shape_functions.mesh.nodes_per_facet; ++f)
					{
						file.Read(heat_flows.entry[c].facet_nodes.entry[f]);
					}
					file.SkipComments(';');
				}

				// Load sources conditions
				int source_count;
				file.Read(section_name, HEAT_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Source}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(source_count);
				file.SkipComments(';');
				sources.Resize(source_count);
				for (int c = 1; c <= source_count; ++c)
				{
					file.Read(sources.entry[c].element_id);
					file.Read(sources.entry[c].user_function_id);
					file.Read(sources.entry[c].heat);
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

				static const char* problem_name[] =
				{
					"(Undefined)",
					"Stationary",
					"Dynamic"
				};

				Log(1, "Problem --------------------------------------------------------------");
				Log(1, "-Problem type:             %s", problem_name[problem_type]);
				Log(1, "-Save mesh:                %s", save_mesh ? "yes" : "no");
				Log(1, "-Save system of equations: %s", save_system_of_equations ? "yes" : "no");
				Log(1, "-Materials used:           %i", materials.size);
				Log(1, "-Temperature conditions:   %i", temperatures.size);
				Log(1, "-Heat flow conditions:     %i", heat_flows.size);
				Log(1, "-Source conditions:        %i", sources.size);
				if (problem_type == time_dependent)
				{
					Log(1, "Time dependent:");
					Log(1, "-Time per step:            %g", time_per_step);
					Log(1, "-Steps:                    %i", steps);
					Log(1, "-Result every steps:       %i", result_every_steps);
					Log(1, "-Time alpha factor:        %g", time_scheme_factor);
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void FillKe(const int element_id, Matrix<T>& Ke) const throw(Memory::Exception)
		{
			try
			{
				int n = shape_functions.mesh.nodes_per_element;
				int d = shape_functions.nodes.dimension;

				Vector<T> N(n);
				Matrix<T> dN(n, d);

				int material_id = material_index.entry[element_id];
				Material& material = materials.entry[material_id];

				Ke.Fill(T(0));
				for (int q = 1; q <= element_integration_rules.size; ++q)
				{
					T det_J;
					shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

					T weight = element_integration_rules.entry[q].weight;
					for (register int i = 1; i <= n; ++i)
					{
						for (register int j = 1; j <= n; ++j)
						{
							register T sum = 0;
							for (register int k = 1; k <= d; ++k)
							{
								sum += dN.entry[i][k]*dN.entry[j][k];
							}
							Ke.entry[i][j] += weight*material.thermal_conductivity*sum*det_J;
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void FillKeCe(const int element_id, Matrix<T>& Ke, Matrix<T>& Ce) const throw(Memory::Exception)
		{
			try
			{
				// R. W. Lewis, P. Nithiarasu, K. N. Seetharamu
				// Fundamentals of the Finite Element Method for Heat and Fluid Flow
				// Wiley, 2004, pp 156-160

				int n = shape_functions.mesh.nodes_per_element;
				int d = shape_functions.nodes.dimension;

				Vector<T> N(n);
				Matrix<T> dN(n, d);

				int material_id = material_index.entry[element_id];
				Material& material = materials.entry[material_id];

				Ke.Fill(T(0));
				Ce.Fill(T(0));
				for (int q = 1; q <= element_integration_rules.size; ++q)
				{
					T det_J;
					shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

					T weight = element_integration_rules.entry[q].weight;
					for (register int i = 1; i <= n; ++i)
					{
						for (register int j = 1; j <= n; ++j)
						{
							register T sum = 0;
							for (register int k = 1; k <= d; ++k)
							{
								sum += dN.entry[i][k]*dN.entry[j][k];
							}
							Ke.entry[i][j] += weight*material.thermal_conductivity*sum*det_J;

							Ce.entry[i][j] += weight*material.mass_density*material.specific_heat_capacity*N.entry[i]*N.entry[j]*det_J;
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

			for (register int i = 1; i <= temperatures.size; ++i)
			{
				Temperature& temperature_node = temperatures.entry[i];

				int n = temperature_node.node_id;
				if (temperature_node.user_function_id == 0)
				{
					u.entry[n] = temperature_node.temperature;
				}
				else
				{
					T x = shape_functions.nodes.coordinate.entry[n][1];
					T y = shape_functions.nodes.coordinate.entry[n][2];
					T z = (shape_functions.nodes.dimension == 3) ? shape_functions.nodes.coordinate.entry[n][3] : 0;
					u.entry[n] = user_function.entry[temperature_node.user_function_id](x, y, z, time);
				}
			}
		}


		void FillFixed(Vector<bool>& fixed) const throw()
		{
			fixed.Fill(false);

			for (register int i = 1; i <= temperatures.size; ++i)
			{
				Temperature& temperature_node = temperatures.entry[i];

				int n = temperature_node.node_id;
				fixed.entry[n] = temperature_node.fixed;
			}
		}


		void FillF(Vector<T>& f) throw(Memory::Exception)
		{
			try
			{
				f.Fill(0);

				// HeatFlow conditions
				Vector<T> facet_N(shape_functions.mesh.nodes_per_facet);

				for (int i = 1; i <= heat_flows.size; ++i)
				{
					HeatFlow& heat_flow = heat_flows.entry[i];
					for (int q = 1; q <= facet_integration_rules.size; ++q)
					{
						T det_J;
						shape_functions.FacetShapeFunctions(heat_flow.facet_nodes, facet_integration_rules.entry[q].point, facet_N, det_J);

						T w = facet_integration_rules.entry[q].weight;
						for (register int j = 1; j <= shape_functions.mesh.nodes_per_facet; ++j)
						{
							register int n = heat_flow.facet_nodes.entry[j];
							if (heat_flow.user_function_id == 0)
							{
								f.entry[n] += w*heat_flow.flux*facet_N.entry[j]*det_J;
							}
							else
							{
								T x = shape_functions.nodes.coordinate.entry[n][1];
								T y = shape_functions.nodes.coordinate.entry[n][2];
								T z = (shape_functions.nodes.dimension == 3) ? shape_functions.nodes.coordinate.entry[n][3] : 0;
								T flux = user_function.entry[heat_flow.user_function_id](x, y, z, time);
								f.entry[n] += w*flux*facet_N.entry[j]*det_J;
							}
						}
					}
				}

				// Source conditions
				const Matrix<int>& connectivity = shape_functions.mesh.connectivity;
				Vector<T> N(shape_functions.mesh.nodes_per_element);
				Matrix<T> dN(shape_functions.mesh.nodes_per_element, shape_functions.nodes.dimension);

				for (int i = 1; i <= sources.size; ++i)
				{
					Source& source = sources.entry[i];
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						int e = source.element_id;
						T det_J;
						shape_functions.ElementShapeFunctions(e, element_integration_rules.entry[q].point, N, dN, det_J);

						T w = element_integration_rules.entry[q].weight;
						for (register int j = 1; j <= shape_functions.mesh.nodes_per_element; ++j)
						{
							register int n = connectivity.entry[e][j];
							if (source.user_function_id == 0)
							{
								f.entry[n] += w*source.heat*N.entry[j]*det_J;
							}
							else
							{
								T x = shape_functions.nodes.coordinate.entry[n][1];
								T y = shape_functions.nodes.coordinate.entry[n][2];
								T z = (shape_functions.nodes.dimension == 3) ? shape_functions.nodes.coordinate.entry[n][3] : 0;
								T heat = user_function.entry[source.user_function_id](x, y, z, time);
								f.entry[n] += w*heat*N.entry[j]*det_J;
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


		void WriteTemperatures(File& file, const Vector<int>& node_index, const Vector<T>& u, int step) const throw(Memory::Exception, File::Exception)
		{
			try
			{
				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("Result \"Temperature\" \"Heat\" ");
				file.Write(step, format_integer);
				file.Write(" Scalar OnNodes\n");
				file.Write("Values\n");
				for (int n = 1; n <= node_index.size; ++n)
				{
					file.Write(node_index.entry[n], format_integer);
					file.Write(" ");
					file.Write(u.entry[n], format_float);
					file.Write("\n");
				}
				file.Write("End Values\n\n");
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void GetFlux(const int element_id, const Vector<T>& u, Vector<Vector<T> >& flux) const throw(Memory::Exception)
		{
			try
			{
				Vector<T> N(shape_functions.mesh.nodes_per_element);
				Matrix<T> dN(shape_functions.mesh.nodes_per_element, shape_functions.nodes.dimension);

				int material_id = material_index.entry[element_id];
				Material& material = materials.entry[material_id];

				for (int q = 1; q <= element_integration_rules.size; ++q)
				{
					T det_J;
					shape_functions.ElementShapeFunctions(element_id, element_integration_rules.entry[q].point, N, dN, det_J);

					for (register int d = 1; d <= shape_functions.nodes.dimension; ++d)
					{
						register T sum = 0;
						for (register int i = 1; i <= shape_functions.mesh.nodes_per_element; ++i)
						{
							int node_id = shape_functions.mesh.connectivity.entry[element_id][i];
							sum += dN.entry[i][d]*u.entry[node_id];
						}
						flux.entry[q].entry[d] = -material.thermal_conductivity*sum;
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void WriteFluxGaussPoints(File& file) const throw(Memory::Exception, File::Exception)
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

				file.Write("GaussPoints \"FluxGaussPoints\" ElemType ");
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


		void WriteFlux(File& file, const Vector<int>& element_index, const Vector<int>& node_index, const Vector<T>& u, int step) const throw(Memory::Exception, File::Exception)
		{
			try
			{
				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("Result \"Flux\" \"Heat\" ");
				file.Write(step, format_integer);
				file.Write(" Vector OnGaussPoints \"FluxGaussPoints\"\n");
				file.Write("Values\n");

				Vector<Vector<T> > flux(element_integration_rules.size);
				for (int i = 1; i <= flux.size; ++i)
				{
					flux.entry[i].Resize(shape_functions.nodes.dimension);
				}

				Vector<T> global_u(shape_functions.mesh.nodes_count);
				for (register int n = 1; n <= node_index.size; ++n)
				{
					register int node_id = node_index.entry[n];
					global_u.entry[node_id] = u.entry[n];
				}

				for (int e = 1; e <= element_index.size; ++e)
				{
					int element_id = element_index.entry[e];
					GetFlux(element_id, global_u, flux);
					file.Write(element_id, format_integer);
					for (int q = 1; q <= element_integration_rules.size; ++q)
					{
						for (int d = 1; d <= shape_functions.nodes.dimension; ++d)
						{
							file.Write(" ");							
							file.Write(flux.entry[q].entry[d], format_float);
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


	private:

		Heat& operator = (const Heat&) throw()
		{
			return *this;
		}
};

#endif
