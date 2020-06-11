// ElectricPotential.h
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

#ifndef _ElectricPotential_h_
#define _ElectricPotential_h_

#include <Basic/Float.h>
#include <Basic/Integer.h>
#include <Basic/Log.h>
#include <Basic/System.h>
#include <Container/CSRMatrix.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>
#include <FiniteElement/Shape.h>
#include <FiniteElement/ShapeFunctions.h>
#include <FiniteElement/ShapeIntegrationRule.h>
#include <File/File.h>
#include <Math/Formula.h>

#include <string.h>


#define ELECTRIC_POTENTIAL_SECTION_NAME_MAX_SIZE 256


enum ProblemType
{
	problem_simple = 1,
	problem_capacitance_matrix = 2,
	problem_sensitivity_analysis = 3
};


template <typename T>
class ElectricPotential
{
	public:

		struct Material
		{
			T permittivity;
			bool sensitivity_use;
		};

		struct Potential
		{
			int node_id;
			int electrode;
			T potential;
		};

		struct ElectricField
		{
			T field_strength;
			Vector<int> facet_nodes;
		};

		struct Source
		{
			int element_id;
			T charge_density;
		};

		// General
		ProblemType problem_type;
		bool save_mesh;

		// Simple
		bool calculate_potential;
		bool calculate_electric_field;
		bool save_system_of_equations;

		// Capacitance matrix
		int electrodes_count;
		T electrodes_voltage;
		int electrode_active;

		// Sensitivity analysis
		int electrode_segments;
		int segments_per_step;
		bool result_on_nodes;

		Vector<Material> materials;
		Vector<Potential> potentials;
		Vector<ElectricField> electric_fields;
		Vector<Source> sources;
		Vector<IntegrationRule<T> > element_integration_rules;
		Vector<IntegrationRule<T> > facet_integration_rules;

		Vector<Vector<int> > electrode_boundary;

		const ShapeFunctions<T>& shape_functions;
		const Vector<int> material_index;


		ElectricPotential(const ElectricPotential& electric_potential) throw(Memory::Exception)
		:	problem_type(electric_potential.problem_type),
			save_mesh(electric_potential.save_mesh),
			calculate_potential(electric_potential.calculate_potential),
			calculate_electric_field(electric_potential.calculate_electric_field),
			save_system_of_equations(electric_potential.save_system_of_equations),
			electrodes_count(electric_potential.electrodes_count),
			electrodes_voltage(electric_potential.electrodes_voltage),
			electrode_active(electric_potential.electrode_active),
			electrode_segments(electric_potential.electrode_segments),
			segments_per_step(electric_potential.segments_per_step),
			result_on_nodes(electric_potential.result_on_nodes),
			materials(electric_potential.materials),
			potentials(electric_potential.potentials),
			electric_fields(electric_potential.electric_fields),
			sources(electric_potential.sources),
			element_integration_rules(electric_potential.element_integration_rules.size),
			facet_integration_rules(electric_potential.facet_integration_rules.size),
			electrode_boundary(electric_potential.electrode_boundary.size),
			shape_functions(electric_potential.shape_functions),
			material_index(electric_potential.material_index)
		{
			for (int q = 1; q <= element_integration_rules.size; ++q)
			{
				element_integration_rules.entry[q].point.Resize(electric_potential.element_integration_rules.entry[q].point.size);
				element_integration_rules.entry[q].point = electric_potential.element_integration_rules.entry[q].point;
				element_integration_rules.entry[q].weight = electric_potential.element_integration_rules.entry[q].weight;
			}

			for (int q = 1; q <= facet_integration_rules.size; ++q)
			{
				facet_integration_rules.entry[q].point.Resize(electric_potential.facet_integration_rules.entry[q].point.size);
				facet_integration_rules.entry[q].point = electric_potential.facet_integration_rules.entry[q].point;
				facet_integration_rules.entry[q].weight = electric_potential.facet_integration_rules.entry[q].weight;
			}

			for (int b = 1; b <= electrode_boundary.size; ++b)
			{
				electrode_boundary.entry[b].Resize(electric_potential.electrode_boundary.entry[b].size);
				electrode_boundary.entry[b] = electric_potential.electrode_boundary.entry[b];
			}
		}


		ElectricPotential(const char* file_name, const ShapeFunctions<T>& shape_functions, const Vector<int>& material_index) throw(Memory::Exception, File::Exception)
		:	problem_type(),
			save_mesh(),
			calculate_potential(),
			calculate_electric_field(),
			save_system_of_equations(),
			electrodes_count(),
			electrodes_voltage(),
			electrode_active(),
			electrode_segments(),
			segments_per_step(),
			result_on_nodes(),
			materials(),
			potentials(),
			electric_fields(),
			sources(),
			element_integration_rules(),
			facet_integration_rules(),
			electrode_boundary(),
			shape_functions(shape_functions),
			material_index(material_index)
		{
			try
			{
				char section_name[ELECTRIC_POTENTIAL_SECTION_NAME_MAX_SIZE];
				int tmp_integer;
				const Mesh& mesh = shape_functions.mesh;

				// Load problem
				File file;
				file.Open(file_name);
				file.SkipComments(';');

				// General
				file.Read(section_name, ELECTRIC_POTENTIAL_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{General}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(tmp_integer);
				problem_type = (ProblemType)tmp_integer;
				file.SkipComments(';');
				file.Read(save_mesh);
				file.SkipComments(';');

				file.Read(section_name, ELECTRIC_POTENTIAL_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Simple}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(calculate_potential);
				file.SkipComments(';');
				file.Read(calculate_electric_field);
				file.SkipComments(';');
				file.Read(save_system_of_equations);
				file.SkipComments(';');

				// Capacitance matrix
				file.Read(section_name, ELECTRIC_POTENTIAL_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{CapacitanceMatrix}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(electrodes_count);
				file.SkipComments(';');
				file.Read(electrodes_voltage);
				file.SkipComments(';');

				// Sensitivity analysis
				file.Read(section_name, ELECTRIC_POTENTIAL_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{SensitivityAnalysis}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(electrode_segments);
				file.SkipComments(';');
				file.Read(segments_per_step);
				file.SkipComments(';');
				file.Read(result_on_nodes);
				file.SkipComments(';');

				// Load materials
				int materials_count;
				file.Read(section_name, ELECTRIC_POTENTIAL_SECTION_NAME_MAX_SIZE);
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
					file.Read(materials.entry[m].sensitivity_use);
					file.SkipComments(';');
					T relative_permittivity;
					file.Read(relative_permittivity);
					materials.entry[m].permittivity = relative_permittivity*8.85418781762038985053656303e-12;
					file.SkipComments(';');
				}

				// Load potentials conditions
				int potential_count;
				file.Read(section_name, ELECTRIC_POTENTIAL_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Potential}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(potential_count);
				file.SkipComments(';');
				potentials.Resize(potential_count);
				for (int c = 1; c <= potential_count; ++c)
				{
					file.Read(potentials.entry[c].node_id);
					file.Read(potentials.entry[c].potential);
					file.Read(potentials.entry[c].electrode);
					file.SkipComments(';');
				}

				// Load electric field conditions
				int electric_field_count;
				file.Read(section_name, ELECTRIC_POTENTIAL_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{ElectricField}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(electric_field_count);
				file.SkipComments(';');
				electric_fields.Resize(electric_field_count);
				for (int c = 1; c <= electric_field_count; ++c)
				{
					file.Read(electric_fields.entry[c].field_strength);
					electric_fields.entry[c].facet_nodes.Resize(mesh.nodes_per_facet);
					for (int f = 1; f <= mesh.nodes_per_facet; ++f)
					{
						file.Read(electric_fields.entry[c].facet_nodes.entry[f]);
					}
					file.SkipComments(';');
				}

				// Load sources conditions
				int source_count;
				file.Read(section_name, ELECTRIC_POTENTIAL_SECTION_NAME_MAX_SIZE);
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
					file.Read(sources.entry[c].charge_density);
					file.SkipComments(';');
				}

				file.Close();

				// Fill integration rules
				int element_integration_points_count;
				int facet_integration_points_count;
				switch (mesh.element_type)
				{
					case shape_triangle:
					{
						element_integration_points_count = (mesh.nodes_per_element == 3) ? 1 : 3;
						facet_integration_points_count = (mesh.nodes_per_facet == 2) ? 1 : 2;
						break;
					}
					case shape_quadrilateral:
					{
						element_integration_points_count = (mesh.nodes_per_element == 4) ? 2 : (mesh.nodes_per_element == 8) ? 2 : 3;
						facet_integration_points_count = (mesh.nodes_per_facet == 2) ? 1 : 2;
						break;
					}
					case shape_tetrahedron:
					{
						element_integration_points_count = (mesh.nodes_per_element == 4) ? 1 : 4;
						facet_integration_points_count = (mesh.nodes_per_facet == 3) ? 1 : 3;
						break;
					}
					case shape_hexahedron:
					{
						element_integration_points_count = (mesh.nodes_per_element == 8) ? 2 : (mesh.nodes_per_element == 20) ? 3 : 3;
						facet_integration_points_count = (mesh.nodes_per_facet == 4) ? 2 : (mesh.nodes_per_facet == 8) ? 2 : 3;
						break;
					}
					default:
					{
						Throw(File::exception_format);
					}
				}

				ShapeIntegrationRule(mesh.element_type, element_integration_points_count, element_integration_rules);
				ShapeIntegrationRule(mesh.facet_type, facet_integration_points_count, facet_integration_rules);

				// Determine node electrode
				Vector<int> node_electrode(mesh.nodes_count);
				node_electrode.Fill(0);
				for (register int i = 1; i <= potentials.size; ++i)
				{
					Potential& potential_node = potentials.entry[i];

					int n = potential_node.node_id;
					if (potential_node.electrode != 0)
					{
						node_electrode.entry[n] = potential_node.electrode;
					}
				}
 
				// Detect nodes on electrode's boundary
				Vector<int> electrode_boundary_count(electrodes_count);
				electrode_boundary_count.Fill(0);
				Vector<bool> node_electrode_in_boundary(mesh.nodes_count);
				node_electrode_in_boundary.Fill(false);
				for (int e = 1; e <= mesh.elements_count; ++e)
				{
					for (int c = 1; c <= mesh.nodes_per_element; ++c)
					{
						int n = mesh.connectivity.entry[e][c];
						int electrode_id = node_electrode.entry[n];
						if (electrode_id != 0)
						{
							if (!node_electrode_in_boundary.entry[n])
							{
								for (int d = 1; d <= mesh.nodes_per_element; ++d)
								{
									int m = mesh.connectivity.entry[e][d];
									if (node_electrode.entry[m] != electrode_id)
									{
										node_electrode_in_boundary.entry[n] = true;
										++electrode_boundary_count.entry[electrode_id];
										break;
									}
								}
							}
						}
					}
				}

				// Store electrode boundaries
				electrode_boundary.Resize(electrodes_count);
				for (int ee = 1; ee <= electrodes_count; ++ee)
				{
					electrode_boundary.entry[ee].Resize(electrode_boundary_count.entry[ee]);
				}
				for (int n = 1; n <= node_electrode_in_boundary.size; ++n)
				{
					if (node_electrode_in_boundary.entry[n])
					{
						int electrode_id = node_electrode.entry[n];
						int& i = electrode_boundary_count.entry[electrode_id];
						electrode_boundary.entry[electrode_id].entry[i] = n;
						--i;
					}
				}

				static const char* problem_name[] =
				{
					"(Undefined)",
					"Simple",
					"Capacitance matrix",
					"Sensitivity analysis"
				};

				Log(1, "Problem --------------------------------------------------------------");
				Log(1, "-Problem type:              %s", problem_name[problem_type]);
				Log(1, "-Save mesh:                 %s", save_mesh ? "yes" : "no");
				Log(1, "-Save system of equations:  %s", save_system_of_equations ? "yes" : "no");
				Log(1, "-Materials used:            %i", materials.size);
				Log(1, "-Potential conditions:      %i", potentials.size);
				Log(1, "-Electric field conditions: %i", electric_fields.size);
				Log(1, "-Source conditions:         %i", sources.size);
				switch (problem_type)
				{
					case problem_simple:
					{
						Log(1, "-Calculate potential:       %s", save_system_of_equations ? "yes" : "no");
						Log(1, "-Calculate electric field:  %s", save_system_of_equations ? "yes" : "no");
						break;
					}
					case problem_capacitance_matrix:
					{
						Log(1, "-Electrodes count:          %i", electrodes_count);
						Log(1, "-Electrodes voltage:        %i", electrodes_voltage);
						break;
					}
					case problem_sensitivity_analysis:
					{
						Log(1, "-Electrode segments:        %i", electrode_segments);
						Log(1, "-Segments per step:         %i", segments_per_step);
						Log(1, "-Result on nodes:           %s", result_on_nodes ? "yes" : "no");
						break;
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void FillKe(const int element_id, Matrix<T>& Ke) throw(Memory::Exception)
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
							Ke.entry[i][j] += weight*material.permittivity*sum*det_J;
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

			if (problem_type != problem_simple)
			{
				for (register int i = 1; i <= potentials.size; ++i)
				{
					Potential& potential_node = potentials.entry[i];

					int n = potential_node.node_id;
					int last_electrode_active = ((electrode_active + electrode_segments - 2) % electrodes_count) + 1;
					if (electrode_active < last_electrode_active)
					{
						u.entry[n] = ((potential_node.electrode >= electrode_active) && (potential_node.electrode <= last_electrode_active)) ? electrodes_voltage : 0;
					}
					else if (electrode_active > last_electrode_active)
					{
						u.entry[n] = ((potential_node.electrode >= electrode_active) || (potential_node.electrode <= last_electrode_active)) ? electrodes_voltage : 0;
					}
					else
					{
						u.entry[n] = (potential_node.electrode == electrode_active) ? electrodes_voltage : 0;;
					}
				}
			}
			else
			{
				for (register int i = 1; i <= potentials.size; ++i)
				{
					Potential& potential_node = potentials.entry[i];

					int n = potential_node.node_id;
					u.entry[n] = potential_node.potential;
				}
			}
		}


		void FillFixed(Vector<bool>& fixed) throw()
		{
			fixed.Fill(false);

			for (register int i = 1; i <= potentials.size; ++i)
			{
				Potential& potential_node = potentials.entry[i];

				int n = potential_node.node_id;
				fixed.entry[n] = true;
			}
		}


		void FillF(Vector<T>& f) throw(Memory::Exception)
		{
			try
			{
				f.Fill(0);

				// ElectricField conditions
				Vector<T> facet_N(shape_functions.mesh.nodes_per_facet);

				for (int i = 1; i <= electric_fields.size; ++i)
				{
					ElectricField& electric_field = electric_fields.entry[i];
					for (int q = 1; q <= facet_integration_rules.size; ++q)
					{
						T det_J;
						shape_functions.FacetShapeFunctions(electric_field.facet_nodes, facet_integration_rules.entry[q].point, facet_N, det_J);

						T w = facet_integration_rules.entry[q].weight;
						for (register int j = 1; j <= shape_functions.mesh.nodes_per_facet; ++j)
						{
							register int n = electric_field.facet_nodes.entry[j];
							f.entry[n] += w*electric_field.field_strength*facet_N.entry[j]*det_J;
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
							f.entry[n] += w*source.charge_density*N.entry[j]*det_J;
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void WritePotentials(File& file, const Vector<int>& node_index, const Vector<T>& u, int step) throw(Memory::Exception, File::Exception)
		{
			try
			{
				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("Result \"Potential\" \"ElectricPotential\" ");
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


		void GetElectricField(const int element_id, const Vector<T>& u, Vector<Vector<T> >& field_strength) throw(Memory::Exception)
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
						field_strength.entry[q].entry[d] = -material.permittivity*sum;
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void WriteFluxGaussPoints(File& file) throw(Memory::Exception, File::Exception)
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


		void WriteElectricField(File& file, const Vector<int>& element_index, const Vector<int>& node_index, const Vector<T>& u, int step) throw(Memory::Exception, File::Exception)
		{
			try
			{
				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("Result \"Flux\" \"ElectricPotential\" ");
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
					int node_id = node_index.entry[n];
					global_u.entry[node_id] = u.entry[n];
				}

				for (int e = 1; e <= element_index.size; ++e)
				{
					int element_id = element_index.entry[e];
					GetElectricField(element_id, global_u, flux);
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


		void GetCapacitance(const Vector<int>& node_index_inverse, const CSRMatrix<T>& Kc, const Vector<T>& u, Vector<T>& capacitance) throw()
		{
			capacitance.Fill(0);
			for (int ee = 1; ee <= electrodes_count; ++ee)
			{
				Vector<int>& boundary = electrode_boundary.entry[ee];
				for (int b = 1; b <= boundary.size; ++b)
				{
					int n = boundary.entry[b];
					int i = node_index_inverse.entry[n];
					T sum = 0;
					int k_max = Kc.Count(i);
					for (int k = 1; k <= k_max; ++k)
					{
						sum += Kc.entry[i][k]*u.entry[Kc.index[i][k]];
					}
					capacitance.entry[ee] += sum;
				}
				capacitance.entry[ee] /= electrodes_voltage;
			}
		}


		void WriteCapacitance(File& capacitances_file, const Vector<int>& node_index_inverse, const CSRMatrix<T>& Kc, const Vector<T>& u) throw(Memory::Exception, File::Exception)
		{
			try
			{
				FormatFloat format_float(false, false, 1, 10, FormatFloat::exponential);

				Vector<T> capacitance(electrodes_count);
				GetCapacitance(node_index_inverse, Kc, u, capacitance);

				for (int electrode_id = 1; electrode_id <= electrodes_count; ++electrode_id)
				{
					capacitances_file.Write(capacitance.entry[electrode_id], format_float);
					if (electrode_id < electrodes_count)
					{
						capacitances_file.Write(", ");
					}
				}
				capacitances_file.Write("\n");
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void WriteSensitivityGaussPoints(File& file) throw(Memory::Exception, File::Exception)
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

				file.Write("GaussPoints \"SensitivityGaussPoints\" ElemType ");
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


		void WriteSensitivity(File& file, const Vector<int>& element_index, const Vector<Matrix<double> >& sensitivity, int electrode_active, int electrode_test) const throw(File::Exception)
		{
			try
			{
				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("Result \"Electrode ");
				file.Write(electrode_active, format_integer);
				file.Write("\" \"Sensitivity\" ");
				file.Write(electrode_test, format_integer);
				file.Write(" Scalar OnGaussPoints \"SensitivityGaussPoints\"\n");
				file.Write("Values\n");

				for (int e = 1; e <= element_index.size; ++e)
				{
					int element_id = element_index.entry[e];
					if ((sensitivity.entry[element_id].rows > 0) && (sensitivity.entry[element_id].columns > 0))
					{
						file.Write(element_id, format_integer);
						file.Write(" ");							
						file.Write(sensitivity.entry[element_id].entry[electrode_active][electrode_test], format_float);
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


		void WriteSensitivityOnNodes(File& file, const Vector<double>& value, const Vector<int>& count, int electrode_active, int electrode_test) const throw(File::Exception)
		{
			try
			{
				FormatFloat format_float(false, false, 1, 5, FormatFloat::exponential);
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

				file.Write("Result \"Electrode ");
				file.Write(electrode_active, format_integer);
				file.Write("\" \"SensitivityOnNodes\" ");
				file.Write(electrode_test, format_integer);
				file.Write(" Scalar OnNodes\n");
				file.Write("Values\n");
				for (int n = 1; n <= value.size; ++n)
				{
					if (count.entry[n])
					{
						file.Write(n, format_integer);
						file.Write(" ");
						file.Write(value.entry[n]/count.entry[n], format_float);
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

		ElectricPotential& operator = (const ElectricPotential&) throw()
		{
			return *this;
		}
};

#endif
