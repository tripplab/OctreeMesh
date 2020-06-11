// Geometry.h
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

#ifndef _Geometry_h_
#define _Geometry_h_

#include <Basic/Float.h>
#include <Basic/Integer.h>
#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <FiniteElement/Mesh.h>
#include <FiniteElement/Nodes.h>
#include <File/File.h>

#include <string.h>


#define GEOMETRY_SECTION_NAME_MAX_SIZE 256


template <typename T>
class Geometry
{
	public:

		Mesh mesh;
		Nodes<T> nodes;
		Vector<int> material_index;


		Geometry()
		:	mesh(),
			nodes(),
			material_index()
		{
		}


		Geometry(const char* file_name) throw(Memory::Exception, File::ExceptionOpen, File::ExceptionEOF, File::ExceptionRead, File::ExceptionFormat)
		:	mesh(),
			nodes(),
			material_index()
		{
			try
			{
				char section_name[GEOMETRY_SECTION_NAME_MAX_SIZE];

				// Load geometry
				File file;
				file.Open(file_name);
				file.SkipComments(';');

				// Load Nodes
				int dimension;
				int nodes_count;

				file.Read(section_name, GEOMETRY_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Nodes}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(dimension);
				file.SkipComments(';');
				file.Read(nodes_count);
				file.SkipComments(';');
				nodes.coordinate.Resize(nodes_count, dimension);
				for (int i = 1; i <= nodes_count; ++i)
				{
					for (int j = 1; j <= dimension; ++j)
					{
						file.Read(nodes.coordinate.entry[i][j]);
					}
				}

				// Load Mesh
				int element_type;
				int nodes_per_element;
				int elements_count;

				file.Read(section_name, GEOMETRY_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Mesh}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(element_type);
				file.SkipComments(';');
				file.Read(nodes_per_element);
				file.SkipComments(';');
				file.Read(elements_count);
				file.SkipComments(';');
				material_index.Resize(elements_count);
				mesh.Resize((ShapeType)element_type, nodes_per_element, elements_count, nodes_count);
				for (int i = 1; i <= elements_count; ++i)
				{
					file.Read(material_index.entry[i]);
					for (int j = 1; j <= nodes_per_element; ++j)
					{
						file.Read(mesh.connectivity.entry[i][j]);
					}
				}
				*(ShapeType*)&mesh.facet_type = FaceTypeMacro(element_type);
				*(int*)&mesh.nodes_per_facet = NodesPerFacetMacro(element_type, nodes_per_element);

				file.Close();
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void PrintInfo(void) throw()
		{
			nodes.PrintInfo();
			mesh.PrintInfo();
		}


		void SaveMesh(const char* file_name, const char* name, const Vector<int>& element_index, const Vector<int>& node_index, int color_id) const throw(File::ExceptionCreate, File::ExceptionWrite)
		{
			try
			{
				FormatInteger format_integer(false, false, 1, FormatInteger::decimal);
				FormatFloat format_float(false, false, 1, 4, FormatFloat::exponential);

				// Write file header
				File file;
				file.Create(file_name);
				file.Write("MESH \"");
				file.Write(name);
				file.Write("\" dimension ");
				file.Write(nodes.dimension, format_integer);
				file.Write(" ElemType ");
				file.Write((mesh.element_type == shape_triangle) ? "Triangle" : (mesh.element_type == shape_quadrilateral) ? "Quadrilateral" : (mesh.element_type == shape_tetrahedron) ? "Tetrahedra" : (mesh.element_type == shape_hexahedron) ? "Hexahedra" : "Undefined");
				file.Write(" Nnode ");
				file.Write(mesh.nodes_per_element, format_integer);
				file.Write("\n\n");

				static const int color[72][3] =
				{
					{0xE0, 0x00, 0x00}, {0x00, 0xE0, 0xC0}, {0xE0, 0x00, 0xE0}, {0xC0, 0xE0, 0x00}, {0x00, 0x00, 0xE0}, {0xE0, 0xC0, 0x00}, {0x00, 0xE0, 0xE0}, {0xC0, 0xC0, 0xC0}, {0x00, 0xE0, 0x00}, {0xC0, 0x00, 0xE0}, {0xE0, 0xE0, 0x00}, {0x00, 0xC0, 0xE0},
					{0xC0, 0x00, 0x00}, {0x00, 0xC0, 0xA0}, {0xC0, 0x00, 0xC0}, {0xA0, 0xC0, 0x00}, {0x00, 0x00, 0xC0}, {0xC0, 0xA0, 0x00}, {0x00, 0xC0, 0xC0}, {0xA0, 0xA0, 0xA0}, {0x00, 0xC0, 0x00}, {0xA0, 0x00, 0xC0}, {0xC0, 0xC0, 0x00}, {0x00, 0xA0, 0xC0},
					{0xA0, 0x00, 0x00}, {0x00, 0xA0, 0x80}, {0xA0, 0x00, 0xA0}, {0x80, 0xA0, 0x00}, {0x00, 0x00, 0xA0}, {0xA0, 0x80, 0x00}, {0x00, 0xA0, 0xA0}, {0x80, 0x80, 0x80}, {0x00, 0xA0, 0x00}, {0x80, 0x00, 0xA0}, {0xA0, 0xA0, 0x00}, {0x00, 0x80, 0xA0},
					{0x80, 0x00, 0x00}, {0x00, 0x80, 0x60}, {0x80, 0x00, 0x80}, {0x60, 0x80, 0x00}, {0x00, 0x00, 0x80}, {0x80, 0x60, 0x00}, {0x00, 0x80, 0x80}, {0x60, 0x60, 0x60}, {0x00, 0x80, 0x00}, {0x60, 0x00, 0x80}, {0x80, 0x80, 0x00}, {0x00, 0x60, 0x80},
					{0x60, 0x00, 0x00}, {0x00, 0x60, 0x40}, {0x60, 0x00, 0x60}, {0x40, 0x60, 0x00}, {0x00, 0x00, 0x60}, {0x60, 0x40, 0x00}, {0x00, 0x60, 0x60}, {0x40, 0x40, 0x40}, {0x00, 0x60, 0x00}, {0x40, 0x00, 0x60}, {0x60, 0x60, 0x00}, {0x00, 0x40, 0x60},
					{0x40, 0x00, 0x00}, {0x00, 0x40, 0x20}, {0x40, 0x00, 0x40}, {0x20, 0x40, 0x00}, {0x00, 0x00, 0x40}, {0x40, 0x20, 0x00}, {0x00, 0x40, 0x40}, {0x20, 0x20, 0x20}, {0x00, 0x40, 0x00}, {0x20, 0x00, 0x40}, {0x40, 0x40, 0x00}, {0x00, 0x20, 0x40}
				};

				// Write color
				color_id = (color_id -1) % 72;
				file.Write("# color ");
				file.Write(color[color_id][0], format_integer);
				file.Write(" ");
				file.Write(color[color_id][1], format_integer);
				file.Write(" ");
				file.Write(color[color_id][2], format_integer);
				file.Write("\n\n");

				// Write coordinates
				file.Write("Coordinates\n");
				int dimension = nodes.dimension;
				int node_id_size = node_index.size;
				for (int i = 1; i <= node_id_size; ++i)
				{
					int n = node_index.entry[i];
					file.Write(n, format_integer);
					for (int c = 1; c <= dimension; ++c)
					{
						file.Write(" ");
						file.Write(nodes.coordinate.entry[n][c], format_float);
					}
					file.Write("\n");
				}
				file.Write("End Coordinates\n");

				// Write elements
				file.Write("Elements\n");
				int nodes_per_element = mesh.nodes_per_element;
				int element_id_size = element_index.size;
				for (int i = 1; i <= element_id_size; ++i)
				{
					int e = element_index.entry[i];
					file.Write(e, format_integer);
					for (int v = 1; v <= nodes_per_element; ++v)
					{
						int n = mesh.connectivity.entry[e][v];
						file.Write(" ");
						file.Write(n, format_integer);
					}
					file.Write("\n");
				}
				file.Write("End Elements\n");

				file.Close();
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}
};

#endif
