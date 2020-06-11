// Mesh.cpp
// Copyright (C) 2015 Miguel Vargas (miguel.vargas@gmail.com)
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

#include <Basic/Debug.h>
#include <Basic/Integer.h>
#include <Basic/Log.h>
#include <Container/Array.h>
#include <Container/Set.h>
#include <Container/Stack.h>
#include <FiniteElement/Mesh.h>

extern "C"
{
	#include <metis.h>
}

#if defined _MSC_VER
	#pragma warning(disable: 4127)
#endif


Mesh::Mesh() throw()
:	element_type(shape_undefined),
	nodes_per_element(connectivity.columns),
	elements_count(connectivity.rows),
	facet_type(shape_undefined),
	nodes_per_facet(),
	nodes_count(),
	connectivity()
{
}


Mesh::Mesh(const ShapeType element_type, const int nodes_per_element, const int elements_count, const int nodes_count) throw(Memory::Exception)
:	element_type(element_type),
	nodes_per_element(connectivity.columns),
	elements_count(connectivity.rows),
	facet_type(FaceTypeMacro(element_type)),
	nodes_per_facet(NodesPerFacetMacro(element_type, nodes_per_element)),
	nodes_count(nodes_count),
	connectivity(elements_count, nodes_per_element)
{
}


Mesh::Mesh(const Mesh& mesh) throw(Memory::Exception)
:	element_type(mesh.element_type),
	nodes_per_element(connectivity.columns),
	elements_count(connectivity.rows),
	facet_type(mesh.facet_type),
	nodes_per_facet(mesh.nodes_per_facet),
	nodes_count(mesh.nodes_count),
	connectivity(mesh.connectivity)
{
}


Mesh& Mesh::operator = (const Mesh& mesh) throw()
{
	if (&mesh != this)
	{
		*(ShapeType*)&element_type = mesh.element_type;
		*(ShapeType*)&facet_type = mesh.facet_type;
		*(int*)&nodes_per_facet = mesh.nodes_per_facet;
		*(int*)&nodes_count = mesh.nodes_count;
		connectivity = mesh.connectivity;
	}
	return *this;
}


void Mesh::PrintInfo(void) throw()
{
	static const char* shape_name[] =
	{
		"(Undefined)",
		"Linear",
		"Triangle",
		"Quadrilateral",
		"Tetrahedron",
		"Hexahedron"
	};

	Log(1, "Mesh -----------------------------------------------------------------");
	Log(1, "-Nodes:             %i", nodes_count);
	Log(1, "-Elements:          %i", elements_count);
	Log(1, "-Element type:      %s", shape_name[element_type]);
	Log(1, "-Nodes per element: %i", nodes_per_element);
	Log(1, "-Facet type:        %s", shape_name[facet_type]);
	Log(1, "-Nodes per facet:   %i", nodes_per_facet);
}


void Mesh::Resize(const ShapeType element_type, const int nodes_per_element, const int elements_count, const int nodes_count) throw(Memory::Exception)
{
	try
	{
		*(ShapeType*)&this->element_type = element_type;
		*(ShapeType*)&facet_type = FaceTypeMacro(element_type);
		*(int*)&nodes_per_facet = NodesPerFacetMacro(element_type, nodes_per_element);
		*(int*)&this->nodes_count = nodes_count;
		connectivity.Resize(elements_count, nodes_per_element);
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


void Mesh::GetElementsNodes(const Vector<int>& element_index, const bool reorder, Vector<int>& node_index) const throw(Memory::Exception)
{
	try
	{
		// Identify partition nodes
		Vector<bool> node_belong(nodes_count);
		node_belong.Fill(false);

		int node_index_size = 0;
		int element_id_size = element_index.size;
		for (int pe = 1; pe <= element_id_size; ++pe)
		{
			register int e = element_index.entry[pe];
			for (register int v = 1; v <= nodes_per_element; ++v)
			{
				register int n = connectivity.entry[e][v];
				if (!node_belong.entry[n])
				{
					node_belong.entry[n] = true;
					++node_index_size;
				}
			}
		}

		// Fill node_index
		node_index.Resize(node_index_size);
		for (register int n = 1, index = 0; n <= nodes_count; ++n)
		{
			if (node_belong.entry[n])
			{
				node_index.entry[++index] = n;
			}
		}

		// Reordering nodes to improve factorizations
		if (reorder)
		{
			// Create global<->local tables
			int node_index_size = node_index.size;
			Vector<int> local_to_global(node_index_size);
			Vector<int> global_to_local(nodes_count);
			for (int n = 1; n <= node_index_size; ++n)
			{
				int gn = node_index.entry[n];
				local_to_global.entry[n] = gn;
				global_to_local.entry[gn] = n;
			}

			// Create adjacency structure of element_index (in local indexes)
			int element_id_size = element_index.size;
			int nodes_per_element = connectivity.columns;
			Vector<Set<int> > adjacency(node_index_size);
			for (int i = 1; i <= element_id_size; ++i)
			{
				int e = element_index.entry[i];
				for (int v1 = 1; v1 <= nodes_per_element; ++v1)
				{
					int n1 = global_to_local.entry[connectivity.entry[e][v1]];
					for (int v2 = v1 + 1; v2 <= nodes_per_element; ++v2)
					{
						int n2 = global_to_local.entry[connectivity.entry[e][v2]];
						adjacency.entry[n1].Append(n2);
						adjacency.entry[n2].Append(n1);
					}
				}
			}

			// Count total number of edges
			int total_edges = 0;
			for (register int n = 1; n <= node_index_size; ++n)
			{
				total_edges += adjacency.entry[n].size;
			}

			// Create adjacency structures for METIS
			Vector<int> adjacency_start(node_index_size + 1);
			Vector<int> adjacency_nodes(total_edges);
			for (int n = 1, ads = 1; n <= node_index_size; ++n)
			{
				adjacency_start.entry[n] = ads;
				for (register SetItem<int>* __restrict item = adjacency.entry[n].first; item; item = item->next, ++ads)
				{
					register int ne = item->value;
					adjacency_nodes.entry[ads] = ne;
				}
			}
			adjacency_start.entry[node_index_size + 1] = total_edges;

			// Call reordering routine
			int numflag = 1; // Used to indicate which numbering scheme is used for the adjacency structure of the graph: 0 C-style, 1 Fortran-style
			int options[8] = {0, 0, 0, 0, 0, 0, 0, 0};
			Vector<int> node_inverse_id(node_index_size);
			METIS_NodeND(&node_index_size, adjacency_start.data, adjacency_nodes.data, &numflag, options, node_index.data, node_inverse_id.data);

			for (int n = 1; n <= node_index_size; ++n)
			{
				int gn = local_to_global.entry[node_index.entry[n]];
				node_index.entry[n] = gn;
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


void Mesh::GetFacetNodes(int element_id, int facet, Vector<int>& facet_nodes) const throw()
{
	Assert(facet_nodes.size == nodes_per_facet);

	--facet;
	switch (element_type)
	{
		case shape_undefined:
		{
			Assert(false); // Unsoported elemety type
			break;
		}
		case shape_linear:
		{
			Assert(false); // Unsoported elemety type
			break;
		}
		case shape_triangle:
		{
			if (nodes_per_facet == 2)
			{
				static const int nodes[3][2] = {{1, 2}, {2, 3}, {3, 1}};
				facet_nodes.entry[1] = connectivity.entry[element_id][nodes[facet][0]];
				facet_nodes.entry[2] = connectivity.entry[element_id][nodes[facet][1]];
			}
			else // nodes_per_facet == 3
			{
				static const int nodes[3][3] = {{1, 2, 4}, {2, 3, 5}, {3, 1, 6}};
				facet_nodes.entry[1] = connectivity.entry[element_id][nodes[facet][0]];
				facet_nodes.entry[2] = connectivity.entry[element_id][nodes[facet][1]];
				facet_nodes.entry[3] = connectivity.entry[element_id][nodes[facet][2]];
			}
			break;
		}
		case shape_quadrilateral:
		{
			if (nodes_per_facet == 2)
			{
				static const int nodes[4][2] = {{1, 2}, {2, 3}, {3, 4}, {4, 1}};
				facet_nodes.entry[1] = connectivity.entry[element_id][nodes[facet][0]];
				facet_nodes.entry[2] = connectivity.entry[element_id][nodes[facet][1]];
			}
			else // nodes_per_facet == 3
			{
				static const int nodes[4][3] = {{1, 2, 5}, {2, 3, 6}, {3, 4, 7}, {4, 1, 8}};
				facet_nodes.entry[1] = connectivity.entry[element_id][nodes[facet][0]];
				facet_nodes.entry[2] = connectivity.entry[element_id][nodes[facet][1]];
				facet_nodes.entry[3] = connectivity.entry[element_id][nodes[facet][2]];
			}
			break;
		}
		case shape_tetrahedron:
		{
			if (nodes_per_facet == 3)
			{
				static const int nodes[4][3] = {{1, 2, 3}, {2, 4, 3}, {3, 4, 1}, {4, 2, 1}};
				facet_nodes.entry[1] = connectivity.entry[element_id][nodes[facet][0]];
				facet_nodes.entry[2] = connectivity.entry[element_id][nodes[facet][1]];
				facet_nodes.entry[3] = connectivity.entry[element_id][nodes[facet][2]];
			}
			else // nodes_per_facet == 6
			{
				static const int nodes[4][6] = {{1, 2, 3, 5, 6, 7}, {2, 4, 3, 9, 10, 6}, {3, 4, 1, 10, 8, 7}, {4, 2, 1, 9, 5, 8}};
				facet_nodes.entry[1] = connectivity.entry[element_id][nodes[facet][0]];
				facet_nodes.entry[2] = connectivity.entry[element_id][nodes[facet][1]];
				facet_nodes.entry[3] = connectivity.entry[element_id][nodes[facet][2]];
				facet_nodes.entry[4] = connectivity.entry[element_id][nodes[facet][3]];
				facet_nodes.entry[5] = connectivity.entry[element_id][nodes[facet][4]];
				facet_nodes.entry[6] = connectivity.entry[element_id][nodes[facet][5]];
			}
			break;
		}
		case shape_hexahedron:
		{
			if (nodes_per_facet == 4)
			{
				static const int nodes[6][4] = {{1, 2, 3, 4}, {1, 4, 8, 5}, {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 8, 4}, {5, 8, 7, 6}};
				facet_nodes.entry[1] = connectivity.entry[element_id][nodes[facet][0]];
				facet_nodes.entry[2] = connectivity.entry[element_id][nodes[facet][1]];
				facet_nodes.entry[3] = connectivity.entry[element_id][nodes[facet][2]];
				facet_nodes.entry[4] = connectivity.entry[element_id][nodes[facet][3]];
			}
			else if (nodes_per_facet == 8)
			{
				static const int nodes[6][8] = {{1, 2, 3, 4, 9, 10, 11, 12}, {1, 4, 8, 5, 12, 16, 20, 13}, {1, 5, 6, 2, 13, 17, 14, 9}, {2, 6, 7, 3, 14, 18, 15, 10}, {3, 7, 8, 4, 15, 19, 16, 11}, {5, 8, 7, 6, 20, 19, 18, 17}};
				facet_nodes.entry[1] = connectivity.entry[element_id][nodes[facet][0]];
				facet_nodes.entry[2] = connectivity.entry[element_id][nodes[facet][1]];
				facet_nodes.entry[3] = connectivity.entry[element_id][nodes[facet][2]];
				facet_nodes.entry[4] = connectivity.entry[element_id][nodes[facet][3]];
				facet_nodes.entry[5] = connectivity.entry[element_id][nodes[facet][4]];
				facet_nodes.entry[6] = connectivity.entry[element_id][nodes[facet][5]];
				facet_nodes.entry[7] = connectivity.entry[element_id][nodes[facet][6]];
				facet_nodes.entry[8] = connectivity.entry[element_id][nodes[facet][7]];
			}
			else // nodes_per_facet == 9
			{
				static const int nodes[6][9] = {{1, 2, 3, 4, 9, 10, 11, 12, 21}, {1, 4, 8, 5, 12, 16, 20, 13, 25}, {1, 5, 6, 2, 13, 17, 14, 9, 22}, {2, 6, 7, 3, 14, 18, 15, 10, 23}, {3, 7, 8, 4, 15, 19, 16, 11, 24}, {5, 8, 7, 6, 20, 19, 18, 17, 26}};
				facet_nodes.entry[1] = connectivity.entry[element_id][nodes[facet][0]];
				facet_nodes.entry[2] = connectivity.entry[element_id][nodes[facet][1]];
				facet_nodes.entry[3] = connectivity.entry[element_id][nodes[facet][2]];
				facet_nodes.entry[4] = connectivity.entry[element_id][nodes[facet][3]];
				facet_nodes.entry[5] = connectivity.entry[element_id][nodes[facet][4]];
				facet_nodes.entry[6] = connectivity.entry[element_id][nodes[facet][5]];
				facet_nodes.entry[7] = connectivity.entry[element_id][nodes[facet][6]];
				facet_nodes.entry[8] = connectivity.entry[element_id][nodes[facet][7]];
				facet_nodes.entry[9] = connectivity.entry[element_id][nodes[facet][8]];
			}
			break;
		}
	}
}


void Mesh::GetAdjacencies(Adjacencies& adjacencies) const throw(Memory::Exception)
{
	try
	{
		// Generate adjacency for dual graph
		int facets_per_element = FacetsPerElement(element_type);
		Matrix<int> adjacency_matrix(elements_count, facets_per_element);
		int total_adjacency_count = 0;
		{
			Vector<Sequence<int, 16> > contact(elements_count);
			{
				int vertex_count = (element_type == shape_triangle) ? 3 : (element_type == shape_tetrahedron) ? 4 : (element_type == shape_hexahedron) ? 8 : 4;

				// Identify elements sharing a node
				Vector<Sequence<int, 16> > node_elements(nodes_count);
				for (int e = 1; e <= elements_count; ++e)
				{
					for (register int v = 1; v <= vertex_count; ++v)
					{
						register int n = connectivity.entry[e][v];
						node_elements.entry[n].Append() = e;
					}
				}

				// List of ejements in contact
				for (int n = 1; n <= nodes_count; ++n)
				{
					for (register SequenceItem<int>* __restrict l1 = node_elements.entry[n].first; l1; l1 = l1->next)
					{
						int e1 = l1->value;
						for (register SequenceItem<int>* __restrict l2 = l1->next; l2; l2 = l2->next)
						{
							int e2 = l2->value;
							contact.entry[e1].Append() = e2;
							contact.entry[e2].Append() = e1;
						}
					}
				}
			}

			// Create temporal adjacecies in a matrix
			int facet_vertex_count = (element_type == shape_triangle) ? 2 : (element_type == shape_tetrahedron) ? 3 : (element_type == shape_hexahedron) ? 4 : 2;
			Vector<int> contact_count(elements_count);
			contact_count.Fill(0);
			for (int e = 1; e <= elements_count; ++e)
			{
				int count = 0;
				for (register SequenceItem<int>* __restrict c = contact.entry[e].first; c; c = c->next)
				{
					++contact_count.entry[c->value];
					if (contact_count.entry[c->value] == facet_vertex_count)
					{
						++count;
					}
				}
				total_adjacency_count += count;

				count = 0;
				for (register SequenceItem<int>* __restrict c = contact.entry[e].first; c; c = c->next)
				{
					if (contact_count.entry[c->value] == facet_vertex_count)
					{
						adjacency_matrix.entry[e][++count] = c->value;
					}
					contact_count.entry[c->value] = 0;
				}
				for (int j = count + 1; j <= facets_per_element; ++j)
				{
					adjacency_matrix.entry[e][j] = 0;
				}
			}
		}

		// Adjacencies vectors
		adjacencies.start.Resize(elements_count + 1);
		adjacencies.data.Resize(total_adjacency_count);
		int start = 1;
		for (int e = 1; e <= elements_count; ++e)
		{
			adjacencies.start.entry[e] = start;
			for (int j = 1; j <= facets_per_element; ++j)
			{
				register int a = adjacency_matrix.entry[e][j];
				if (a)
				{
					adjacencies.data.entry[start++] = a;
				}
			}
		}
		adjacencies.start.entry[elements_count + 1] = start;
	}
	catch (Exception&)
	{
		ReThrow();
	}
}



void Mesh::ElementColoring(Vector<Vector<int> >& element_index) const throw(Memory::Exception)
{
	// Get elements adjacent by a side or by a node
	Vector<Stack<int, 16> > adjacencies(elements_count);
	{
		int vertex_count = (element_type == shape_triangle) ? 3 : (element_type == shape_tetrahedron) ? 4 : (element_type == shape_hexahedron) ? 8 : 4;

		// Identify elements sharing a node
		Vector<Stack<int, 16> > node_elements(nodes_count);
		for (int e = 1; e <= elements_count; ++e)
		{
			for (register int v = 1; v <= vertex_count; ++v)
			{
				register int n = connectivity.entry[e][v];
				node_elements.entry[n].Push(e);
			}
		}

		// List of ejements in adjacencies
		for (int n = 1; n <= nodes_count; ++n)
		{
			Array<int, 2048> contact;
			Stack<int, 16>& elements = node_elements.entry[n];
			int size = elements.size;
			for (register int i = 1; i <= size; ++i)
			{
				int e;
				elements.Pop(e);
				contact.entry[i] = e;
			}
			for (register int i = 1; i <= size; ++i)
			{
				int e1 = contact.entry[i];
				for (register int j = i + 1; j <= size; ++j)
				{
					int e2 = contact.entry[j];
					adjacencies.entry[e1].Push(e2);
					adjacencies.entry[e2].Push(e1);
				}
			}
		}
	}

	Array<int, 256> color_count; // 256 colors should be enough (even for 3D meshes with tetrahedra)
	color_count.Fill(0);
	Array<bool, 256> color_used;
	color_used.Fill(false);

	Vector<short> elements_color(elements_count);
	elements_color.Fill(0);

	short number_of_colors = 0;
	for (int e = 1; e <= elements_count; ++e)
	{
		Stack<int, 16>& adjacency = adjacencies.entry[e];
		for (register int i = adjacency.size; i; --i)
		{
			int a;
			adjacency.Pop(a);
			register int c = elements_color.entry[a];
			if (c)
			{
				color_used.entry[c] = true;
			}
		}

		short color = 0;
		for (register short d = number_of_colors; d; --d)
		{
			register short c = ((e + d - 2) % number_of_colors) + 1;
			if (!color_used.entry[c])
			{
				color = c;
				break;
			}
		}
		if (!color)
		{
			++number_of_colors;
			color = number_of_colors;
		}

		elements_color.entry[e] = color;
		++color_count.entry[color];

		for (register short c = 1; c <= number_of_colors; ++c)
		{
			color_used.entry[c] = false;
		}
	}

	element_index.Resize(number_of_colors);
	for (int c = 1; c <= element_index.size; ++c)
	{
		Vector<int>& index = element_index.entry[c];
		index.Resize(color_count.entry[c]);
	}
	color_count.Fill(0);
	for (register int e = 1; e <= elements_count; ++e)
	{
		register int c = elements_color.entry[e];
		element_index.entry[c].entry[++color_count.entry[c]] = e;
	}
}


void Mesh::GetAdjacencies(Matrix<int>& adjacencies) const throw(Memory::Exception)
{
	try
	{
		int vertex_count = (element_type == shape_triangle) ? 3 : (element_type == shape_tetrahedron) ? 4 : (element_type == shape_hexahedron) ? 8 : 4;

		// Generate adjacencies for dual graph
		Vector<Sequence<int, 16> > contact(elements_count);
		{
			// Identify elements sharing a node
			Vector<Sequence<int, 16> > node_elements(nodes_count);
			for (int e = 1; e <= elements_count; ++e)
			{
				for (register int v = 1; v <= vertex_count; ++v)
				{
					register int n = connectivity.entry[e][v];
					node_elements.entry[n].Append() = e;
				}
			}

			// List of elements in contact
			for (int n = 1; n <= nodes_count; ++n)
			{
				for (register SequenceItem<int>* __restrict l1 = node_elements.entry[n].first; l1; l1 = l1->next)
				{
					int e1 = l1->value;
					for (register SequenceItem<int>* __restrict l2 = l1->next; l2; l2 = l2->next)
					{
						int e2 = l2->value;
						contact.entry[e1].Append() = e2;
						contact.entry[e2].Append() = e1;
					}
				}
			}
		}

		// Create adjacecies in a matrix
		int facets_per_element = FacetsPerElement(element_type);
		int facet_vertex_count = (element_type == shape_triangle) ? 2 : (element_type == shape_tetrahedron) ? 3 : (element_type == shape_hexahedron) ? 4 : 2;

		adjacencies.Resize(elements_count, facets_per_element);
		adjacencies.Fill(0);

		Vector<int> common_nodes(nodes_count);
		common_nodes.Fill(0);

		Vector<int> contact_count(elements_count);
		contact_count.Fill(0);
		for (int e = 1; e <= elements_count; ++e)
		{
			// Store all nodes in element (using local numeration)
			for (register int l = 1; l <= vertex_count; ++l)
			{
				register int n = connectivity.entry[e][l];
				common_nodes.entry[n] = l;
			}

			// Identify how many nodes are shared with each element
			for (register SequenceItem<int>* __restrict node = contact.entry[e].first; node; node = node->next)
			{
				int ce = node->value;
				++contact_count.entry[ce];
				if (contact_count.entry[ce] == facet_vertex_count)
				{
					// Identify the facet that is adjacent to the element ce
					unsigned int facet_mask = 0;
					for (register int v = 1; v <= vertex_count; ++v)
					{
						register int n = connectivity.entry[ce][v];
						register int l = common_nodes.entry[n];
						if (l)
						{
							facet_mask |= 1 << l;
						}
					}

					static const unsigned int triangle_masks[] = {0, 6, 12, 10}; // 1-2, 2-3, 3-1
					static const unsigned int quadrilateral_masks[] = {0, 6, 12, 24, 18}; // 1-2, 2-3, 3-4, 4-1
					static const unsigned int tetrahedra_masks[] = {0, 14, 28, 26, 22}; // 1-2-3, 2-4-3, 3-4-1, 4-2-1
					static const unsigned int hexahedra_masks[] = {0, 30, 306, 102, 204, 408, 480}; // 1-2-3-4, 1-4-8-5, 1-5-6-2, 2-6-7-3, 3-7-8-4, 5-8-7-6
					static const unsigned int* all_masks[] = {(unsigned int*)0, (unsigned int*)0, triangle_masks, quadrilateral_masks, tetrahedra_masks, hexahedra_masks};

					// Add element adjacencies for the corresponding facet
					const unsigned int* mask = all_masks[element_type];
					for (register int f = 1; f <= facets_per_element; ++f)
					{
						if (facet_mask == mask[f])
						{
							adjacencies.entry[e][f] = ce;
							break;
						}
					}
				}
			}

			// Clear contacts
			for (register SequenceItem<int>* __restrict node = contact.entry[e].first; node; node = node->next)
			{
				contact_count.entry[node->value] = 0;
			}

			// Clear nodes
			for (register int v = 1; v <= vertex_count; ++v)
			{
				register int n = connectivity.entry[e][v];
				common_nodes.entry[n] = 0;
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


void Mesh::ElementColoring(const Matrix<int>& adjacencies, Vector<Vector<int> >& element_index) const throw(Memory::Exception)
{
	Array<int, 256> color_count; // 256 colors should be enough (even for 3D meshes with tetrahedra)
	color_count.Fill(0);
	Array<bool, 256> color_used;
	color_used.Fill(false);

	Vector<short> elements_color(elements_count);
	elements_color.Fill(0);

	short number_of_colors = 0;
	for (int e = 1; e <= elements_count; ++e)
	{
		for (int c = 1; c <= adjacencies.columns; ++c)
		{
			int a = adjacencies.entry[e][c];
			if (a)
			{
				register int c = elements_color.entry[a];
				if (c)
				{
					color_used.entry[c] = true;
				}
			}
		}

		short color = 0;
		for (register short d = number_of_colors; d; --d)
		{
			register short c = ((e + d - 2) % number_of_colors) + 1;
			if (!color_used.entry[c])
			{
				color = c;
				break;
			}
		}
		if (!color)
		{
			++number_of_colors;
			color = number_of_colors;
		}

		elements_color.entry[e] = color;
		++color_count.entry[color];

		for (register short c = 1; c <= number_of_colors; ++c)
		{
			color_used.entry[c] = false;
		}
	}

	element_index.Resize(number_of_colors);
	for (int c = 1; c <= element_index.size; ++c)
	{
		Vector<int>& index = element_index.entry[c];
		index.Resize(color_count.entry[c]);
	}
	color_count.Fill(0);
	for (register int e = 1; e <= elements_count; ++e)
	{
		register int c = elements_color.entry[e];
		element_index.entry[c].entry[++color_count.entry[c]] = e;
	}
}


bool Mesh::CheckNodes(Sequence<int, 32>& error_nodes) const throw(Memory::Exception)
{
	Vector<bool> nodes(nodes_count);
	nodes.Fill(false);
	for (int e = 1; e <= connectivity.rows; ++e)
	{
		for (int v = 1; v <= connectivity.columns; ++v)
		{
			int n = connectivity.entry[e][v];
			nodes.entry[n] = true;
		}
	}

	error_nodes.Clear();
	for (int n = 1; n <= nodes_count; ++n)
	{
		if (!nodes.entry[n])
		{
			error_nodes.Append() = n;
		}
	}

	return error_nodes.size > 0 ? false : true;
}


void Mesh::CheckElements(Vector<Sequence<int, 32> >& group) const throw(Memory::Exception)
{
	Matrix<int> adjacencies;
	GetAdjacencies(adjacencies);

	Vector<int> element_group(connectivity.rows);
	element_group.Fill(0);

	int group_id = 1;
	int element_min = 2;
	Stack<int, 512> to_check;
	to_check.Push(1);
	do
	{
		int e = to_check.Pop();
		element_group.entry[e] = group_id;
		for (int v = 1; v <= adjacencies.columns; ++v)
		{
			int a = adjacencies.entry[e][v];
			if (a > 0)
			{
				if (element_group.entry[a] == 0)
				{
					to_check.Push(a);
					if (a == element_min)
					{
						++element_min;
					}
				}
			}
		}
		if (to_check.size == 0)
		{
			for (int e = element_min; e <= element_group.size; ++e)
			{
				if (element_group.entry[e] == 0)
				{
					++group_id;
					to_check.Push(e);
					break;
				}
			}
		}
	} while (to_check.size > 0);

	group.Resize(group_id);
	for (int e = 1; e <= element_group.size; ++e)
	{
		int g = element_group.entry[e];
		group.entry[g].Append() = e;
	}
}
