// Partition.cpp
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

#include <Basic/Debug.h>
#include <Container/List.h>
#include <Container/Set.h>
#include <FiniteElement/Adjacencies.h>
#include <FiniteElement/Partition.h>

extern "C"
{
	#include <metis.h>
}


void MetisPartitionDualGraph(const Vector<idxtype>& adjacency_start, const Vector<idxtype>& adjacencies, const int partitions_count, Vector<idxtype>& element_partition_id) throw(Memory::Exception)
{
	Assert(partitions_count >= 2);

	try
	{
		int elements_count = adjacency_start.size - 1;

		// Call METIS partitioning routine for dual graph
		int wgtflag = 0; // Used to indicate if the graph is weighted
		int numflag = 1; // Used to indicate which numbering scheme is used for the adjacency structure of the graph
		int options[5] = {0, 0, 0, 0, 0};
		int edgecut; // Upon successful completion, this variable stores the number of edges that are cut by the partition.
		element_partition_id.Resize(elements_count); // Stores the partition number assigned to each element
		METIS_PartGraphKway((int*)&elements_count, adjacency_start.data, adjacencies.data, (idxtype*)0, (idxtype*)0, &wgtflag, &numflag, (int*)&partitions_count, options, &edgecut, element_partition_id.data);
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


void SimplePartitioning(const Mesh& mesh, const int partitions_count, Vector<Partition>& partitions, bool reorder) throw(Memory::Exception)
{
	Assert(partitions_count >= 2);

	try
	{
		int elements_count = mesh.elements_count;

		// Partition dual graph
		Vector<idxtype> element_partition_id;
		{
			// Generate adjacency for dual graph
			Adjacencies adjacencies;
			mesh.GetAdjacencies(adjacencies);

			// Partition dual graph
			MetisPartitionDualGraph(adjacencies.start, adjacencies.data, partitions_count, element_partition_id);
		}

		// Count elements per partition
		partitions.Resize(partitions_count);
		Vector<int> partition_elements_count(partitions_count);
		partition_elements_count.Fill(0);
		for (int e = 1; e <= elements_count; ++e)
		{
			int p = element_partition_id.entry[e];
			++partition_elements_count.entry[p];
		}

		// Fill partition element_index and node_index
		for (int p = 1; p <= partitions_count; ++p)
		{
			Vector<int>& element_index = partitions.entry[p].element_index;
			Vector<int>& node_index = partitions.entry[p].node_index;

			element_index.Resize(partition_elements_count.entry[p]);
			for (register int e = 1, index = 0; e <= elements_count; ++e)
			{
				if (element_partition_id.entry[e] == p)
				{
					element_index.entry[++index] = e;
				}
			}
			mesh.GetElementsNodes(element_index, reorder, node_index);
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


void StructurePartitioning(const Mesh& mesh, const int partitions_count, Vector<Substructure>& substructures, Boundary& boundary) throw(Memory::Exception)
{
	Assert(partitions_count >= 2);

	try
	{
		int nodes_count = mesh.nodes_count;

		{
			int elements_count = mesh.elements_count;

			// Partition dual graph
			Vector<idxtype> element_partition_id;
			{
				// Generate adjacency for dual graph
				Adjacencies adjacencies;
				mesh.GetAdjacencies(adjacencies);

				// Partition dual graph
				MetisPartitionDualGraph(adjacencies.start, adjacencies.data, partitions_count, element_partition_id);
			}

			// Count elements per partition
			substructures.Resize(partitions_count);
			Vector<int> partition_elements_count(partitions_count);
			partition_elements_count.Fill(0);
			for (int e = 1; e <= elements_count; ++e)
			{
				int p = element_partition_id.entry[e];
				++partition_elements_count.entry[p];
			}

			// Fill partition element_index and node_index
			for (int p = 1; p <= partitions_count; ++p)
			{
				Vector<int>& element_index = substructures.entry[p].element_index;

				element_index.Resize(partition_elements_count.entry[p]);
				for (register int e = 1, index = 0; e <= elements_count; ++e)
				{
					if (element_partition_id.entry[e] == p)
					{
						element_index.entry[++index] = e;
					}
				}
				mesh.GetElementsNodes(substructures.entry[p].element_index, true, substructures.entry[p].node_index);
			}
		}

		// Identify boundary nodes
		Vector<int> used_node(nodes_count);
		used_node.Fill(0);
		for (int p = 1; p <= partitions_count; ++p)
		{
			Vector<int>& node_index = substructures.entry[p].node_index;
			for (int i = 1; i <= node_index.size; ++i)
			{
				int n = node_index.entry[i];
				++used_node.entry[n];
			}
		}

		// Fill boundary.node_index and boundary.inverse_node_index
		int boundary_node_count = 0;
		for (int n = 1; n <= nodes_count; ++n)
		{
			if (used_node.entry[n] > 1)
			{
				++boundary_node_count;
			}
		}
		boundary.node_index.Resize(boundary_node_count);
		boundary.inverse_node_index.Resize(nodes_count);
		boundary.inverse_node_index.Fill(0x8FFFFFFF); // To identify invalid indexes
		for (int n = 1, b = 0; n <= nodes_count; ++n)
		{
			if (used_node.entry[n] > 1)
			{
				boundary.node_index.entry[++b] = n;
				boundary.inverse_node_index.entry[n] = b;
			}
		}

		// Fill boundary.element_index
		{
			List<int> boundary_elements;
			const Matrix<int>& connectivity = mesh.connectivity;
			for (int e = 1; e <= connectivity.rows; ++e)
			{
				for (int c = 1; c <= connectivity.columns; ++c)
				{
					int n = connectivity.entry[e][c];
					if (used_node.entry[n] > 1)
					{
						boundary_elements.AppendLast(e);
						break;
					}
				}
			}
			boundary.element_index.Resize(boundary_elements.size);
			int b = 0;
			for (ListItem<int>* item = boundary_elements.first; item; item = item->next)
			{
				boundary.element_index.entry[++b] = item->value;
			}
		}

		// Move boundary nodes to the end of node_index
		for (int p = 1; p <= partitions_count; ++p)
		{
			Vector<int>& node_index = substructures.entry[p].node_index;
			Vector<int> boundary_node_index(node_index.size);
			Vector<int> boundary(node_index.size);
			int boundary_count = 0;
			
			for (int i = 1; i <= node_index.size; ++i)
			{
				int n = node_index.entry[i];
				if (used_node.entry[n] > 1)
				{
					++boundary_count;
					boundary.entry[boundary_count] = i;
					boundary_node_index.entry[boundary_count] = n;
				}
			}
			if (boundary_count > 0)
			{
				for (int i = boundary.entry[1], j = boundary.entry[1] + 1, b = 2; j <= node_index.size; ++j)
				{
					if ((b <= boundary_count) && (j == boundary.entry[b]))
					{
						++b;
					}
					else
					{
						node_index.entry[i] = node_index.entry[j];
						++i;
					}
				}
				for (int b = 1; b <= boundary_count; ++b)
				{
					node_index.entry[node_index.size - boundary_count + b] = boundary_node_index.entry[b];
				}
			}
			substructures.entry[p].boundary_count = boundary_count;
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


void StructurePartitioning(const CSRMatrix<double>& A, const int partitions_count, Vector<Substructure>& substructures, Boundary& boundary) throw(Memory::Exception)
{
	Assert(partitions_count >= 2);

	try
	{
		int nodes_count = A.rows;
		{
			Vector<idxtype> node_partition_id;
			Vector<int> partition_nodes_count(partitions_count);
			Vector<Set<int, 32> > boundary(partitions_count);
			{
				// Create adjacency structures for METIS
				int total_edges = A.NonZero() - nodes_count;
				Vector<int> adjacency_start(nodes_count + 1);
				Vector<int> adjacency_nodes(total_edges);
				for (int n = 1, a = 1; n <= nodes_count; ++n)
				{
					adjacency_start.entry[n] = a;
					int k_max = A.Count(n);
					for (int k = 1; k <= k_max; ++k)
					{
						int ne = A.index[n][k];
						if (ne != n)
						{			
							adjacency_nodes.entry[a] = ne;
							++a;
						}
					}
				}
				adjacency_start.entry[nodes_count + 1] = total_edges + 1;

				// Call reordering routine
				int wgtflag = 0; // 0 No weights (vwgts and adjwgt are NULL)
				int numflag = 1; // Used to indicate which numbering scheme is used for the adjacency structure of the graph: 0 C-style, 1 Fortran-style
				int nparts = partitions_count;
				int options[8] = {0, 0, 0, 0, 0, 0, 0, 0};
				int edgecut;
				node_partition_id.Resize(nodes_count);
				METIS_PartGraphKway(&nodes_count, adjacency_start.data, adjacency_nodes.data, (idxtype*)0, (idxtype*)0, &wgtflag, &numflag, &nparts, options, &edgecut, node_partition_id.data);

				// Count nodes per partition and add boundary nodes
				substructures.Resize(partitions_count);
				partition_nodes_count.Fill(0);
				for (register int n = 1; n <= nodes_count; ++n)
				{
					int p = node_partition_id.entry[n];
					++partition_nodes_count.entry[p];

					int start = adjacency_start.entry[n];
					int size = adjacency_start.entry[n + 1] - start;
					for (int i = 0; i < size; ++i)
					{
						int a = adjacency_nodes.entry[start + i];
						int q = node_partition_id.entry[a];
						if (p < q)
						{
							boundary.entry[p].Append(a);
						}
					}
				}
			}

			// Fill partition node_index
			for (int p = 1; p <= partitions_count; ++p)
			{
				partition_nodes_count.entry[p] += boundary.entry[p].size;
				substructures.entry[p].node_index.Resize(partition_nodes_count.entry[p]);
			}
			partition_nodes_count.Fill(0);
			for (register int n = 1; n <= nodes_count; ++n)
			{
				int p = node_partition_id.entry[n];
				substructures.entry[p].node_index.entry[++partition_nodes_count.entry[p]] = n;
			}
			for (int p = 1; p <= partitions_count; ++p)
			{
				for (SetItem<int>* item = boundary.entry[p].first; item; item = item->next)
				{
					substructures.entry[p].node_index.entry[++partition_nodes_count.entry[p]] = item->value;
				}
			}
		}

		// Identify boundary nodes
		Vector<int> used_node(nodes_count);
		used_node.Fill(0);
		for (int p = 1; p <= partitions_count; ++p)
		{
			Vector<int>& node_index = substructures.entry[p].node_index;
			for (int i = 1; i <= node_index.size; ++i)
			{
				int n = node_index.entry[i];
				++used_node.entry[n];
			}
		}

		// Fill boundary.node_index and boundary.inverse_node_index
		int boundary_node_count = 0;
		for (int n = 1; n <= nodes_count; ++n)
		{
			if (used_node.entry[n] > 1)
			{
				++boundary_node_count;
			}
		}
		boundary.node_index.Resize(boundary_node_count);
		boundary.inverse_node_index.Resize(nodes_count);
		boundary.inverse_node_index.Fill(0x8FFFFFFF); // To identify invalid indexes
		for (int n = 1, b = 0; n <= nodes_count; ++n)
		{
			if (used_node.entry[n] > 1)
			{
				boundary.node_index.entry[++b] = n;
				boundary.inverse_node_index.entry[n] = b;
			}
		}

		// Move boundary nodes to the end of node_index
		for (int p = 1; p <= partitions_count; ++p)
		{
			Vector<int>& node_index = substructures.entry[p].node_index;
			Vector<int> boundary_node_index(node_index.size);
			Vector<int> boundary(node_index.size);
			int boundary_count = 0;
			
			for (int i = 1; i <= node_index.size; ++i)
			{
				int n = node_index.entry[i];
				if (used_node.entry[n] > 1)
				{
					++boundary_count;
					boundary.entry[boundary_count] = i;
					boundary_node_index.entry[boundary_count] = n;
				}
			}
			if (boundary_count > 0)
			{
				for (int i = boundary.entry[1], j = boundary.entry[1] + 1, b = 2; j <= node_index.size; ++j)
				{
					if ((b <= boundary_count) && (j == boundary.entry[b]))
					{
						++b;
					}
					else
					{
						node_index.entry[i] = node_index.entry[j];
						++i;
					}
				}
				for (int b = 1; b <= boundary_count; ++b)
				{
					node_index.entry[node_index.size - boundary_count + b] = boundary_node_index.entry[b];
				}
			}
			substructures.entry[p].boundary_count = boundary_count;
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}