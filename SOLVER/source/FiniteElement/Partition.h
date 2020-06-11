// Partition.h
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

#ifndef _Partition_h_
#define _Partition_h_

#include <Container/CSRMatrix.h>
#include <Container/Vector.h>
#include <FiniteElement/Mesh.h>


struct Partition
{
	Vector<int> element_index;
	Vector<int> node_index;
};


void SimplePartitioning(const Mesh& mesh, const int partitions_count, Vector<Partition>& partitions, bool reorder) throw(Memory::Exception);


struct Substructure
{
	Vector<int> element_index;
	Vector<int> node_index;
	int boundary_count;
};


struct Boundary
{
	Vector<int> element_index;
	Vector<int> node_index;
	Vector<int> inverse_node_index;
};


void StructurePartitioning(const Mesh& mesh, const int partitions_count, Vector<Substructure>& substructures, Boundary& boundary) throw(Memory::Exception);


void StructurePartitioning(const CSRMatrix<double>& A, const int partitions_count, Vector<Substructure>& substructures, Boundary& boundary) throw(Memory::Exception);

#endif
