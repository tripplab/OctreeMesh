// Mesh.h
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

#ifndef _Mesh_h_
#define _Mesh_h_

#include <Basic/Memory.h>
#include <Container/Matrix.h>
#include <Container/Sequence.h>
#include <Container/Vector.h>
#include <FiniteElement/Adjacencies.h>
#include <FiniteElement/Shape.h>


class Mesh
{
	public:

		const ShapeType element_type;
		const int& nodes_per_element;
		const int& elements_count;

		const ShapeType facet_type;
		const int nodes_per_facet;

		const int nodes_count;

		Matrix<int> connectivity;

		Mesh() throw();

		Mesh(const ShapeType element_type, const int nodes_per_element, const int elements_count, const int nodes_count) throw(Memory::Exception);

		Mesh(const Mesh& mesh) throw(Memory::Exception);

		Mesh& operator = (const Mesh& mesh) throw();

		void PrintInfo(void) throw();

		void Resize(const ShapeType element_type, const int nodes_per_element, const int elements_count, const int nodes_count) throw(Memory::Exception);

		void GetElementsNodes(const Vector<int>& element_index, const bool reorder, Vector<int>& node_index) const throw(Memory::Exception);

		void GetFacetNodes(int element_id, int facet, Vector<int>& facet_nodes) const throw();

		void GetAdjacencies(Adjacencies& adjacencies) const throw(Memory::Exception);

		void ElementColoring(Vector<Vector<int> >& element_index) const throw(Memory::Exception);

		void GetAdjacencies(Matrix<int>& adjacencies) const throw(Memory::Exception);

		void ElementColoring(const Matrix<int>& adjacencies, Vector<Vector<int> >& element_index) const throw(Memory::Exception);

		bool CheckNodes(Sequence<int, 32>& error_nodes) const throw(Memory::Exception);

		void CheckElements(Vector<Sequence<int, 32> >& group) const throw(Memory::Exception);
};

#endif
