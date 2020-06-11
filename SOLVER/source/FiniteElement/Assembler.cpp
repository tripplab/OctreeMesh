// Assembler.cpp
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

#include <FiniteElement/Assembler.h>


Assembler::Assembler(const Mesh& mesh, const Vector<int>& element_index, const Vector<int>& node_index, const int degrees_of_freedom) throw(Memory::Exception)
:	mesh(mesh),
	element_index(element_index),
	node_index(node_index),
	degrees_of_freedom(degrees_of_freedom),
	global_to_local(mesh.nodes_count)
{
	// Generate global to local node conversion vector
	global_to_local.Fill(0);
	for (register int i = 1; i <= node_index.size; ++i)
	{
		register int node_id = node_index.entry[i];
		global_to_local.entry[node_id] = i;
	}
}


Assembler& Assembler::operator = (const Assembler&) throw()
{
	return *this;
}
